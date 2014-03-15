
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// integrators/dipolesubsurface.cpp*
#include "stdafx.h"
#include "integrators/dipolesubsurface.h"
#include "scene.h"
#include "montecarlo.h"
#include "sampler.h"
#include "progressreporter.h"
#include "intersection.h"
#include "paramset.h"
#include "reflection.h"
#include "octree.h"
#include "camera.h"
#include "floatfile.h"
#include "diffusionutil.h"

// DipoleSubsurfaceIntegrator Local Declarations


// DipoleSubsurfaceIntegrator Method Definitions
DipoleSubsurfaceIntegrator::~DipoleSubsurfaceIntegrator() {
    delete[] lightSampleOffsets;
    delete[] bsdfSampleOffsets;
}


void DipoleSubsurfaceIntegrator::RequestSamples(Sampler *sampler, Sample *sample,
        const Scene *scene) {
    // Allocate and request samples for sampling all lights
    uint32_t nLights = scene->lights.size();
    lightSampleOffsets = new LightSampleOffsets[nLights];
    bsdfSampleOffsets = new BSDFSampleOffsets[nLights];
    for (uint32_t i = 0; i < nLights; ++i) {
        const Light *light = scene->lights[i];
        int nSamples = light->nSamples;
        if (sampler) nSamples = sampler->RoundSize(nSamples);
        lightSampleOffsets[i] = LightSampleOffsets(nSamples, sample);
        bsdfSampleOffsets[i] = BSDFSampleOffsets(nSamples, sample);
    }
}


void DipoleSubsurfaceIntegrator::Preprocess(const Scene *scene,
        const Camera *camera, const Renderer *renderer) {
    if (scene->lights.size() == 0) return;
    vector<SurfacePoint> pts;
    // Get _SurfacePoint_s for translucent objects in scene
    if (filename != "") {
        // Initialize _SurfacePoint_s from file
        vector<float> fpts;
        if (ReadFloatFile(filename.c_str(), &fpts)) {
            if ((fpts.size() % 8) != 0)
                Error("Excess values (%d) in points file \"%s\"", int(fpts.size() % 8),
                      filename.c_str());
            for (u_int i = 0; i < fpts.size(); i += 8)
                pts.push_back(SurfacePoint(Point(fpts[i], fpts[i+1], fpts[i+2]),
                                           Normal(fpts[i+3], fpts[i+4], fpts[i+5]),
                                           fpts[i+6], fpts[i+7]));
        }
    }
    if (pts.size() == 0) {
        Point pCamera = camera->CameraToWorld(camera->shutterOpen,
                                              Point(0, 0, 0));
        FindPoissonPointDistribution(pCamera, camera->shutterOpen,
                                     minSampleDist, scene, &pts);
    }

    // Compute irradiance values at sample points
    RNG rng;
    MemoryArena arena;
    PBRT_SUBSURFACE_STARTED_COMPUTING_IRRADIANCE_VALUES();
    ProgressReporter progress(pts.size(), "Computing Irradiances");
    for (uint32_t i = 0; i < pts.size(); ++i) {
        SurfacePoint &sp = pts[i];
        Spectrum E(0.f);
        for (uint32_t j = 0; j < scene->lights.size(); ++j) {
            // Add irradiance from light at point
            const Light *light = scene->lights[j];
            Spectrum Elight = 0.f;
            int nSamples = RoundUpPow2(light->nSamples);
            uint32_t scramble[2] = { rng.RandomUInt(), rng.RandomUInt() };
            uint32_t compScramble = rng.RandomUInt();
            for (int s = 0; s < nSamples; ++s) {
                float lpos[2];
                Sample02(s, scramble, lpos);
                float lcomp = VanDerCorput(s, compScramble);
                LightSample ls(lpos[0], lpos[1], lcomp);
                Vector wi;
                float lightPdf;
                VisibilityTester visibility;
                Spectrum Li = light->Sample_L(sp.p, sp.rayEpsilon,
                    ls, camera->shutterOpen, &wi, &lightPdf, &visibility);
                if (Dot(wi, sp.n) <= 0.) continue;
                if (Li.IsBlack() || lightPdf == 0.f) continue;
                Li *= visibility.Transmittance(scene, renderer, NULL, rng, arena);
                if (visibility.Unoccluded(scene))
                    Elight += Li * AbsDot(wi, sp.n) / lightPdf;
            }
            E += Elight / nSamples;
        }
        irradiancePoints.push_back(IrradiancePoint(sp, E));
        PBRT_SUBSURFACE_COMPUTED_IRRADIANCE_AT_POINT(&sp, &E);
        arena.FreeAll();
        progress.Update();
    }
    progress.Done();
    PBRT_SUBSURFACE_FINISHED_COMPUTING_IRRADIANCE_VALUES();

    // Create octree of clustered irradiance samples
    octree = octreeArena.Alloc<SubsurfaceOctreeNode>();
    for (uint32_t i = 0; i < irradiancePoints.size(); ++i)
        octreeBounds = Union(octreeBounds, irradiancePoints[i].p);
    for (uint32_t i = 0; i < irradiancePoints.size(); ++i)
        octree->Insert(octreeBounds, &irradiancePoints[i], octreeArena);
    octree->InitHierarchy();
}


Spectrum DipoleSubsurfaceIntegrator::Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena) const {
    Spectrum L(0.);
    Vector wo = -ray.d;
    // Compute emitted light if ray hit an area light source
    L += isect.Le(wo);

    // Evaluate BSDF at hit point
    BSDF *bsdf = isect.GetBSDF(ray, arena);
    const Point &p = bsdf->dgShading.p;
    const Normal &n = bsdf->dgShading.nn;
    // Evaluate BSSRDF and possibly compute subsurface scattering
    BSSRDF *bssrdf = isect.GetBSSRDF(ray, arena);
    if (bssrdf && octree) {
        Spectrum sigma_a  = bssrdf->sigma_a();
        Spectrum sigmap_s = bssrdf->sigma_prime_s();
        Spectrum sigmap_t = sigmap_s + sigma_a;
        if (!sigmap_t.IsBlack()) {
            // Use hierarchical integration to evaluate reflection from dipole model
            PBRT_SUBSURFACE_STARTED_OCTREE_LOOKUP(const_cast<Point *>(&p));
            DiffusionReflectance Rd(sigma_a, sigmap_s, bssrdf->eta());
            Spectrum Mo = octree->Mo(octreeBounds, p, Rd, maxError);
            FresnelDielectric fresnel(1.f, bssrdf->eta());
            Spectrum Ft = Spectrum(1.f) - fresnel.Evaluate(AbsDot(wo, n));
            float Fdt = 1.f - Fdr(bssrdf->eta());
            L += (INV_PI * Ft) * (Fdt * Mo);
            PBRT_SUBSURFACE_FINISHED_OCTREE_LOOKUP();
        }
    }
    L += UniformSampleAllLights(scene, renderer, arena, p, n,
        wo, isect.rayEpsilon, ray.time, bsdf, sample, rng, lightSampleOffsets,
        bsdfSampleOffsets);
    if (ray.depth < maxSpecularDepth) {
        // Trace rays for specular reflection and refraction
        L += SpecularReflect(ray, bsdf, rng, isect, renderer, scene, sample,
                             arena);
        L += SpecularTransmit(ray, bsdf, rng, isect, renderer, scene, sample,
                              arena);
    }
    return L;
}


DipoleSubsurfaceIntegrator *CreateDipoleSubsurfaceIntegrator(const ParamSet &params) {
    int maxDepth = params.FindOneInt("maxdepth", 5);
    float maxError = params.FindOneFloat("maxerror", .05f);
    float minDist = params.FindOneFloat("minsampledistance", .25f);
    string pointsfile = params.FindOneString("pointsfile", "");
    if (PbrtOptions.quickRender) { maxError *= 4.f; minDist *= 4.f; }
    return new DipoleSubsurfaceIntegrator(maxDepth, maxError, minDist, pointsfile);
}


