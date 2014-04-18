
/*
    Copyright(c) 2013-2014 Yifan Wu.

    This file is part of fork of pbrt (pbrt-v2-skin).

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

#include "StdAfx.h"

#include "multipolesubsurface.h"
#include "paramset.h"
#include "floatfile.h"
#include "progressreporter.h"
#include "diffusionutil.h"
#include "parallel.h"

// IrradianceTask Declarations
class IrradianceTask : public Task {
public:
	IrradianceTask(const Scene* sc, const Renderer* ren,
		const Camera* c, ProgressReporter& pr,
		const vector<SurfacePoint>& pts,
		vector<IrradiancePoint>& irpts, Mutex& mt,
		int tn, int tc)
		: reporter(pr), surfacePoints(pts), irradiancePoints(irpts),
		  mutex(mt)
	{
		scene = sc; renderer = ren; camera = c;
		taskNum = tn; taskCount = tc;
	}
	void Run() override;
private:
	const Scene* scene;
	const Renderer* renderer;
	const Camera* camera;
	ProgressReporter& reporter;
	const vector<SurfacePoint>& surfacePoints;
	vector<IrradiancePoint>& irradiancePoints;
	Mutex& mutex;
	int taskNum;
	int taskCount;
};


void IrradianceTask::Run() {
	size_t idxStart = (int)((uint64_t)taskNum * surfacePoints.size() / taskCount);
	size_t idxEnd = (int)((uint64_t)(taskNum + 1) * surfacePoints.size() / taskCount);

	if (idxStart == idxEnd) {
		reporter.Update();
		return;
	}

    RNG rng(taskNum * 47);
    MemoryArena arena;

	// Store locally and transfer to the integrator at once
	vector<IrradiancePoint> localIrradiancePoints;
	localIrradiancePoints.reserve(idxEnd - idxStart);

    for (size_t i = idxStart; i < idxEnd; ++i) {
        const SurfacePoint &sp = surfacePoints[i];
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
        localIrradiancePoints.push_back(IrradiancePoint(sp, E));
        PBRT_SUBSURFACE_COMPUTED_IRRADIANCE_AT_POINT(&sp, &E);
        arena.FreeAll();
    }
	{
		MutexLock lock(mutex);
		irradiancePoints.insert(irradiancePoints.end(),
			localIrradiancePoints.begin(), localIrradiancePoints.end());
	}
	reporter.Update();
}

void MultipoleSubsurfaceIntegrator::RequestSamples(Sampler *sampler, Sample *sample,
	const Scene *scene)
{
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

void MultipoleSubsurfaceIntegrator::Preprocess(const Scene *scene, const Camera *camera,
	const Renderer *renderer)
{
	if (scene->lights.size() == 0) return;
	vector<SurfacePoint> pts;
	// Get _SurfacePoint_s for translucent objects in scene
	if (filename != "") {
        // Initialize _SurfacePoint_s from file
		ReadBinaryFile(filename.c_str(), &pts);
		Info("%d surface points read from file \"%s\"",
			pts.size(), filename.c_str());
	}
	if (pts.size() == 0) {
		GetSurfacePointsThroughTessellation(camera->shutterOpen,
			minSampleDist, scene, &originalPrimitives, &pts);
	}

    // Compute irradiance values at sample points
    // Create and launch _IrradianceTask_s for rendering image
    PBRT_SUBSURFACE_STARTED_COMPUTING_IRRADIANCE_VALUES();
	int nPoints = pts.size();
	irradiancePoints.reserve(nPoints);
	int nTasks = max(32 * NumSystemCores(), nPoints / 1024);
	nTasks = RoundUpPow2(nTasks);
    ProgressReporter reporter(nTasks, "Computing Irradiances");
	Mutex* mutex = Mutex::Create();
	vector<Task*> renderTasks;
	renderTasks.reserve(nTasks);
	for (int i = 0; i < nTasks; ++i) {
		renderTasks.push_back(new IrradianceTask(scene, renderer,
		camera, reporter, pts, irradiancePoints, *mutex, i, nTasks));
	}
	EnqueueTasks(renderTasks);
	WaitForAllTasks();
	for (Task* task : renderTasks)
		delete task;
	Mutex::Destroy(mutex);
	reporter.Done();
    PBRT_SUBSURFACE_FINISHED_COMPUTING_IRRADIANCE_VALUES();

    // Create octree of clustered irradiance samples
    octree = octreeArena.Alloc<SubsurfaceOctreeNode>();
    for (uint32_t i = 0; i < irradiancePoints.size(); ++i)
        octreeBounds = Union(octreeBounds, irradiancePoints[i].p);
    for (uint32_t i = 0; i < irradiancePoints.size(); ++i)
        octree->Insert(octreeBounds, &irradiancePoints[i], octreeArena);
    octree->InitHierarchy();
}


class MultipoleReflectance {
public:
	MultipoleReflectance(const MultipoleBSSRDF* bssrdf) {
		this->bssrdf = bssrdf;
	}
	Spectrum operator()(float d2) const {
		return bssrdf->reflectance(d2);
	}
private:
	const MultipoleBSSRDF* bssrdf;
};

Spectrum MultipoleSubsurfaceIntegrator::Li(const Scene *scene, const Renderer *renderer,
	const RayDifferential &ray, const Intersection &isect,
	const Sample *sample, RNG &rng, MemoryArena &arena) const
{
    Spectrum L(0.);
    Vector wo = -ray.d;
    // Compute emitted light if ray hit an area light source
    L += isect.Le(wo);

    // Evaluate BSDF at hit point
    BSDF *bsdf = isect.GetBSDF(ray, arena);
    const Point &p = bsdf->dgShading.p;
    const Normal &n = bsdf->dgShading.nn;
    // Evaluate MultipoleBSSRDF and possibly compute subsurface scattering
    const MultipoleBSSRDF *bssrdf = isect.GetMultipoleBSSRDF(ray, arena);
    if (bssrdf && octree) {
        // Use hierarchical integration to evaluate reflection from dipole model
        PBRT_SUBSURFACE_STARTED_OCTREE_LOOKUP(const_cast<Point *>(&p));
		MultipoleReflectance mr(bssrdf);
        Spectrum Mo = octree->Mo(octreeBounds, p, n, mr, maxError);
        FresnelDielectric fresnel(1.f, bssrdf->eta(0));
        Spectrum Ft = Spectrum(1.f) - fresnel.Evaluate(AbsDot(wo, n));
		// Bypass outgoing fresnel term if it's a Monte-Carlo profile,
		// because the term is already included in the profile
        float Fdt = bssrdf->IsMonteCarlo() ? 1.f : (1.f - Fdr(bssrdf->eta(0)));
		L += ((INV_PI * Ft) * (Fdt * Mo) * bssrdf->albedo()).Clamp(0.f);
        PBRT_SUBSURFACE_FINISHED_OCTREE_LOOKUP();
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

MultipoleSubsurfaceIntegrator *CreateMultipoleSubsurfaceIntegrator(const ParamSet &params,
	const vector<Reference<Primitive> >* originalPrimitives)
{
    int maxDepth = params.FindOneInt("maxdepth", 5);
    float maxError = params.FindOneFloat("maxerror", .05f);
    float minDist = params.FindOneFloat("minsampledistance", .25f);
    string pointsfile = params.FindOneString("pointsfile", "");
    if (PbrtOptions.quickRender) { maxError *= 4.f; minDist *= 4.f; }
    return new MultipoleSubsurfaceIntegrator(maxDepth, maxError, minDist, pointsfile,
		originalPrimitives);
}
