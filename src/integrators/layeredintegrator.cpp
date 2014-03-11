
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

#include "stdafx.h"

#include "layeredintegrator.h"
#include "paramset.h"

void LayeredIntegrator::Preprocess(const Scene *scene, const Camera *camera,
								   const Renderer *renderer)
{
	// TODO
}

void LayeredIntegrator::RequestSamples(Sampler *sampler, Sample *sample,
									   const Scene *scene)
{
	for (uint32_t wlIndex = 0; wlIndex < nSpectralSamples; wlIndex++)
		for (int i = 0; i < SAMPLE_DEPTH; ++i) {
			lightSampleOffsets[wlIndex][i] = LightSampleOffsets(1, sample);
			lightNumOffset[wlIndex][i] = sample->Add1D(1);
			bsdfSampleOffsets[wlIndex][i] = BSDFSampleOffsets(1, sample);
			pathSampleOffsets[wlIndex][i] = BSDFSampleOffsets(1, sample);
		}
}

Spectrum LayeredIntegrator::Li(const Scene *scene, const Renderer *renderer,
							   const RayDifferential &r, const Intersection &isect,
							   const Sample *sample, RNG &rng, MemoryArena &arena) const
{
	Spectrum pathThroughput = 1., L = 0.;
	RayDifferential ray(r);
	bool specularBounce = false;
    Intersection localIsect;
    const Intersection *isectp = &isect;
	for (int bounces = 0; ; ++bounces) {
		// Possibly add emitted light at path vertex
		if (bounces == 0 || specularBounce)
			L += pathThroughput * isectp->Le(-ray.d);

		// Sample illumination from lights to find path contribution
        BSDF *bsdf = isectp->GetBSDF(ray, arena);
        const Point &p = bsdf->dgShading.p;
        const Normal &n = bsdf->dgShading.nn;
        Vector wo = -ray.d;
        if (bounces < SAMPLE_DEPTH)
            L += pathThroughput *
                 UniformSampleOneLight(scene, renderer, arena, p, n, wo,
                     isectp->rayEpsilon, ray.time, bsdf, sample, rng,
                     lightNumOffset[0][bounces], &lightSampleOffsets[0][bounces],
                     &bsdfSampleOffsets[0][bounces]);
        else
            L += pathThroughput *
                 UniformSampleOneLight(scene, renderer, arena, p, n, wo,
                     isectp->rayEpsilon, ray.time, bsdf, sample, rng);

        if (bounces == maxDepth)
            break;

        // Sample BSDF to get new path direction

        // Get _outgoingBSDFSample_ for sampling new path direction
        BSDFSample outgoingBSDFSample;
        if (bounces < SAMPLE_DEPTH)
            outgoingBSDFSample = BSDFSample(sample, pathSampleOffsets[0][bounces],
                                            0);
        else
            outgoingBSDFSample = BSDFSample(rng);
        Vector wi;
        float pdf;
        BxDFType flags;
        Spectrum f = bsdf->Sample_f(wo, &wi, outgoingBSDFSample, &pdf,
                                    BSDF_ALL, &flags);
        if (f.IsBlack() || pdf == 0.)
            break;
        specularBounce = (flags & BSDF_SPECULAR) != 0;
        pathThroughput *= f * (AbsDot(wi, n) / pdf);

        // Possibly terminate the path
        if (bounces > TERM_DEPTH) {
            float continueProbability = min(.5f, pathThroughput.y());
            if (rng.RandomFloat() > continueProbability)
                break;
            pathThroughput /= continueProbability;
        }

        ray = RayDifferential(p, wi, ray, 0.f);
		if (!specularBounce && Dot(ray.d, isectp->dg.nn) < 0) {
			// Ray enters the material. See if it's layered
			const LayeredGeometricPrimitive* lprim = isectp->primitive->ToLayered();
			if (lprim) {
				SampledSpectrum sL = L.ToSampledSpectrum();
				SampledSpectrum sPathThroughput = pathThroughput.ToSampledSpectrum();
				// Switch to spectral shading
				for (int i = 0; i < SampledSpectrum::NumComponents(); i++) {
					sL[i] += LiSpectral(i, sPathThroughput[i], bounces, specularBounce,
						scene, renderer, ray, *isectp, sample, rng, arena);
				}
				return Spectrum::FromSampledSpectrum(sL);
			}
		}

		// Find next vertex of path
		if (!scene->IntersectExcept(ray, &localIsect, isectp->primitiveId)) {
            if (specularBounce)
                for (uint32_t i = 0; i < scene->lights.size(); ++i)
                   L += pathThroughput * scene->lights[i]->Le(ray);
            break;
        }
        pathThroughput *= renderer->Transmittance(scene, ray, NULL, rng, arena);
        isectp = &localIsect;
	}
	return L;
}

float LayeredIntegrator::LiSpectral(uint32_t wlIndex, float pathThroughput, int bounces, bool specularBounce,
									const Scene *scene,	const Renderer *renderer,
									const RayDifferential &r, const Intersection &isect,
									const Sample *sample, RNG &rng, MemoryArena &arena) const
{
	static bool t = false;
	if (!t) {
		t = true;
		Warning("Spectral path tracing being used. Rendering time may increase dramatically.");
	}
	float L = 0.;
	RayDifferential ray(r);
	Intersection localIsect;
	const Intersection *isectp = &isect;
	while (true) {
		if (!specularBounce && Dot(ray.d, isectp->dg.nn) < 0) {
			// Ray enters the material. See if it's layered
			const LayeredGeometricPrimitive* lprim = isectp->primitive->ToLayered();
			if (lprim) {
				// Modify pathThroughput & ray intersection to reflect subsurface scattering
				RayDifferential outRay;
				uint32_t outPrimitiveId;
				pathThroughput = RandomWalk(wlIndex, pathThroughput,
					lprim, ray, isectp->primitiveId, rng, arena, &outRay, &outPrimitiveId);
				if (pathThroughput == 0.f)
					break;

				ray = outRay;
				localIsect.primitiveId = outPrimitiveId;
				isectp = &localIsect;
			}
		}

		// Find next vertex of path
		if (!scene->IntersectExcept(ray, &localIsect, isectp->primitiveId)) {
			if (specularBounce)
				for (uint32_t i = 0; i < scene->lights.size(); ++i)
					L += pathThroughput * scene->lights[i]->Le(ray)[wlIndex];
			break;
		}
		pathThroughput *= renderer->Transmittance(scene, ray, NULL, rng, arena)[wlIndex];
		isectp = &localIsect;

		bounces++;
		// End of original loop ================================

		// Possibly add emitted light at path vertex
		if (specularBounce)
			L += pathThroughput * isectp->Le(-ray.d)[wlIndex];

		// Sample illumination from lights to find path contribution
		BSDF *bsdf = isectp->GetBSDF(ray, arena);
		const Point &p = bsdf->dgShading.p;
		const Normal &n = bsdf->dgShading.nn;
		Vector wo = -ray.d;
        if (bounces < SAMPLE_DEPTH)
            L += pathThroughput *
                 UniformSampleOneLight(scene, renderer, arena, p, n, wo,
                     isectp->rayEpsilon, ray.time, bsdf, sample, rng,
                     lightNumOffset[wlIndex][bounces], &lightSampleOffsets[wlIndex][bounces],
                     &bsdfSampleOffsets[wlIndex][bounces])[wlIndex];
        else
            L += pathThroughput *
                 UniformSampleOneLight(scene, renderer, arena, p, n, wo,
                     isectp->rayEpsilon, ray.time, bsdf, sample, rng)[wlIndex];

		if (bounces == maxDepth)
			break;

		// Sample BSDF to get new path direction

		// Get _outgoingBSDFSample_ for sampling new path direction
		BSDFSample outgoingBSDFSample;
		if (bounces < SAMPLE_DEPTH)
			outgoingBSDFSample = BSDFSample(sample, pathSampleOffsets[wlIndex][bounces],
											0);
		else
			outgoingBSDFSample = BSDFSample(rng);
		Vector wi;
		float pdf;
		BxDFType flags;
		float f = bsdf->Sample_f(wo, &wi, outgoingBSDFSample, &pdf,
								 BSDF_ALL, &flags)[wlIndex];
		if (f == 0. || pdf == 0.)
			break;
		specularBounce = (flags & BSDF_SPECULAR) != 0;
		pathThroughput *= f * (AbsDot(wi, n) / pdf);

		// Possibly terminate the path
		if (bounces > TERM_DEPTH) {
			float continueProbability = min(.5f, pathThroughput);
			if (rng.RandomFloat() > continueProbability)
				break;
			pathThroughput /= continueProbability;
		}

        ray = RayDifferential(p, wi, ray, 0.f);
	}
	return L;
}


float LayeredIntegrator::RandomWalk(uint32_t wlIndex, float pathThroughput,
									const LayeredGeometricPrimitive* lprim,
									const RayDifferential& r, uint32_t primitiveId,
									RNG& rng, MemoryArena& arena,
									RayDifferential* outRay, uint32_t* outPrimitiveId) const
{
	float lambda = SampledSpectrum::WaveLength(wlIndex);

	int currentLayer = 0;
	int targetLayer = 0;
	RayDifferential ray(r);
	float numMFPa = 0.f;
	Intersection isect;
	uint32_t lastPrimitiveId = primitiveId;
	const LayeredMaterial* pmat = static_cast<const LayeredMaterial*>(lprim->GetMaterial());

	while (currentLayer >= 0) {
		// Fetch params for current layer
		LayerParam lp = pmat->GetLayerParam(currentLayer);
		float mua = lp.mua.getValueForWL(lambda);
		float musp = lp.musp.getValueForWL(lambda);
		float ga = lp.ga;

		while (currentLayer == targetLayer) {
			// Sample mfp(s) according to the exponential distribution: musp * exp(-mfp * musp)
			float mfps = -logf(rng.RandomFloat()) / musp;
		
			// Hit test
			ray.maxt = mfps;

			int layerIndex;
			if (lprim->IntersectInternal(ray, lastPrimitiveId, &isect, &layerIndex)) {
				// hit something, sample BSDF to get outgoing direction
				BSDF* bsdf = isect.GetBSDF(ray, arena);
				const Point& p = bsdf->dgShading.p;
				const Normal& n = bsdf->dgShading.nn;
				Vector wo = -ray.d;
				Vector wi;
				BSDFSample outgoingBSDFSample(rng);
				float pdf;
				BxDFType flags;
				float f = bsdf->Sample_f(wo, &wi, outgoingBSDFSample, &pdf,
										 BSDF_ALL, &flags)[wlIndex];
				if (f == 0. || pdf == 0.)
					return 0.f;

				// Absorption
				float distance = (p - ray.o).Length();
				pathThroughput *= exp(-mua * distance);
				numMFPa += mua * distance;

				// Determine the target traveling layer
				targetLayer = (Dot(ray.d, isect.dg.nn) < 0.f) ? layerIndex : layerIndex - 1;

				if (targetLayer < 0) {
					break;
				}

				// Update results
				pathThroughput *= f * (AbsDot(wi, n) / pdf);
				ray = RayDifferential(p, wi, ray, 0.f);
				lastPrimitiveId = isect.primitiveId;
			} else {
				// Absorption
				float distance = mfps;
				pathThroughput *= exp(-mua * distance);
				numMFPa += mua * distance;

				// Doesn't hit anything, scatter into a new direction
				Point p = ray(mfps);
				Vector wi = SampleHG(ray.d, ga, rng.RandomFloat(), rng.RandomFloat());
				ray = RayDifferential(p, wi, ray, 0.f);

				// Allow re-intersecting the same primitive due to scattering
				lastPrimitiveId = *outPrimitiveId = 0;
			}

			// Possibly terminate the path
			if (numMFPa > 3.f) {
				float continueProbability = min(.5f, pathThroughput);
				if (rng.RandomFloat() > continueProbability)
					return 0.f;
				pathThroughput /= continueProbability;
			}
		}

		currentLayer = targetLayer;
	}

	*outRay = ray;
	outRay->maxt = FLT_MAX;
	*outPrimitiveId = lastPrimitiveId;
	return pathThroughput;
}

LayeredIntegrator *CreateLayeredSurfaceIntegrator(const ParamSet &params) {
    int maxDepth = params.FindOneInt("maxdepth", 5);
    return new LayeredIntegrator(maxDepth);
}
