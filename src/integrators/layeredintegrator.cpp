
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

#define USE_RANDOMWALK_PROBES

#ifdef USE_RANDOMWALK_PROBES
#include <fstream>
#include <iomanip>
#include "parallel.h"
using namespace std;
void volatile_memset(volatile void* dst, int32_t val, size_t size) {
	for (size_t i = 0; i < size / sizeof(val); i++) {
		((volatile decltype(val)*)dst)[i] = val;
	}
}
#endif

class RandomWalkProbes {
public:
	RandomWalkProbes() {
#ifdef USE_RANDOMWALK_PROBES
#define CLEAR_ARRAY(arr) volatile_memset(arr, 0, sizeof(arr))
		CLEAR_ARRAY(spectralRayCount);
		CLEAR_ARRAY(discardedRayCount);
		CLEAR_ARRAY(contributedRayCount);
		volatile_memset(minBounces, 0x7fffffff, sizeof(minBounces));
		volatile_memset(maxBounces, 0x80000000, sizeof(minBounces));
		CLEAR_ARRAY(bouncesDiscarded);
		CLEAR_ARRAY(bouncesContributed);
		CLEAR_ARRAY(mfpDiscarded);
		CLEAR_ARRAY(mfpContributed);
		CLEAR_ARRAY(discardedRayCountInLayers);
#endif
	}

	~RandomWalkProbes() {
#ifdef USE_RANDOMWALK_PROBES
		bool hasProbes = false;
		for (int i = 0; i < nSpectralSamples; i++) {
			if (spectralRayCount[i]) {
				hasProbes = true;
				break;
			}
		}
		if (!hasProbes)
			return;

		ofstream out(L"randomwalk-probes.txt", ios::out | ios::app);
		time_t tm = time(NULL);
		out << ctime(&tm) << endl;
		out << "Wave Length\tTotal\tContributed\tDiscarded\tMin Bounces\tMax Bounces\tBounces discarded\tBounces Contributed"
			"\tMFP discarded\tMFP contributed\tLayer 0\tLayer 1\tLayer 2" << endl;
		for (int i = 0; i < nSpectralSamples; i++) {
			out << SampledSpectrum::WaveLength(i);
			out << "\t" << spectralRayCount[i];
			out << "\t" << contributedRayCount[i];
			out << "\t" << discardedRayCount[i];
			out << "\t" << minBounces[i];
			out << "\t" << maxBounces[i];
			out << "\t" << bouncesDiscarded[i];
			out << "\t" << bouncesContributed[i];
			out << "\t" << mfpDiscarded[i];
			out << "\t" << mfpContributed[i];
			out << "\t" << discardedRayCountInLayers[i][0];
			out << "\t" << discardedRayCountInLayers[i][1];
			out << "\t" << discardedRayCountInLayers[i][2];
			out << endl;
		}
		out.close();
#endif
	}

	void spectralRayTraced(int wlIndex) {
#ifdef USE_RANDOMWALK_PROBES
		AtomicAdd(spectralRayCount + wlIndex, 1);
#endif
	}

	void spectralRayOut(int wlIndex, int bounces, int layer, float L, float mfp) {
#ifdef USE_RANDOMWALK_PROBES
		if (layer < 0 || L > 0.f) {
			AtomicAdd(contributedRayCount + wlIndex, 1);
			AtomicAdd(bouncesContributed + wlIndex, bounces);
			AtomicAdd(mfpContributed + wlIndex, mfp);
		} else {
			AtomicAdd(discardedRayCount + wlIndex, 1);
			AtomicAdd(bouncesDiscarded + wlIndex, bounces);
			AtomicAdd(mfpDiscarded + wlIndex, mfp);
			int l = max(0, min(layer, 2));
			AtomicAdd(discardedRayCountInLayers[wlIndex] + l, 1);
		}
		AtomicMin(minBounces + wlIndex, bounces);
		AtomicMax(maxBounces + wlIndex, bounces);
#endif
	}
private:
#ifdef USE_RANDOMWALK_PROBES
	AtomicInt64 spectralRayCount[nSpectralSamples];
	AtomicInt64 discardedRayCount[nSpectralSamples];
	AtomicInt64 contributedRayCount[nSpectralSamples];

	AtomicInt32 minBounces[nSpectralSamples], maxBounces[nSpectralSamples];
	AtomicInt64 bouncesContributed[nSpectralSamples], bouncesDiscarded[nSpectralSamples];

	volatile double mfpDiscarded[nSpectralSamples], mfpContributed[nSpectralSamples];

	AtomicInt64 discardedRayCountInLayers[nSpectralSamples][3];
#endif
};


static RandomWalkProbes probes;


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

		ray = RayDifferential(p, wi, ray, isectp->rayEpsilon);
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
		bool randomWalked = false;
		if (!specularBounce && Dot(ray.d, isectp->dg.nn) < 0) {
			// Ray enters the material. See if it's layered
			const LayeredGeometricPrimitive* lprim = isectp->primitive->ToLayered();
			if (lprim) {
				// Modify pathThroughput & ray intersection to reflect subsurface scattering
				RayDifferential outRay;
				Intersection outIsect;
				L += RandomWalk(wlIndex, pathThroughput, lprim, scene, renderer,
								ray, isectp->primitiveId, sample, rng, arena, &outRay, &outIsect);
				if (pathThroughput == 0.f)
					break;

				ray = outRay;
				localIsect = outIsect;
				isectp = &localIsect;
				randomWalked = true;
			}
		}

		if (!randomWalked) {
			// Find next vertex of path
			if (!scene->IntersectExcept(ray, &localIsect, isectp->primitiveId)) {
				if (specularBounce)
					for (uint32_t i = 0; i < scene->lights.size(); ++i)
						L += pathThroughput * scene->lights[i]->Le(ray)[wlIndex];
				break;
			}
			pathThroughput *= renderer->Transmittance(scene, ray, NULL, rng, arena)[wlIndex];
			isectp = &localIsect;
		}

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

		ray = RayDifferential(p, wi, ray, isectp->rayEpsilon);
	}
	return L;
}


float LayeredIntegrator::RandomWalk(uint32_t wlIndex, float& pathThroughput,
									const LayeredGeometricPrimitive* lprim,
									const Scene *scene,	const Renderer *renderer,
									const RayDifferential& r, uint32_t primitiveId,
									const Sample *sample, RNG& rng, MemoryArena& arena,
									RayDifferential* outRay, Intersection* outIsect) const
{
	probes.spectralRayTraced(wlIndex);

	float lambda = SampledSpectrum::WaveLength(wlIndex);
	float L = 0.f;
	int currentLayer = 0;
	int targetLayer = 0;
	RayDifferential ray(r);
	float numMFPa = 0.f;
	Intersection isect;
	int bounces = 0;
	bool hit = true;
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
		
			// Hit test, using rayEpsilon for double ensurance
			ray.mint = min(ray.mint, min(lp.thickness, mfps) / 10.f);
			if (hit)
				ray.mint = max(ray.mint, lp.thickness / 100.f);
			ray.maxt = mfps;

			int layerIndex;
			if (hit = lprim->IntersectInternal(ray, lastPrimitiveId, &isect, &layerIndex)) {
				BSDF* bsdf = isect.GetBSDF(ray, arena);
				const Point& p = bsdf->dgShading.p;
				const Normal& n = bsdf->dgShading.nn;
				// Absorption
				float distance = (p - ray.o).Length();
				pathThroughput *= expf(-mua * distance);
				numMFPa += mua * distance;

				Vector wo = -ray.d;
				if (layerIndex == 0) {
					// Sample illumination from lights to find path contribution
					L += pathThroughput *
						 UniformSampleOneLight(scene, renderer, arena, p, n, wo,
							isect.rayEpsilon, ray.time, bsdf, sample, rng)[wlIndex];
				}

				// Sample BSDF to get outgoing direction
				Vector wi;
				BSDFSample outgoingBSDFSample(rng);
				float pdf;
				BxDFType flags;
				float f = bsdf->Sample_f(wo, &wi, outgoingBSDFSample, &pdf,
										 BSDF_ALL, &flags)[wlIndex];
				if (f == 0. || pdf == 0.) {
					probes.spectralRayOut(wlIndex, bounces + 1, currentLayer, L, numMFPa);
					pathThroughput = 0.f;
					return L;
				}


				// Update results
				pathThroughput *= f * (AbsDot(wi, n) / pdf);
				ray = RayDifferential(p, wi, ray, isect.rayEpsilon);
				lastPrimitiveId = isect.primitiveId;

				// Determine the target traveling layer
				targetLayer = (Dot(wi, isect.dg.nn) < 0.f) ? layerIndex : layerIndex - 1;
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
				lastPrimitiveId = 0;
			}

			// Possibly terminate the path
			bounces++;
			if (numMFPa > 3.f || bounces >= 5) {
				float continueProbability = min(.5f, pathThroughput * 10.f);
				if (rng.RandomFloat() > continueProbability) {
					probes.spectralRayOut(wlIndex, bounces, currentLayer, L, numMFPa);
					pathThroughput = 0.f;
					return L;
				}
				pathThroughput /= continueProbability;
			}
		} // Inner while loop

		currentLayer = targetLayer;
	} // Outer while loop

	probes.spectralRayOut(wlIndex, bounces, currentLayer, L, numMFPa);
	*outRay = ray;
	*outIsect = isect;
	return L;
}

LayeredIntegrator *CreateLayeredSurfaceIntegrator(const ParamSet &params) {
    int maxDepth = params.FindOneInt("maxdepth", 5);
    return new LayeredIntegrator(maxDepth);
}
