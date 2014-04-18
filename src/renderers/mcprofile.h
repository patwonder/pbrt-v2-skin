
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

#pragma once

#include "pbrt.h"
#include "renderer.h"
#include "layer.h"

struct MCProfileEntry {
	MCProfileEntry() {
		reflectance = transmittance = 0.;
	}
	MCProfileEntry(double r, double t)
		: reflectance(r), transmittance(t)
	{
	}
	double reflectance;
	double transmittance;
};

typedef vector<MCProfileEntry> MCProfile;

struct MCProfileResult {
	double totalMCReflectance;
	double totalMCTransmittance;
	double totalNoLerpReflectance;
	double totalNoLerpTransmittance;
	double totalLerpReflectance;
	double totalLerpTransmittance;

	MCProfileResult() {
		totalMCReflectance = 0.;
		totalMCTransmittance = 0.;
		totalNoLerpReflectance = 0.;
		totalNoLerpTransmittance = 0.;
		totalLerpReflectance = 0.;
		totalLerpTransmittance = 0.;
	}
};

class MonteCarloProfileRenderer : public Renderer {
public:
	// MonteCarloProfileRenderer Public Methods
	MonteCarloProfileRenderer(const Layer* layers, int nLayers, float mfpRange,
		int segments, uint64_t photons, string filename,
		bool noClearCache = false, bool silent = false, bool noCompare = false)
		: layers(layers, layers + nLayers), mfpRange(mfpRange), nSegments(segments),
		  nPhotons(photons), filename(filename), noClearCache(noClearCache), silent(silent),
		  noCompare(noCompare)
	{
	}
	void Render(const Scene* scene) override;
	Spectrum Li(const Scene* scene, const RayDifferential& ray,
		const Sample* sample, RNG& rng, MemoryArena& arena,
		Intersection* isect, Spectrum* T) const override;
	Spectrum Transmittance(const Scene* scene, const RayDifferential& ray,
		const Sample* sample, RNG& rng, MemoryArena& arena) const override;
	const MCProfileResult& GetResult() const { return result; }
	const MCProfile& GetProfile() const { return resultProfile; }

	static void ClearCache();
private:
	// MonteCarloProfileRenderer Private Data
	vector<Layer> layers;
	float mfpRange;
	int nSegments;
	string filename;
	uint64_t nPhotons;
	bool noClearCache;
	bool silent;
	bool noCompare;
	
	MCProfileResult result;
	MCProfile resultProfile;

	// MonteCarloProfileRenderer Private Methods
};

MonteCarloProfileRenderer* CreateMonteCarloProfileRenderer(const ParamSet& params);
