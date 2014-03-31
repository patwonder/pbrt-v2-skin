
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

#include "mcprofile.h"
#include "paramset.h"
#include "parallel.h"
#include "progressreporter.h"


void MonteCarloProfileRenderer::Render(const Scene* scene) {
	Info("Monte Carlo Profile Renderer");
	for (int i = 0; i < layers.size(); i++) {
		const Layer& l = layers[i];
		Info("Layer %d: mua=%f, musp=%f, ior=%f, thickness=%f",
			i, l.mua, l.musp, l.ior, l.thickness);
	}
	Info("Gather samples within %f mean-free-paths, divided into %d segments.",
		mfpRange, nSegments);
	Info("Output filename: %s", filename.c_str());
}

Spectrum MonteCarloProfileRenderer::Li(const Scene* scene, const RayDifferential& ray,
	const Sample* sample, RNG& rng, MemoryArena& arena,
	Intersection* isect, Spectrum* T) const
{
	return Spectrum(0.f);
}

Spectrum MonteCarloProfileRenderer::Transmittance(const Scene* scene, const RayDifferential& ray,
	const Sample* sample, RNG& rng, MemoryArena& arena) const
{
	return Spectrum(0.f);
}

MonteCarloProfileRenderer* CreateMonteCarloProfileRenderer(const ParamSet& params) {
	int nLayers;
	const Layer* layers = params.FindLayer("layers", &nLayers);
	if (!layers || nLayers < 1) {
		Error("No layers param set for MCProfileRenderer.");
		abort();
	}
	float mfpRange = params.FindOneFloat("mfprange", 16.f);
	int segments = params.FindOneInt("segments", 1024);
	string filename = params.FindOneString("filename", "");

	return new MonteCarloProfileRenderer(layers, nLayers, mfpRange, segments, filename);
}
