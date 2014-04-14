
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

#include "batchmcprofile.h"
#include "mcprofile.h"
#include "paramset.h"
#include "progressreporter.h"
#include <fstream>
#include <algorithm>

using namespace std;

void BatchMCProfileRenderer::Render(const Scene* scene) {
	Info("Batch Monte-Carlo Profile Renderer");

	vector<pair<float, MCProfileResult> > results;

	ProgressReporter reporter(nThicknessSegments + 1, "Batch MC Profile");

	vector<Layer> layers = this->layers;
	// No need to parallelize since MCProfileRenderer already does a good job on that
	for (int i = 0; i <= nThicknessSegments; i++) {
		int ii = (i & 1) ? (nThicknessSegments - i / 2) : i / 2;
		float thickness = thicknessMin + (float)ii * (thicknessMax - thicknessMin) / nThicknessSegments;
		layers[0].thickness = thickness;

		MonteCarloProfileRenderer mcprofilerenderer(&layers[0], (int)layers.size(),
			mfpRange, nSegments, nPhotons, "", true, true);

		mcprofilerenderer.Render(scene);
		results.push_back(make_pair(thickness, mcprofilerenderer.GetResult()));

		reporter.Update();
	}
	sort(results.begin(), results.end(),
		[] (const pair<float, MCProfileResult>& p1, const pair<float, MCProfileResult>& p2) {
			return p1.first < p2.first;
		}
	);

	MonteCarloProfileRenderer::ClearCache();
	reporter.Done();

	Info("Writing to file...");

	// Output result
	if (filename != "") {
		ofstream out(filename, ios::out | ios::trunc);
		if (!out) {
			Error("Failed to open output file: %s", filename.c_str());
			abort();
		}

		// Header
		out << "Thickness"
			   "\tMonte-Carlo Reflectance\tMultipole Reflectance\tLerped Reflectance"
			   "\tMonte-Carlo Transmittance\tMultipole Transmittance\tLerped Transmittance" << endl;
		// Body
		for (const auto& pair : results) {
			const auto& result = pair.second;
			out << pair.first
				<< "\t" << result.totalMCReflectance << "\t" << result.totalNoLerpReflectance
				<< "\t" << result.totalLerpReflectance
				<< "\t" << result.totalMCTransmittance << "\t" << result.totalNoLerpTransmittance
				<< "\t" << result.totalLerpTransmittance << endl;
		}

		out.close();
	}
}

Spectrum BatchMCProfileRenderer::Li(const Scene* scene, const RayDifferential& ray,
	const Sample* sample, RNG& rng, MemoryArena& arena,
	Intersection* isect, Spectrum* T) const
{
	return Spectrum(0.f);
}

Spectrum BatchMCProfileRenderer::Transmittance(const Scene* scene, const RayDifferential& ray,
	const Sample* sample, RNG& rng, MemoryArena& arena) const
{
	return Spectrum(0.f);
}

BatchMCProfileRenderer* CreateBatchMCProfileRenderer(const ParamSet& params) {
	int nLayers;
	const Layer* layers = params.FindLayer("layers", &nLayers);
	if (!layers || nLayers < 1) {
		Error("No layers param set for BatchMCProfileRenderer.");
		abort();
	}
	float mfpRange = params.FindOneFloat("mfprange", 16.f);
	int segments = params.FindOneInt("segments", 1024);

	int nThicknessRange;
	const float* thicknessRange = params.FindFloat("thicknessrange", &nThicknessRange);
	if (!thicknessRange || nThicknessRange != 2 || thicknessRange[0] < 0.f || thicknessRange[1] < thicknessRange[0]) {
		Error("Thickness range incorrectly set for BatchMCProfileRenderer.");
		abort();
	}
	float thicknessMin = thicknessRange[0];
	float thicknessMax = thicknessRange[1];
	int nThicknessSegments = params.FindOneInt("thicknesssegments", 256);

	string strPhotons = params.FindOneString("photons", "100");
	uint64_t photons = _strtoui64(strPhotons.c_str(), NULL, 10);
	string filename = params.FindOneFilename("filename", "");

	return new BatchMCProfileRenderer(layers, nLayers, mfpRange, segments,
		photons, thicknessMin, thicknessMax, nThicknessSegments, filename);
}
