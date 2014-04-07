
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
#include "skinlayer.h"

typedef vector<float> ParamRange;

class SkinCoefficients;
struct VariableParams;
class GaussianFitTask;
struct SpectralGaussianCoeffs {
	vector<SampledSpectrum> coeffs;
	SampledSpectrum error;
};

class MultipoleProfileFitRenderer : public Renderer {
public:
	// MultipoleProfileFitRenderer Public Methods
	MultipoleProfileFitRenderer(const vector<SkinLayer>& layers,
		ParamRange f_mel, ParamRange f_eu, ParamRange f_blood,
		ParamRange f_ohg, uint32_t desiredLength, uint32_t segments, const string& filename,
		int idRange[2], int splitTasks, const string& pbrtFilePath);
	void Render(const Scene* scene) override;
	Spectrum Li(const Scene* scene, const RayDifferential& ray,
		const Sample* sample, RNG& rng, MemoryArena& arena,
		Intersection* isect, Spectrum* T) const override;
	Spectrum Transmittance(const Scene* scene, const RayDifferential& ray,
		const Sample* sample, RNG& rng, MemoryArena& arena) const override;
private:
	// MultipoleProfileFitRenderer Private Data
	vector<SkinLayer> layers;
	ParamRange pr_f_mel;
	ParamRange pr_f_eu;
	ParamRange pr_f_blood;
	ParamRange pr_f_ohg;
	uint32_t desiredLength;
	uint32_t nSegments;
	string filename;
	int idRange[2];
	int splitTasks;
	string pbrtFilePath;

	// MultipoleProfileFitRenderer Private Methods
	template <class Processor>
	void ForRanges(const Processor& p) const;

	static float NextParamFromId(uint32_t& id, ParamRange pr);
	VariableParams ParamsFromId(uint32_t id) const;

	vector<Task*> CreateGaussianFitTasks(const SkinCoefficients& coeffs,
		const vector<float>& sigmas, SpectralGaussianCoeffs& sgc,
		ProgressReporter& reporter, uint32_t id) const;

	int GetTotalNumTasks() const;

	void DoProfileFit() const;
	void DoPbrtFileSplit() const;
};

MultipoleProfileFitRenderer* CreateMultipoleProfileFitRenderer(const ParamSet& params, const string& pbrtFilePath);
