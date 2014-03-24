
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

#include "multipole.h"
#include "multipole/MultipoleProfileCalculator/MultipoleProfileCalculator.h"
#include "parallel.h"
#include "progressreporter.h"

struct MultipoleProfileDataEntry {
	float distanceSquared;
	float reflectance;
	float transmittance;
};

struct Profile {
	static const uint32_t LUT_SIZE = 1024;

	uint32_t lut[LUT_SIZE];
	float rcpChunkLength;
	vector<MultipoleProfileDataEntry> data;
};


static void generateLUT(Profile& profile) {
	// divide extent of the profile into equal-sized chunks,
	// saving a few binary search steps and reducing cache misses
	uint32_t current = 0;
	const auto& data = profile.data;
	float chunkLength = data.back().distanceSquared / Profile::LUT_SIZE;
	profile.rcpChunkLength = 1.f / chunkLength;
	for (uint32_t i = 0; i < Profile::LUT_SIZE; i++) {
		float target = chunkLength * (i + 1);
		// find first entry out of the chunk
		while (current < data.size() - 1 && data[current].distanceSquared < target)
			current++;

		profile.lut[i] = current;
	}
}


static float sampleProfile(const Profile& profile, float distanceSquared, float MultipoleProfileDataEntry::* dataField) {
	const auto& data = profile.data;
	float extent = data.back().distanceSquared;
	if (distanceSquared > extent)
		return 0.f; // ensure integral convergence

	uint32_t chunkId = min((uint32_t)(distanceSquared * profile.rcpChunkLength), Profile::LUT_SIZE - 1);
	uint32_t lo = chunkId ? profile.lut[chunkId - 1] : 0;
	uint32_t hi = profile.lut[chunkId];

	// binary search to find the interval to interpolate
	while (lo < hi) {
		uint32_t mid = (lo + hi) / 2;
		if (distanceSquared > data[mid].distanceSquared)
			lo = mid + 1;
		else
			hi = mid;
	}
	// lo points to the entry with minimum dsq no less than distanceSquared
	if (lo) {
		float lerpAmount = Clamp((distanceSquared - data[lo - 1].distanceSquared) /
			(data[lo].distanceSquared - data[lo - 1].distanceSquared), 0.f, 1.f);
		if (isnan(lerpAmount))
			lerpAmount = 0.5;
		return Lerp(lerpAmount, data[lo - 1].*dataField, data[lo].*dataField);
	} else {
		return data[0].*dataField;
	}
}

static SampledSpectrum sampleSpectralProfile(const Profile spectralProfile[nSpectralSamples], float distanceSquared,
	float MultipoleProfileDataEntry::* dataField)
{
	SampledSpectrum result;
	for (int i = 0; i < nSpectralSamples; i++) {
		result[i] = sampleProfile(spectralProfile[i], distanceSquared, dataField);
	}
	return result;
}

struct MultipoleProfileData {
	Profile spectralProfile[nSpectralSamples];

	SampledSpectrum reflectance(float distanceSquared) const {
		return sampleSpectralProfile(spectralProfile, distanceSquared, &MultipoleProfileDataEntry::reflectance);
	}
	SampledSpectrum transmittance(float distanceSquared) const {
		return sampleSpectralProfile(spectralProfile, distanceSquared, &MultipoleProfileDataEntry::transmittance);
	}
};


Spectrum MultipoleBSSRDFData::reflectance(float distanceSquared) const {
	return Spectrum::FromSampledSpectrum(pData->reflectance(distanceSquared));
}


Spectrum MultipoleBSSRDFData::transmittance(float distanceSquared) const {
	return Spectrum::FromSampledSpectrum(pData->transmittance(distanceSquared));
}


float integrateProfile(const Profile& profile, float MultipoleProfileDataEntry::* dataField) {
	// Computes the integral {r,0,+inf}2*pi*r*Rd(r)dr = {r^2,0,+inf}pi*Rd(r)d(r^2) as total reflectance
	// The integrals does not include delta distributions
	const auto& data = profile.data;
	float lastD2 = 0.f;
	float integral = 0.f;
	for (size_t i = 0; i < data.size(); i++) {
		integral += M_PI * data[i].*dataField * (data[i].distanceSquared - lastD2);
		lastD2 = data[i].distanceSquared;
	}
	return integral;
}


SampledSpectrum integrateSpectralProfile(const Profile profile[nSpectralSamples], float MultipoleProfileDataEntry::* dataField) {
	Spectrum result;
	for (int i = 0; i < nSpectralSamples; i++) {
		result[i] = integrateProfile(profile[i], dataField);
	}
	return result;
}


Spectrum MultipoleBSSRDFData::totalReflectance() const {
	return Spectrum::FromSampledSpectrum(integrateSpectralProfile(pData->spectralProfile,
		&MultipoleProfileDataEntry::reflectance));
}


Spectrum MultipoleBSSRDFData::totalTransmittance() const {
	return Spectrum::FromSampledSpectrum(integrateSpectralProfile(pData->spectralProfile,
		&MultipoleProfileDataEntry::transmittance));
}


class MultipoleProfileTask : public Task {
public:
	MultipoleProfileTask(int spectralChannel, int layers, const SampledSpectrum mua[],
		const SampledSpectrum musp[], float et[], float thickness[],
		MultipoleProfileData* pData, ProgressReporter& pr) : reporter(pr)
	{
		sc = spectralChannel;
		this->layers = layers;
		this->mua = mua;
		this->musp = musp;
		this->et = et;
		this->thickness = thickness;
		this->pData = pData;
	}

	void Run() override;
private:
	int sc;
	int layers;
	const SampledSpectrum* mua;
	const SampledSpectrum* musp;
	float* et;
	float* thickness;
	MultipoleProfileData* pData;
	ProgressReporter& reporter;
};

void MultipoleProfileTask::Run() {
	MPC_LayerSpec* pLayerSpecs = new MPC_LayerSpec[layers];
	MPC_Options options;
	options.desiredLength = 512;

	// Compute mfp
	float mfpTotal = 0.f;
	for (int layer = 0; layer < layers; layer++) {
		mfpTotal += 1.f / (mua[layer][sc] + musp[layer][sc]);
	}
	float mfp = mfpTotal / (float)layers;
	
	// Fill in layer information
	for (int layer = 0; layer < layers; layer++) {
		MPC_LayerSpec& ls = pLayerSpecs[layer];
		ls.g_HG = 0.f;
		ls.ior = et[layer];
		ls.mua = mua[layer][sc];
		ls.musp = musp[layer][sc];
		ls.thickness = thickness[layer];
	}
	// Compute desired step size
	options.desiredStepSize = 16.f * mfp / (float)options.desiredLength;

	// Do the computation
	MPC_Output* pOutput;
	MPC_ComputeDiffusionProfile(layers, pLayerSpecs, &options, &pOutput);

	// Collect result
	Profile& profile = pData->spectralProfile[sc];
	auto& data = profile.data;
	data.clear();
	data.reserve(pOutput->length);

	for (uint32 ri = 0; ri < pOutput->length; ri++) {
		MultipoleProfileDataEntry entry;
		entry.distanceSquared = pOutput->pDistanceSquared[ri];
		entry.reflectance = pOutput->pReflectance[ri];
		entry.transmittance = pOutput->pTransmittance[ri];
		data.push_back(entry);
	}

	MPC_FreeOutput(pOutput);
	delete [] pLayerSpecs;

	// Generate lut for profile to accelerate lookup
	generateLUT(profile);

	reporter.Update();
}


void ComputeMultipoleProfile(int layers, const SampledSpectrum mua[], const SampledSpectrum musp[], float et[], float thickness[],
	MultipoleProfileData** oppData)
{
	if (!oppData) {
		Error("Cannot output multipole profile.");
	}

	MultipoleProfileData* pData = *oppData = new MultipoleProfileData;

	int nTasks = nSpectralSamples;
    ProgressReporter reporter(nTasks, "Profile");
	vector<Task*> profileTasks;
	profileTasks.reserve(nTasks);
	for (int i = 0; i < nTasks; ++i) {
		profileTasks.push_back(new MultipoleProfileTask(i, layers, mua,
			musp, et, thickness, pData, reporter));
	}
	EnqueueTasks(profileTasks);
	WaitForAllTasks();
	for (Task* task : profileTasks)
		delete task;
	reporter.Done();

	Spectrum tr = Spectrum::FromSampledSpectrum(integrateSpectralProfile(pData->spectralProfile,
		&MultipoleProfileDataEntry::reflectance));
	Spectrum tt = Spectrum::FromSampledSpectrum(integrateSpectralProfile(pData->spectralProfile,
		&MultipoleProfileDataEntry::transmittance));
	Info("Total Reflectance: %s", tr.ToString().c_str());
	Info("Total Transmittance: %s", tt.ToString().c_str());
}

void ReleaseMultipoleProfile(MultipoleProfileData* pData) {
	delete pData;
}
