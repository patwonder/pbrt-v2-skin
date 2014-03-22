
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

typedef vector<MultipoleProfileDataEntry> Profile;

static float sampleProfile(const Profile& profile, float distanceSquared, float MultipoleProfileDataEntry::* dataField) {
	if (distanceSquared <= profile[0].distanceSquared)
		return profile[0].*dataField;
	if (distanceSquared > profile.back().distanceSquared)
		return 0.f; // ensure integral convergence
	size_t lo = 1, hi = profile.size() - 1;

	// binary search to find the interval to interpolate
	while (lo < hi) {
		size_t mid = (lo + hi) / 2;
		if (distanceSquared > profile[mid].distanceSquared)
			lo = mid + 1;
		else
			hi = mid;
	}
	// lo points to the entry with minimum dsq no less than distanceSquared
	float lerpAmount = Clamp((distanceSquared - profile[lo - 1].distanceSquared) /
		(profile[lo].distanceSquared - profile[lo - 1].distanceSquared), 0.f, 1.f);
	if (isnan(lerpAmount))
		lerpAmount = 0.5;
	return Lerp(lerpAmount, profile[lo - 1].*dataField, profile[lo].*dataField);
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


Spectrum MultipoleBSSRDF::reflectance(float distanceSquared) const {
	return Spectrum::FromSampledSpectrum(pData->reflectance(distanceSquared));
}


Spectrum MultipoleBSSRDF::transmittance(float distanceSquared) const {
	return Spectrum::FromSampledSpectrum(pData->transmittance(distanceSquared));
}


float integrateProfile(const Profile& profile, float MultipoleProfileDataEntry::* dataField) {
	// Computes the integral {r,0,+inf}2*pi*r*Rd(r)dr = {r^2,0,+inf}pi*Rd(r)d(r^2) as total reflectance
	// The integrals does not include delta distributions
	float lastD2 = 0.f;
	float integral = 0.f;
	for (size_t i = 0; i < profile.size(); i++) {
		integral += M_PI * profile[i].*dataField * (profile[i].distanceSquared - lastD2);
		lastD2 = profile[i].distanceSquared;
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


Spectrum MultipoleBSSRDF::totalReflectance() const {
	return Spectrum::FromSampledSpectrum(integrateSpectralProfile(pData->spectralProfile,
		&MultipoleProfileDataEntry::reflectance));
}


Spectrum MultipoleBSSRDF::totalTransmittance() const {
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

	// Compute mfps as spectral data
	SampledSpectrum mfpTotal = Spectrum(0.f);
	for (int layer = 0; layer < layers; layer++) {
		mfpTotal += Spectrum(1.f) / (mua[layer] + musp[layer]);
	}
	SampledSpectrum mfp = mfpTotal / (float)layers;
	
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
	options.desiredStepSize = 16.f * mfp[sc] / (float)options.desiredLength;

	// Do the computation
	MPC_Output* pOutput;
	MPC_ComputeDiffusionProfile(layers, pLayerSpecs, &options, &pOutput);

	// Collect result
	Profile& profile = pData->spectralProfile[sc];
	profile.clear();
	profile.reserve(pOutput->length);

	for (uint32 ri = 0; ri < pOutput->length; ri++) {
		MultipoleProfileDataEntry entry;
		entry.distanceSquared = pOutput->pDistanceSquared[ri];
		entry.reflectance = pOutput->pReflectance[ri];
		entry.transmittance = pOutput->pTransmittance[ri];
		profile.push_back(entry);
	}

	MPC_FreeOutput(pOutput);
	delete [] pLayerSpecs;

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
}

void ReleaseMultipoleProfile(MultipoleProfileData* pData) {
	delete pData;
}
