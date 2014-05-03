
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
#include "renderers/mcprofile.h"
#include "kahansum.h"


struct MultipoleProfileDataEntry {
	// omitted, assuming uniform distribution of dsq
	//float distanceSquared;
	float reflectance;
#ifdef SAMPLE_TRANSMITTANCE
	float transmittance;
#endif
};

struct Profile {
	vector<MultipoleProfileDataEntry> data;
	float dsqSpacing;
	float rcpDsqSpacing;
	float totalReflectance;
	float totalTransmittance;
};


static float sampleProfile(const Profile& profile, float distanceSquared, float MultipoleProfileDataEntry::* dataField) {
	const auto& data = profile.data;

	double fSegId = distanceSquared * profile.rcpDsqSpacing;
	// Use double to compare - prevent overflowing the uint32_t
	if (fSegId >= (double)(data.size() - 1))
		return 0.f; // ensure integral convergence

	uint32_t segId = (uint32_t)fSegId;

	// segId points to the segment containing dsq
	float lerpAmount = (float)(fSegId - (double)segId);
	return Lerp(lerpAmount, data[segId].*dataField, data[segId + 1].*dataField);
}

static SampledSpectrum sampleSpectralProfile(const Profile* spectralProfile, float distanceSquared,
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
#ifdef SAMPLE_TRANSMITTANCE
		return sampleSpectralProfile(spectralProfile, distanceSquared, &MultipoleProfileDataEntry::transmittance);
#else
		Severe("Transmittance profile sampling not implemented.");
		return 0.f;
#endif
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
#if 0
	const auto& data = profile.data;
	float integral = 0.f;
	for (size_t i = 0; i < data.size(); i++) {
		integral += M_PI * data[i].*dataField;
	}
	return integral * profile.dsqSpacing;
#else
	// The easy way...
	if (dataField == &MultipoleProfileDataEntry::reflectance) {
		return profile.totalReflectance;
	}
	return profile.totalTransmittance;
#endif
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
#ifdef SAMPLE_TRANSMITTANCE
	return Spectrum::FromSampledSpectrum(integrateSpectralProfile(pData->spectralProfile,
		&MultipoleProfileDataEntry::transmittance));
#else
	Severe("Transmittance profile sampling not implemented.");
	return 0.f;
#endif
}


class MultipoleProfileTask : public Task {
public:
	MultipoleProfileTask(int spectralChannel, int layers, const SampledSpectrum mua[],
		const SampledSpectrum musp[], float et[], float thickness[],
		MultipoleProfileData* pData, ProgressReporter& pr, bool lerpOnThinSlab)
		: reporter(pr)
	{
		sc = spectralChannel;
		this->layers = layers;
		this->mua = mua;
		this->musp = musp;
		this->et = et;
		this->thickness = thickness;
		this->pData = pData;
		this->lerpOnThinSlab = lerpOnThinSlab;
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
	bool lerpOnThinSlab;
};

void MultipoleProfileTask::Run() {
	MPC_LayerSpec* pLayerSpecs = new MPC_LayerSpec[layers];
	MPC_Options options;
	options.desiredLength = 512;
	options.lerpOnThinSlab = lerpOnThinSlab;

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
	options.desiredStepSize = 12.f * mfp / (float)options.desiredLength;

	// Do the computation
	MPC_Output* pOutput;
	MPC_ComputeDiffusionProfile(layers, pLayerSpecs, &options, &pOutput);
	MPC_ResampleForUniformDistanceSquaredDistribution(pOutput);

	// Collect result
	Profile& profile = pData->spectralProfile[sc];
	uint32 length = pOutput->length;
	profile.dsqSpacing = pOutput->pDistanceSquared[length - 1] / (float)(length - 1);
	profile.rcpDsqSpacing = (float)(length - 1) / pOutput->pDistanceSquared[length - 1];
	profile.totalReflectance = pOutput->totalReflectance;
	profile.totalTransmittance = pOutput->totalTransmittance;
	auto& data = profile.data;
	data.clear();
	data.reserve(length);

	for (uint32 ri = 0; ri < length; ri++) {
		MultipoleProfileDataEntry entry;
		entry.reflectance = pOutput->pReflectance[ri];
#ifdef SAMPLE_TRANSMITTANCE
		entry.transmittance = pOutput->pTransmittance[ri];
#endif
		data.push_back(entry);
	}

	MPC_FreeOutput(pOutput);
	delete [] pLayerSpecs;

	reporter.Update();
}


static void ComputeMonteCarloProfile(int nLayers, const SampledSpectrum mua[], const SampledSpectrum musp[], float et[], float thickness[],
	uint64_t nPhotons, MultipoleProfileData* pData)
{
	ProgressReporter reporter(nSpectralSamples, "MC Profile");

	// Prepare layers
	vector<Layer> layers(nLayers);
	for (int i = 0; i < nLayers; i++) {
		layers[i].ior = et[i];
		layers[i].thickness = thickness[i];
	}

	// Compute spectral MC profile
	for (int sc = 0; sc < nSpectralSamples; sc++) {
		// Compute mfp
		double mfpTotal = 0.;
		for (int i = 0; i < nLayers; i++) {
			layers[i].mua = mua[i][sc];
			layers[i].musp = musp[i][sc];
			mfpTotal += 1. / (mua[i][sc] + musp[i][sc]);
		}
		double mfp = mfpTotal / (double)layers.size();

		// Do random walk
		const uint32 nSegments = 4096;
		const double mfpRange = 12.0f;
		double extent = mfpRange * mfp;
		MonteCarloProfileRenderer renderer(&layers[0], nLayers, mfpRange, nSegments, nPhotons, "", true, true, true);
		renderer.Render(NULL);
		MCProfile mcprofile = renderer.GetProfile();

		// Convert to MPC_Output*
		MPC_Output* pOutput;
		MPC_CreateOutput(nSegments, &pOutput);
		for (int i = 0; i < nSegments; i++) {
			pOutput->pDistanceSquared[i] = (float)(((i + 0.5) * extent / (double)nSegments)
										 * ((i + 0.5) * extent / (double)nSegments));
			pOutput->pReflectance[i] = mcprofile[i].reflectance;
			pOutput->pTransmittance[i] = mcprofile[i].transmittance;
		}
		pOutput->totalReflectance = renderer.GetResult().totalMCReflectance;
		pOutput->totalTransmittance = renderer.GetResult().totalMCTransmittance;
		MPC_ResampleForUniformDistanceSquaredDistribution(pOutput, nSegments * 16);

		// Collect result
		Profile& profile = pData->spectralProfile[sc];
		uint32 length = pOutput->length;
		profile.dsqSpacing = pOutput->pDistanceSquared[length - 1] / (float)(length - 1);
		profile.rcpDsqSpacing = (float)(length - 1) / pOutput->pDistanceSquared[length - 1];
		profile.totalReflectance = pOutput->totalReflectance;
		profile.totalTransmittance = pOutput->totalTransmittance;
		auto& data = profile.data;
		data.clear();
		data.reserve(length);

		for (uint32 ri = 0; ri < length; ri++) {
			MultipoleProfileDataEntry entry;
			entry.reflectance = pOutput->pReflectance[ri];
#ifdef SAMPLE_TRANSMITTANCE
			entry.transmittance = pOutput->pTransmittance[ri];
#endif
			data.push_back(entry);
		}

		MPC_FreeOutput(pOutput);
		reporter.Update();
	}

	MonteCarloProfileRenderer::ClearCache();
	reporter.Done();
}


void ComputeMultipoleProfile(int layers, const SampledSpectrum mua[], const SampledSpectrum musp[], float et[], float thickness[],
	MultipoleProfileData** oppData, bool useMonteCarlo, bool lerpOnThinSlab, uint64_t nPhotons)
{
	if (!oppData) {
		Error("Cannot output multipole profile.");
	}

	MultipoleProfileData* pData = *oppData = new MultipoleProfileData;

	if (useMonteCarlo) {
		ComputeMonteCarloProfile(layers, mua, musp, et, thickness, nPhotons, pData);
	} else {
		int nTasks = nSpectralSamples;
		ProgressReporter reporter(nTasks, "Profile");
		vector<Task*> profileTasks;
		profileTasks.reserve(nTasks);
		for (int i = 0; i < nTasks; ++i) {
			profileTasks.push_back(new MultipoleProfileTask(i, layers, mua,
				musp, et, thickness, pData, reporter, lerpOnThinSlab));
		}
		EnqueueTasks(profileTasks);
		WaitForAllTasks();
		for (Task* task : profileTasks)
			delete task;
		MPC_ClearCache();
		reporter.Done();
	}

	Spectrum tr = Spectrum::FromSampledSpectrum(integrateSpectralProfile(pData->spectralProfile,
		&MultipoleProfileDataEntry::reflectance));
	Info("Total Reflectance: %s", tr.ToString().c_str());

#ifdef SAMPLE_TRANSMITTANCE
	Spectrum tt = Spectrum::FromSampledSpectrum(integrateSpectralProfile(pData->spectralProfile,
		&MultipoleProfileDataEntry::transmittance));
	Info("Total Transmittance: %s", tt.ToString().c_str());
#endif
}

void ReleaseMultipoleProfile(MultipoleProfileData* pData) {
	delete pData;
}


Spectrum MultipoleBSSRDFData::rho(float costheta) const {
	float fid = costheta * (rhoData.hd.size() - 1);
	int id = Clamp((int)fid, 0, (int)(rhoData.hd.size() - 2));
	float lerp = fid - (float)id;
	return Lerp(lerp, rhoData.hd[id], rhoData.hd[id + 1]);
}


Spectrum ComputeRhoHHFromBxDF(const BxDF* bxdf) {
	RNG rng(6428263 * 3 * 7);

	int sqrtSamples = 256;
    int nSamples = sqrtSamples * sqrtSamples;
    float* s1 = (float*)malloc(sizeof(float) * 2 * nSamples);
    StratifiedSample2D(s1, sqrtSamples, sqrtSamples, rng);
    float* s2 = (float*)malloc(sizeof(float) * 2 * nSamples);
    StratifiedSample2D(s2, sqrtSamples, sqrtSamples, rng);
    Spectrum rhohh = bxdf->rho(nSamples, s1, s2);
	free(s2);
	free(s1);

	return rhohh;
}


Spectrum MultipoleBSSRDFData::rho() const {
	return rhoData.hh;
}


class RhoTask : public Task {
public:
	RhoTask(Spectrum& dest, const BxDF* bxdf, float costheta, ProgressReporter& pr, int id)
		: destination(dest), bxdf(bxdf), costheta(costheta), reporter(pr), id(id)
	{

	}

	void Run() override;
private:
	Spectrum& destination;
	const BxDF* bxdf;
	float costheta;
	ProgressReporter& reporter;
	int id;
};


void RhoTask::Run() {
	RNG rng(6428263 * id);

	int sqrtSamples = 256;
	int nSamples = sqrtSamples * sqrtSamples;
	float* s1 = (float*)malloc(sizeof(float) * 2 * nSamples);
	StratifiedSample2D(s1, sqrtSamples, sqrtSamples, rng);
	Vector wo = SphericalDirection<float>(sqrtf(1 - costheta * costheta), costheta, 0.f);
	destination = bxdf->rho(wo, nSamples, s1);
	free(s1);

	reporter.Update();
}


RhoData ComputeRhoDataFromBxDF(const BxDF* bxdf) {
	RhoData result;

	const int rhoDataSize = 1024;
	result.hd.resize(rhoDataSize + 1, Spectrum(0.f));

	int nTasks = rhoDataSize + 1;
	vector<Task*> tasks;
	tasks.reserve(nTasks);

	ProgressReporter reporter(nTasks + 1, "Rho Data");
	for (int taskId = 0; taskId < nTasks; taskId++) {
		float costheta = (float)taskId / rhoDataSize;
		if (costheta == 0.f)
			costheta = 0.01f / rhoDataSize;
		tasks.push_back(new RhoTask(result.hd[taskId], bxdf,
			costheta, reporter, taskId));
	}
	EnqueueTasks(tasks);
	WaitForAllTasks();
	for (Task* task : tasks)
		delete task;

	result.hh = ComputeRhoHHFromBxDF(bxdf);
	reporter.Update();
	reporter.Done();

	return result;
}

void ComputeIrradiancePointsProfile(MultipoleProfileData** oppData, float radius) {
	if (!oppData) return;

	MultipoleProfileData* pData = *oppData = new MultipoleProfileData;
	SampledSpectrum fullSpectrum(1.f);
	for (int i = 0; i < SampledSpectrum::nComponents; i++) {
		Profile& profile = pData->spectralProfile[i];
		float area = M_PI * radius * radius;
		MultipoleProfileDataEntry entry = { fullSpectrum[i] / area };
		profile.data.push_back(entry);
		profile.data.push_back(entry);
		profile.dsqSpacing = radius * radius;
		profile.rcpDsqSpacing = 1.f / profile.dsqSpacing;
		profile.totalReflectance = fullSpectrum[i];
		profile.totalTransmittance = 0.f;
	}
}

RhoData ComputeRoughRhoData() {
	RhoData result = { vector<Spectrum>(2, Spectrum(0.f)), Spectrum(0.f) };
	return result;
}
