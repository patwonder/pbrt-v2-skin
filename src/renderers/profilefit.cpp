
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

#include "profilefit.h"
#include "paramset.h"
#include "materials/skincoeffs.h"
#include "parallel.h"
#include "progressreporter.h"
#include "multipole/MultipoleProfileCalculator/MultipoleProfileCalculator.h"
#include "multipole/MultipoleProfileCalculator/gaussianfit.h"
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace std;

struct VariableParams {
	VariableParams(float f_mel, float f_eu, float f_blood, float f_ohg)
		: f_mel(f_mel), f_eu(f_eu), f_blood(f_blood), f_ohg(f_ohg) { }
	float f_mel, f_eu, f_blood, f_ohg;
};

static ostream& operator<<(ostream& os, const VariableParams& vps) {
	return os << vps.f_mel << "\t" << vps.f_eu << "\t" << vps.f_blood << "\t" << vps.f_ohg;
}

static ostream& operator<<(ostream& os, const ParamRange& pr) {
	os << pr.size();
	for (float value : pr) {
		os << "\t" << value;
	}
	return os;
}

template <class Processor>
void ForParamRange(const ParamRange& pr, const Processor& processor) {
	for (float value : pr) {
		processor(value);
	}
}

template <class Processor>
void MultipoleProfileFitRenderer::ForRanges(const Processor& p) const {
	uint32_t id = 0;
	ForParamRange(pr_f_mel, [&](float f_mel) {
		ForParamRange(pr_f_eu, [&](float f_eu) {
			ForParamRange(pr_f_blood, [&](float f_blood) {
				ForParamRange(pr_f_ohg, [&](float f_ohg) {
					p(VariableParams(f_mel, f_eu, f_blood, f_ohg), id);
					id++;
				});
			});
		});
	});
}


int MultipoleProfileFitRenderer::GetTotalNumTasks() const {
	return (int)(pr_f_mel.size() * pr_f_eu.size() * pr_f_blood.size() * pr_f_ohg.size());
}


float MultipoleProfileFitRenderer::NextParamFromId(uint32_t& id, ParamRange pr) {
	uint32_t remainder = id % (uint32_t)pr.size();
	id /= (uint32_t)pr.size();
	return pr[remainder];
}

VariableParams MultipoleProfileFitRenderer::ParamsFromId(uint32_t id) const {
	float f_ohg = NextParamFromId(id, pr_f_ohg);
	float f_blood = NextParamFromId(id, pr_f_blood);
	float f_eu = NextParamFromId(id, pr_f_eu);
	float f_mel = NextParamFromId(id, pr_f_mel);
	return VariableParams(f_mel, f_eu, f_blood, f_ohg);
}

MultipoleProfileFitRenderer::MultipoleProfileFitRenderer(const vector<SkinLayer>& layers,
	ParamRange f_mel, ParamRange f_eu, ParamRange f_blood, ParamRange f_ohg,
	uint32_t desiredLength, uint32_t segments, const string& filename,
	int idRange[2], int splitTasks, const string& pbrtFilePath)
	: layers(layers), pr_f_mel(f_mel), pr_f_eu(f_eu), pr_f_blood(f_blood), pr_f_ohg(f_ohg),
	  desiredLength(desiredLength), nSegments(segments), filename(filename),
	  splitTasks(splitTasks), pbrtFilePath(pbrtFilePath)
{
	this->idRange[0] = idRange[0];
	this->idRange[1] = idRange[1];
}


class GaussianFitTask : public Task {
public:
	static const int nLayers = 2;

	GaussianFitTask(const SampledSpectrum* mua, const SampledSpectrum* musp,
		const float* et, const float* thickness, const vector<float>& sigmas,
		SpectralGaussianCoeffs& coeffs, ProgressReporter& pr, uint32_t id,
		int sc, uint32_t desiredLength)
		: sigmas(sigmas), coeffs(coeffs), reporter(pr), id(id), sc(sc),
		  desiredLength(desiredLength)
	{
		for (int i = 0; i < nLayers; i++) {
			this->mua[i] = mua[i][sc];
			this->musp[i] = musp[i][sc];
			this->et[i] = et[i];
			this->thickness[i] = thickness[i];
		}
	}

private:
	float mua[nLayers];
	float musp[nLayers];
	float et[nLayers];
	float thickness[nLayers];
	const vector<float>& sigmas;
	SpectralGaussianCoeffs& coeffs;
	uint32_t id;
	uint32_t desiredLength;
	int sc;

	ProgressReporter& reporter;

	void Run() override;
	void DoGaussianFit(const MPC_Output* pOutput, float mfpMin, float mfpMax);
};

void GaussianFitTask::Run() {
	MPC_LayerSpec* pLayerSpecs = new MPC_LayerSpec[nLayers];
	MPC_Options options;
	options.desiredLength = desiredLength;
	options.lerpOnThinSlab = true;

	// Compute mfp
	float mfpMin = INFINITY;
	float mfpMax = 0.f;
	float mfpTotal = 0.f;
	for (int layer = 0; layer < nLayers; layer++) {
		float mfp = 1.f / (mua[layer] + musp[layer]);
		mfpTotal += mfp;
		mfpMin = min(mfpMin, mfp);
		mfpMax = max(mfpMax, mfp);
	}
	float mfp = mfpTotal / (float)nLayers;
	
	// Fill in layer information
	for (int layer = 0; layer < nLayers; layer++) {
		MPC_LayerSpec& ls = pLayerSpecs[layer];
		ls.g_HG = 0.f;
		ls.ior = et[layer];
		ls.mua = mua[layer];
		ls.musp = musp[layer];
		ls.thickness = thickness[layer];
	}
	// Compute desired step size
	options.desiredStepSize = 12.f * mfp / (float)options.desiredLength;

	// Do the computation
	MPC_Output* pOutput;
	MPC_ComputeDiffusionProfile(nLayers, pLayerSpecs, &options, &pOutput);
	MPC_ResampleForUniformDistanceSquaredDistribution(pOutput, pOutput->length);

	// Do gaussian fit
	DoGaussianFit(pOutput, mfpMin, mfpMax);
		
	MPC_FreeOutput(pOutput);

	delete [] pLayerSpecs;

	reporter.Update();
}

void GaussianFitTask::DoGaussianFit(const MPC_Output* pOutput, float mfpMin, float mfpMax) {
	int nSigmas = sigmas.size();
	const int nTargetSigmas = 6;
	const int sigmaNegExtent = nTargetSigmas / 2;
	const int sigmaPosExtent = nTargetSigmas - sigmaNegExtent - 1;

	// Compute the best fit of nTargetSigmas gaussians
	uint32 length = pOutput->length;
	vector<float> distArray(length);
	for (uint32 i = 0; i < length; i++) {
		float dsq = pOutput->pDistanceSquared[i];
		distArray[i] = sqrt(dsq);
	}
	GF_Output* pGFOutput;

	float minError = INFINITY;
	int minErrorSigmaCenter = 0;
	std::map<int, GF_Output*> oCache;
	// Try sigmas between mfpMin / 2 and mfpMax * 2
	for (int iSigmaCenter = sigmaNegExtent; iSigmaCenter < nSigmas - sigmaPosExtent; iSigmaCenter++) {
		float sigma = sigmas[iSigmaCenter];
		if (sigma < mfpMin / 2.f) continue;
		if (sigma > mfpMax * 2.f) break;
		GF_FitSumGaussians(length, &distArray[0], pOutput->pReflectance,
			nTargetSigmas, &sigmas[iSigmaCenter - sigmaNegExtent], &pGFOutput);
		if (pGFOutput->overallError < minError) {
			minError = pGFOutput->overallError;
			minErrorSigmaCenter = iSigmaCenter;
		}
		oCache[iSigmaCenter] = pGFOutput;
		pGFOutput = NULL;
	}

	// Use minErrorSigmaCenter to do fitting
	//GF_FitSumGaussians(length, &distArray[0], pOutput->pReflectance,
	//	nTargetSigmas, &sigmas[minErrorSigmaCenter - sigmaNegExtent], &pGFOutput);
	pGFOutput = oCache[minErrorSigmaCenter];

	for (int iSigma = minErrorSigmaCenter - sigmaNegExtent;
		iSigma <= minErrorSigmaCenter + sigmaPosExtent; iSigma++)
	{
		coeffs.coeffs[iSigma][sc] =
			pGFOutput->pNormalizedCoeffs[iSigma - (minErrorSigmaCenter - sigmaNegExtent)];
	}
	coeffs.error[sc] = pGFOutput->overallError;
	if (pGFOutput->overallError > 0.05f) {
		Warning("ID %d wavelength %f fit has error %.2f%%", id,
			SampledSpectrum::WaveLength(sc), pGFOutput->overallError * 100.f);
	}

	for (auto& pair : oCache)
		GF_FreeOutput(pair.second);
}

vector<Task*> MultipoleProfileFitRenderer::CreateGaussianFitTasks(const SkinCoefficients& coeffs,
	const vector<float>& sigmas, SpectralGaussianCoeffs& sgc, ProgressReporter& reporter, uint32_t id) const
{
	sgc.coeffs.resize(sigmas.size(), SampledSpectrum(0.f));

	SampledSpectrum mua[2] = {
		coeffs.mua_epi().toSampledSpectrum() / 10.f,
		coeffs.mua_derm().toSampledSpectrum() / 10.f
	};
	SampledSpectrum musp[2] = {
		coeffs.musp_epi().toSampledSpectrum() / 10.f,
		coeffs.musp_derm().toSampledSpectrum() / 10.f
	};
	float et[2] = { layers[0].ior, layers[1].ior };
	float thickness[2] = {
		layers[0].thickness / 1e6f,
		layers[1].thickness / 1e6f
	};

	vector<Task*> tasks;
	for (int sc = 0; sc < SampledSpectrum::nComponents; sc++) {
		tasks.push_back(new GaussianFitTask(mua, musp, et, thickness,
			sigmas, sgc, reporter, id, sc, desiredLength));
	}
	return tasks;
}

void MultipoleProfileFitRenderer::DoProfileFit() const {
	// All the units here are based on mm (millimeter)
	// Calculate min/max mfps
	WLDValue mfp_min(INFINITY), mfp_max(0.f);
	int nTotalTasks = 0;
	ForRanges([&] (const VariableParams& vps, uint32_t id) {
		SkinCoefficients coeffs(vps.f_mel, vps.f_eu, vps.f_blood, vps.f_ohg, 0.f, 0.f, 0.f);
		WLDValue mutp_epi = coeffs.mua_epi() + coeffs.musp_epi();
		WLDValue mutp_derm = coeffs.mua_derm() + coeffs.musp_derm();
		WLDValue mfp_epi = WLDValue(10) / mutp_epi;
		WLDValue mfp_derm = WLDValue(10) / mutp_derm;
		mfp_min = Min(Min(mfp_min, mfp_epi), mfp_derm);
		mfp_max = Max(Max(mfp_max, mfp_epi), mfp_derm);
		nTotalTasks++;
	});
	float minmfp = mfp_min.Min();
	float maxmfp = mfp_max.Max();
	Info("Min and max mfp, in millimeters: %f %f", minmfp, maxmfp);

	// Generate sigmas
	vector<float> sigmas;
	for (float sigma = minmfp / 8.f; sigma <= maxmfp * 8.f; sigma *= 2.f) {
		sigmas.push_back(sigma);
	}
	// Compute profiles and do fitting
	uint32_t idMin = (uint32_t)max(0, idRange[0]);
	uint32_t idMax = (uint32_t)(idRange[1] < 0 ? nTotalTasks : min(nTotalTasks, idRange[1]));
	if (idMax <= idMin) return;
	int nTasks = idMax - idMin;

	vector<SpectralGaussianCoeffs> gaussianCoeffs(nTasks);
	ProgressReporter reporter(nTasks * SampledSpectrum::nComponents, "Gaussian Fit");
	vector<Task*> fitTasks;
	fitTasks.reserve(nTasks * SampledSpectrum::nComponents);
	ForRanges([&] (const VariableParams& vps, uint32_t id) {
		if (id < idMin || id >= idMax) return;

		SkinCoefficients coeffs(vps.f_mel, vps.f_eu, vps.f_blood, vps.f_ohg, 0.f, 0.f, 0.f);
		vector<Task*> newTasks = CreateGaussianFitTasks(coeffs, sigmas, gaussianCoeffs[id - idMin], reporter, id);
		fitTasks.insert(fitTasks.end(), newTasks.begin(), newTasks.end());
	});
	EnqueueTasks(fitTasks);
	WaitForAllTasks();
	for (Task* task : fitTasks) {
		delete task;
	}
	MPC_ClearCache();
	reporter.Done();

	if (filename != "") {
		string filename = AbsolutePath(ResolveFilename(this->filename));
		Info("Writing to file \"%s\"", filename.c_str());
		ofstream out(filename, ios::out | ios::trunc);
		if (!out) {
			Error("Cannot open %s", filename.c_str());
			abort();
		}
		// Only print the header when range contains id 0.
		// This is to allow easier concatenation of results.
		if (idMin == 0) {
			// Output _ParamRange_s
			out << "Param\tNum Samples\tSample Points" << endl;
			out << "f_mel\t" << pr_f_mel << endl;
			out << "f_eu\t" << pr_f_eu << endl;
			out << "f_blood\t" << pr_f_blood << endl;
			out << "f_ohg\t" << pr_f_ohg << endl;
			// Output sigmas
			out << "ID\tf_mel\tf_eu\tf_blood\tf_ohg\tRGB";
			for (float sigma : sigmas) {
				out << "\t" << sigma;
			}
			out << "\tError" << endl;
		}
		for (uint32_t id = idMin; id < idMax; id++) {
			const SpectralGaussianCoeffs& coeffs = gaussianCoeffs[id - idMin];
			VariableParams vps = ParamsFromId(id);
			for (int sc = 0; sc < SampledSpectrum::nComponents; sc++) {
				out << id << "\t" << vps << "\t" << SampledSpectrum::WaveLength(sc);
				for (int i = 0; i < sigmas.size(); i++)
					out << "\t" << coeffs.coeffs[i][sc];
				out << "\t" << coeffs.error[sc] << endl;
			}
			vector<RGBSpectrum> rgbCoeffs;
			for (const SampledSpectrum& spectrum : coeffs.coeffs) {
				rgbCoeffs.push_back(spectrum.ToRGBSpectrum());
			}
			for (int rgb = 0; rgb < RGBSpectrum::nComponents; rgb++) {
				out << id << "\t" << vps << "\t" << "RGB"[rgb];
				for (int i = 0; i < sigmas.size(); i++)
					out << "\t" << rgbCoeffs[i][rgb];
				out << endl;
			}
		}
		out.close();
	}
}


static string CombineFileName(const string& prefix, int id, const string& suffix) {
	stringstream ss;
	ss << prefix << id << suffix;
	return ss.str();
}


static string GetParamRangeStringInPbrtFile(const ParamRange& pr) {
	if (pr.size() == 0) return "";

	stringstream ss;
	ss << pr[0];
	for (size_t i = 1; i < pr.size(); i++) {
		ss << " " << pr[i];
	}
	return ss.str();
}


static ostream& operator<<(ostream& os, const SkinLayer& layer) {
	return os << layer.thickness << " " << layer.ior;
}

static string GetSkinLayerStringInPbrtFile(const vector<SkinLayer>& layers) {
	if (layers.size() == 0) return "";
	
	stringstream ss;
	ss << layers[0];
	for (size_t i = 1; i < layers.size(); i++) {
		ss << " " << layers[i];
	}
	return ss.str();
}


void MultipoleProfileFitRenderer::DoPbrtFileSplit() const {
	int nTotalTasks = GetTotalNumTasks();
	int nSubPbrtFiles = (nTotalTasks + splitTasks - 1) / splitTasks;

	if (nSubPbrtFiles <= 1) {
		Warning("No need to split!");
		return;
	}

	// Avoid doing stupid things...
	if (nSubPbrtFiles > 20) {
		char ch;
		Warning("%d files will be generated. Are you sure? [Y/n]", nSubPbrtFiles);
		cin >> ch;
		if (ch != 'Y' && ch != 'y') {
			Warning("Aborted.");
			return;
		}
	}

	string pbrtFilePrefix = AbsolutePath(ResolveFilename(FileNameWithoutExt(pbrtFilePath))) + "-";
	string outFilePrefix = FileNameWithoutExt(filename) + "-";
	string pbrtFileSuffix = ".pbrt";
	string outFileSuffix = ".txt";

	for (int i = 0; i < nSubPbrtFiles; i++) {
		int idMin = i * splitTasks;
		int idMax = min((i + 1) * splitTasks, nTotalTasks);

		string pbrtFilename = CombineFileName(pbrtFilePrefix, i, pbrtFileSuffix);
		string outFilename = CombineFileName(outFilePrefix, i, outFileSuffix);
		
		Info("%d - %d -> %s", idMin, idMax, pbrtFilename.c_str());

		ofstream out(pbrtFilename, ios::out | ios::trunc);
		if (!out) {
			Error("Failed to open pbrt output file: %s", pbrtFilename.c_str());
			continue;
		}

		// Print out parameters for the profile fit
		out << "Renderer \"multipoleprofilefit\"" << endl
			<< "         \"float f_mel\" [ " << GetParamRangeStringInPbrtFile(pr_f_mel) << " ]" << endl
			<< "         \"float f_eu\" [ " << GetParamRangeStringInPbrtFile(pr_f_eu) << " ]" << endl
			<< "         \"float f_blood\" [ " << GetParamRangeStringInPbrtFile(pr_f_blood) << " ]" << endl
			<< "         \"float f_ohg\" [ " << GetParamRangeStringInPbrtFile(pr_f_ohg) << " ]" << endl
			<< "         \"skinlayer layers\" [ " << GetSkinLayerStringInPbrtFile(layers) << " ]" << endl
			<< "         \"string filename\" \"" << outFilename << "\"" << endl
			<< "         \"integer desiredProfileLength\" " << desiredLength << endl
			<< "         \"integer idrange\" [ " << idMin << " " << idMax << " ]" << endl
			<< "WorldBegin" << endl
			<< "WorldEnd" << endl;
		out.close();
	}

	Info("Generated %d files.", nSubPbrtFiles);
}


void MultipoleProfileFitRenderer::Render(const Scene* scene) {
	Info("Computing using MultipoleProfileFitRenderer.");

	if (splitTasks) {
		DoPbrtFileSplit();
	} else {
		DoProfileFit();
	}
}

Spectrum MultipoleProfileFitRenderer::Li(const Scene* scene, const RayDifferential& ray,
	const Sample* sample, RNG& rng, MemoryArena& arena,
	Intersection* isect, Spectrum* T) const
{
	return Spectrum(0.f);
}

Spectrum MultipoleProfileFitRenderer::Transmittance(const Scene* scene, const RayDifferential& ray,
	const Sample* sample, RNG& rng, MemoryArena& arena) const
{
	return Spectrum(0.f);
}

static ParamRange FindParamRange(const ParamSet& params, const string& name) {
	int nItems;
	const float* pData = params.FindFloat(name, &nItems);
	if (!pData) {
		Error("Param %s not specified for MultipoleProfileFitRenderer.", name);
		abort();
	}
	if (nItems < 1) {
		Error("Param %s doesn't specify any float.", name);
		abort();
	}
	ParamRange pr(pData, pData + nItems);
	sort(pr.begin(), pr.end());
	return pr;
}

MultipoleProfileFitRenderer* CreateMultipoleProfileFitRenderer(const ParamSet& params, const string& pbrtFilePath) {
	int nLayers;
	const SkinLayer* layers = params.FindSkinLayer("layers", &nLayers);
	if (!layers) {
		Error("No layers param set for MultipoleProfileFitRenderer.");
		abort();
	}
	if (nLayers < 2) {
		Error("Not enough layers specified for MultipoleProfileFitRenderer.");
		abort();
	}
	ParamRange pr_f_mel = FindParamRange(params, "f_mel");
	ParamRange pr_f_eu = FindParamRange(params, "f_eu");
	ParamRange pr_f_blood = FindParamRange(params, "f_blood");
	ParamRange pr_f_ohg = FindParamRange(params, "f_ohg");
	uint32_t segments = params.FindOneInt("segments", 10);
    string filename = params.FindOneString("filename", "");
	uint32_t desiredLength = params.FindOneInt("desiredProfileLength", 256);
	int nRanges;
	const int* idRanges = params.FindInt("idrange", &nRanges);
	int idRangesArray[2];
	if (!idRanges || nRanges != 2) {
		if (idRanges)
			Warning("Param \"idrange\" ignored: incorrect number of items, 2 expected.");
		idRangesArray[0] = 0;
		idRangesArray[1] = -1;
	} else {
		idRangesArray[0] = idRanges[0];
		idRangesArray[1] = idRanges[1];
	}
	int splitTasks = params.FindOneInt("splittasks", 0);
	return new MultipoleProfileFitRenderer(vector<SkinLayer>(layers, layers + nLayers),
		pr_f_mel, pr_f_eu, pr_f_blood, pr_f_ohg, desiredLength, segments, filename,
		idRangesArray, splitTasks, pbrtFilePath);
}
