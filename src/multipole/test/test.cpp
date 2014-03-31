// test.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>
#include "MultipoleProfileCalculator.h"
#include "gaussianfit.h"
#include <vector>
#include <iomanip>

using namespace std;

void computeGaussianFit(const MPC_Output* pOutput) {
	uint32 length = pOutput->length;
	float mfp = sqrt(pOutput->pDistanceSquared[length - 1]) / 12.f;
	vector<float> distArray(length);
	for (uint32 i = 0; i < length; i++) {
		float dsq = pOutput->pDistanceSquared[i];
		distArray[i] = sqrt(dsq);
	}
	const uint32 nSigmas = 4;
	float sigmaArray[nSigmas];
	for (uint32 i = 0; i < nSigmas; i++) {
		sigmaArray[i] = (float)(mfp * pow(2, -((int)nSigmas / 2) + (int)i));
	}

	GF_Output* pGFOutput;
	GF_FitSumGaussians(length, &distArray[0], pOutput->pReflectance, nSigmas, sigmaArray, &pGFOutput);

	cout << "Gaussian Fit: Error = " 
		 << setiosflags(ios::fixed) << setprecision(4) << pGFOutput->overallError * 100.f << "%" << endl;
	cout.copyfmt(ios(NULL));
	cout << "{ ";
	for (uint32 i = 0; i < pGFOutput->nCoeffs; i++) {
		if (i) cout << ", ";
		cout << pGFOutput->pNormalizedCoeffs[i];
	}
	cout << " }" << endl;

	GF_FreeOutput(pGFOutput);
}

void computeConfiguration(uint32 numLayers, const MPC_LayerSpec* pLayerSpecs,
	const MPC_Options* pOptions) {
	MPC_Output* pOutput;
	MPC_ComputeDiffusionProfile(numLayers, pLayerSpecs, pOptions, &pOutput);
	MPC_ResampleForUniformDistanceSquaredDistribution(pOutput);

	// Computes the integral {r,0,+inf}2*pi*r*Rd(r)dr = {r^2,0,+inf}pi*Rd(r)d(r^2) as total reflectance
	// The integrals does not include delta distributions
	cout << "Reflectance(" << pOutput->length << "):";
	float integralR = 0.f;
	float prevDistanceSquared = 0.f;
	for (uint32 i = 0; i < pOutput->length; i++) {
		if (i < 32)
			cout << " " << pOutput->pReflectance[i];
		integralR += 3.141592654f * pOutput->pReflectance[i] * (pOutput->pDistanceSquared[i] - prevDistanceSquared);
		prevDistanceSquared = pOutput->pDistanceSquared[i];
	}
	if (pOutput->length > 32)
		cout << " ...";
	cout << endl << "Transmittance(" << pOutput->length << "):";
	float integralT = 0;
	prevDistanceSquared = 0.f;
	for (uint32 i = 0; i < pOutput->length; i++) {
		if (i < 32)
			cout << " " << pOutput->pTransmittance[i];
		integralT += 3.141592654f * pOutput->pTransmittance[i] * (pOutput->pDistanceSquared[i] - prevDistanceSquared);
		prevDistanceSquared = pOutput->pDistanceSquared[i];
	}
	if (pOutput->length > 32)
		cout << " ...";
	cout << endl;
	cout << "Reflectance Integral: " << integralR << endl;
	cout << "Transmittance Integral: " << integralT << endl;

	computeGaussianFit(pOutput);

	MPC_FreeOutput(pOutput);
}

int _tmain(int argc, _TCHAR* argv[])
{
	MPC_LayerSpec specs[2];
	specs[0].mua = 13.8629f;
	specs[0].musp = 19.4879f;
	specs[0].thickness = 0.025f;
	specs[0].ior = 1.4f;
	specs[0].g_HG = 0.90f;
	specs[1].mua = 1.6f;
	specs[1].musp = 19.4879f;
	specs[1].thickness = 0.2f;
	specs[1].ior = 1.45f;
	specs[1].g_HG = 0.80f;
	//specs[1].mua = 0.1f;
	//specs[1].musp = 0.9f;
	//specs[1].thickness = 2.f;
	//specs[1].ior = 1.f;
	//specs[1].g_HG = 0.80f;
	MPC_Options options;
	options.desiredLength = 512;
	float mfp0 = 1.f / (specs[0].mua + specs[0].musp);
	float mfp1 = 1.f / (specs[1].mua + specs[1].musp);
	
	cout << "First layer: " << endl;
	options.desiredStepSize = 12.f * mfp0 / (float)options.desiredLength;
	computeConfiguration(1, &specs[0], &options);

	cout << endl << "Second Layer: " << endl;
	options.desiredStepSize = 12.f * mfp1 / (float)options.desiredLength;
	computeConfiguration(1, &specs[1], &options);

	cout << endl << "Combined Layer: " << endl;
	options.desiredStepSize = 12.f * min(mfp0, mfp1) / (float)options.desiredLength;
	computeConfiguration(2, &specs[0], &options);

	MPC_ClearCache();

	system("PAUSE");

	return 0;
}

