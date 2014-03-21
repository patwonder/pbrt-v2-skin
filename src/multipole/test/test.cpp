// test.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>
#include "MultipoleProfileCalculator.h"

using namespace std;

void computeConfiguration(uint32 numLayers, const MPC_LayerSpec* pLayerSpecs,
	const MPC_Options* pOptions) {
	MPC_Output* pOutput;
	MPC_ComputeDiffusionProfile(numLayers, pLayerSpecs, pOptions, &pOutput);
	
	// The integrals does not include delta distributions
	cout << "Reflectance(" << pOutput->length << "):";
	float integralR = 0.f;
	float prevDistance = 0.f;
	for (uint32 i = 0; i < pOutput->length; i++) {
		if (i < 32)
			cout << " " << pOutput->pReflectance[i];
		integralR += 2.f * 3.141592654f * pOutput->pDistance[i] * pOutput->pReflectance[i] * (pOutput->pDistance[i] - prevDistance);
		prevDistance = pOutput->pDistance[i];
	}
	if (pOutput->length > 32)
		cout << " ...";
	cout << endl << "Transmittance(" << pOutput->length << "):";
	float integralT = 0;
	prevDistance = 0.f;
	for (uint32 i = 0; i < pOutput->length; i++) {
		if (i < 32)
			cout << " " << pOutput->pTransmittance[i];
		integralT += 2.f * 3.141592654f * pOutput->pDistance[i] * pOutput->pTransmittance[i] * (pOutput->pDistance[i] - prevDistance);
		prevDistance = pOutput->pDistance[i];
	}
	if (pOutput->length > 32)
		cout << " ...";
	cout << endl;
	cout << "Reflectance Integral: " << integralR << endl;
	cout << "Transmittance Integral: " << integralT << endl;
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
	MPC_Options options;
	options.desiredLength = 512;
	float mfp0 = 1.f / (specs[0].mua + specs[0].musp);
	float mfp1 = 1.f / (specs[1].mua + specs[1].musp);
	options.desiredStepSize = 8.f * (mfp0 + mfp1) / (float)options.desiredLength;
	
	cout << "First layer: " << endl;
	computeConfiguration(1, &specs[0], &options);
	cout << endl << "Second Layer: " << endl;
	computeConfiguration(1, &specs[1], &options);
	cout << endl << "Combined Layer: " << endl;
	computeConfiguration(2, &specs[0], &options);

	system("PAUSE");

	return 0;
}

