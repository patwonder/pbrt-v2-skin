// test.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>
#include "MultipoleProfileCalculator.h"

using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{
	MPC_LayerSpec specs[1];
	//specs[0].mua = 13.8629f;
	//specs[0].musp = 19.4879f;
	specs[0].mua = 0;
	specs[0].musp = 1;
	specs[0].thickness = 1 / (specs[0].mua + specs[0].musp);
	specs[0].ior = 1.0f;
	specs[0].g_HG = 0.90f;
	MPC_Options options;
	options.desiredLength = 1024;
	float mfp0 = 1.f / (specs[0].mua + specs[0].musp);
	options.desiredStepSize = 16.f * mfp0 / (float)options.desiredLength;
	
	MPC_Output* pOutput;
	MPC_ComputeDiffusionProfile(1, specs, &options, &pOutput);
	
	cout << "Reflectance:";
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
	cout << endl << "Transmittance:";
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

	return 0;
}

