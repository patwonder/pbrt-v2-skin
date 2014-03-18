// test.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>
#include "MultipoleProfileCalculator.h"

using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{
	MPC_LayerSpec specs[1];
	specs[0].mua = 13.8629f;
	specs[0].musp = 19.4879f;
	specs[0].thickness = 0.035f;
	specs[0].ior = 1.4f;
	specs[0].g_HG = 0.90f;
	MPC_Options options;
	options.desiredLength = 64;
	float mfp0 = 1.f / (specs[0].mua + specs[0].musp);
	options.desiredStepSize = 10.f * mfp0 / (float)options.desiredLength;
	
	MPC_Output* pOutput;
	MPC_ComputeDiffusionProfile(1, specs, &options, &pOutput);
	
	cout << "Step size: " << pOutput->stepSize << endl;
	cout << "Reflectance:";
	for (uint32 i = 0; i < pOutput->length; i++) {
		cout << " " << pOutput->pReflectance[i];
	}
	cout << endl << "Transmittance:";
	for (uint32 i = 0; i < pOutput->length; i++) {
		cout << " " << pOutput->pTransmittance[i];
	}
	cout << endl;
	MPC_FreeOutput(pOutput);

	return 0;
}

