
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
#include "MultipoleProfileCalculator.h"
#include <vector>
#include "DipoleCalculator.h"
#include "numutil.h"
#include "tools/kiss_fftndr.h"

using namespace std;

inline kiss_fft_cpx& operator+=(kiss_fft_cpx& cpx, const kiss_fft_cpx& other) {
	cpx.r += other.r;
	cpx.i += other.i;
	return cpx;
}

inline kiss_fft_cpx& operator*=(kiss_fft_cpx& cpx, const kiss_fft_cpx& other) {
	float r = cpx.r * other.r - cpx.i * other.i;
	float i = cpx.r * other.i + cpx.i * other.r;
	cpx.r = r; cpx.i = i;
	return cpx;
}

inline kiss_fft_cpx& operator/=(kiss_fft_cpx& cpx, const kiss_fft_cpx& other) {
	float divisor = other.r * other.r + other.i * other.i;
	float acbd = cpx.r * other.r + cpx.i * other.i;
	float bcad = cpx.i * other.r - cpx.r * other.i;
	cpx.r = acbd / divisor;
	cpx.i = bcad / divisor;
	return cpx;
}

struct MatrixProfile {
	MatrixProfile(uint32 nRows, uint32 nCols)
		: reflectance(nRows, nCols), transmittance(nRows, nCols) { }
	Matrix<float> reflectance;
	Matrix<float> transmittance;
	uint32 GetNumRows() { return reflectance.GetNumRows(); }
	uint32 GetNumCols() { return reflectance.GetNumCols(); }
	void Clear() { reflectance.Clear(); transmittance.Clear(); }
};


void FFT(const Matrix<float>& profile, Matrix<kiss_fft_cpx>& out) {
	int dims[2] = { profile.GetNumRows(), profile.GetNumCols() };
	kiss_fftndr_cfg cfg = kiss_fftndr_alloc(dims, 2, 0, NULL, NULL);
	kiss_fftndr(cfg, profile.GetData(), out.GetData());
	kiss_fft_free(cfg);
}


void IFFT(const Matrix<kiss_fft_cpx>& profile, Matrix<float>& out) {
	int dims[2] = { profile.GetNumRows(), profile.GetNumCols() };
	kiss_fftndr_cfg cfg = kiss_fftndr_alloc(dims, 2, 1, NULL, NULL);
	kiss_fftndri(cfg, profile.GetData(), out.GetData());
	kiss_fft_free(cfg);
}


static inline uint32 RoundUpPow2(uint32 v) {
    v--;
    v |= v >> 1;    v |= v >> 2;
    v |= v >> 4;    v |= v >> 8;
    v |= v >> 16;
    return v+1;
}


void ComputeLayerProfile(const MPC_LayerSpec& spec, float iorUpper, float iorLower,
	float stepSize, MatrixProfile& profile)
{
	profile.Clear();

	uint32 numRows = profile.GetNumRows(), numCols = profile.GetNumCols();
	uint32 centerRow = (numRows - 1) / 2, centerCol = (numCols - 1) / 2;
	uint32 extentRow = centerRow, extentCol = centerCol;
	const uint32 numDipolePairs = 15;
	for (int32 pair = -((int32)numDipolePairs - 1) / 2; pair <= (int)(numDipolePairs - 1) / 2; pair++) {
		DipoleCalculator dc(iorUpper, iorLower, spec.thickness, spec.mua, spec.musp, pair);
		for (uint32 sampleRow = 0; sampleRow <= extentRow; sampleRow++) {
			for (uint32 sampleCol = 0; sampleCol <= extentCol; sampleCol++) {
				float r2 = (sampleRow * sampleRow + sampleCol * sampleCol) * (stepSize * stepSize);
				float rd = dc.Rd(r2), td = dc.Td(r2);
				profile.reflectance[centerRow + sampleRow][centerCol + sampleCol] += rd;
				profile.transmittance[centerRow + sampleRow][centerCol + sampleCol] += td;
			}
		}
	}
	for (uint32 sampleRow = 0; sampleRow <= extentRow; sampleRow++) {
		for (uint32 sampleCol = 0; sampleCol <= extentCol; sampleCol++) {
			float rd = profile.reflectance[centerRow + sampleRow][centerCol + sampleCol];
			profile.reflectance[centerRow - sampleRow][centerCol + sampleCol] = rd;
			profile.reflectance[centerRow + sampleRow][centerCol - sampleCol] = rd;
			profile.reflectance[centerRow - sampleRow][centerCol - sampleCol] = rd;
			float td = profile.transmittance[centerRow + sampleRow][centerCol + sampleCol];
			profile.transmittance[centerRow - sampleRow][centerCol + sampleCol] = td;
			profile.transmittance[centerRow + sampleRow][centerCol - sampleCol] = td;
			profile.transmittance[centerRow - sampleRow][centerCol - sampleCol] = td;
		}
	}
}


void CombineLayerProfiles(const MatrixProfile& layer1, const MatrixProfile& layer2,
	MatrixProfile& combined)
{
}


MULTIPOLEPROFILECALCULATOR_API void MPC_ComputeDiffusionProfile(uint32 numLayers, const MPC_LayerSpec* pLayerSpecs,
	const MPC_Options* pOptions, MPC_Output** oppOutput)
{
	if (!oppOutput || !numLayers) return;

	uint32 length = RoundUpPow2(pOptions->desiredLength);
	float stepSize = pOptions->desiredStepSize;

	MatrixProfile mp0(length * 2, length * 2);
	float iorLower = numLayers > 1 ? pLayerSpecs[0].ior / pLayerSpecs[1].ior : pLayerSpecs[0].ior;
	ComputeLayerProfile(pLayerSpecs[0], pLayerSpecs[0].ior, iorLower, stepSize, mp0);
	for (uint32 i = 1; i < numLayers; i++) {
		iorLower = numLayers > i ? pLayerSpecs[i].ior / pLayerSpecs[i + 1].ior : pLayerSpecs[i].ior;
		MatrixProfile mp1(length * 2, length * 2);
		ComputeLayerProfile(pLayerSpecs[i], pLayerSpecs[i].ior / pLayerSpecs[i - 1].ior, iorLower, stepSize, mp1);
		MatrixProfile combined(length * 2, length * 2);
		CombineLayerProfiles(mp0, mp1, combined);
		mp0.reflectance = std::move(combined.reflectance);
		mp0.transmittance = std::move(combined.transmittance);
	}

	MPC_Output* pOut = *oppOutput = new MPC_Output;

	pOut->stepSize = stepSize;
	pOut->length = length;

	pOut->pReflectance = new float[length];
	pOut->pTransmittance = new float[pOut->length];

	uint32 center = length - 1;
	uint32 extent = center;
	for (uint32 i = 0; i <= extent; i++) {
		pOut->pReflectance[i] = mp0.reflectance[center][center + i];
		pOut->pTransmittance[i] = mp0.transmittance[center][center + i];
	}
}

MULTIPOLEPROFILECALCULATOR_API void MPC_FreeOutput(MPC_Output* pOutput) {
	delete [] pOutput->pReflectance;
	delete [] pOutput->pTransmittance;
	delete pOutput;
}
