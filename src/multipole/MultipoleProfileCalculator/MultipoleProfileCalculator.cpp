
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

inline kiss_fft_cpx& operator-=(kiss_fft_cpx& cpx, const kiss_fft_cpx& other) {
	cpx.r -= other.r;
	cpx.i -= other.i;
	return cpx;
}

inline kiss_fft_cpx& operator*=(kiss_fft_cpx& cpx, const kiss_fft_cpx& other) {
	kiss_fft_scalar r = cpx.r * other.r - cpx.i * other.i;
	kiss_fft_scalar i = cpx.r * other.i + cpx.i * other.r;
	cpx.r = r; cpx.i = i;
	return cpx;
}

inline kiss_fft_cpx& operator/=(kiss_fft_cpx& cpx, const kiss_fft_cpx& other) {
	kiss_fft_scalar divisor = other.r * other.r + other.i * other.i;
	kiss_fft_scalar acbd = cpx.r * other.r + cpx.i * other.i;
	kiss_fft_scalar bcad = cpx.i * other.r - cpx.r * other.i;
	cpx.r = acbd / divisor;
	cpx.i = bcad / divisor;
	return cpx;
}

inline const kiss_fft_cpx operator-(const kiss_fft_cpx& cpx) {
	kiss_fft_cpx res;
	res.r = -cpx.r;
	res.i = -cpx.i;
	return res;
}

inline const kiss_fft_cpx operator+(const kiss_fft_cpx& cpx, const kiss_fft_cpx& other) {
	return kiss_fft_cpx(cpx) += other;
}

inline const kiss_fft_cpx operator-(const kiss_fft_cpx& cpx, const kiss_fft_cpx& other) {
	return kiss_fft_cpx(cpx) -= other;
}

inline const kiss_fft_cpx operator*(const kiss_fft_cpx& cpx, const kiss_fft_cpx& other) {
	return kiss_fft_cpx(cpx) *= other;
}

inline const kiss_fft_cpx operator/(const kiss_fft_cpx& cpx, const kiss_fft_cpx& other) {
	return kiss_fft_cpx(cpx) /= other;
}

kiss_fft_cpx makecpx(kiss_fft_scalar r, kiss_fft_scalar i = 0) {
	kiss_fft_cpx ret;
	ret.r = r;
	ret.i = i;
	return ret;
}

template <>
class MTX_Traits<kiss_fft_cpx> {
public:
	static const kiss_fft_cpx one;
	static const kiss_fft_cpx zero;
};

const kiss_fft_cpx MTX_Traits<kiss_fft_cpx>::one = makecpx(1);
const kiss_fft_cpx MTX_Traits<kiss_fft_cpx>::zero = makecpx(0);

struct MatrixProfile {
	MatrixProfile(uint32 length)
		: reflectance(length, length), transmittance(length, length) { }
	Matrix<kiss_fft_scalar> reflectance;
	Matrix<kiss_fft_scalar> transmittance;
	uint32 GetLength() const { return reflectance.GetNumRows(); }
	void Clear() { reflectance.Clear(); transmittance.Clear(); }
};


void FFT(const Matrix<kiss_fft_scalar>& profile, Matrix<kiss_fft_cpx>& out) {
	int dims[2] = { profile.GetNumRows(), profile.GetNumCols() };
	kiss_fftndr_cfg cfg = kiss_fftndr_alloc(dims, 2, 0, NULL, NULL);
	kiss_fftndr(cfg, profile.GetData(), out.GetData());
	kiss_fft_free(cfg);
}


void IFFT(const Matrix<kiss_fft_cpx>& profile, Matrix<kiss_fft_scalar>& out) {
	int dims[2] = { out.GetNumRows(), out.GetNumCols() };
	kiss_fftndr_cfg cfg = kiss_fftndr_alloc(dims, 2, 1, NULL, NULL);
	kiss_fftndri(cfg, profile.GetData(), out.GetData());
	kiss_fft_free(cfg);
	out *= (kiss_fft_scalar)1 / (dims[0] * dims[1]);
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

	float thickness = spec.thickness;
	kiss_fft_scalar mfp2 = 2. / (spec.mua + spec.musp);
	// For slabs thicker than 2mfps, multipole approximation is quite accurate.
	// For slabs thinner than that, light becomes less diffuse in the slab.
	// Thus we lerp the reflectance to 0 and transmittance to delta distribution
	// for slabs thinner than 2mfps.
	kiss_fft_scalar lerp = (thickness < mfp2) ? thickness / mfp2 : 1.;
	thickness = max(thickness, (float)mfp2);

	uint32 length = profile.GetLength();
	uint32 center = (length - 1) / 2;
	uint32 extent = center;
	// Using 11 dipole pairs is enough for all normal circumstances,
	// even in the extreme configuration that scattering/absorption ratio is 99.
	// (Planar integral within 20 mfps yields less than 0.01% error)
	const uint32 numDipolePairs = 11;
	kiss_fft_scalar normalizeFactor = stepSize * stepSize;
	for (int32 pair = -((int32)numDipolePairs - 1) / 2; pair <= (int)(numDipolePairs - 1) / 2; pair++) {
		DipoleCalculator dc(iorUpper, iorLower, thickness, spec.mua, spec.musp, pair);
		for (uint32 sampleRow = 0; sampleRow <= extent; sampleRow++) {
			for (uint32 sampleCol = sampleRow; sampleCol <= extent; sampleCol++) {
				kiss_fft_scalar drow2 = sampleRow * sampleRow;
				kiss_fft_scalar dcol2 = sampleCol * sampleCol;
				kiss_fft_scalar r2 = (drow2 + dcol2) * (stepSize * stepSize);
				kiss_fft_scalar rd = dc.Rd((float)r2) * normalizeFactor;
				kiss_fft_scalar td = dc.Td((float)r2) * normalizeFactor;
				profile.reflectance[center + sampleRow][center + sampleCol] += rd;
				profile.transmittance[center + sampleRow][center + sampleCol] += td;
			}
		}
	}
	if (lerp < 1.) {
		for (uint32 sampleRow = 0; sampleRow <= extent; sampleRow++) {
			for (uint32 sampleCol = sampleRow; sampleCol <= extent; sampleCol++) {
				profile.reflectance[center + sampleRow][center + sampleCol] *= lerp;
				profile.transmittance[center + sampleRow][center + sampleCol] *= lerp;
			}
		}
		// Lerp to (approximately) delta distribution
		// This is for discrete convolution to work properly in 2D
		profile.transmittance[center][center] += 1. - lerp;
	}
	// Fill the rest of the matrix through symmetricity
	for (uint32 sampleRow = 1; sampleRow <= extent; sampleRow++) {
		for (uint32 sampleCol = 0; sampleCol < sampleRow; sampleCol++) {
			kiss_fft_scalar rd = profile.reflectance[center + sampleCol][center + sampleRow];
			profile.reflectance[center + sampleRow][center + sampleCol] = rd;
			kiss_fft_scalar td = profile.transmittance[center + sampleCol][center + sampleRow];
			profile.transmittance[center + sampleRow][center + sampleCol] = td;
		}
	}
	for (uint32 sampleRow = 0; sampleRow <= extent; sampleRow++) {
		for (uint32 sampleCol = 0; sampleCol <= extent; sampleCol++) {
			kiss_fft_scalar rd = profile.reflectance[center + sampleRow][center + sampleCol];
			profile.reflectance[center - sampleRow][center + sampleCol] = rd;
			profile.reflectance[center + sampleRow][center - sampleCol] = rd;
			profile.reflectance[center - sampleRow][center - sampleCol] = rd;
			kiss_fft_scalar td = profile.transmittance[center + sampleRow][center + sampleCol];
			profile.transmittance[center - sampleRow][center + sampleCol] = td;
			profile.transmittance[center + sampleRow][center - sampleCol] = td;
			profile.transmittance[center - sampleRow][center - sampleCol] = td;
		}
	}
}


Matrix<kiss_fft_cpx> ToFrequencyDomain(const Matrix<kiss_fft_scalar>& profile) {
	uint32 length = profile.GetNumRows();
	uint32 cl = length * 2; // convolution length
	uint32 center = (length - 1) / 2;
	Matrix<kiss_fft_cpx> out(cl, cl / 2 + 1);
	FFT(profile.ScaleAndShift(cl, cl, center, center), out);
	return out;
}


Matrix<kiss_fft_scalar> ToTimeDomain(const Matrix<kiss_fft_cpx>& freq) {
	uint32 cl = freq.GetNumRows(); // convolution length
	uint32 length = cl / 2;
	uint32 center = (length - 1) / 2;
	Matrix<kiss_fft_scalar> cOut(cl, cl);
	IFFT(freq, cOut);
	return cOut.ScaleAndShiftReversed(length, length, center, center);
}


Matrix<kiss_fft_scalar>& Convolution(Matrix<kiss_fft_scalar>& out, const Matrix<kiss_fft_scalar>& profile) {
	uint32 length = out.GetNumRows();
	uint32 cl = length * 2; // convolution length

	Matrix<kiss_fft_scalar> cOut = out.ChangeSize(cl, cl, 0, 0, 0, 0);
	Matrix<kiss_fft_scalar> cPro = profile.ChangeSize(cl, cl, 0, 0, 0, 0);
	Matrix<kiss_fft_cpx> fout(cl, cl / 2 + 1); FFT(cOut, fout);
	Matrix<kiss_fft_cpx> fpro(cl, cl / 2 + 1); FFT(cPro, fpro);
	fout *= fpro;
	IFFT(fout, cOut);
	uint32 center = (length - 1) / 2;
	out.ChangeSizeFrom(cOut, 0, 0, center, center);

	return out;
}


void CombineLayerProfiles(const MatrixProfile& layer1, const MatrixProfile& layer2,
	MatrixProfile& combined)
{
	// T12 = T1*T2/(1-R2*R1)
	// R12 = R1+T1*R2*T1/(1-R2*R1)
	Matrix<kiss_fft_cpx> fR1 = ToFrequencyDomain(layer1.reflectance);
	Matrix<kiss_fft_cpx> fR2 = ToFrequencyDomain(layer2.reflectance);
	Matrix<kiss_fft_cpx> fT1 = ToFrequencyDomain(layer1.transmittance);
	Matrix<kiss_fft_cpx> fT2 = ToFrequencyDomain(layer2.transmittance);

	Matrix<kiss_fft_cpx> fOneR2R1 = fR2;
	fOneR2R1 *= fR1;
	fOneR2R1.OneMinusSelf();

	Matrix<kiss_fft_cpx> fR12 = fT1;
	fR12 *= fR2;
	fR12 *= fT1;
	fR12 /= fOneR2R1;
	fR12 += fR1;

	combined.reflectance = ToTimeDomain(fR12);

	Matrix<kiss_fft_cpx> fT12 = fT1;
	fT12 *= fT2;
	fT12 /= fOneR2R1;

	combined.transmittance = ToTimeDomain(fT12);
}


MULTIPOLEPROFILECALCULATOR_API void MPC_ComputeDiffusionProfile(uint32 numLayers, const MPC_LayerSpec* pLayerSpecs,
	const MPC_Options* pOptions, MPC_Output** oppOutput)
{
	if (!oppOutput || !numLayers) return;

	uint32 length = RoundUpPow2(pOptions->desiredLength);
	float stepSize = pOptions->desiredStepSize;

	MatrixProfile mp0(length * 2);
	float iorLower = numLayers > 1 ? pLayerSpecs[0].ior / pLayerSpecs[1].ior : pLayerSpecs[0].ior;
	ComputeLayerProfile(pLayerSpecs[0], pLayerSpecs[0].ior, iorLower, stepSize, mp0);
	for (uint32 i = 1; i < numLayers; i++) {
		iorLower = numLayers > (i + 1) ? pLayerSpecs[i].ior / pLayerSpecs[i + 1].ior : pLayerSpecs[i].ior;
		MatrixProfile mp1(length * 2);
		ComputeLayerProfile(pLayerSpecs[i], pLayerSpecs[i].ior / pLayerSpecs[i - 1].ior, iorLower, stepSize, mp1);
		MatrixProfile combined(length * 2);
		CombineLayerProfiles(mp0, mp1, combined);
		mp0.reflectance = std::move(combined.reflectance);
		mp0.transmittance = std::move(combined.transmittance);
	}

	MPC_Output* pOut = *oppOutput = new MPC_Output;

	pOut->length = length;
	
	pOut->pDistance = new float[length];
	pOut->pReflectance = new float[length];
	pOut->pTransmittance = new float[pOut->length];

	uint32 center = length - 1;
	uint32 extent = center;
	float denormalizeFactor = 1.f / (stepSize * stepSize);
	for (uint32 i = 0; i <= extent; i++) {
		pOut->pDistance[i] = i * stepSize;
		pOut->pReflectance[i] = (float)mp0.reflectance[center][center + i] * denormalizeFactor;
		pOut->pTransmittance[i] = (float)mp0.transmittance[center][center + i] * denormalizeFactor;
	}
}

MULTIPOLEPROFILECALCULATOR_API void MPC_FreeOutput(MPC_Output* pOutput) {
	delete [] pOutput->pDistance;
	delete [] pOutput->pReflectance;
	delete [] pOutput->pTransmittance;
	delete pOutput;
}
