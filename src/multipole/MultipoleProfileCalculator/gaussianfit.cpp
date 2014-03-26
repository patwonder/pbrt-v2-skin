
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
#include "gaussianfit.h"
#include "numutil.h"
#include <vector>

using namespace std;

typedef Matrix<double> Mat;
typedef vector<double> Vec;

static double GaussianFunc(double r, double sigma) {
	// We are not normalizing the Gaussians here.
	// The normalization could be applied after calculating the coefficients of the basis functions.
	return exp(-(r * r) / (2 * sigma * sigma));
}

static Vec LinearSolve(const Mat& a, const Vec& b) {
	assert(a.GetNumRows() == a.GetNumCols() && a.GetNumRows() == b.size());

	Mat aa = a;
	Vec bb = b;
	// Gaussian elimination
	uint32 length = a.GetNumRows();
	for (uint32 row = 0; row < length - 1; row++) {
		// Find max abs value
		uint32 pivotRow = row;
		double maxAbs = abs(a[row][row]);
		for (uint32 row2 = row + 1; row2 < length; row2++) {
			double absval = abs(a[row2][row]);
			if (absval > maxAbs) {
				maxAbs = absval;
				pivotRow = row2;
			}
		}
		if (maxAbs > 0) {
			if (pivotRow != row) {
				aa.swapRow(row, pivotRow, row);
				swap(bb[row], bb[pivotRow]);
			}
			// elimination
			double div = aa[row][row];
			aa[row][row] = 1;
			for (uint32 col = row + 1; col < length; col++) {
				aa[row][col] /= div;
			}
			bb[row] /= div;
			for (uint32 row2 = row + 1; row2 < length; row2++) {
				for (uint32 col2 = row + 1; col2 < length; col2++) {
					aa[row2][col2] -= aa[row2][row] * aa[row][col2];
				}
				bb[row2] -= aa[row2][row] * bb[row];
			}
		}
	}
	// Solve
	Vec res(length);
	for (uint32 row = length - 1; row != (uint32)(-1); row--) {
		double remainder = bb[row];
		for (uint32 col = length - 1; col != row; col--) {
			remainder -= aa[row][col] * res[col];
		}
		res[row] = remainder / aa[row][row];
	}
	return res;
}

MULTIPOLEPROFILECALCULATOR_API void GF_FitSumGaussians(uint32 length, const float* pDistance, const float* pData,
	uint32 nSigmas, const float* pSigmas, GF_Output** oppOutput)
{
	if (!oppOutput) return;

	// Linear least squares fit of pData using Gaussian bases (defined in pSigmas)

	// Prepare the matrix and vector
	Mat X(nSigmas, nSigmas);
	X.Clear();
	Vec Y(nSigmas, 0);

	for (uint32 i = 0; i < length; i++) {
		double xi = pDistance[i];
		double yi = pData[i];
		for (uint32 j = 0; j < nSigmas; j++) {
			double sigma_j = pSigmas[j];
			double Gj = GaussianFunc(xi, sigma_j);
			Y[j] += xi * xi * Gj * yi;
			for (uint32 k = 0; k < nSigmas; k++) {
				double sigma_k = pSigmas[k];
				X[j][k] += xi * xi * Gj * GaussianFunc(xi, sigma_k);
			}
		}
	}

	// Solve XB==Y
	Vec B = LinearSolve(X, Y);

	// Compute error term
	double error = 0;
	double total = 0;
	double prevDistance = 0;
	for (uint32 i = 0; i < length; i++) {
		double xi = pDistance[i];
		double yi = pData[i];
		double fitted = 0;
		for (uint32 j = 0; j < nSigmas; j++) {
			fitted += B[j] * GaussianFunc(xi, pSigmas[j]);
		}
		error += xi * (fitted - yi) * (fitted - yi) * (xi - prevDistance);
		total += xi * yi * yi * (xi - prevDistance);
		prevDistance = pDistance[i];
	}

	// Return result
	GF_Output* pOutput = *oppOutput = new GF_Output;
	pOutput->nCoeffs = nSigmas;
	pOutput->overallError = error / total;
	pOutput->pNormalizedCoeffs = new float[nSigmas];

	for (uint32 k = 0; k < nSigmas; k++) {
		double sigma = pSigmas[k];
		// normalize coeffs
		pOutput->pNormalizedCoeffs[k] = (float)(B[k] * 2 * PI * sigma * sigma);
	}
}

MULTIPOLEPROFILECALCULATOR_API void GF_FreeOutput(GF_Output* pOutput) {
	delete [] pOutput->pNormalizedCoeffs;
	delete pOutput;
}
