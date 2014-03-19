
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

#pragma once

#include <cfloat>
#include <cassert>

template <class T>
inline T clamp(T value, T low, T high) {
	if (value <= low) return low;
	if (value >= high) return high;
	return value;
}

static const float PI = 3.141592654f;
static const float INV_FOURPI = 0.25f / PI;

template <class Type>
class Matrix {
public:
	Matrix(uint32 rows, uint32 cols) {
		data = new Type[rows * cols];
		nRows = rows; nCols = cols;
	}
	~Matrix() {
		delete [] data;
	}
	Matrix(const Matrix& other) {
		nRows = other.nRows;
		nCols = other.nCols;
		data = new Type[nRows * nCols];
		memcpy(data, other.data, sizeof(Type) * nRows * nCols);
	}
	Matrix(Matrix&& other) {
		nRows = other.nRows;
		nCols = other.nCols;
		data = other.data;
		other.data = NULL;
	}
	Matrix& operator=(const Matrix& other) {
		if (this == &other) return *this;

		nRows = other.nRows;
		nCols = other.nCols;
		delete[] data;
		data = new Type[nRows * nCols];
		memcpy(data, other.data, sizeof(Type) * nRows * nCols);

		return *this;
	}

	Matrix& operator=(Matrix&& other) {
		if (this == &other) return *this;

		nRows = other.nRows;
		nCols = other.nCols;
		Type* p = data;
		data = other.data;
		other.data = p;

		return *this;
	}

	Type* operator[](uint32 row) {
		return data + row * nCols;
	}
	const Type* operator[](uint32 row) const {
		return data + row * nCols;
	}
	Matrix& operator+=(const Matrix& other) {
		assert(nRows == other.nRows && nCols == other.nCols);
		for (uint32 i = 0; i < nRows * nCols; i++) {
			data[i] += other.data[i];
		}
	}
	Matrix operator+(const Matrix& other) {
		return Matrix(*this) += other;
	}
	// Component-wise multiplication
	Matrix& operator*=(const Matrix& other) {
		assert(nRows == other.nRows && nCols == other.nCols);
		for (uint32 i = 0; i < nRows * nCols; i++) {
			data[i] *= other.data[i];
		}
	}
	Matrix operator*(const Matrix& other) {
		return Matrix(*this) *= other;
	}
	// Component-wise division
	Matrix& operator/=(const Matrix& other) {
		assert(nRows == other.nRows && nCols == other.nCols);
		for (uint32 i = 0; i < nRows * nCols; i++) {
			data[i] /= other.data[i];
		}
	}
	Matrix operator/(const Matrix& other) {
		return Matrix(*this) /= other;
	}

	Matrix ChangeSize(uint32 newRows, uint32 newCols, uint32 baseRow, uint32 baseCol, uint32 startRow, uint32 startCol) {
		Matrix ret(newRows, newCols);
		uint32 minRows = min(newRows - baseRow, nRows), minCols = min(newCols - baseCol, nCols);
		for (uint32 i = 0; i < newRows * newCols; i++)
			ret.data[i] = Type();
		for (uint32 i = startRow; i < minRows; i++) {
			for (uint32 j = startCol; j < minCols; j++) {
				ret[baseRow + i][baseCol + j] = (*this)[i][j];
			}
		}
	}

	Type* GetData() { return data; }
	const Type* GetData() const { return data; }
	uint32 GetNumRows() const { return nRows; }
	uint32 GetNumCols() const { return nCols; }
	void Clear() {
		for (uint32 i = 0; i < nRows * nCols; i++)
			data[i] = Type();
	}
private:
	Type* data;
	uint32 nRows, nCols;
};
