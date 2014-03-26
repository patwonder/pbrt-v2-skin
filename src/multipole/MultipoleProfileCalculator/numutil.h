
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

template<class Type>
class MTX_Traits {
public:
	static const Type zero;
	static const Type one;
};

template<class Type>
const Type MTX_Traits<Type>::zero = 0;
template<class Type>
const Type MTX_Traits<Type>::one = 1;


template <class Type, class Type_traits=MTX_Traits<Type> >
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
		return *this;
	}
	const Matrix operator+(const Matrix& other) const {
		return Matrix(*this) += other;
	}
	Matrix& operator-=(const Matrix& other) {
		assert(nRows == other.nRows && nCols == other.nCols);
		for (uint32 i = 0; i < nRows * nCols; i++) {
			data[i] -= other.data[i];
		}
		return *this;
	}
	const Matrix operator-(const Matrix& other) const {
		return Matrix(*this) -= other;
	}
	// Component-wise multiplication
	Matrix& operator*=(const Matrix& other) {
		assert(nRows == other.nRows && nCols == other.nCols);
		for (uint32 i = 0; i < nRows * nCols; i++) {
			data[i] *= other.data[i];
		}
		return *this;
	}
	const Matrix operator*(const Matrix& other) const {
		return Matrix(*this) *= other;
	}
	Matrix& operator*=(Type mul) {
		for (uint32 i = 0; i < nRows * nCols; i++) {
			data[i] *= mul;
		}
		return *this;
	}
	const Matrix operator*(Type mul) const {
		return Matrix(*this) *= mul;
	}
	// Component-wise division
	Matrix& operator/=(const Matrix& other) {
		assert(nRows == other.nRows && nCols == other.nCols);
		for (uint32 i = 0; i < nRows * nCols; i++) {
			data[i] /= other.data[i];
		}
		return *this;
	}
	const Matrix operator/(const Matrix& other) const {
		return Matrix(*this) /= other;
	}
	Matrix& operator/=(Type div) {
		for (uint32 i = 0; i < nRows * nCols; i++) {
			data[i] /= div;
		}
		return *this;
	}
	const Matrix operator/(Type div) const {
		return Matrix(*this) /= div;
	}
	Matrix& OneMinusSelf() {
		for (uint32 i = 0; i < nRows * nCols; i++) {
			data[i] = Type_traits::one - data[i];
		}
		return *this;
	}

	Matrix ChangeSize(uint32 newRows, uint32 newCols, uint32 baseRow, uint32 baseCol, uint32 startRow, uint32 startCol) const {
		Matrix ret(newRows, newCols);
		uint32 minRows = min(newRows - baseRow, nRows - startRow), minCols = min(newCols - baseCol, nCols - startCol);
		ret.Clear();
		for (uint32 i = 0; i < minRows; i++) {
			for (uint32 j = 0; j < minCols; j++) {
				ret[baseRow + i][baseCol + j] = (*this)[startRow + i][startCol + j];
			}
		}
		return ret;
	}
	Matrix& ChangeSizeFrom(const Matrix& from, uint32 baseRow, uint32 baseCol, uint32 startRow, uint32 startCol) {
		Clear();
		uint32 minRows = min(from.nRows - startRow, nRows - baseRow), minCols = min(from.nCols - startCol, nCols - baseCol);
		for (uint32 i = 0; i < minRows; i++) {
			for (uint32 j = 0; j < minCols; j++) {
				(*this)[baseRow + i][baseCol + j] = from[startRow + i][startCol + j];
			}
		}
		return *this;
	}
	Matrix ScaleAndShift(uint32 newRows, uint32 newCols, uint32 shRow, uint32 shCol) const {
		Matrix ret(newRows, newCols);
		uint32 minRows = min(nRows, newRows), minCols = min(nCols, newCols);
		ret.Clear();
		for (uint32 i = 0; i < minRows; i++) {
			uint32 ii = newRows - shRow + i; if (ii >= newRows) ii -= newRows;
			for (uint32 j = 0; j < minCols; j++) {
				uint32 jj = newCols - shCol + j; if (jj >= newCols) jj -= newCols;
				ret[ii][jj] = (*this)[i][j];
			}
		}
		return ret;
	}
	Matrix ScaleAndShiftReversed(uint32 newRows, uint32 newCols, uint32 shRow, uint32 shCol) const {
		Matrix ret(newRows, newCols);
		uint32 minRows = min(nRows, newRows), minCols = min(nCols, newCols);
		ret.Clear();
		for (uint32 i = 0; i < minRows; i++) {
			uint32 ii = nRows - shRow + i; if (ii >= nRows) ii -= nRows;
			for (uint32 j = 0; j < minCols; j++) {
				uint32 jj = nRows - shCol + j; if (jj >= nCols) jj -= nCols;
				ret[i][j] = (*this)[ii][jj];
			}
		}
		return ret;
	}
	void swapRow(uint32 nRow1, uint32 nRow2, uint32 startCol = 0) {
		for (uint32 i = startCol; i < nCols; i++) {
			std::swap((*this)[nRow1][i], (*this)[nRow2][i]);
		}
	}
	void swapCol(uint32 nCol1, uint32 nCol2, uint32 startRow = 0) {
		for (uint32 i = startRow; i < nRows; i++) {
			swap((*this)[i][nCol1], (*this)[i][nCol2]);
		}
	}

	Type* GetData() { return data; }
	const Type* GetData() const { return data; }
	uint32 GetNumRows() const { return nRows; }
	uint32 GetNumCols() const { return nCols; }
	void Clear() {
		for (uint32 i = 0; i < nRows * nCols; i++)
			data[i] = Type_traits::zero;
	}
private:
	Type* data;
	uint32 nRows, nCols;
};
