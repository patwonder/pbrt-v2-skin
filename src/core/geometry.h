
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

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

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CORE_GEOMETRY_H
#define PBRT_CORE_GEOMETRY_H

// core/geometry.h*
#include "pbrt.h"

// Geometry Declarations
template <class scalar>
class VectorBase {
protected:
	typedef ScalarTraits<scalar> Traits;
public:
    // VectorBase Public Methods
    VectorBase() { x = y = z = Traits::zero(); }
    VectorBase(scalar xx, scalar yy, scalar zz)
        : x(xx), y(yy), z(zz) {
        Assert(!HasNaNs());
    }
    bool HasNaNs() const { return Traits::isNaN(x) || Traits::isNaN(y) || Traits::isNaN(z); }
    explicit VectorBase(const PointBase<scalar> &p);
#ifndef NDEBUG
    // The default versions of these are fine for release builds; for debug
    // we define them so that we can add the Assert checks.
    VectorBase(const VectorBase &v) {
        Assert(!v.HasNaNs());
        x = v.x; y = v.y; z = v.z;
    }
    
    VectorBase &operator=(const VectorBase &v) {
        Assert(!v.HasNaNs());
        x = v.x; y = v.y; z = v.z;
        return *this;
    }
#endif // !NDEBUG
    VectorBase operator+(const VectorBase &v) const {
        Assert(!v.HasNaNs());
        return VectorBase(x + v.x, y + v.y, z + v.z);
    }
    
    VectorBase& operator+=(const VectorBase &v) {
        Assert(!v.HasNaNs());
        x += v.x; y += v.y; z += v.z;
        return *this;
    }
    VectorBase operator-(const VectorBase &v) const {
        Assert(!v.HasNaNs());
        return VectorBase(x - v.x, y - v.y, z - v.z);
    }
    
    VectorBase& operator-=(const VectorBase &v) {
        Assert(!v.HasNaNs());
        x -= v.x; y -= v.y; z -= v.z;
        return *this;
    }
	template <class scalartype>
    VectorBase operator*(scalartype f) const {
		return VectorBase(Traits::value(f)*x, Traits::value(f)*y, Traits::value(f)*z);
	}
    
	template <class scalartype>
    VectorBase &operator*=(scalartype f) {
        Assert(!Traits::isNaN(Traits::value(f)));
        x *= Traits::value(f); y *= Traits::value(f); z *= Traits::value(f);
        return *this;
    }
	template <class scalartype>
    VectorBase operator/(scalartype f) const {
        Assert(Traits::value(f) != Traits::zero());
        scalar inv = Traits::one() / Traits::value(f);
        return VectorBase(x * inv, y * inv, z * inv);
    }
    
	template <class scalartype>
    VectorBase &operator/=(scalartype f) {
        Assert(Traits::value(f) != Traits::zero());
        scalar inv = Traits::one() / Traits::value(f);
        x *= inv; y *= inv; z *= inv;
        return *this;
    }
    VectorBase operator-() const { return VectorBase(-x, -y, -z); }
    scalar operator[](int i) const {
        Assert(i >= 0 && i <= 2);
        return (&x)[i];
    }
    
    scalar &operator[](int i) {
        Assert(i >= 0 && i <= 2);
        return (&x)[i];
    }
    scalar LengthSquared() const { return x*x + y*y + z*z; }
    scalar Length() const { return sqrt(LengthSquared()); }
    explicit VectorBase(const NormalBase<scalar> &n);

    bool operator==(const VectorBase &v) const {
        return x == v.x && y == v.y && z == v.z;
    }
    bool operator!=(const VectorBase &v) const {
        return x != v.x || y != v.y || z != v.z;
    }

    // VectorBase Public Data
    scalar x, y, z;

	// VectorBase Public Static Data
	static const VectorBase Zero;
};

template<class scalar>
const VectorBase<scalar> VectorBase<scalar>::Zero(Traits::zero(), Traits::zero(), Traits::zero());

template <class scalar>
class PointBase {
protected:
	typedef ScalarTraits<scalar> Traits;
public:
    // PointBase Public Methods
    PointBase() { x = y = z = Traits::zero(); }
    PointBase(scalar xx, scalar yy, scalar zz)
        : x(xx), y(yy), z(zz) {
        Assert(!HasNaNs());
    }
#ifndef NDEBUG
    PointBase(const PointBase &p) {
        Assert(!p.HasNaNs());
        x = p.x; y = p.y; z = p.z;
    }
    
    PointBase &operator=(const PointBase &p) {
        Assert(!p.HasNaNs());
        x = p.x; y = p.y; z = p.z;
        return *this;
    }
#endif // !NDEBUG
    PointBase operator+(const VectorBase<scalar> &v) const {
        Assert(!v.HasNaNs());
        return PointBase(x + v.x, y + v.y, z + v.z);
    }
    
    PointBase &operator+=(const VectorBase<scalar> &v) {
        Assert(!v.HasNaNs());
        x += v.x; y += v.y; z += v.z;
        return *this;
    }
    VectorBase<scalar> operator-(const PointBase &p) const {
        Assert(!p.HasNaNs());
        return VectorBase<scalar>(x - p.x, y - p.y, z - p.z);
    }
    
    PointBase operator-(const VectorBase<scalar> &v) const {
        Assert(!v.HasNaNs());
        return PointBase(x - v.x, y - v.y, z - v.z);
    }
    
    PointBase &operator-=(const VectorBase<scalar> &v) {
        Assert(!v.HasNaNs());
        x -= v.x; y -= v.y; z -= v.z;
        return *this;
    }
    PointBase &operator+=(const PointBase &p) {
        Assert(!p.HasNaNs());
        x += p.x; y += p.y; z += p.z;
        return *this;
    }
    PointBase operator+(const PointBase &p) const {
        Assert(!p.HasNaNs());
        return PointBase(x + p.x, y + p.y, z + p.z);
    }
	template <class scalartype>
    PointBase operator* (scalartype f) const {
        return PointBase(Traits::value(f)*x, Traits::value(f)*y, Traits::value(f)*z);
    }
	template <class scalartype>
    PointBase &operator*=(scalartype f) {
        x *= Traits::value(f); y *= Traits::value(f); z *= Traits::value(f);
        return *this;
    }
	template <class scalartype>
    PointBase operator/ (scalartype f) const {
        scalar inv = Traits::one()/Traits::value(f);
        return PointBase(inv*x, inv*y, inv*z);
    }
	template <class scalartype>
    PointBase &operator/=(scalartype f) {
        scalar inv = Traits::one()/Traits::value(f);
        x *= inv; y *= inv; z *= inv;
        return *this;
    }
    scalar operator[](int i) const {
        Assert(i >= 0 && i <= 2);
        return (&x)[i];
    }
    
    scalar &operator[](int i) {
        Assert(i >= 0 && i <= 2);
        return (&x)[i];
    }
    bool HasNaNs() const {
        return Traits::isNaN(x) || Traits::isNaN(y) || Traits::isNaN(z);
    }

    bool operator==(const PointBase &p) const {
        return x == p.x && y == p.y && z == p.z;
    }
    bool operator!=(const PointBase &p) const {
        return x != p.x || y != p.y || z != p.z;
    }

    // PointBase Public Data
    scalar x, y, z;

	// PointBase Public Static Data
	static const PointBase Zero;
};

template<class scalar>
const PointBase<scalar> PointBase<scalar>::Zero(Traits::zero(), Traits::zero(), Traits::zero());

template <class scalar>
class NormalBase {
protected:
	typedef ScalarTraits<scalar> Traits;
public:
    // NormalBase Public Methods
    NormalBase() { x = y = z = Traits::zero(); }
    NormalBase(scalar xx, scalar yy, scalar zz)
        : x(xx), y(yy), z(zz) {
        Assert(!HasNaNs());
    }
    NormalBase operator-() const {
        return NormalBase(-x, -y, -z);
    }
    NormalBase operator+ (const NormalBase &n) const {
        Assert(!n.HasNaNs());
        return NormalBase(x + n.x, y + n.y, z + n.z);
    }
    
    NormalBase& operator+=(const NormalBase &n) {
        Assert(!n.HasNaNs());
        x += n.x; y += n.y; z += n.z;
        return *this;
    }
    NormalBase operator- (const NormalBase &n) const {
        Assert(!n.HasNaNs());
        return NormalBase(x - n.x, y - n.y, z - n.z);
    }
    
    NormalBase& operator-=(const NormalBase &n) {
        Assert(!n.HasNaNs());
        x -= n.x; y -= n.y; z -= n.z;
        return *this;
    }
    bool HasNaNs() const {
        return Traits::isNaN(x) || Traits::isNaN(y) || Traits::isNaN(z);
    }
	template <class scalartype>
    NormalBase operator*(scalartype f) const {
        return NormalBase(Traits::value(f)*x, Traits::value(f)*y, Traits::value(f)*z);
    }
	template <class scalartype>
    NormalBase &operator*=(scalartype f) {
        x *= Traits::value(f); y *= Traits::value(f); z *= Traits::value(f);
        return *this;
    }
	template <class scalartype>
    NormalBase operator/(scalartype f) const {
        Assert(Traits::value(f) != Traits::zero());
        scalar inv = Traits::one()/Traits::value(f);
        return NormalBase(x * inv, y * inv, z * inv);
    }
	template <class scalartype>
    NormalBase &operator/=(scalartype f) {
        Assert(Traits::value(f) != Traits::zero());
        scalar inv = Traits::one()/Traits::value(f);
        x *= inv; y *= inv; z *= inv;
        return *this;
    }
    scalar LengthSquared() const { return x*x + y*y + z*z; }
    scalar Length() const        { return sqrt(LengthSquared()); }
    
#ifndef NDEBUG
    NormalBase(const NormalBase &n) {
        Assert(!n.HasNaNs());
        x = n.x; y = n.y; z = n.z;
    }
    
    NormalBase &operator=(const NormalBase &n) {
        Assert(!n.HasNaNs());
        x = n.x; y = n.y; z = n.z;
        return *this;
    }
#endif // !NDEBUG
    explicit NormalBase(const VectorBase<scalar> &v)
      : x(v.x), y(v.y), z(v.z) {
        Assert(!v.HasNaNs());
    }
    scalar operator[](int i) const {
        Assert(i >= 0 && i <= 2);
        return (&x)[i];
    }
    
    scalar &operator[](int i) {
        Assert(i >= 0 && i <= 2);
        return (&x)[i];
    }

    bool operator==(const NormalBase &n) const {
        return x == n.x && y == n.y && z == n.z;
    }
    bool operator!=(const NormalBase &n) const {
        return x != n.x || y != n.y || z != n.z;
    }

    // NormalBase Public Data
    scalar x, y, z;

	// NormalBase Public Static Data
	static const NormalBase Zero;
};

template<class scalar>
const NormalBase<scalar> NormalBase<scalar>::Zero(Traits::zero(), Traits::zero(), Traits::zero());

template <class scalar>
class RayBase {
protected:
	typedef ScalarTraits<scalar> Traits;
public:
    // RayBase Public Methods
    RayBase() : mint(ScalarTraits<scalar>::zero()), maxt(Traits::max()), time(ScalarTraits<scalar>::zero()), depth(0) { }
    RayBase(const PointBase<scalar> &origin, const VectorBase<scalar> &direction,
        scalar start, scalar end = Traits::max(), scalar t = ScalarTraits<scalar>::zero(), int d = 0)
        : o(origin), d(direction), mint(start), maxt(end), time(t), depth(d) { }
    RayBase(const PointBase<scalar> &origin, const VectorBase<scalar> &direction, const RayBase &parent,
        scalar start, scalar end = Traits::max())
        : o(origin), d(direction), mint(start), maxt(end),
          time(parent.time), depth(parent.depth+1) { }
	template <class scalartype>
    PointBase<scalar> operator()(scalartype t) const { return o + d * Traits::value(t); }
    bool HasNaNs() const {
        return (o.HasNaNs() || d.HasNaNs() ||
                Traits::isNaN(mint) || Traits::isNaN(maxt));
    }

    // RayBase Public Data
    PointBase<scalar> o;
    VectorBase<scalar> d;
    mutable scalar mint, maxt;
    scalar time;
    int depth;
};


template <class scalar>
class RayDifferentialBase : public RayBase<scalar> {
public:
    // RayDifferentialBase Public Methods
    RayDifferentialBase() { hasDifferentials = false; }
    RayDifferentialBase(const PointBase<scalar> &org, const VectorBase<scalar> &dir, scalar start,
        scalar end = Traits::max(), scalar t = ScalarTraits<scalar>::zero(), int d = 0)
            : RayBase<scalar>(org, dir, start, end, t, d) {
        hasDifferentials = false;
    }
    RayDifferentialBase(const PointBase<scalar> &org, const VectorBase<scalar> &dir, const RayBase<scalar> &parent,
        scalar start, scalar end = Traits::max())
            : RayBase<scalar>(org, dir, start, end, parent.time, parent.depth+1) {
        hasDifferentials = false;
    }
    explicit RayDifferentialBase(const RayBase<scalar> &ray) : RayBase<scalar>(ray) {
        hasDifferentials = false;
    }
    bool HasNaNs() const {
        return RayBase<scalar>::HasNaNs() ||
           (hasDifferentials && (rxOrigin.HasNaNs() || ryOrigin.HasNaNs() ||
                                 rxDirection.HasNaNs() || ryDirection.HasNaNs()));
    }
	template <class scalartype>
    void ScaleDifferentials(scalartype s) {
        rxOrigin = o + (rxOrigin - o) * Traits::value(s);
        ryOrigin = o + (ryOrigin - o) * Traits::value(s);
        rxDirection = d + (rxDirection - d) * Traits::value(s);
        ryDirection = d + (ryDirection - d) * Traits::value(s);
    }

    // RayDifferentialBase Public Data
    bool hasDifferentials;
    PointBase<scalar> rxOrigin, ryOrigin;
    VectorBase<scalar> rxDirection, ryDirection;
};


template <class scalar>
class BBoxBase {
private:
	typedef ScalarTraits<scalar> Traits;
public:
    // BBoxBase Public Methods
    BBoxBase() {
        pMin = PointBase<scalar>(Traits::max(), Traits::max(), Traits::max());
        pMax = PointBase<scalar>(Traits::negmax(), Traits::negmax(), Traits::negmax());
    }
    BBoxBase(const PointBase<scalar> &p) : pMin(p), pMax(p) { }
    BBoxBase(const PointBase<scalar> &p1, const PointBase<scalar> &p2) {
        pMin = PointBase<scalar>(min(p1.x, p2.x), min(p1.y, p2.y), min(p1.z, p2.z));
        pMax = PointBase<scalar>(max(p1.x, p2.x), max(p1.y, p2.y), max(p1.z, p2.z));
    }
    friend BBoxBase Union(const BBoxBase &b, const PointBase<scalar> &p);
    friend BBoxBase Union(const BBoxBase &b, const BBoxBase &b2);
    bool Overlaps(const BBoxBase &b) const {
        bool x = (pMax.x >= b.pMin.x) && (pMin.x <= b.pMax.x);
        bool y = (pMax.y >= b.pMin.y) && (pMin.y <= b.pMax.y);
        bool z = (pMax.z >= b.pMin.z) && (pMin.z <= b.pMax.z);
        return (x && y && z);
    }
    bool Inside(const PointBase<scalar> &pt) const {
        return (pt.x >= pMin.x && pt.x <= pMax.x &&
                pt.y >= pMin.y && pt.y <= pMax.y &&
                pt.z >= pMin.z && pt.z <= pMax.z);
    }
	template <class scalartype>
    void Expand(scalartype delta) {
        pMin -= VectorBase<scalar>(Traits::value(delta), Traits::value(delta), Traits::value(delta));
        pMax += VectorBase<scalar>(Traits::value(delta), Traits::value(delta), Traits::value(delta));
    }
    scalar SurfaceArea() const {
        VectorBase<scalar> d = pMax - pMin;
        return Traits::value(2) * (d.x * d.y + d.x * d.z + d.y * d.z);
    }
    scalar Volume() const {
        VectorBase<scalar> d = pMax - pMin;
        return d.x * d.y * d.z;
    }
    int MaximumExtent() const {
        VectorBase<scalar> diag = pMax - pMin;
        if (diag.x > diag.y && diag.x > diag.z)
            return 0;
        else if (diag.y > diag.z)
            return 1;
        else
            return 2;
    }
    const PointBase<scalar> &operator[](int i) const;
    PointBase<scalar> &operator[](int i);
	template <class scalar1, class scalar2, class scalar3>
	PointBase<scalar> Lerp(scalar1 tx, scalar2 ty, scalar3 tz) const {
        return PointBase<scalar>(::Lerp(Traits::value(tx), pMin.x, pMax.x), ::Lerp(Traits::value(ty), pMin.y, pMax.y),
                     ::Lerp(Traits::value(tz), pMin.z, pMax.z));
    }
    VectorBase<scalar> Offset(const PointBase<scalar> &p) const {
        return VectorBase<scalar>((p.x - pMin.x) / (pMax.x - pMin.x),
                      (p.y - pMin.y) / (pMax.y - pMin.y),
                      (p.z - pMin.z) / (pMax.z - pMin.z));
    }
    void BoundingSphere(PointBase<scalar> *c, scalar *rad) const;
    bool IntersectP(const RayBase<scalar> &ray, scalar *hitt0 = NULL, scalar *hitt1 = NULL) const;

    bool operator==(const BBoxBase &b) const {
        return b.pMin == pMin && b.pMax == pMax;
    }
    bool operator!=(const BBoxBase &b) const {
        return b.pMin != pMin || b.pMax != pMax;
    }

    // BBoxBase Public Data
    PointBase<scalar> pMin, pMax;
};



// Geometry Inline Functions
template <class scalar>
inline VectorBase<scalar>::VectorBase(const PointBase<scalar> &p)
    : x(p.x), y(p.y), z(p.z) {
    Assert(!HasNaNs());
}


template <class scalar, class scalartype>
inline VectorBase<scalar> operator*(scalartype f, const VectorBase<scalar> &v) { return v*f; }
template <class scalar>
inline float Dot(const VectorBase<scalar> &v1, const VectorBase<scalar> &v2) {
    Assert(!v1.HasNaNs() && !v2.HasNaNs());
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}


template <class scalar>
inline float AbsDot(const VectorBase<scalar> &v1, const VectorBase<scalar> &v2) {
    Assert(!v1.HasNaNs() && !v2.HasNaNs());
    return abs(Dot(v1, v2));
}


template <class scalar>
inline VectorBase<scalar> Cross(const VectorBase<scalar> &v1, const VectorBase<scalar> &v2) {
    Assert(!v1.HasNaNs() && !v2.HasNaNs());
    double v1x = v1.x, v1y = v1.y, v1z = v1.z;
    double v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return VectorBase<scalar>((v1y * v2z) - (v1z * v2y),
                  (v1z * v2x) - (v1x * v2z),
                  (v1x * v2y) - (v1y * v2x));
}


template <class scalar>
inline VectorBase<scalar> Cross(const VectorBase<scalar> &v1, const NormalBase<scalar> &v2) {
    Assert(!v1.HasNaNs() && !v2.HasNaNs());
    double v1x = v1.x, v1y = v1.y, v1z = v1.z;
    double v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return VectorBase<scalar>((v1y * v2z) - (v1z * v2y),
                  (v1z * v2x) - (v1x * v2z),
                  (v1x * v2y) - (v1y * v2x));
}


template <class scalar>
inline VectorBase<scalar> Cross(const NormalBase<scalar> &v1, const VectorBase<scalar> &v2) {
    Assert(!v1.HasNaNs() && !v2.HasNaNs());
    double v1x = v1.x, v1y = v1.y, v1z = v1.z;
    double v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return VectorBase<scalar>((v1y * v2z) - (v1z * v2y),
                  (v1z * v2x) - (v1x * v2z),
                  (v1x * v2y) - (v1y * v2x));
}


template <class scalar>
inline VectorBase<scalar> Normalize(const VectorBase<scalar> &v) { return v / v.Length(); }
template <class scalar>
inline void CoordinateSystem(const VectorBase<scalar> &v1, VectorBase<scalar> *v2, VectorBase<scalar> *v3) {
    if (abs(v1.x) > abs(v1.y)) {
        scalar invLen = ScalarTraits<scalar>::one() / sqrt(v1.x*v1.x + v1.z*v1.z);
        *v2 = VectorBase<scalar>(-v1.z * invLen, ScalarTraits<scalar>::zero(), v1.x * invLen);
    }
    else {
        scalar invLen = ScalarTraits<scalar>::one() / sqrt(v1.y*v1.y + v1.z*v1.z);
        *v2 = VectorBase<scalar>(ScalarTraits<scalar>::zero(), v1.z * invLen, -v1.y * invLen);
    }
    *v3 = Cross(v1, *v2);
}


template <class scalar>
inline scalar Distance(const PointBase<scalar> &p1, const PointBase<scalar> &p2) {
    return (p1 - p2).Length();
}


template <class scalar>
inline scalar DistanceSquared(const PointBase<scalar> &p1, const PointBase<scalar> &p2) {
    return (p1 - p2).LengthSquared();
}


template <class scalar, class scalartype>
inline PointBase<scalar> operator*(scalartype f, const PointBase<scalar> &p) {
    Assert(!p.HasNaNs());
    return p*f;
}


template <class scalar, class scalartype>
inline NormalBase<scalar> operator*(scalartype f, const NormalBase<scalar> &n) {
    return NormalBase<scalar>(ScaVal(f)*n.x,
		ScaVal(f)*n.y, ScaVal(f)*n.z);
}


template <class scalar>
inline NormalBase<scalar> Normalize(const NormalBase<scalar> &n) {
    return n / n.Length();
}


template <class scalar>
inline VectorBase<scalar>::VectorBase(const NormalBase<scalar> &n)
  : x(n.x), y(n.y), z(n.z) {
    Assert(!n.HasNaNs());
}


template <class scalar>
inline scalar Dot(const NormalBase<scalar> &n1, const VectorBase<scalar> &v2) {
    Assert(!n1.HasNaNs() && !v2.HasNaNs());
    return n1.x * v2.x + n1.y * v2.y + n1.z * v2.z;
}


template <class scalar>
inline scalar Dot(const VectorBase<scalar> &v1, const NormalBase<scalar> &n2) {
    Assert(!v1.HasNaNs() && !n2.HasNaNs());
    return v1.x * n2.x + v1.y * n2.y + v1.z * n2.z;
}


template <class scalar>
inline scalar Dot(const NormalBase<scalar> &n1, const NormalBase<scalar> &n2) {
    Assert(!n1.HasNaNs() && !n2.HasNaNs());
    return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
}


template <class scalar>
inline scalar AbsDot(const NormalBase<scalar> &n1, const VectorBase<scalar> &v2) {
    Assert(!n1.HasNaNs() && !v2.HasNaNs());
    return abs(n1.x * v2.x + n1.y * v2.y + n1.z * v2.z);
}


template <class scalar>
inline scalar AbsDot(const VectorBase<scalar> &v1, const NormalBase<scalar> &n2) {
    Assert(!v1.HasNaNs() && !n2.HasNaNs());
    return abs(v1.x * n2.x + v1.y * n2.y + v1.z * n2.z);
}


template <class scalar>
inline scalar AbsDot(const NormalBase<scalar> &n1, const NormalBase<scalar> &n2) {
    Assert(!n1.HasNaNs() && !n2.HasNaNs());
    return abs(n1.x * n2.x + n1.y * n2.y + n1.z * n2.z);
}


template <class scalar>
inline NormalBase<scalar> Faceforward(const NormalBase<scalar> &n, const VectorBase<scalar> &v) {
    return (Dot(n, v) < ScalarTraits<scalar>::zero()) ? -n : n;
}


template <class scalar>
inline NormalBase<scalar> Faceforward(const NormalBase<scalar> &n, const NormalBase<scalar> &n2) {
    return (Dot(n, n2) < ScalarTraits<scalar>::zero()) ? -n : n;
}



template <class scalar>
inline VectorBase<scalar> Faceforward(const VectorBase<scalar> &v, const VectorBase<scalar> &v2) {
    return (Dot(v, v2) < ScalarTraits<scalar>::zero()) ? -v : v;
}



template <class scalar>
inline VectorBase<scalar> Faceforward(const VectorBase<scalar> &v, const NormalBase<scalar> &n2) {
    return (Dot(v, n2) < ScalarTraits<scalar>::zero()) ? -v : v;
}


template <class scalar>
inline const PointBase<scalar> &BBoxBase<scalar>::operator[](int i) const {
    Assert(i == 0 || i == 1);
    return (&pMin)[i];
}



template <class scalar>
inline PointBase<scalar> &BBoxBase<scalar>::operator[](int i) {
    Assert(i == 0 || i == 1);
    return (&pMin)[i];
}


template <class scalar, class scalar1, class scalar2, class scalar3>
inline VectorBase<scalar> SphericalDirection(scalar1 sintheta,
                                 scalar2 costheta, scalar3 phi) {
    return VectorBase<scalar>(ScaVal(sintheta) * cos(ScaVal(phi)),
                  ScaVal(sintheta) * sin(ScaVal(phi)),
                  ScaVal(costheta));
}


template <class scalar, class scalar1, class scalar2, class scalar3>
inline VectorBase<scalar> SphericalDirection(scalar1 sintheta, scalar2 costheta,
                                 scalar3 phi, const VectorBase<scalar> &x,
                                 const VectorBase<scalar> &y, const VectorBase<scalar> &z) {
    return ScaVal(sintheta) * cos(ScaVal(phi)) * x +
           ScaVal(sintheta) * sin(ScaVal(phi)) * y +
		   ScaVal(costheta) * z;
}


template <class scalar>
inline scalar SphericalTheta(const VectorBase<scalar> &v) {
    return acos(Clamp(v.z, -ScalarTraits<scalar>::one(), ScalarTraits<scalar>::one()));
}


template <class scalar>
inline scalar SphericalPhi(const VectorBase<scalar> &v) {
    scalar p = atan2(v.y, v.x);
    return (p < ScalarTraits<scalar>::zero()) ? p + ScaVal(2*M_PI) : p;
}



#endif // PBRT_CORE_GEOMETRY_H
