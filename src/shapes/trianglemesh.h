
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

#ifndef PBRT_SHAPES_TRIANGLEMESH_H
#define PBRT_SHAPES_TRIANGLEMESH_H

// shapes/trianglemesh.h*
#include "shape.h"
#include "renderers/surfacepoints.h"
#include <map>
#include <functional>
using std::map;

struct BarycentricCoordinate;
class BumpMapping;

// TriangleMesh Declarations
class TriangleMesh : public ShrinkableShape, public Tessellatable {
public:
    // TriangleMesh Public Methods
    TriangleMesh(const Transform *o2w, const Transform *w2o, bool ro,
                 int ntris, int nverts, const int *vptr,
                 const Point *P, const Normal *N, const Vector *S,
                 const float *uv, const Reference<Texture<float> > &atex);
    ~TriangleMesh();
    BBox ObjectBound() const;
    BBox WorldBound() const;
    bool CanIntersect() const { return false; }
    void Refine(vector<Reference<Shape> > &refined) const;
	Reference<ShrinkableShape> Shrink(float_type distance) const override;
	void TessellateSurfacePoints(float minDist, const BumpMapping& bump,
		vector<SurfacePoint>& points, ProgressReporter* pr = NULL) const override;
	int GetTessellationWork() const override;
	template<class MeshReferenceType>
	friend class TriangleBase;
    template <typename T> friend class VertexTexture;
protected:
    // TriangleMesh Protected Data
    int ntris, nverts;
    int *vertexIndex;
    Point *p;
    Normal *n;
    Vector *s;
    float *uvs;
    Reference<Texture<float> > alphaTexture;
private:
	// TriangleMesh Private Methods
	// Just obtain a copy of the base class subobject
	TriangleMesh(const TriangleMesh&);
	// Fixed-function tessellator similar to DX11's
	static void tessellator(float tfe0, float tfe1, float tfe2, float tfc,
		const std::function<void (BarycentricCoordinate bv0,
		BarycentricCoordinate bv1, BarycentricCoordinate bv2)>& domainShader);
	static void matching(BarycentricCoordinate b0Inner, BarycentricCoordinate b1Inner, int segsInner,
		BarycentricCoordinate b0Outer, BarycentricCoordinate b1Outer, int segsOuter,
		const std::function<void (BarycentricCoordinate bv0, BarycentricCoordinate bv1,
		BarycentricCoordinate bv2)>& domainShader);
};


template<class MeshReferenceType>
class TriangleBase : public Shape {
public:
    // Triangle Public Methods
    TriangleBase(const Transform *o2w, const Transform *w2o, bool ro,
                 const TriangleMesh *m, int n)
        : Shape(o2w, w2o, ro) {
        mesh = m;
        v = &mesh->vertexIndex[3*n];
        PBRT_CREATED_TRIANGLE(this);
    }
    BBox ObjectBound() const;
    BBox WorldBound() const;
    bool Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
                   DifferentialGeometry *dg) const;
    bool IntersectP(const Ray &ray) const;
    void GetUVs(float uv[3][2]) const {
        if (mesh->uvs) {
            uv[0][0] = mesh->uvs[2*v[0]];
            uv[0][1] = mesh->uvs[2*v[0]+1];
            uv[1][0] = mesh->uvs[2*v[1]];
            uv[1][1] = mesh->uvs[2*v[1]+1];
            uv[2][0] = mesh->uvs[2*v[2]];
            uv[2][1] = mesh->uvs[2*v[2]+1];
        }
        else {
            uv[0][0] = 0.; uv[0][1] = 0.;
            uv[1][0] = 1.; uv[1][1] = 0.;
            uv[2][0] = 1.; uv[2][1] = 1.;
        }
    }
    float Area() const;
    virtual void GetShadingGeometry(const Transform &obj2world,
            const DifferentialGeometry &dg,
            DifferentialGeometry *dgShading) const;
	void GetDifferentialGeometries(const BarycentricCoordinate& bc,
		DifferentialGeometry* dgGeom, DifferentialGeometry* dgShading) const;
    Point Sample(float u1, float u2, Normal *Ns) const;
private:
    // Triangle Private Data
    MeshReferenceType mesh;
    int *v;
};


typedef TriangleBase<Reference<const TriangleMesh> > Triangle;
typedef TriangleBase<const TriangleMesh*> TempTriangle;

#include "trianglemesh.inl"


TriangleMesh *CreateTriangleMeshShape(const Transform *o2w, const Transform *w2o,
    bool reverseOrientation, const ParamSet &params,
    map<string, Reference<Texture<float> > > *floatTextures = NULL);

#endif // PBRT_SHAPES_TRIANGLEMESH_H
