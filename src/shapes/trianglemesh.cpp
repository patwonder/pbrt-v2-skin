
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


// shapes/trianglemesh.cpp*
#include "stdafx.h"
#include "shapes/trianglemesh.h"
#include "texture.h"
#include "textures/constant.h"
#include "paramset.h"
#include "montecarlo.h"


// TriangleMesh Method Definitions
TriangleMesh::TriangleMesh(const Transform *o2w, const Transform *w2o,
        bool ro, int nt, int nv, const int *vi, const Point *P,
        const Normal *N, const Vector *S, const float *uv,
        const Reference<Texture<float> > &atex)
    : ShrinkableShape(o2w, w2o, ro), alphaTexture(atex) {
    ntris = nt;
    nverts = nv;
    vertexIndex = new int[3 * ntris];
    memcpy(vertexIndex, vi, 3 * ntris * sizeof(int));
    // Copy _uv_, _N_, and _S_ vertex data, if present
    if (uv) {
        uvs = new float[2*nverts];
        memcpy(uvs, uv, 2*nverts*sizeof(float));
    }
    else uvs = NULL;
    p = new Point[nverts];
    if (N) {
        n = new Normal[nverts];
        memcpy(n, N, nverts*sizeof(Normal));
    }
    else n = NULL;
    if (S) {
        s = new Vector[nverts];
        memcpy(s, S, nverts*sizeof(Vector));
    }
    else s = NULL;

    // Transform mesh vertices to world space
    for (int i = 0; i < nverts; ++i)
        p[i] = (*ObjectToWorld)(P[i]);
}


TriangleMesh::~TriangleMesh() {
    delete[] vertexIndex;
    delete[] p;
    delete[] s;
    delete[] n;
    delete[] uvs;
}


BBox TriangleMesh::ObjectBound() const {
    BBox objectBounds;
    for (int i = 0; i < nverts; i++)
        objectBounds = Union(objectBounds, (*WorldToObject)(p[i]));
    return objectBounds;
}


BBox TriangleMesh::WorldBound() const {
    BBox worldBounds;
    for (int i = 0; i < nverts; i++)
        worldBounds = Union(worldBounds, p[i]);
    return worldBounds;
}


void TriangleMesh::Refine(vector<Reference<Shape> > &refined) const {
    for (int i = 0; i < ntris; ++i)
        refined.push_back(new Triangle(ObjectToWorld,
                          WorldToObject, ReverseOrientation,
                          (TriangleMesh *)this, i));
}

TriangleMesh::TriangleMesh(const TriangleMesh& mesh)
	: ShrinkableShape(mesh) { }

Reference<ShrinkableShape> TriangleMesh::Shrink(float_type distance) const {
	// Copy Shape subobject
	TriangleMesh* pNewMesh = new TriangleMesh(*this);

	// Copy data members
	pNewMesh->ntris = ntris;
	pNewMesh->nverts = nverts;
	pNewMesh->vertexIndex = new int[3 * ntris];
	memcpy(pNewMesh->vertexIndex, vertexIndex, 3 * ntris * sizeof(int));
	if (n) {
		pNewMesh->n = new Normal[nverts];
		memcpy(pNewMesh->n, n, nverts * sizeof(Normal));
	} else {
		pNewMesh->n = NULL;
	}
	if (s) {
		pNewMesh->s = new Vector[nverts];
		memcpy(pNewMesh->s, s, nverts * sizeof(Vector));
	} else {
		pNewMesh->s = NULL;
	}
	if (uvs) {
		pNewMesh->uvs = new float[2 * nverts];
		memcpy(pNewMesh->uvs, uvs, 2 * nverts * sizeof(float));
	} else {
		pNewMesh->uvs = NULL;
	}
	pNewMesh->alphaTexture = alphaTexture;

	// Copy vertex positions, shrink the model if possible
	pNewMesh->p = new Point[nverts];
	if (n) {
		// Apply shrinking for each vertex
		for (int i = 0; i < nverts; i++) {
			Normal worldNormal = Normalize((*ObjectToWorld)(n[i]));
			pNewMesh->p[i] = p[i] - Vector(worldNormal) * distance;
		}
	} else {
		memcpy(pNewMesh->p, p, nverts * sizeof(Point));
	}

	return pNewMesh;
}


struct BarycentricCoordinate {
	BarycentricCoordinate() {};
	BarycentricCoordinate(float b0, float b1, float b2)
		: b0(b0), b1(b1), b2(b2) {}
	float b0, b1, b2;
	Point Evaluate(const Point& v0, const Point& v1, const Point& v2) const {
		return b0 * v0 + b1 * v1 + b2 * v2;
	}
	Normal Evaluate(const Normal& n0, const Normal& n1, const Normal& n2) const {
		return Normalize(b0 * n0 + b1 * n1 + b2 * n2);
	}
	BarycentricCoordinate Evaluate(const BarycentricCoordinate& b0, const BarycentricCoordinate& b1, const BarycentricCoordinate& b2) const {
		return BarycentricCoordinate(
			this->b0 * b0.b0 + this->b1 * b1.b0 + this->b2 * b2.b0,
			this->b0 * b0.b1 + this->b1 * b1.b1 + this->b2 * b2.b1,
			this->b0 * b0.b2 + this->b1 * b1.b2 + this->b2 * b2.b2
		);
	}
	static BarycentricCoordinate Lerp(float value, const BarycentricCoordinate& b0, const BarycentricCoordinate& b1) {
		return BarycentricCoordinate(
			::Lerp(value, b0.b0, b1.b0),
			::Lerp(value, b0.b1, b1.b1),
			::Lerp(value, b0.b2, b1.b2)
		);
	}
	static const BarycentricCoordinate baryCentric;
};

const BarycentricCoordinate BarycentricCoordinate::baryCentric(1.f / 3.f, 1.f / 3.f, 1.f / 3.f);

void TriangleMesh::TessellateSurfacePoints(float minDist, vector<SurfacePoint>& points) const {
	// Assume tessellator splits a triangle into tf^2 pieces,
	// That's roughtly 1 point per 1/2*(l/tf)^2 area (assume right triangles)
	// Our goal is to achieve 1 point per pi * (minDist/2)^2 area
	// hence tf ~= l / (sqrt(pi * 2) * minDist / 2)
    //          ~= l / minDist * 0.8
	// Actually, precision is not quite needed here, a rough tf will do pretty well.
	for (int itri = 0; itri < ntris; itri++) {
		// Compute tessellation factors
		int vi0 = vertexIndex[itri * 3];
		int vi1 = vertexIndex[itri * 3 + 1];
		int vi2 = vertexIndex[itri * 3 + 2];
		const Point& v0 = p[vi0];
		const Point& v1 = p[vi1];
		const Point& v2 = p[vi2];
		const Normal& n0 = (*ObjectToWorld)(n[vi0]);
		const Normal& n1 = (*ObjectToWorld)(n[vi1]);
		const Normal& n2 = (*ObjectToWorld)(n[vi2]);
		float le0 = (v1 - v2).Length();
		float le1 = (v2 - v0).Length();
		float le2 = (v0 - v1).Length();
		float tfe0 = le0 / minDist * 0.8f;
		float tfe1 = le1 / minDist * 0.8f;
		float tfe2 = le2 / minDist * 0.8f;
		float tfc = floorf((tfe0 + tfe1 + tfe2) / 3.f + .5f);
		tfe0 = floorf(tfe0 + .5f);
		tfe1 = floorf(tfe1 + .5f);
		tfe2 = floorf(tfe2 + .5f);
		// Let tessellator do the work
		tessellator(tfe0, tfe1, tfe2, tfc, [&] (BarycentricCoordinate bv0, BarycentricCoordinate bv1, BarycentricCoordinate bv2) {
			// Compute information about barycentric point
			BarycentricCoordinate bc = BarycentricCoordinate::baryCentric.Evaluate(bv0, bv1, bv2);
			SurfacePoint sp;
			sp.p = bc.Evaluate(v0, v1, v2);
			sp.n = bc.Evaluate(n0, n1, n2);
			sp.area = M_PI * (minDist / 2.f) * (minDist / 2.f);
			sp.rayEpsilon = minDist / 100.f;
			// Add barycentric point of the subdivided triangle to the collection of surface points
			points.push_back(sp);
		});
	}
}


void TriangleMesh::tessellator(float tfe0, float tfe1, float tfe2, float tfc,
		const std::function<void (BarycentricCoordinate bv0,
		BarycentricCoordinate bv1, BarycentricCoordinate bv2)>& domainShader)
{
	typedef BarycentricCoordinate BC;

	// http://fgiesen.wordpress.com/2011/09/06/a-trip-through-the-graphics-pipeline-2011-part-12/
	int itfe0 = max(Ceil2Int(tfe0), 1),
		itfe1 = max(Ceil2Int(tfe1), 1),
		itfe2 = max(Ceil2Int(tfe2), 1),
		itfc = max(Ceil2Int(tfc), 1);
	// Should at least tessellate a little in the center if edges need tessellation
	if (itfe0 > 1 || itfe1 > 1 || itfe2 > 1)
		itfc = max(itfc, 2);

	BC b0(1.f, 0.f, 0.f), b1(0.f, 1.f, 0.f), b2(0.f, 0.f, 1.f), bc(BC::baryCentric);
	// Generate inside tessellated triangles
	int rings = (itfc + 1) / 2;
	for (int r = 0; r < rings - 1; r++) {
		// Central edge factor
		int edgeInner = itfc - (rings - r) * 2;
		if (edgeInner >= 0) {
			int edgeOuter = edgeInner + 2;
			BC b0Inner = BC::Lerp((float)r / rings, bc, b0);
			BC b1Inner = BC::Lerp((float)r / rings, bc, b1);
			BC b2Inner = BC::Lerp((float)r / rings, bc, b2);
			BC b0Outer = BC::Lerp((float)(r + 1) / rings, bc, b0);
			BC b1Outer = BC::Lerp((float)(r + 1) / rings, bc, b1);
			BC b2Outer = BC::Lerp((float)(r + 1) / rings, bc, b2);
			matching(b0Inner, b1Inner, edgeInner, b0Outer, b1Outer, edgeOuter, domainShader);
			matching(b1Inner, b2Inner, edgeInner, b1Outer, b2Outer, edgeOuter, domainShader);
			matching(b2Inner, b0Inner, edgeInner, b2Outer, b0Outer, edgeOuter, domainShader);
		} else {
			// Generate a single triangle at center
			BC b0Outer = BC::Lerp((float)(r + 1) / rings, bc, b0);
			BC b1Outer = BC::Lerp((float)(r + 1) / rings, bc, b1);
			BC b2Outer = BC::Lerp((float)(r + 1) / rings, bc, b2);
			domainShader(b0Outer, b1Outer, b2Outer);
		}
	}
	// Generate the outermost ring of triangles
	int edgeInner = itfc - 2;
	if (edgeInner >= 0) {
		BC b0Inner = BC::Lerp((float)(rings - 1) / rings, bc, b0);
		BC b1Inner = BC::Lerp((float)(rings - 1) / rings, bc, b1);
		BC b2Inner = BC::Lerp((float)(rings - 1) / rings, bc, b2);
		matching(b0Inner, b1Inner, edgeInner, b0, b1, itfe2, domainShader);
		matching(b1Inner, b2Inner, edgeInner, b1, b2, itfe0, domainShader);
		matching(b2Inner, b0Inner, edgeInner, b2, b0, itfe1, domainShader);
	} else {
		// Send the original triangle back - no tessellation needed here
		domainShader(b0, b1, b2);
	}
}


void TriangleMesh::matching(BarycentricCoordinate b0Inner, BarycentricCoordinate b1Inner, int segsInner,
	BarycentricCoordinate b0Outer, BarycentricCoordinate b1Outer, int segsOuter,
	const std::function<void (BarycentricCoordinate bv0, BarycentricCoordinate bv1,
	BarycentricCoordinate bv2)>& domainShader)
{
	typedef BarycentricCoordinate BC;

	int innerPos = 0, outerPos = 0;
	while (innerPos < segsInner && outerPos < segsOuter) {
		BC bInner = BC::Lerp((float)(innerPos) / segsInner, b0Inner, b1Inner);
		BC bOuter = BC::Lerp((float)(outerPos) / segsOuter, b0Outer, b1Outer);
		// Choose inner/outer edge based on "Least Slope Criteria"
		float slopeInner = (innerPos < segsInner) ? 
			fabsf((float)(innerPos + 1) + 1.f - (float)outerPos / segsOuter * (segsInner + 2)) :
			INFINITY;
		float slopeOuter = (outerPos < segsOuter) ?
			fabsf((float)innerPos + 1.f - (float)(outerPos + 1) / segsOuter * (segsInner + 2)) :
			INFINITY;
		if (slopeInner < slopeOuter) {
			// Advance inner edge
			BC bInnerNew = BC::Lerp((float)(innerPos + 1) / segsInner, b0Inner, b1Inner);
			domainShader(bInnerNew, bInner, bOuter);
			innerPos++;
		} else {
			// Advance outer edge
			BC bOuterNew = BC::Lerp((float)(outerPos + 1) / segsOuter, b0Outer, b1Outer);
			domainShader(bInner, bOuter, bOuterNew);
			outerPos++;
		}
	}
}


BBox Triangle::ObjectBound() const {
    // Get triangle vertices in _p1_, _p2_, and _p3_
    const Point &p1 = mesh->p[v[0]];
    const Point &p2 = mesh->p[v[1]];
    const Point &p3 = mesh->p[v[2]];
    return Union(BBox((*WorldToObject)(p1), (*WorldToObject)(p2)),
                 (*WorldToObject)(p3));
}


BBox Triangle::WorldBound() const {
    // Get triangle vertices in _p1_, _p2_, and _p3_
    const Point &p1 = mesh->p[v[0]];
    const Point &p2 = mesh->p[v[1]];
    const Point &p3 = mesh->p[v[2]];
    return Union(BBox(p1, p2), p3);
}


bool Triangle::Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
                         DifferentialGeometry *dg) const {
    PBRT_RAY_TRIANGLE_INTERSECTION_TEST(const_cast<Ray *>(&ray), const_cast<Triangle *>(this));
    // Compute $\VEC{s}_1$

    // Get triangle vertices in _p1_, _p2_, and _p3_
    const Point &p1 = mesh->p[v[0]];
    const Point &p2 = mesh->p[v[1]];
    const Point &p3 = mesh->p[v[2]];
    Vector e1 = p2 - p1;
    Vector e2 = p3 - p1;
    Vector s1 = Cross(ray.d, e2);
    float divisor = Dot(s1, e1);
    
    if (divisor == 0.)
        return false;
    float invDivisor = 1.f / divisor;

    // Compute first barycentric coordinate
    Vector d = ray.o - p1;
    float b1 = Dot(d, s1) * invDivisor;
    if (b1 < 0. || b1 > 1.)
        return false;

    // Compute second barycentric coordinate
    Vector s2 = Cross(d, e1);
    float b2 = Dot(ray.d, s2) * invDivisor;
    if (b2 < 0. || b1 + b2 > 1.)
        return false;

    // Compute _t_ to intersection point
    float t = Dot(e2, s2) * invDivisor;
    if (t < ray.mint || t > ray.maxt)
        return false;

    // Compute triangle partial derivatives
    Vector dpdu, dpdv;
    float uvs[3][2];
    GetUVs(uvs);

    // Compute deltas for triangle partial derivatives
    float du1 = uvs[0][0] - uvs[2][0];
    float du2 = uvs[1][0] - uvs[2][0];
    float dv1 = uvs[0][1] - uvs[2][1];
    float dv2 = uvs[1][1] - uvs[2][1];
    Vector dp1 = p1 - p3, dp2 = p2 - p3;
    float determinant = du1 * dv2 - dv1 * du2;
    if (determinant == 0.f) {
        // Handle zero determinant for triangle partial derivative matrix
        CoordinateSystem(Normalize(Cross(e2, e1)), &dpdu, &dpdv);
    }
    else {
        float invdet = 1.f / determinant;
        dpdu = ( dv2 * dp1 - dv1 * dp2) * invdet;
        dpdv = (-du2 * dp1 + du1 * dp2) * invdet;
    }

    // Interpolate $(u,v)$ triangle parametric coordinates
    float b0 = 1 - b1 - b2;
    float tu = b0*uvs[0][0] + b1*uvs[1][0] + b2*uvs[2][0];
    float tv = b0*uvs[0][1] + b1*uvs[1][1] + b2*uvs[2][1];

    // Test intersection against alpha texture, if present
    if (ray.depth != -1) {
    if (mesh->alphaTexture) {
        DifferentialGeometry dgLocal(ray(t), dpdu, dpdv,
                                     Normal(0,0,0), Normal(0,0,0),
                                     tu, tv, this);
        if (mesh->alphaTexture->Evaluate(dgLocal) == 0.f)
            return false;
    }
    }

    // Fill in _DifferentialGeometry_ from triangle hit
    *dg = DifferentialGeometry(ray(t), dpdu, dpdv,
                               Normal(0,0,0), Normal(0,0,0),
                               tu, tv, this);
    *tHit = t;
    *rayEpsilon = 1e-3f * *tHit;
    PBRT_RAY_TRIANGLE_INTERSECTION_HIT(const_cast<Ray *>(&ray), t);
    return true;
}


bool Triangle::IntersectP(const Ray &ray) const {
    PBRT_RAY_TRIANGLE_INTERSECTIONP_TEST(const_cast<Ray *>(&ray), const_cast<Triangle *>(this));
    // Compute $\VEC{s}_1$

    // Get triangle vertices in _p1_, _p2_, and _p3_
    const Point &p1 = mesh->p[v[0]];
    const Point &p2 = mesh->p[v[1]];
    const Point &p3 = mesh->p[v[2]];
    Vector e1 = p2 - p1;
    Vector e2 = p3 - p1;
    Vector s1 = Cross(ray.d, e2);
    float divisor = Dot(s1, e1);
    
    if (divisor == 0.)
        return false;
    float invDivisor = 1.f / divisor;

    // Compute first barycentric coordinate
    Vector d = ray.o - p1;
    float b1 = Dot(d, s1) * invDivisor;
    if (b1 < 0. || b1 > 1.)
        return false;

    // Compute second barycentric coordinate
    Vector s2 = Cross(d, e1);
    float b2 = Dot(ray.d, s2) * invDivisor;
    if (b2 < 0. || b1 + b2 > 1.)
        return false;

    // Compute _t_ to intersection point
    float t = Dot(e2, s2) * invDivisor;
    if (t < ray.mint || t > ray.maxt)
        return false;

    // Test shadow ray intersection against alpha texture, if present
    if (ray.depth != -1 && mesh->alphaTexture) {
        // Compute triangle partial derivatives
        Vector dpdu, dpdv;
        float uvs[3][2];
        GetUVs(uvs);

        // Compute deltas for triangle partial derivatives
        float du1 = uvs[0][0] - uvs[2][0];
        float du2 = uvs[1][0] - uvs[2][0];
        float dv1 = uvs[0][1] - uvs[2][1];
        float dv2 = uvs[1][1] - uvs[2][1];
        Vector dp1 = p1 - p3, dp2 = p2 - p3;
        float determinant = du1 * dv2 - dv1 * du2;
        if (determinant == 0.f) {
            // Handle zero determinant for triangle partial derivative matrix
            CoordinateSystem(Normalize(Cross(e2, e1)), &dpdu, &dpdv);
        }
        else {
            float invdet = 1.f / determinant;
            dpdu = ( dv2 * dp1 - dv1 * dp2) * invdet;
            dpdv = (-du2 * dp1 + du1 * dp2) * invdet;
        }

        // Interpolate $(u,v)$ triangle parametric coordinates
        float b0 = 1 - b1 - b2;
        float tu = b0*uvs[0][0] + b1*uvs[1][0] + b2*uvs[2][0];
        float tv = b0*uvs[0][1] + b1*uvs[1][1] + b2*uvs[2][1];
        DifferentialGeometry dgLocal(ray(t), dpdu, dpdv,
                                     Normal(0,0,0), Normal(0,0,0),
                                     tu, tv, this);
        if (mesh->alphaTexture->Evaluate(dgLocal) == 0.f)
            return false;
    }
    PBRT_RAY_TRIANGLE_INTERSECTIONP_HIT(const_cast<Ray *>(&ray), t);
    return true;
}


float Triangle::Area() const {
    // Get triangle vertices in _p1_, _p2_, and _p3_
    const Point &p1 = mesh->p[v[0]];
    const Point &p2 = mesh->p[v[1]];
    const Point &p3 = mesh->p[v[2]];
    return 0.5f * Cross(p2-p1, p3-p1).Length();
}


void Triangle::GetShadingGeometry(const Transform &obj2world,
        const DifferentialGeometry &dg,
        DifferentialGeometry *dgShading) const {
    if (!mesh->n && !mesh->s) {
        *dgShading = dg;
        return;
    }
    // Initialize _Triangle_ shading geometry with _n_ and _s_

    // Compute barycentric coordinates for point
    float b[3];

    // Initialize _A_ and _C_ matrices for barycentrics
    float uv[3][2];
    GetUVs(uv);
    float A[2][2] =
        { { uv[1][0] - uv[0][0], uv[2][0] - uv[0][0] },
          { uv[1][1] - uv[0][1], uv[2][1] - uv[0][1] } };
    float C[2] = { dg.u - uv[0][0], dg.v - uv[0][1] };
    if (!SolveLinearSystem2x2(A, C, &b[1], &b[2])) {
        // Handle degenerate parametric mapping
        b[0] = b[1] = b[2] = 1.f/3.f;
    }
    else
        b[0] = 1.f - b[1] - b[2];

    // Use _n_ and _s_ to compute shading tangents for triangle, _ss_ and _ts_
    Normal ns;
    Vector ss, ts;
    if (mesh->n) ns = Normalize(obj2world(b[0] * mesh->n[v[0]] +
                                          b[1] * mesh->n[v[1]] +
                                          b[2] * mesh->n[v[2]]));
    else   ns = dg.nn;
    if (mesh->s) ss = Normalize(obj2world(b[0] * mesh->s[v[0]] +
                                          b[1] * mesh->s[v[1]] +
                                          b[2] * mesh->s[v[2]]));
    else   ss = Normalize(dg.dpdu);
    
    ts = Cross(ss, ns);
    if (ts.LengthSquared() > 0.f) {
        ts = Normalize(ts);
        ss = Cross(ts, ns);
    }
    else
        CoordinateSystem((Vector)ns, &ss, &ts);
    Normal dndu, dndv;

    // Compute $\dndu$ and $\dndv$ for triangle shading geometry
    if (mesh->n) {
        float uvs[3][2];
        GetUVs(uvs);
        // Compute deltas for triangle partial derivatives of normal
        float du1 = uvs[0][0] - uvs[2][0];
        float du2 = uvs[1][0] - uvs[2][0];
        float dv1 = uvs[0][1] - uvs[2][1];
        float dv2 = uvs[1][1] - uvs[2][1];
        Normal dn1 = mesh->n[v[0]] - mesh->n[v[2]];
        Normal dn2 = mesh->n[v[1]] - mesh->n[v[2]];
        float determinant = du1 * dv2 - dv1 * du2;
        if (determinant == 0.f)
            dndu = dndv = Normal(0,0,0);
        else {
            float invdet = 1.f / determinant;
            dndu = ( dv2 * dn1 - dv1 * dn2) * invdet;
            dndv = (-du2 * dn1 + du1 * dn2) * invdet;
        }
    }
    else
        dndu = dndv = Normal(0,0,0);
    *dgShading = DifferentialGeometry(dg.p, ss, ts,
        (*ObjectToWorld)(dndu), (*ObjectToWorld)(dndv),
        dg.u, dg.v, dg.shape);
    dgShading->dudx = dg.dudx;  dgShading->dvdx = dg.dvdx;
    dgShading->dudy = dg.dudy;  dgShading->dvdy = dg.dvdy;
    dgShading->dpdx = dg.dpdx;  dgShading->dpdy = dg.dpdy;
}


TriangleMesh *CreateTriangleMeshShape(const Transform *o2w, const Transform *w2o,
        bool reverseOrientation, const ParamSet &params,
        map<string, Reference<Texture<float> > > *floatTextures) {
    int nvi, npi, nuvi, nsi, nni;
    const int *vi = params.FindInt("indices", &nvi);
    const Point *P = params.FindPoint("P", &npi);
    const float *uvs = params.FindFloat("uv", &nuvi);
    if (!uvs) uvs = params.FindFloat("st", &nuvi);
    bool discardDegnerateUVs = params.FindOneBool("discarddegenerateUVs", false);
    // XXX should complain if uvs aren't an array of 2...
    if (uvs) {
        if (nuvi < 2 * npi) {
            Error("Not enough of \"uv\"s for triangle mesh.  Expencted %d, "
                  "found %d.  Discarding.", 2*npi, nuvi);
            uvs = NULL;
        }
        else if (nuvi > 2 * npi)
            Warning("More \"uv\"s provided than will be used for triangle "
                    "mesh.  (%d expcted, %d found)", 2*npi, nuvi);
    }
    if (!vi || !P) return NULL;
    const Vector *S = params.FindVector("S", &nsi);
    if (S && nsi != npi) {
        Error("Number of \"S\"s for triangle mesh must match \"P\"s");
        S = NULL;
    }
    const Normal *N = params.FindNormal("N", &nni);
    if (N && nni != npi) {
        Error("Number of \"N\"s for triangle mesh must match \"P\"s");
        N = NULL;
    }
    if (discardDegnerateUVs && uvs && N) {
        // if there are normals, check for bad uv's that
        // give degenerate mappings; discard them if so
        const int *vp = vi;
        for (int i = 0; i < nvi; i += 3, vp += 3) {
            float area = .5f * Cross(P[vp[0]]-P[vp[1]], P[vp[2]]-P[vp[1]]).Length();
            if (area < 1e-7) continue; // ignore degenerate tris.
            if ((uvs[2*vp[0]] == uvs[2*vp[1]] &&
                uvs[2*vp[0]+1] == uvs[2*vp[1]+1]) ||
                (uvs[2*vp[1]] == uvs[2*vp[2]] &&
                uvs[2*vp[1]+1] == uvs[2*vp[2]+1]) ||
                (uvs[2*vp[2]] == uvs[2*vp[0]] &&
                uvs[2*vp[2]+1] == uvs[2*vp[0]+1])) {
                Warning("Degenerate uv coordinates in triangle mesh.  Discarding all uvs.");
                uvs = NULL;
                break;
            }
        }
    }
    for (int i = 0; i < nvi; ++i)
        if (vi[i] >= npi) {
            Error("trianglemesh has out of-bounds vertex index %d (%d \"P\" values were given",
                vi[i], npi);
            return NULL;
        }

    Reference<Texture<float> > alphaTex = NULL;
    string alphaTexName = params.FindTexture("alpha");
    if (alphaTexName != "") {
        if (floatTextures->find(alphaTexName) != floatTextures->end())
            alphaTex = (*floatTextures)[alphaTexName];
        else
            Error("Couldn't find float texture \"%s\" for \"alpha\" parameter",
                  alphaTexName.c_str());
    }
    else if (params.FindOneFloat("alpha", 1.f) == 0.f)
        alphaTex = new ConstantTexture<float>(0.f);
    return new TriangleMesh(o2w, w2o, reverseOrientation, nvi/3, npi, vi, P,
        N, S, uvs, alphaTex);
}


Point Triangle::Sample(float u1, float u2, Normal *Ns) const {
    float b1, b2;
    UniformSampleTriangle(u1, u2, &b1, &b2);
    // Get triangle vertices in _p1_, _p2_, and _p3_
    const Point &p1 = mesh->p[v[0]];
    const Point &p2 = mesh->p[v[1]];
    const Point &p3 = mesh->p[v[2]];
    Point p = b1 * p1 + b2 * p2 + (1.f - b1 - b2) * p3;
    Normal n = Normal(Cross(p2-p1, p3-p1));
    *Ns = Normalize(n);
    if (ReverseOrientation) *Ns *= -1.f;
    return p;
}


