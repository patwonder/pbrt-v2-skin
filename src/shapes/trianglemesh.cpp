
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
#include "progressreporter.h"


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
                          this, i));
}

TriangleMesh::TriangleMesh(const TriangleMesh& mesh)
	: ShrinkableShape(mesh) { }

Reference<ShrinkableShape> TriangleMesh::Shrink(float distance) const {
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

void TriangleMesh::TessellateSurfacePoints(float minDist, const BumpMapping& bump,
	uint32_t materialId, vector<SurfacePoint>& points, ProgressReporter* pr) const
{
	// Assume tessellator splits a triangle into tf^2 pieces,
	// That's roughtly 1 point per 1/2*(l/tf)^2 area (assume right triangles)
	// Our goal is to achieve 1 point per pi * (minDist/2)^2 area
	// hence tf ~= l / (sqrt(pi * 2) * minDist / 2)
    //          ~= l / minDist * 0.8
	// Actually, precision is not quite needed here, a rough tf will do pretty well.
	int block = max(ntris / 100, 100);
	for (int itri = 0; itri < ntris; itri++) {
		// Compute tessellation factors
		int vi0 = vertexIndex[itri * 3];
		int vi1 = vertexIndex[itri * 3 + 1];
		int vi2 = vertexIndex[itri * 3 + 2];
		const Point& v0 = p[vi0];
		const Point& v1 = p[vi1];
		const Point& v2 = p[vi2];
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
		TempTriangle triangle(ObjectToWorld, WorldToObject, ReverseOrientation, this, itri);
		// Let tessellator do the work
		tessellator(tfe0, tfe1, tfe2, tfc, [&] (BarycentricCoordinate bv0, BarycentricCoordinate bv1, BarycentricCoordinate bv2) {
			Point sv0 = bv0.Evaluate(v0, v1, v2);
			Point sv1 = bv1.Evaluate(v0, v1, v2);
			Point sv2 = bv2.Evaluate(v0, v1, v2);
#if 1
			// Compute information about barycentric point
			BarycentricCoordinate bc = BarycentricCoordinate::baryCentric.Evaluate(bv0, bv1, bv2);
#else
			// Compute information about incenter point
			float l0 = (sv1 - sv2).Length();
			float l1 = (sv2 - sv0).Length();
			float l2 = (sv0 - sv1).Length();
			BarycentricCoordinate bic(l0 / (l0 + l1 + l2), l1 / (l0 + l1 + l2), l2 / (l0 + l1 + l2));
			BarycentricCoordinate bc = bic.Evaluate(bv0, bv1, bv2);
#endif
			SurfacePoint sp;
			sp.p = bc.Evaluate(v0, v1, v2);
			DifferentialGeometry dgGeom, dgShading, dgBump;
			triangle.GetDifferentialGeometries(bc, &dgGeom, &dgShading);
#if 1
			bump.Bump(dgGeom, dgShading, &dgBump);
#else
			dgBump = dgShading;
#endif
			sp.n = dgBump.nn;
			sp.u = dgBump.u;
			sp.v = dgBump.v;
			sp.materialId = materialId;
			// Compute triangle area
			sp.area = .5f * Cross(sv1 - sv0, sv2 - sv0).Length();
			sp.rayEpsilon = minDist / 10.f;
			// Add barycentric point of the subdivided triangle to the collection of surface points
			points.push_back(sp);
		});
		if (pr && (itri + 1) % block == 0)
			pr->Update(block);
	}
	if (pr)
		pr->Update(ntris % block);
}


int TriangleMesh::GetTessellationWork() const {
	return ntris;
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
	while (innerPos < segsInner || outerPos < segsOuter) {
		BC bInner = segsInner ? BC::Lerp((float)(innerPos) / segsInner, b0Inner, b1Inner) : b0Inner;
		BC bOuter = BC::Lerp((float)(outerPos) / segsOuter, b0Outer, b1Outer); // segsOuter should always > 0
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

