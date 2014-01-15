
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

#include "StdAfx.h"
#include "layeredskin.h"
#include "texture.h"
#include "textures/constant.h"
#include "paramset.h"

LayeredSkin::LayeredSkin(const Transform *o2w, const Transform *w2o, bool ro,
	int nt, int nv, const int *vi,
	const Point *P, const Normal *N, const Vector *S,
	const float *uv, const Reference<Texture<float> > &atex,
	int nSkinLayers, const SkinLayer layers[])
 : TriangleMesh(o2w, w2o, ro, nt, nv, vi, P, N, S, uv, atex),
   nSkinLayers(nSkinLayers)
{
	this->layers = new SkinLayer[nSkinLayers];
	memcpy(this->layers, layers, nSkinLayers * sizeof(SkinLayer));
}

LayeredSkin::~LayeredSkin() {
	delete [] layers;
}

void LayeredSkin::Refine(vector<Reference<Shape> >& refined) const {
	TriangleMesh::Refine(refined);
	float totalThickness = 0;
	for (int i = 0; i < nSkinLayers; i++) {
		totalThickness += layers[i].thickness;
		Shrink(totalThickness)->Refine(refined);
	}
}

LayeredSkin* CreateLayeredSkinShape(const Transform *o2w, const Transform *w2o,
	bool reverseOrientation, const ParamSet &params,
	map<string, Reference<Texture<float> > > *floatTextures)
{
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

	int nLayers;
	const SkinLayer* layers = params.FindSkinLayer("layers", &nLayers);

	return new LayeredSkin(o2w, w2o, reverseOrientation, nvi/3, npi, vi, P,
		N, S, uvs, alphaTex, nLayers, layers);
}
