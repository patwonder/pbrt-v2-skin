
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

#include "trianglemesh.h"

struct SkinLayer;
class LayeredSkin : public TriangleMesh {
public:
	// LoopSubdiv Public Methods
	LayeredSkin(const Transform *o2w, const Transform *w2o, bool ro,
		int nt, int nv, const int *vi,
		const Point *P, const Normal *N, const Vector *S,
		const float *uv, const Reference<Texture<float> > &atex,
		int nSkinLayers, const SkinLayer layers[]);
	~LayeredSkin();

	void Refine(vector<Reference<Shape> >& refined) const override;
private:
	// LoopSubdiv Private Data
	int nSkinLayers;
	SkinLayer* layers;
};

// Creator function
LayeredSkin* CreateLayeredSkinShape(const Transform *o2w, const Transform *w2o,
	bool reverseOrientation, const ParamSet &params,
	map<string, Reference<Texture<float> > > *floatTextures);
