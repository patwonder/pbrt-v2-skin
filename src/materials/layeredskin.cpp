
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

#include "layeredskin.h"
#include "paramset.h"

LayeredSkin::LayeredSkin(const vector<SkinLayer>& layers)
	: layers(layers)
{

}

vector<float_type> LayeredSkin::GetLayerThickness() const {
	vector<float_type> ret;
	for (const auto& layer : layers) {
		ret.push_back(layer.thickness);
	}
	return ret;
}

BSDF* LayeredSkin::GetBSDF(const DifferentialGeometry &dgGeom,
	const DifferentialGeometry &dgShading,
	MemoryArena &arena) const
{
	float_type ior = layers[0].ior;
	BSDF* bsdf = BSDF_ALLOC(arena, BSDF)(dgShading, dgGeom.nn, ior);
	float rough = 0.4f;
	Fresnel *fresnel = BSDF_ALLOC(arena, FresnelDielectric)(1.f, ior);
	bsdf->Add(BSDF_ALLOC(arena, Microfacet)(Spectrum(1.), fresnel,
		BSDF_ALLOC(arena, Blinn)(1.f / rough)));
	bsdf->Add(BSDF_ALLOC(arena, MicrofacetTransmission)(Spectrum(1.), fresnel,
		BSDF_ALLOC(arena, Blinn)(1.f / rough), ior));
	//bsdf->Add(BSDF_ALLOC(arena, BRDFToBTDF)(BSDF_ALLOC(arena, Microfacet)(Spectrum(1.), fresnel,
	//	BSDF_ALLOC(arena, Blinn)(1.f / rough))));
	//bsdf->Add(BSDF_ALLOC(arena, BRDFToBTDF)(BSDF_ALLOC(arena, MicrofacetTransmission)(Spectrum(1.), fresnel,
	//	BSDF_ALLOC(arena, Blinn)(1.f / rough), ior)));
	return bsdf;
}


BSDF* LayeredSkin::GetLayeredBSDF(int layerIndex,
	const DifferentialGeometry &dgGeom,
	const DifferentialGeometry &dgShading,
	MemoryArena &arena) const
{
	return NULL;
}

BSSRDF* LayeredSkin::GetLayeredBSSRDF(int layerIndex,
	const DifferentialGeometry &dgGeom,
	const DifferentialGeometry &dgShading,
	MemoryArena &arena) const
{
	return NULL;
}

LayeredSkin* CreateLayeredSkinMaterial(const ParamSet& ps)
{
	int nLayers;
	const SkinLayer* layers = ps.FindSkinLayer("layers", &nLayers);
	if (!layers)
		Error("No layers param set for LayeredSkin material.");
	return new LayeredSkin(vector<SkinLayer>(layers, layers + nLayers));
}
