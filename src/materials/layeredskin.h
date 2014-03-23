
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

#include "material.h"
#include "skinlayer.h"

class SkinCoefficients;
struct MultipoleProfileData;
class LayeredSkin : public LayeredMaterial {
public:
    // LayeredSkin Public Methods
	LayeredSkin(const vector<SkinLayer>& layers, float roughness,
		float_type nmperunit, const SkinCoefficients& coeff,
		Reference<Texture<Spectrum> > Kr, Reference<Texture<Spectrum> > Kt,
		Reference<Texture<float> > bumpMap, Reference<Texture<Spectrum> > albedo);
	~LayeredSkin();

	vector<float_type> GetLayerThickness() const override;
	BSDF* GetBSDF(const DifferentialGeometry &dgGeom,
		const DifferentialGeometry &dgShading,
		MemoryArena &arena) const override;
	BSSRDF* GetBSSRDF(const DifferentialGeometry &dgGeom,
		const DifferentialGeometry &dgShading,
		MemoryArena &arena) const override;
	const MultipoleBSSRDF* GetMultipoleBSSRDF(
		const DifferentialGeometry &dgGeom,
		const DifferentialGeometry &dgShading,
		MemoryArena &arena) const override;
	bool HasSubsurfaceScattering() const override { return true; }
	BumpMapping GetBumpMapping() const override {
		return BumpMapping(bumpMap);
	}

	BSDF* GetLayeredBSDF(int layerIndex,
		const DifferentialGeometry &dgGeom,
		const DifferentialGeometry &dgShading,
		MemoryArena &arena) const override;
	BSSRDF* GetLayeredBSSRDF(int layerIndex,
		const DifferentialGeometry &dgGeom,
		const DifferentialGeometry &dgShading,
		MemoryArena &arena) const override;
	LayerParam GetLayerParam(int index) const override; 
private:
    // LayeredSkin Private Data
	vector<SkinLayer> layers;
	float roughness;
	float_type nmperunit;
	SkinCoefficients* pcoeff;
	LayerParam lps[3];
	Reference<Texture<Spectrum> > Kr;
	Reference<Texture<Spectrum> > Kt;
	Reference<Texture<float> > bumpMap;
	Reference<Texture<Spectrum> > albedo;
	MultipoleBSSRDFData* preparedBSSRDFData;
	MultipoleProfileData* profileData;
};

// Creator function
LayeredSkin* CreateLayeredSkinMaterial(const ParamSet& ps, const TextureParams& mp);
