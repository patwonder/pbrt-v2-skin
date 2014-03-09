
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

#ifndef PBRT_CORE_MATERIAL_H
#define PBRT_CORE_MATERIAL_H

// core/material.h*
#include "pbrt.h"
#include "memory.h"
#include "spectrum.h"
#include <vector>
using std::vector;

// Material Declarations
class Material : public ReferenceCounted {
public:
    // Material Interface
    virtual BSDF *GetBSDF(const DifferentialGeometry &dgGeom,
                          const DifferentialGeometry &dgShading,
                          MemoryArena &arena) const = 0;
    virtual BSSRDF *GetBSSRDF(const DifferentialGeometry &dgGeom,
                              const DifferentialGeometry &dgShading,
                              MemoryArena &arena) const {
        return NULL;
    }
    virtual ~Material();
    static void Bump(const Reference<Texture<float> > &d, const DifferentialGeometry &dgGeom,
        const DifferentialGeometry &dgShading, DifferentialGeometry *dgBump);
};


// Wavelength dependent value
const int WLD_nSamples = 31;
extern const float WLD_lambdas[WLD_nSamples];
class WLDValue : public CoefficientSpectrum<WLD_nSamples> {
public:
	typedef CoefficientSpectrum<WLD_nSamples> Base;

	Spectrum toSpectrum() const {
		return Spectrum::FromSampled(WLD_lambdas, c, WLD_nSamples);
	}
	float& operator[](uint32_t index) {
		return c[index];
	}
	const float& operator[](uint32_t index) const {
		return c[index];
	}

	WLDValue(float v = 0.f) : Base(v) { }
    WLDValue(const Base& v) : Base(v) { }

	static WLDValue FromSampled(const float* lambdas, const float* vals, int n) {
		WLDValue res;

		// Use linear interpolation to get value at each wavelength of interest
		int idxSampledLambda = 0;
		float lambda0 = 0.f, lambda1 = lambdas[0];
		float val0 = vals[0], val1 = vals[0];
		for (int i = 0; i < WLD_nSamples; i++) {
			float lambda = WLD_lambdas[i];
			while (lambda > lambda1 && idxSampledLambda < n - 1) {
				++idxSampledLambda;
				lambda0 = lambda1;
				lambda1 = lambdas[idxSampledLambda];
				val0 = val1;
				val1 = vals[idxSampledLambda];
			}
			if (lambda <= lambda1) {
				res[i] = Lerp((lambda - lambda0) / (lambda1 - lambda0), val0, val1);
			} else
				res[i] = val1;
		}

		return res;
	}
};

struct LayerParam {
	LayerParam() {
		mua = musp = WLDValue(0.f);
		ga = b = 0.f;
		isotropicHGPF = false;
	}

	// Thickness
	float thickness;
	// Absorption and (reduced) scattering coefficient 
	WLDValue mua, musp;
	// Anisotropy and isotropy associated with HGPF
	float ga, b;
	// Whether to use modified version of HGPF that accounts for isotropic scattering
	bool isotropicHGPF;
};


class LayeredMaterial : public Material {
public:
    // LayeredMaterial Interface
	virtual vector<float_type> GetLayerThickness() const = 0;
	virtual BSDF* GetLayeredBSDF(int layerIndex,
		const DifferentialGeometry &dgGeom,
		const DifferentialGeometry &dgShading,
		MemoryArena &arena) const = 0;
	virtual BSSRDF* GetLayeredBSSRDF(int layerIndex,
		const DifferentialGeometry &dgGeom,
		const DifferentialGeometry &dgShading,
		MemoryArena &arena) const
	{
		return NULL;
	}
	virtual LayerParam GetLayerParam(int layerIndex) const = 0;
};

class LayeredMaterialWrapper : public Material {
public:
	LayeredMaterialWrapper(const Reference<LayeredMaterial>& layeredMaterial,
		int layerIndex);
	int GetLayerIndex() const;
	BSDF* GetBSDF(const DifferentialGeometry &dgGeom,
		const DifferentialGeometry &dgShading,
		MemoryArena &arena) const override;
	BSSRDF* GetBSSRDF(const DifferentialGeometry &dgGeom,
		const DifferentialGeometry &dgShading,
		MemoryArena &arena) const override;
private:
	int layerIndex;
	Reference<LayeredMaterial> layeredMaterial;
};


#endif // PBRT_CORE_MATERIAL_H
