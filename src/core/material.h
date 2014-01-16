
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
