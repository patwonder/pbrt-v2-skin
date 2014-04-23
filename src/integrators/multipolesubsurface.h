
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

#include "integrator.h"
#include "irradiancepoint.h"
struct SubsurfaceOctreeNode;

class MultipoleSubsurfaceIntegrator : public SurfaceIntegrator {
public:
	// MultipoleSubsurfaceIntegrator Public Methods
	void Preprocess(const Scene *scene, const Camera *camera,
		const Renderer *renderer) override;
	void RequestSamples(Sampler *sampler, Sample *sample,
		const Scene *scene) override;
	Spectrum Li(const Scene *scene, const Renderer *renderer,
		const RayDifferential &ray, const Intersection &isect,
		const Sample *sample, RNG &rng, MemoryArena &arena) const override;
	MultipoleSubsurfaceIntegrator(int mdepth, float merror, float mindist,
		const string &fn, const vector<Reference<Primitive> >* ops, float mix,
		bool showIrradiancePoints, bool usePoissonPointFinder)
		: originalPrimitives(*ops), mix(mix), showIrradiancePoints(showIrradiancePoints),
		  usePoissonPointFinder(usePoissonPointFinder)
	{
		maxSpecularDepth = mdepth;
		maxError = merror;
		minSampleDist = mindist;
		filename = fn;
		octree = NULL;
		lightSampleOffsets = NULL;
		bsdfSampleOffsets = NULL;
	}
	~MultipoleSubsurfaceIntegrator() {
		delete [] lightSampleOffsets;
		delete [] bsdfSampleOffsets;
	}
private:
	// MultipoleSubsurfaceIntegrator Private Data
	int maxSpecularDepth;
	float maxError, minSampleDist;
	float mix;
	bool showIrradiancePoints, usePoissonPointFinder;
	string filename;
	vector<IrradiancePoint> irradiancePoints;
	BBox octreeBounds;
	SubsurfaceOctreeNode *octree;
	MemoryArena octreeArena;
	const vector<Reference<Primitive> >& originalPrimitives;

	// Declare sample parameters for light source sampling
	LightSampleOffsets *lightSampleOffsets;
	BSDFSampleOffsets *bsdfSampleOffsets;
};


MultipoleSubsurfaceIntegrator *CreateMultipoleSubsurfaceIntegrator(const ParamSet &params,
	const vector<Reference<Primitive> >* originalPrimitives);
