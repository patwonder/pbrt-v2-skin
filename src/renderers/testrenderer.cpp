
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

#include "testrenderer.h"
#include "materials/skincoeffs.h"
#include <fstream>
#include "multipole.h"

using namespace std;

void TestRenderer::Render(const Scene* scene) {
	// Do whatever you want here!
#if 0
	WLDValue mua_ohg = SkinCoefficients::mua_ohg();
	WLDValue mua_dhg = SkinCoefficients::mua_dhg();
	ofstream out("hg.txt", ios::trunc);

	for (int i = 0; i < WLD_nSamples; i++) {
		float wl = WLD_lambdas[i];
		out << wl << "\t" << mua_ohg[i] << "\t" << mua_dhg[i] << endl;
	}

	out.close();
#elif 0
	WLDValue musp_mie = SkinCoefficients::musp_Mie_fibers() / 2;
	WLDValue musp_rayleigh = SkinCoefficients::musp_Rayleigh() / 2;

	ofstream out("musp.txt", ios::trunc);

	for (int i = 0; i < WLD_nSamples; i++) {
		float wl = WLD_lambdas[i];
		out << wl << "\t" << musp_mie[i] << "\t" << musp_rayleigh[i] << endl;
	}

	out.close();
#else
	MemoryArena arena;
	Fresnel *fresnel = BSDF_ALLOC(arena, FixedFresnelDielectric)(1.f, 1.4f);
	BxDF* bxdf = BSDF_ALLOC(arena, Microfacet)(Spectrum(1.f), fresnel,
		BSDF_ALLOC(arena, Beckmann)(0.3));
	auto rhoData = ComputeRhoDataFromBxDF(bxdf);
	arena.FreeAll();

	ofstream out("rho.txt", ios::trunc);

	for (int i = 0; i < rhoData.hd.size(); i++) {
		float rho = rhoData.hd[i].y();
		float costh = (float)i / (rhoData.hd.size() - 1);
		out << costh << "\t" << rho << endl;
	}

	out.close();
#endif
}


Spectrum TestRenderer::Li(const Scene* scene, const RayDifferential& ray,
		const Sample* sample, RNG& rng, MemoryArena& arena,
		Intersection* isect, Spectrum* T) const
{
	return Spectrum(0.f);
}


Spectrum TestRenderer::Transmittance(const Scene* scene, const RayDifferential& ray,
		const Sample* sample, RNG& rng, MemoryArena& arena) const
{
	return Spectrum(0.f);
}

TestRenderer* CreateTestRenderer(const ParamSet& params) {
	return new TestRenderer();
}
