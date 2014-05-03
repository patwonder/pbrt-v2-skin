
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
#include "multipole.h"
#include "skincoeffs.h"

LayeredSkin::LayeredSkin(const vector<SkinLayer>& layers, float r, float npu,
	const SkinCoefficients& coeff, Reference<Texture<Spectrum> > Kr, Reference<Texture<Spectrum> > Kt,
	Reference<Texture<float> > bumpMap, Reference<Texture<Spectrum> > albedo, bool doubleRefSSLF,
	bool generateProfile, bool useMonteCarloProfile, uint64_t nPhotons, bool lerpOnThinSlab,
	bool showIrradiancePoints, float irradiancePointSize)
	: layers(layers), roughness(r), nmperunit(npu), pcoeff(new SkinCoefficients(coeff)), Kr(Kr), Kt(Kt),
	  bumpMap(bumpMap), albedo(albedo), doubleRefSSLF(doubleRefSSLF)
{
	const float NM_PER_CM = 1e7f;
	// Calculate layer params
	// Epidermis
	lps[0].thickness = layers[0].thickness / nmperunit;
	lps[0].mua = pcoeff->mua_epi() * (nmperunit / NM_PER_CM);
	lps[0].musp = pcoeff->musp_epi() * (nmperunit / NM_PER_CM);
	lps[0].ga = pcoeff->ga_epi;
	lps[0].isotropicHGPF = false;
	// Dermis
	lps[1].thickness = layers[1].thickness / nmperunit;
	lps[1].mua = pcoeff->mua_derm() * (nmperunit / NM_PER_CM);
	lps[1].musp = pcoeff->musp_derm() * (nmperunit / NM_PER_CM);
	lps[1].ga = pcoeff->ga_derm;
	lps[1].b = pcoeff->b_derm;
	lps[1].isotropicHGPF = true;
	// The thing below...
	lps[2].thickness = 2e6f / nmperunit; // 2mm
	lps[2].mua = lps[1].mua;
	lps[2].musp = lps[1].musp;
	lps[2].ga = 1.f;
	lps[2].isotropicHGPF = false;

	// Prepare BSSRDF for better performance
	if (generateProfile) {
		Spectrum mua[2];
		Spectrum musp[2];
		SampledSpectrum smua[2];
		SampledSpectrum smusp[2];

		float eta[2];
		float thickness[2];
		for (int i = 0; i < 2; i++) {
			mua[i] = lps[i].mua.toSpectrum();
			musp[i] = lps[i].musp.toSpectrum();
			smua[i] = lps[i].mua.toSampledSpectrum();
			smusp[i] = lps[i].musp.toSampledSpectrum();
			eta[i] = layers[i].ior;
			thickness[i] = lps[i].thickness;
		}

		RhoData rhoData;

		if (showIrradiancePoints) {
			ComputeIrradiancePointsProfile(&profileData, irradiancePointSize);
			rhoData = ComputeRoughRhoData();
		} else {
			ComputeMultipoleProfile(2, smua, smusp, eta, thickness, &profileData, useMonteCarloProfile, lerpOnThinSlab, nPhotons);
			MemoryArena arena;
			Fresnel *fresnel = doubleRefSSLF ? BSDF_ALLOC(arena, FixedFresnelDielectric)(1.f, layers[0].ior)
											 : BSDF_ALLOC(arena, FresnelDielectric)(1.f, layers[0].ior);
			BxDF* bxdf = BSDF_ALLOC(arena, Microfacet)(Spectrum(1.f), fresnel,
				BSDF_ALLOC(arena, Beckmann)(roughness));
			rhoData = ComputeRhoDataFromBxDF(bxdf);
			arena.FreeAll();
		}

		preparedBSSRDFData = new MultipoleBSSRDFData(2, mua, musp, eta, thickness, profileData, rhoData, useMonteCarloProfile);

		//for (int i = 0; i < 11; i++) {
		//	float cost = (float)i / 10;
		//	float val = preparedBSSRDFData->rho(SphericalDirection<float>(sqrt(1 - cost * cost), cost, 0.f))[0];
		//	Info("Rho[%f]=%f", cost, val);
		//}
	} else {
		preparedBSSRDFData = NULL;
	}
}

LayeredSkin::~LayeredSkin() {
	delete pcoeff;
	ReleaseMultipoleProfile(profileData);
	delete preparedBSSRDFData;
}

vector<float> LayeredSkin::GetLayerThickness() const {
	vector<float> ret;
	for (const auto& layer : layers) {
		ret.push_back(layer.thickness / nmperunit);
	}
	return ret;
}

BSDF* LayeredSkin::GetBSDF(const DifferentialGeometry &dgGeom,
	const DifferentialGeometry &dgShading,
	MemoryArena &arena) const
{
	DifferentialGeometry dgs;
	if (bumpMap)
		Bump(bumpMap, dgGeom, dgShading, &dgs);
	else
		dgs = dgShading;

	float ior = layers[0].ior;
	BSDF* bsdf = BSDF_ALLOC(arena, BSDF)(dgs, dgGeom.nn, ior);
	float rough = roughness;
	Fresnel *fresnel = doubleRefSSLF ? BSDF_ALLOC(arena, FixedFresnelDielectric)(1.f, ior)
									 : BSDF_ALLOC(arena, FresnelDielectric)(1.f, ior);
	Spectrum R = Kr->Evaluate(dgs);
	Spectrum T = Kt->Evaluate(dgs);
	if (!R.IsBlack())
		bsdf->Add(BSDF_ALLOC(arena, Microfacet)(R, fresnel,
			BSDF_ALLOC(arena, Beckmann)(rough)));
	if (!T.IsBlack())
		bsdf->Add(BSDF_ALLOC(arena, MicrofacetTransmission)(T, fresnel,
			BSDF_ALLOC(arena, Beckmann)(rough), ior));
	//bsdf->Add(BSDF_ALLOC(arena, BRDFToBTDF)(BSDF_ALLOC(arena, Microfacet)(Spectrum(1.), fresnel,
	//	BSDF_ALLOC(arena, Blinn)(1.f / rough))));
	//bsdf->Add(BSDF_ALLOC(arena, BRDFToBTDF)(BSDF_ALLOC(arena, MicrofacetTransmission)(Spectrum(1.), fresnel,
	//	BSDF_ALLOC(arena, Blinn)(1.f / rough), ior)));
	return bsdf;
}


BSSRDF* LayeredSkin::GetBSSRDF(const DifferentialGeometry &dgGeom,
	const DifferentialGeometry &dgShading, MemoryArena &arena) const
{
	// Trick the tracing SurfacePointsRenderer to think I'm translucent!
	const LayerParam& lp = lps[1];
	return BSDF_ALLOC(arena, BSSRDF)(lp.mua.toSpectrum(),
		lp.musp.toSpectrum(), layers[1].ior);
}


const MultipoleBSSRDF* LayeredSkin::GetMultipoleBSSRDF(const DifferentialGeometry &dgGeom,
	const DifferentialGeometry &dgShading, MemoryArena &arena) const
{
	Spectrum al = albedo->Evaluate(dgShading);
	return BSDF_ALLOC(arena, MultipoleBSSRDF)(preparedBSSRDFData, al);
}


BSDF* LayeredSkin::GetLayeredBSDF(int layerIndex,
	const DifferentialGeometry &dgGeom,
	const DifferentialGeometry &dgShading,
	MemoryArena &arena) const
{
	if (layerIndex == 0)
		return GetBSDF(dgGeom, dgShading, arena);

	float ior;
	if ((size_t)layerIndex == layers.size())
		ior = 1.f;
	else
		ior = layers[layerIndex].ior / layers[layerIndex - 1].ior;
	BSDF* bsdf = BSDF_ALLOC(arena, BSDF)(dgShading, dgGeom.nn, ior);
	Fresnel *fresnel = BSDF_ALLOC(arena, FresnelDielectric)(1.f, ior);
	bsdf->Add(BSDF_ALLOC(arena, SpecularReflection)(Spectrum(1.f), fresnel));
	bsdf->Add(BSDF_ALLOC(arena, SpecularTransmission)(Spectrum(1.f), 1.f, ior));
	return bsdf;
}

BSSRDF* LayeredSkin::GetLayeredBSSRDF(int layerIndex,
	const DifferentialGeometry &dgGeom,
	const DifferentialGeometry &dgShading,
	MemoryArena &arena) const
{
	return NULL;
}

LayerParam LayeredSkin::GetLayerParam(int index) const {
	if (index > 2)
		index = 2;
	return lps[index];
}

LayeredSkin* CreateLayeredSkinMaterial(const ParamSet& ps, const TextureParams& mp)
{
	int nLayers;
	const SkinLayer* layers = ps.FindSkinLayer("layers", &nLayers);
	if (!layers) {
		Error("No layers param set for LayeredSkin material.");
		abort();
	}
	if (nLayers < 2) {
		Error("Not enough layers specified for LayeredSkin material.");
		abort();
	}
	float roughness = ps.FindOneFloat("roughness", 0.4f);
	float nmperunit = ps.FindOneFloat("nmperunit", 100e6f);
	float f_mel = ps.FindOneFloat("f_mel", 0.15f);
	float f_eu = ps.FindOneFloat("f_eu", 1.f);
	float f_blood = ps.FindOneFloat("f_blood", 0.002f);
	float f_ohg = ps.FindOneFloat("f_ohg", 0.3f);
	float ga_epi = ps.FindOneFloat("ga_epi", 0.9f);
	float ga_derm = ps.FindOneFloat("ga_derm", 0.8f);
	float b_derm = ps.FindOneFloat("b_derm", 0.4f);
	SkinCoefficients coeff(f_mel, f_eu, f_blood, f_ohg, ga_epi, ga_derm, b_derm);
    Reference<Texture<Spectrum> > Kr = mp.GetSpectrumTexture("Kr", Spectrum(1.f));
    Reference<Texture<Spectrum> > Kt = mp.GetSpectrumTexture("Kt", Spectrum(1.f));
    Reference<Texture<float> > bumpMap = mp.GetFloatTextureOrNull("bumpmap");
	Reference<Texture<Spectrum> > albedo = mp.GetSpectrumTexture("albedo", Spectrum(1.f));
	bool doubleRefSSLF = ps.FindOneBool("doublerefsslf", false);
	bool generateProfile = ps.FindOneBool("genprofile", true);
	bool useMonteCarloProfile = ps.FindOneBool("usemontecarlo", false);
	string strPhotons = ps.FindOneString("photons", "10000000");
	uint64_t photons = _strtoui64(strPhotons.c_str(), NULL, 10);
	bool lerpOnThinSlab = ps.FindOneBool("lerponthinslab", true);
	bool showIrradiancePoints = ps.FindOneBool("showirradiancepoints", false);
	float irradiancePointSize = ps.FindOneFloat("irradiancepointsize", 0.002f);
	return new LayeredSkin(vector<SkinLayer>(layers, layers + nLayers),
		roughness, nmperunit, coeff, Kr, Kt, bumpMap, albedo, doubleRefSSLF,
		generateProfile, useMonteCarloProfile, photons, lerpOnThinSlab,
		showIrradiancePoints, irradiancePointSize);
}
