
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
#include <functional>


// Calculate various absorption/scattering coefficients of human skin
// http://omlc.ogi.edu/news/jan98/skinoptics.html
class SkinCoefficients {
public:
	static const float ohg_lambdas[];
	static const float ohg_vals[];
	static const int ohg_n;
	static const float dhg_lambdas[];
	static const float dhg_vals[];
	static const int dhg_n;

	SkinCoefficients(float f_mel, float f_blood, float f_ohg,
		float ga_epi, float ga_derm, float b_derm)
	: f_mel(f_mel), f_blood(f_blood), f_ohg(f_ohg),
	  ga_epi(ga_epi), ga_derm(ga_derm), b_derm(b_derm) {}

	// Baseline absorption coefficient, mua.skinbaseline
	static WLDValue mua_skinbaseline() {
		return calculate([] (float wl) {
			return 0.244f + 85.3f * expf(-(wl - 154.f) / 66.2f);
		});
	}

	// == Epidermis ============================================
	// Absorption coefficient of a single melanosome, mua.mel
	static WLDValue mua_mel() {
		return calculate([] (float wl) {
			return 6.6e11f * powf(wl, -3.33f);
		});
	}
	// Volume fraction of melanosomes in epidermis
	float f_mel;
	// Net epidermal absorption coefficient, mua.epi
	WLDValue mua_epi() const {
		return f_mel * mua_mel() + (1 - f_mel) * mua_skinbaseline();
	}
	//// Scattering coefficient of the epidermis, mus.epi
	//Spectrum mus_epi() const;
	//// Anisotropy of the epidermis, g.epi
	//float g_epi;
	// Reduced scattering coefficient of the epidermis, musp.epi
	static WLDValue musp_epi() {
		return musp_derm();
	}

	// == Dermis ============================================
	// Volume fraction of oxyhemoglobin, f_ohg
	float f_ohg;
	// Absorption coefficient of whole blood, mua.blood
	// Data from: http://www.npsg.uwaterloo.ca/data/blood.php
	WLDValue mua_blood() const {
		const float ln10 = 2.303f;
		const float molarWeight = 64500.f; // g/mole
		const float concentration = 150.f; // g/L
		return ln10 / molarWeight * concentration *
			(f_ohg * WLDValue::FromSampled(
			ohg_lambdas, ohg_vals, ohg_n) +
			(1.f - f_ohg) * WLDValue::FromSampled(
			dhg_lambdas, dhg_vals, dhg_n));
	}
	// Average volume fraction of blood, f.blood
	float f_blood;
	// Absorption coefficient of dermis perfused with blood, mua.derm
	WLDValue mua_derm() const {
		return f_blood * mua_blood() + (1 - f_blood) * mua_skinbaseline();
	}
	// (Reduced) Mie scattering coefficient of collagen fibers, musp_Mie.fibers
	static WLDValue musp_Mie_fibers() {
		return calculate([] (float wl) {
			return 2e5f * powf(wl, -1.5f);
		});
	}
	// (Reduced) Rayleigh scattering coefficient of the dermis, musp_Rayleigh
	static WLDValue musp_Rayleigh() {
		return calculate([] (float wl) {
			return 2e12f * pow(wl, -4.f);
		});
	}
	// (Reduced) scattering coefficient of dermis, musp.derm
	static WLDValue musp_derm() {
		return musp_Rayleigh() + musp_Mie_fibers();
	}
	// Anisotropy of the epidermis, ga.epi
	float ga_epi;
	// Anisotropy of the epidermis, ga.derm
	float ga_derm;
	// Isotropic scattering coefficient of the dermis, b.derm;
	float b_derm;
private:
	// Calculated SampledSpectrum from wavelength->value mapping
	static WLDValue calculate(std::function<float(float wl)> mapping) {
		WLDValue res;
		for (int i = 0; i < WLD_nSamples; i++) {
			res[i] = mapping(WLD_lambdas[i]);
		}
		return res;
	}
};


const float SkinCoefficients::ohg_lambdas[] = {
	400, 405, 410, 415, 420, 425, 430, 435, 440, 445,
	450, 455, 460, 465, 470, 475, 480, 485, 490, 495,
	500, 505, 510, 515, 520, 525, 530, 535, 540, 545,
	550, 555, 560, 565, 570, 575, 580, 585, 590, 595,
	600, 605, 610, 615, 620, 625, 630, 635, 640, 645,
	650, 655, 660, 665, 670, 675, 680, 685, 690, 695,
	700
};
const float SkinCoefficients::ohg_vals[] = {
	266200, 331450, 466800, 523100, 480400, 351100, 246100, 149050,
	102600, 78880, 62820, 51525, 44480, 38440, 33210, 29480, 26630,
	24925, 23680, 22155, 20930, 20185, 20040, 20715, 24200, 30885,
	39960, 48335, 53240, 50985, 43020, 35650, 32610, 35210, 44500,
	54425, 50100, 30620, 14400, 6681.5, 3200, 1958.5, 1506, 1166.5,
	942, 740.8, 610, 495.6, 442, 397.7, 368, 340.3, 319.6, 305.6,
	294, 283.8, 277.6, 273.6, 276, 280.6, 290
};
const int SkinCoefficients::ohg_n = ARRAYSIZE(ohg_lambdas);

const float SkinCoefficients::dhg_lambdas[] = {
	400, 405, 410, 415, 420, 425, 430, 435, 440, 445,
	450, 455, 460, 465, 470, 475, 480, 485, 490, 495,
	500, 505, 510, 515, 520, 525, 530, 535, 540, 545,
	550, 555, 560, 565, 570, 575, 580, 585, 590, 595,
	600, 605, 610, 615, 620, 625, 630, 635, 640, 645,
	650, 655, 660, 665, 670, 675, 680, 685, 690, 695,
	700
};
const float SkinCoefficients::dhg_vals[] = {
	223300, 261950, 304000, 353200, 407600, 471500, 528600, 549600,
	413300, 259950, 103300, 33435, 23390, 18700, 16160, 14920, 14550,
	15375, 16680, 18650, 20860, 23285, 25770, 28680, 31590, 35170,
	39040, 42840, 46590, 50490, 53410, 54530, 53790, 49700, 45070,
	40905, 37020, 33590, 28320, 21185, 14680, 12040, 9444, 7553.5,
	6510, 5763.5, 5149, 4666.5, 4345, 4026.5, 3750, 3481.5, 3227,
	3011, 2795, 2591, 2408, 2224.5, 2052, 1923.5, 1794
};
const int SkinCoefficients::dhg_n = ARRAYSIZE(dhg_lambdas);


LayeredSkin::LayeredSkin(const vector<SkinLayer>& layers, float r, float npu,
	const SkinCoefficients& coeff, Reference<Texture<Spectrum> > Kr, Reference<Texture<Spectrum> > Kt,
	Reference<Texture<float> > bumpMap)
	: layers(layers), roughness(r), nmperunit(npu), pcoeff(new SkinCoefficients(coeff)), Kr(Kr), Kt(Kt),
	bumpMap(bumpMap)
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
	lps[2].mua = lps[1].mua * 3.f;
	lps[2].musp = lps[1].musp * 0.333f;
	lps[2].ga = 1.f;
	lps[2].isotropicHGPF = false;
}

LayeredSkin::~LayeredSkin() {
	delete pcoeff;
}

vector<float_type> LayeredSkin::GetLayerThickness() const {
	vector<float_type> ret;
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
	Fresnel *fresnel = BSDF_ALLOC(arena, FresnelDielectric)(1.f, ior);
	Spectrum R = Kr->Evaluate(dgs);
	Spectrum T = Kt->Evaluate(dgs);
	if (!R.IsBlack())
		bsdf->Add(BSDF_ALLOC(arena, Microfacet)(R, fresnel,
			BSDF_ALLOC(arena, Blinn)(1.f / rough)));
	if (!T.IsBlack())
		bsdf->Add(BSDF_ALLOC(arena, MicrofacetTransmission)(T, fresnel,
			BSDF_ALLOC(arena, Blinn)(1.f / rough), ior));
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


MultipoleBSSRDF* LayeredSkin::GetMultipoleBSSRDF(const DifferentialGeometry &dgGeom,
	const DifferentialGeometry &dgShading, MemoryArena &arena) const
{
	const LayerParam& lp = lps[1];
	return BSDF_ALLOC(arena, MultipoleBSSRDF)(lp.mua.toSpectrum(),
		lp.musp.toSpectrum(), layers[1].ior);
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
		ior = 1.5f;
	else
		ior = layers[layerIndex].ior / layers[layerIndex - 1].ior;
	BSDF* bsdf = BSDF_ALLOC(arena, BSDF)(dgShading, dgGeom.nn, ior);
	Fresnel *fresnel = BSDF_ALLOC(arena, FresnelDielectric)(1.f, ior);
	bsdf->Add(BSDF_ALLOC(arena, SpecularReflection)(Spectrum(1.f), fresnel));
	if ((size_t)layerIndex != layers.size())
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
	if (!layers)
		Error("No layers param set for LayeredSkin material.");
	float roughness = ps.FindOneFloat("roughness", 0.4f);
	float nmperunit = ps.FindOneFloat("nmperunit", 100e6f);
	float f_mel = ps.FindOneFloat("f_mel", 0.15f);
	float f_blood = ps.FindOneFloat("f_blood", 0.002f);
	float f_ohg = ps.FindOneFloat("f_ohg", 0.3f);
	float ga_epi = ps.FindOneFloat("ga_epi", 0.9f);
	float ga_derm = ps.FindOneFloat("ga_derm", 0.8f);
	float b_derm = ps.FindOneFloat("b_derm", 0.4f);
	SkinCoefficients coeff(f_mel, f_blood, f_ohg, ga_epi, ga_derm, b_derm);
    Reference<Texture<Spectrum> > Kr = mp.GetSpectrumTexture("Kr", Spectrum(1.f));
    Reference<Texture<Spectrum> > Kt = mp.GetSpectrumTexture("Kt", Spectrum(1.f));
    Reference<Texture<float> > bumpMap = mp.GetFloatTextureOrNull("bumpmap");
	return new LayeredSkin(vector<SkinLayer>(layers, layers + nLayers),
		roughness, nmperunit, coeff, Kr, Kt, bumpMap);
}
