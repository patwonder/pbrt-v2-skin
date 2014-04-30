
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

// Calculate various absorption/scattering coefficients of human skin
// http://omlc.ogi.edu/news/jan98/skinoptics.html
class SkinCoefficients {
public:
	SkinCoefficients(float f_mel, float f_eu, float f_blood, float f_ohg,
		float ga_epi, float ga_derm, float b_derm)
	: f_mel(f_mel), f_eu(f_eu), f_blood(f_blood), f_ohg(f_ohg),
	  ga_epi(ga_epi), ga_derm(ga_derm), b_derm(b_derm) {}

	// Baseline absorption coefficient, mua.skinbaseline
	static WLDValue mua_skinbaseline() {
		return calculate([] (float wl) {
			return 0.244f + 85.3f * expf(-(wl - 154.f) / 66.2f);
		});
	}

	// == Epidermis ============================================
	// Absorption coefficient of a single eumelanosome, mua.eumel
	static WLDValue mua_eumel() {
		return calculate([] (float wl) {
			return 6.6e11f * powf(wl, -3.33f);
		});
	}
	// Absorption coefficient of a single pheomelanosome, mua.pheomel
	static WLDValue mua_pheomel() {
		return calculate([] (float wl) {
			return 2.9e15f * powf(wl, -4.75f);
		});
	}
	// Volume fraction of melanosomes in epidermis
	float f_mel;
	// Fraction of eumelanin in melanosomes
	float f_eu;
	// Net epidermal absorption coefficient, mua.epi
	WLDValue mua_epi() const {
		return f_mel * (f_eu * mua_eumel() + (1 - f_eu) * mua_pheomel()) +
			(1 - f_mel) * mua_skinbaseline();
	}
	//// Scattering coefficient of the epidermis, mus.epi
	//Spectrum mus_epi() const;
	//// Anisotropy of the epidermis, g.epi
	//float g_epi;
	// Reduced scattering coefficient of the epidermis, musp.epi
	static WLDValue musp_epi() {
		return musp_Rayleigh() + musp_Mie_fibers();
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
	static WLDValue mua_ohg() {
		const float ln10 = 2.303f;
		const float molarWeight = 64500.f; // g/mole
		const float concentration = 150.f; // g/L
		return ln10 / molarWeight * concentration *
			WLDValue::FromSampled(ohg_lambdas, ohg_vals, ohg_n);
	}
	static WLDValue mua_dhg() {
		const float ln10 = 2.303f;
		const float molarWeight = 64500.f; // g/mole
		const float concentration = 150.f; // g/L
		return ln10 / molarWeight * concentration *
			WLDValue::FromSampled(dhg_lambdas, dhg_vals, dhg_n);
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
			//return 2e5f * powf(wl, -1.5f);
			return 147.4f * powf(wl, -0.22);
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
		// scale coeff by 50% as the dermis is more translucent
		return musp_epi() * 0.5f;
	}
	// Anisotropy of the epidermis, ga.epi
	float ga_epi;
	// Anisotropy of the epidermis, ga.derm
	float ga_derm;
	// Isotropic scattering coefficient of the dermis, b.derm;
	float b_derm;
private:
	// Calculated SampledSpectrum from wavelength->value mapping
	template <class MappingFunction>
	static WLDValue calculate(const MappingFunction& mapping) {
		WLDValue res;
		for (int i = 0; i < WLD_nSamples; i++) {
			res[i] = mapping(WLD_lambdas[i]);
		}
		return res;
	}

	static const float ohg_lambdas[];
	static const float ohg_vals[];
	static const int ohg_n;
	static const float dhg_lambdas[];
	static const float dhg_vals[];
	static const int dhg_n;
};

