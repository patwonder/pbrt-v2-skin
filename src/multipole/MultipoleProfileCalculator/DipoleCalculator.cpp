
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

#include "DipoleCalculator.h"
#include <cmath>
#include "numutil.h"

static float FresnelDiffuseReflectance(float eta) {
	if (eta >= 1.f) {
		return -1.4399f / (eta * eta) + 0.7099f / eta + 0.6681f + 0.0636f * eta;
	} else {
		float eta2 = eta * eta;
		return -0.4399f + 0.7099f / eta - 0.3319f / eta2 + 0.0636f / (eta2 * eta);
	}
}

static float A(float Fdr) {
	return (1.f + Fdr) / (1.f - Fdr);
}

DipoleCalculator::DipoleCalculator(float eta_0, float eta_d, float d,
	float sigma_a, float sigmap_s, int32 zi, bool lerpOnThinSlab)
{
	this->d = d;
	float sigmap_t = sigma_a + sigmap_s;
	sigma_tr = sqrt(3 * sigma_a * sigmap_t);
	alphap = sigmap_s / sigmap_t;
	float Fdr0 = FresnelDiffuseReflectance(eta_0);
	float Fdrd = FresnelDiffuseReflectance(eta_d);
	float A_0 = A(Fdr0), A_d = A(Fdrd);
	float D = 1.f / (3.f * sigmap_t);
	float zb_0 = 2.f * A_0 * D, zb_d = 2.f * A_d * D;
	float l = 1.f / sigmap_t;
	if (lerpOnThinSlab && l > d * .5f) {
		l = d * .5f;
	}
	zpos = 2.f * (float)zi * (d + zb_0 + zb_d) + l;
	zneg = zpos - 2.f * (l + zb_0);
}

float DipoleCalculator::Rd(float dsq) const {
	float dpos = sqrt(dsq + zpos * zpos);
	float dneg = sqrt(dsq + zneg * zneg);
	float dpos3 = dpos * dpos * dpos;
	float dneg3 = dneg * dneg * dneg;
	float rd = alphap * INV_FOURPI * (
		zpos * (1 + sigma_tr * dpos) * exp(-sigma_tr * dpos) / dpos3 -
		zneg * (1 + sigma_tr * dneg) * exp(-sigma_tr * dneg) / dneg3);
	return rd;
}

float DipoleCalculator::Td(float dsq) const {
	float dpos = sqrt(dsq + (d - zpos) * (d - zpos));
	float dneg = sqrt(dsq + (d - zneg) * (d - zneg));
	float dpos3 = dpos * dpos * dpos;
	float dneg3 = dneg * dneg * dneg;
	float td = alphap * INV_FOURPI * (
		(d - zpos) * (1 + sigma_tr * dpos) * exp(-sigma_tr * dpos) / dpos3 -
		(d - zneg) * (1 + sigma_tr * dneg) * exp(-sigma_tr * dneg) / dneg3);
	return td;
}
