
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

struct MultipoleProfileData;

class MultipoleBSSRDFData {
public:
	// MultipoleBSSRDF Public Methods
    MultipoleBSSRDFData(int layers, const Spectrum mua[], const Spectrum musp[], float et[], float thickness[],
		const MultipoleProfileData* pData, const vector<Spectrum>& rhoData, bool isMonteCarlo = false) {
		if (layers > MAX_LAYERS)
			layers = MAX_LAYERS;
		nLayers = layers;
		for (int i = 0; i < layers; i++) {
			sig_a[i] = mua[i];
			sigp_s[i] = musp[i];
			e[i] = et[i];
			d[i] = thickness[i];
		}
		this->pData = pData;
		this->isMonteCarlo = isMonteCarlo;
		this->rhoData = rhoData;
		computeRhoIntegral();
	}
	int numLayers() const { return nLayers; }
	float thickness(int layer) const { return d[layer]; }
    float eta(int layer) const { return e[layer]; }
    Spectrum sigma_a(int layer) const { return sig_a[layer]; }
    Spectrum sigma_prime_s(int layer) const { return sigp_s[layer]; }
	bool IsMonteCarlo() const { return isMonteCarlo; }
	Spectrum rho(const Vector& wo) const;
	Spectrum rho() const;

	Spectrum reflectance(float distanceSquared) const;
	Spectrum transmittance(float distanceSquared) const;
	Spectrum totalReflectance() const;
	Spectrum totalTransmittance() const;

	// MultipoleBSSRDFData Public Data
	static const int MAX_LAYERS = 4;
private:
    // MultipoleBSSRDFData Private Data
	int nLayers;
	float d[MAX_LAYERS];
    float e[MAX_LAYERS];
    Spectrum sig_a[MAX_LAYERS];
	Spectrum sigp_s[MAX_LAYERS];
	const MultipoleProfileData* pData;
	bool isMonteCarlo;

	vector<Spectrum> rhoData;
	Spectrum rhoIntegral;

	// MultipoleBSSRDF Private Methods
	void computeRhoIntegral();
};

void ComputeMultipoleProfile(int layers, const SampledSpectrum mua[], const SampledSpectrum musp[], float et[], float thickness[],
							 MultipoleProfileData** oppData, bool useMonteCarlo = false, bool lerpOnThinSlab = true);
void ComputeIrradiancePointsProfile(MultipoleProfileData** oppData, float radius);
void ReleaseMultipoleProfile(MultipoleProfileData* pData);

vector<Spectrum> ComputeRhoDataFromBxDF(const BxDF* bxdf);
vector<Spectrum> ComputeRoughRhoData();
