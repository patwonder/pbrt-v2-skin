
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

#include "mpc-types.h"

typedef struct MPC_LayerSpec {
	// Index of refraction
	float ior;
	// Layer thickness
	float thickness;
	// Absorption coefficient
	float mua;
	// Reduced scattering coefficient
	float musp;
	// Anisotropy for the HGPF (probably not needed)
	float g_HG;

} MPC_LayerSpec;

typedef struct MPC_Options {
	// distance between adjacent samples
	float desiredStepSize;
	// number of samples in one direction
	uint32 desiredLength;
} MPC_Options;

typedef struct MPC_Output {
	// number of samples
	uint32 length;
	// Profile data
	float* pDistanceSquared;
	float* pReflectance;
	float* pTransmittance;
} MPC_Output;

// Profile calculation
MULTIPOLEPROFILECALCULATOR_API void MPC_ComputeDiffusionProfile(uint32 numLayers, const MPC_LayerSpec* pLayerSpecs,
	const MPC_Options* pOptions, MPC_Output** oppOutput);
MULTIPOLEPROFILECALCULATOR_API void MPC_ResampleForUniformDistanceSquaredDistribution(MPC_Output* pOutput);
MULTIPOLEPROFILECALCULATOR_API void MPC_FreeOutput(MPC_Output* pOutput);
