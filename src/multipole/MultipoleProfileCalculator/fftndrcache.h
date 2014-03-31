
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

#include "tools/kiss_fftndr.h"
#include "cachedallocator.h"

struct NDRParams {
	int dims[2];
	bool inverse;
	NDRParams(int dim0, int dim1, bool inv) {
		dims[0] = dim0;
		dims[1] = dim1;
		inverse = inv;
	}
	bool operator==(const NDRParams& p) const {
		return dims[0] == p.dims[0] && dims[1] == p.dims[1] && inverse == p.inverse;
	}
};

namespace std {

	template<>
	class hash<NDRParams> {
	public:
		size_t operator()(const NDRParams& key) const {
			return inthasher(key.dims[0]) + inthasher(key.dims[1]) + boolhasher(key.inverse);
		}
	private:
		static hash<int> inthasher;
		static hash<bool> boolhasher;
	};

}

class NDRAllocator {
public:
	kiss_fftndr_cfg alloc(const NDRParams& params) {
		return kiss_fftndr_alloc(params.dims, 2, params.inverse, NULL, NULL);
	}
	void free(kiss_fftndr_cfg cfg) {
		::kiss_fft_free(cfg);
	}
};

typedef CachedAllocator<NDRParams, kiss_fftndr_cfg, NDRAllocator> NDRCache;
