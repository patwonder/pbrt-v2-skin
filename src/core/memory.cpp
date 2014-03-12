
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


// core/memory.cpp*
#include "stdafx.h"
#include "memory.h"

// Memory Allocation Functions
void *AllocAligned(size_t size) {
#if defined(PBRT_IS_WINDOWS)
    return _aligned_malloc(size, PBRT_L1_CACHE_LINE_SIZE);
#elif defined (PBRT_IS_OPENBSD) || defined(PBRT_IS_APPLE)
    // Allocate excess memory to ensure an aligned pointer can be returned
    void *mem = malloc(size + (PBRT_L1_CACHE_LINE_SIZE-1) + sizeof(void*));
    char *amem = ((char*)mem) + sizeof(void*);
#if (PBRT_POINTER_SIZE == 8)
    amem += PBRT_L1_CACHE_LINE_SIZE - (reinterpret_cast<uint64_t>(amem) &
                                       (PBRT_L1_CACHE_LINE_SIZE - 1));
#else
    amem += PBRT_L1_CACHE_LINE_SIZE - (reinterpret_cast<uint32_t>(amem) &
                                       (PBRT_L1_CACHE_LINE_SIZE - 1));
#endif
    ((void**)amem)[-1] = mem;
    return amem;
#else
    return memalign(PBRT_L1_CACHE_LINE_SIZE, size);
#endif
}


void FreeAligned(void *ptr) {
    if (!ptr) return;
#if defined(PBRT_IS_WINDOWS)
    _aligned_free(ptr);
#elif defined (PBRT_IS_OPENBSD) || defined(PBRT_IS_APPLE)
    free(((void**)ptr)[-1]);
#else
    free(ptr);
#endif
}


#include <unordered_map>


class SelfishAllocator {
public:
	SelfishAllocator(uint32_t sz)
		: blockSize(sz) { }
	void* AllocBlock(uint32_t minSize = 0) {
		void* block;
		if (bucket.size() && minSize <= blockSize) {
			block = bucket.back();
			bucket.pop_back();
		} else {
			block = AllocAligned(max(minSize, blockSize));
		}
		return block;
	}
	void FreeBlock(void* block) {
		// Selfish Allocator never returns memory to OS
		// It eats the block itself
		bucket.push_back(block);
	}
	~SelfishAllocator() {
		// Actually, not that selfish - -
		for (void* block : bucket)
			FreeAligned(block);
	}
private:
	uint32_t blockSize;
	vector<void*> bucket;
};


// Good name!
class SelfishAllocatorAllocator {
public:
	SelfishAllocator& getAllocator(uint32_t blockSize) {
		SelfishAllocator* pAllocator;
		auto& iter = mapSelfishAllocators.find(blockSize);
		if (iter == mapSelfishAllocators.end()) {
			pAllocator = new SelfishAllocator(blockSize);
			mapSelfishAllocators[blockSize] = pAllocator;
		} else {
			pAllocator = iter->second;
		}
		return *pAllocator;
	}
	~SelfishAllocatorAllocator() {
		for (auto& pair : mapSelfishAllocators)
			delete pair.second;
	}
private:
	std::unordered_map<uint32_t, SelfishAllocator*> mapSelfishAllocators;
};


__declspec(thread) static SelfishAllocatorAllocator* pAA = NULL;


void InitTLSAllocator() {
	pAA = new SelfishAllocatorAllocator();
}


void CleanupTLSAllocator() {
	delete pAA;
	pAA = NULL;
}


static SelfishAllocator& GetAllocatorForBlockSize(uint32_t blockSize) {
	if (!pAA) // Main thread
		InitTLSAllocator();
	return pAA->getAllocator(blockSize);
}


MemoryArena::MemoryArena(uint32_t bs)
	: allocator(GetAllocatorForBlockSize(bs))
{
	blockSize = bs;
	curBlockPos = 0;
	currentBlock = (char*)allocator.AllocBlock();
}


MemoryArena::~MemoryArena() {
	allocator.FreeBlock(currentBlock);
	for (uint32_t i = 0; i < usedBlocks.size(); ++i)
		allocator.FreeBlock(usedBlocks[i]);
	for (uint32_t i = 0; i < availableBlocks.size(); ++i)
		allocator.FreeBlock(availableBlocks[i]);
}


void* MemoryArena::Alloc(uint32_t sz) {
	// Round up _sz_ to minimum machine alignment
	sz = ((sz + 15) & (~15));
	if (curBlockPos + sz > blockSize) {
		// Get new block of memory for _MemoryArena_
		usedBlocks.push_back(currentBlock);
		if (availableBlocks.size() && sz <= blockSize) {
			currentBlock = availableBlocks.back();
			availableBlocks.pop_back();
		}
		else
			currentBlock = (char*)allocator.AllocBlock(sz);
		curBlockPos = 0;
	}
	void *ret = currentBlock + curBlockPos;
	curBlockPos += sz;
	return ret;
}


void MemoryArena::FreeAll() {
	curBlockPos = 0;
	while (usedBlocks.size()) {
#ifndef NDEBUG
		memset(usedBlocks.back(), 0xfa, blockSize);
#endif
		availableBlocks.push_back(usedBlocks.back());
		usedBlocks.pop_back();
	}
}

