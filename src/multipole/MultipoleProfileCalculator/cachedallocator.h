
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

#include <unordered_map>
#include <vector>
#include <mutex>
#include <memory.h>

template <class Param, class Type>
class CPPAllocator {
public:
	Type alloc(const Param& param) {
		return ::operator new(param);
	}
	void free(const Type& p) {
		::operator delete(p);
	}
};

template <class Param, class Type, class TrueAllocator = CPPAllocator<Param, Type>, class Mutex = std::mutex>
class CachedAllocator {
public:
	~CachedAllocator() {
		clear();
	}
	Type alloc(const Param& param) {
		{
			std::lock_guard<Mutex> lock(mtx);
			auto iter = cache.find(param);
			if (iter != cache.end()) {
				auto& blocks = iter->second;
				if (blocks.size()) {
					Type block = blocks.back();
					blocks.pop_back();
					return block;
				}
			}
		}
		Type block = allocator.alloc(param);
		{
			std::lock_guard<Mutex> lock(mtx);
			paramMap.insert(make_pair(block, param));
		}
		return block;
	}
	void free(const Type& p) {
		std::lock_guard<Mutex> lock(mtx);
		auto iter = paramMap.find(p);
		if (iter != paramMap.end()) {
			cache[iter->second].push_back(p);
		} else {
			allocator.free(p);
		}
	}
	void safeClear() {
		std::lock_guard<Mutex> lock(mtx);
		clear();
	}
	void clear() {
		for (auto& pair : cache)
			for (auto& block : pair.second)
				allocator.free(block);
		cache.clear();
		paramMap.clear();
	}
private:
	std::unordered_map<Param, std::vector<Type> > cache;
	std::unordered_map<Type, Param> paramMap;
	Mutex mtx;
	static TrueAllocator allocator;
};

template <class Param, class Type, class TrueAllocator, class Mutex>
TrueAllocator CachedAllocator<Param, Type, TrueAllocator, Mutex>::allocator;

template <class Param, class Type, class TrueAllocator, class Mutex>
void* operator new(size_t size, CachedAllocator<Param, Type, TrueAllocator, Mutex>& alloc) {
	return alloc.alloc(size);
}

template <class Param, class Type, class TrueAllocator, class Mutex>
void operator delete(void* p, CachedAllocator<Param, Type, TrueAllocator, Mutex>& alloc) {
	alloc.free(p);
}
