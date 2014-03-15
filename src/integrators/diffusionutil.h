
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

#pragma once

#include "octree.h"
#include "renderers/surfacepoints.h"
#include "irradiancepoint.h"

struct DiffusionReflectance {
    // DiffusionReflectance Public Methods
    DiffusionReflectance(const Spectrum &sigma_a, const Spectrum &sigmap_s,
                         float eta) {
        A = (1.f + Fdr(eta)) / (1.f - Fdr(eta));
        sigmap_t = sigma_a + sigmap_s;
        sigma_tr = Sqrt(3.f * sigma_a * sigmap_t);
        alphap = sigmap_s / sigmap_t;
        zpos = Spectrum(1.f) / sigmap_t;
        zneg = -zpos * (1.f + (4.f/3.f) * A);
    }
    Spectrum operator()(float d2) const {
        Spectrum dpos = Sqrt(Spectrum(d2) + zpos * zpos);
        Spectrum dneg = Sqrt(Spectrum(d2) + zneg * zneg);
        Spectrum Rd = (alphap / (4.f * M_PI)) *
            ((zpos * (dpos * sigma_tr + Spectrum(1.f)) *
              Exp(-sigma_tr * dpos)) / (dpos * dpos * dpos) -
             (zneg * (dneg * sigma_tr + Spectrum(1.f)) *
              Exp(-sigma_tr * dneg)) / (dneg * dneg * dneg));
        return Rd.Clamp();
    }

    // DiffusionReflectance Data
    Spectrum zpos, zneg, sigmap_t, sigma_tr, alphap;
    float A;
};


struct SubsurfaceOctreeNode {
    // SubsurfaceOctreeNode Methods
    SubsurfaceOctreeNode() {
        isLeaf = true;
        sumArea = 0.f;
        for (int i = 0; i < 8; ++i)
            ips[i] = NULL;
    }
    void Insert(const BBox &nodeBound, IrradiancePoint *ip,
                MemoryArena &arena) {
        Point pMid = .5f * nodeBound.pMin + .5f * nodeBound.pMax;
        if (isLeaf) {
            // Add _IrradiancePoint_ to leaf octree node
            for (int i = 0; i < 8; ++i) {
                if (!ips[i]) {
                    ips[i] = ip;
                    return;
                }
            }

            // Convert leaf node to interior node, redistribute points
            isLeaf = false;
            IrradiancePoint *localIps[8];
            for (int i = 0; i < 8; ++i) {
                localIps[i] = ips[i];
                children[i] = NULL;
            }
            for (int i = 0; i < 8; ++i)  {
                IrradiancePoint *ip = localIps[i];
                // Add _IrradiancePoint_ _ip_ to interior octree node
                int child = (ip->p.x > pMid.x ? 4 : 0) +
                    (ip->p.y > pMid.y ? 2 : 0) + (ip->p.z > pMid.z ? 1 : 0);
                if (!children[child])
                    children[child] = arena.Alloc<SubsurfaceOctreeNode>();
                BBox childBound = octreeChildBound(child, nodeBound, pMid);
                children[child]->Insert(childBound, ip, arena);
            }
            /* fall through to interior case to insert the new point... */
        }
        // Add _IrradiancePoint_ _ip_ to interior octree node
        int child = (ip->p.x > pMid.x ? 4 : 0) +
            (ip->p.y > pMid.y ? 2 : 0) + (ip->p.z > pMid.z ? 1 : 0);
        if (!children[child])
            children[child] = arena.Alloc<SubsurfaceOctreeNode>();
        BBox childBound = octreeChildBound(child, nodeBound, pMid);
        children[child]->Insert(childBound, ip, arena);
    }
    void InitHierarchy() {
        if (isLeaf) {
            // Init _SubsurfaceOctreeNode_ leaf from _IrradiancePoint_s
            float sumWt = 0.f;
            uint32_t i;
            for (i = 0; i < 8; ++i) {
                if (!ips[i]) break;
                float wt = ips[i]->E.y();
                E += ips[i]->E;
                p += wt * ips[i]->p;
                sumWt += wt;
                sumArea += ips[i]->area;
            }
            if (sumWt > 0.f) p /= sumWt;
            E /= i;
        }
        else {
            // Init interior _SubsurfaceOctreeNode_
            float sumWt = 0.f;
            uint32_t nChildren = 0;
            for (uint32_t i = 0; i < 8; ++i) {
                if (!children[i]) continue;
                ++nChildren;
                children[i]->InitHierarchy();
                float wt = children[i]->E.y();
                E += children[i]->E;
                p += wt * children[i]->p;
                sumWt += wt;
                sumArea += children[i]->sumArea;
            }
            if (sumWt > 0.f) p /= sumWt;
            E /= nChildren;
        }
    }
    Spectrum Mo(const BBox &nodeBound, const Point &p, const DiffusionReflectance &Rd,
                float maxError);

    // SubsurfaceOctreeNode Public Data
    Point p;
    bool isLeaf;
    Spectrum E;
    float sumArea;
    union {
        SubsurfaceOctreeNode *children[8];
        IrradiancePoint *ips[8];
    };
};
