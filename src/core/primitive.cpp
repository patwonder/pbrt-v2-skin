
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


// core/primitive.cpp*
#include "stdafx.h"
#include "primitive.h"
#include "light.h"
#include "intersection.h"
#include "accelerators/bvh.h"
#include "paramset.h"

// Primitive Method Definitions
uint32_t Primitive::nextprimitiveId = 1;
Primitive::~Primitive() { }

bool Primitive::CanIntersect() const {
    return true;
}



void Primitive::Refine(vector<Reference<Primitive> > &refined) const {
    Severe("Unimplemented Primitive::Refine() method called!");
}


void
Primitive::FullyRefine(vector<Reference<Primitive> > &refined) const {
    vector<Reference<Primitive> > todo;
    todo.push_back(const_cast<Primitive *>(this));
    while (todo.size()) {
        // Refine last primitive in todo list
        Reference<Primitive> prim = todo.back();
        todo.pop_back();
        if (prim->CanIntersect())
            refined.push_back(prim);
        else
            prim->Refine(todo);
    }
}


const AreaLight *Aggregate::GetAreaLight() const {
    Severe("Aggregate::GetAreaLight() method"
         "called; should have gone to GeometricPrimitive");
    return NULL;
}


BSDF *Aggregate::GetBSDF(const DifferentialGeometry &,
        const Transform &, MemoryArena &) const {
    Severe("Aggregate::GetBSDF() method"
        "called; should have gone to GeometricPrimitive");
    return NULL;
}


BSSRDF *Aggregate::GetBSSRDF(const DifferentialGeometry &,
        const Transform &, MemoryArena &) const {
    Severe("Aggregate::GetBSSRDF() method"
        "called; should have gone to GeometricPrimitive");
    return NULL;
}


const MultipoleBSSRDF* Aggregate::GetMultipoleBSSRDF(const DifferentialGeometry &dg,
	const Transform &ObjectToWorld, MemoryArena &arena) const
{
    Severe("Aggregate::GetMultipoleBSSRDF() method"
        "called; should have gone to GeometricPrimitive");
    return NULL;
}


// TransformedPrimitive Method Definitions
bool TransformedPrimitive::Intersect(const Ray &r,
                                     Intersection *isect) const {
    Transform w2p;
    WorldToPrimitive.Interpolate(r.time, &w2p);
    Ray ray = w2p(r);
    if (!primitive->Intersect(ray, isect))
        return false;
    r.maxt = ray.maxt;
    isect->primitiveId = primitiveId;
    if (!w2p.IsIdentity()) {
        // Compute world-to-object transformation for instance
        isect->WorldToObject = isect->WorldToObject * w2p;
        isect->ObjectToWorld = Inverse(isect->WorldToObject);

        // Transform instance's differential geometry to world space
        Transform PrimitiveToWorld = Inverse(w2p);
        isect->dg.p = PrimitiveToWorld(isect->dg.p);
        isect->dg.nn = Normalize(PrimitiveToWorld(isect->dg.nn));
        isect->dg.dpdu = PrimitiveToWorld(isect->dg.dpdu);
        isect->dg.dpdv = PrimitiveToWorld(isect->dg.dpdv);
        isect->dg.dndu = PrimitiveToWorld(isect->dg.dndu);
        isect->dg.dndv = PrimitiveToWorld(isect->dg.dndv);
    }
    return true;
}


bool TransformedPrimitive::IntersectP(const Ray &r) const {
    return primitive->IntersectP(WorldToPrimitive(r));
}


// GeometricPrimitive Method Definitions
BBox GeometricPrimitive::WorldBound() const {
    return shape->WorldBound();
}


bool GeometricPrimitive::IntersectP(const Ray &r) const {
    return shape->IntersectP(r);
}


bool GeometricPrimitive::CanIntersect() const {
    return shape->CanIntersect();
}


void GeometricPrimitive::
        Refine(vector<Reference<Primitive> > &refined)
        const {
    vector<Reference<Shape> > r;
    shape->Refine(r);
    for (uint32_t i = 0; i < r.size(); ++i) {
        GeometricPrimitive *gp = new GeometricPrimitive(r[i],
               material, areaLight);
        refined.push_back(gp);
    }
}


GeometricPrimitive::GeometricPrimitive(const Reference<Shape> &s,
        const Reference<Material> &m, const AreaLight *a)
    : shape(s), material(m), areaLight(a) {
}


bool GeometricPrimitive::Intersect(const Ray &r,
                                   Intersection *isect) const {
    float thit, rayEpsilon;
    if (!shape->Intersect(r, &thit, &rayEpsilon, &isect->dg))
        return false;
    isect->primitive = this;
    isect->WorldToObject = *shape->WorldToObject;
    isect->ObjectToWorld = *shape->ObjectToWorld;
    isect->shapeId = shape->shapeId;
    isect->primitiveId = primitiveId;
    isect->rayEpsilon = rayEpsilon;
    r.maxt = thit;
    return true;
}


const AreaLight *GeometricPrimitive::GetAreaLight() const {
    return areaLight;
}


BSDF *GeometricPrimitive::GetBSDF(const DifferentialGeometry &dg,
                                  const Transform &ObjectToWorld,
                                  MemoryArena &arena) const {
    DifferentialGeometry dgs;
    shape->GetShadingGeometry(ObjectToWorld, dg, &dgs);
    return material->GetBSDF(dg, dgs, arena);
}


BSSRDF *GeometricPrimitive::GetBSSRDF(const DifferentialGeometry &dg,
                                  const Transform &ObjectToWorld,
                                  MemoryArena &arena) const {
    DifferentialGeometry dgs;
    shape->GetShadingGeometry(ObjectToWorld, dg, &dgs);
    return material->GetBSSRDF(dg, dgs, arena);
}


const MultipoleBSSRDF* GeometricPrimitive::GetMultipoleBSSRDF(const DifferentialGeometry &dg,
	const Transform &ObjectToWorld, MemoryArena &arena) const
{
    DifferentialGeometry dgs;
    shape->GetShadingGeometry(ObjectToWorld, dg, &dgs);
    return material->GetMultipoleBSSRDF(dg, dgs, arena);
}


const Shape* GeometricPrimitive::GetShape() const {
	return shape.GetPtr();
}


const Material* GeometricPrimitive::GetMaterial() const {
	return material.GetPtr();
}


class LayeredGeometricPrimitiveMain : public LayeredGeometricPrimitive {
public:
    // LayeredGeometricPrimitiveMain Public Methods
	LayeredGeometricPrimitiveMain(const Reference<ShrinkableShape>& s,
		const Reference<LayeredMaterial>& m, const AreaLight* a);
 	const ShrinkableShape* GetShrinkableShape() const;
	const LayeredMaterial* GetLayeredMaterial() const;
	void Refine(vector<Reference<Primitive> > &refined) const override;
	bool IntersectInternal(const Ray& r, uint32_t primitiveId,
		Intersection* isect, int* layerIndex) const override;
	BSDF *GetLayeredBSDF(int layerIndex,
		const DifferentialGeometry &dg,
		const Transform &ObjectToWorld, MemoryArena &arena) const override;
	BSSRDF *GetLayeredBSSRDF(int layerIndex,
		const DifferentialGeometry &dg,
		const Transform &ObjectToWorld, MemoryArena &arena) const override;
private:
	friend class LayeredGeometricPrimitivePatch;
    // LayeredGeometricPrimitiveMain Private Data
	vector<Reference<ShrinkableShape> > internalShapes;
	Reference<Aggregate> internalAggregate;
	vector<Reference<Primitive> > fullyRefined;
};


class LayeredGeometricPrimitivePatch : public LayeredGeometricPrimitive {
public:
    // LayeredGeometricPrimitivePatch Public Methods
	LayeredGeometricPrimitivePatch(const Reference<Shape>& s,
		const Reference<const LayeredGeometricPrimitiveMain>& main);
    void Refine(vector<Reference<Primitive> > &refined) const override;
	bool IntersectInternal(const Ray& r, uint32_t primitiveId,
		Intersection* isect, int* layerIndex) const override;
	BSDF *GetLayeredBSDF(int layerIndex,
		const DifferentialGeometry &dg,
		const Transform &ObjectToWorld, MemoryArena &arena) const override;
	BSSRDF *GetLayeredBSSRDF(int layerIndex,
		const DifferentialGeometry &dg,
		const Transform &ObjectToWorld, MemoryArena &arena) const override;
private:
    // LayeredGeometricPrimitivePatch Private Data
	Reference<const LayeredGeometricPrimitiveMain> main;
};


LayeredGeometricPrimitive::LayeredGeometricPrimitive(const Reference<Shape>& s,
	const Reference<LayeredMaterial>& m, const AreaLight *a)
	: GeometricPrimitive(s, m, a)
{
}

LayeredGeometricPrimitiveMain::LayeredGeometricPrimitiveMain(
	const Reference<ShrinkableShape>& s,
	const Reference<LayeredMaterial>& m, const AreaLight* a)
	: LayeredGeometricPrimitive(s, m, a)
{
	vector<float_type> thicknesses = m->GetLayerThickness();
	// Add shrinked primitives for internal intersection test
	this->FullyRefine(fullyRefined);
	vector<Reference<Primitive> > internalPrimitives(fullyRefined);
	internalShapes.push_back(s);

	float_type thickness = 0;
	for (size_t i = 0; i < thicknesses.size(); i++) {
		thickness += thicknesses[i];
		Reference<ShrinkableShape> shrinked = s->Shrink(thickness);
		internalShapes.push_back(shrinked);
		internalPrimitives.push_back(new GeometricPrimitive(
			shrinked, new LayeredMaterialWrapper(m, i + 1), a
		));
	}
	internalAggregate = CreateBVHAccelerator(internalPrimitives,
		ParamSet());
}

void LayeredGeometricPrimitiveMain::Refine(vector<Reference<Primitive> > &refined) const {
	if (fullyRefined.size()) {
		refined.insert(refined.end(), fullyRefined.begin(), fullyRefined.end());
	} else {
		vector<Reference<Shape> > r;
		shape->Refine(r);
		for (uint32_t i = 0; i < r.size(); ++i) {
			LayeredGeometricPrimitivePatch* patch =
				new LayeredGeometricPrimitivePatch(r[i], this);
			refined.push_back(patch);
		}
	}
}

const ShrinkableShape* LayeredGeometricPrimitiveMain::GetShrinkableShape() const {
	return static_cast<const ShrinkableShape*>(shape.GetPtr());
}

const LayeredMaterial* LayeredGeometricPrimitiveMain::GetLayeredMaterial() const {
	return static_cast<const LayeredMaterial*>(material.GetPtr());
}

bool LayeredGeometricPrimitiveMain::IntersectInternal(const Ray& r,
	uint32_t primitiveId, Intersection* isect, int* layerIndex) const
{
	if (!internalAggregate->IntersectExcept(r, isect, primitiveId))
		return false;
	
	if (isect->primitive->ToLayered()) {
		// outer layer (shares LayeredGeometricPrimitivePatch s with the scene)
		*layerIndex = 0;
	} else {
		const LayeredMaterialWrapper* material =
			static_cast<const LayeredMaterialWrapper*>(
			static_cast<const GeometricPrimitive*>(isect->primitive)->GetMaterial());

		*layerIndex = material->GetLayerIndex();
	}
	return true;
}

BSDF* LayeredGeometricPrimitiveMain::GetLayeredBSDF(int layerIndex,
	const DifferentialGeometry &dg,
	const Transform &ObjectToWorld, MemoryArena &arena) const
{
	DifferentialGeometry dgs;
	internalShapes[layerIndex]->GetShadingGeometry(ObjectToWorld, dg, &dgs);
	return GetLayeredMaterial()->GetLayeredBSDF(layerIndex, dg, dgs, arena);
}

BSSRDF* LayeredGeometricPrimitiveMain::GetLayeredBSSRDF(int layerIndex,
	const DifferentialGeometry &dg,
	const Transform &ObjectToWorld, MemoryArena &arena) const
{
	DifferentialGeometry dgs;
	internalShapes[layerIndex]->GetShadingGeometry(ObjectToWorld, dg, &dgs);
	return GetLayeredMaterial()->GetLayeredBSSRDF(layerIndex, dg, dgs, arena);
}

LayeredGeometricPrimitivePatch::LayeredGeometricPrimitivePatch(
	const Reference<Shape>& s,
	const Reference<const LayeredGeometricPrimitiveMain>& main)
	: LayeredGeometricPrimitive(s, static_cast<LayeredMaterial*>(main->material.GetPtr()),
		main->GetAreaLight()), main(main)
{
}

void LayeredGeometricPrimitivePatch::Refine(vector<Reference<Primitive> > &refined) const {
    vector<Reference<Shape> > r;
    shape->Refine(r);
    for (uint32_t i = 0; i < r.size(); ++i) {
        LayeredGeometricPrimitivePatch* patch =
			new LayeredGeometricPrimitivePatch(r[i], main);
        refined.push_back(patch);
    }
}

bool LayeredGeometricPrimitivePatch::IntersectInternal(const Ray& r,
	uint32_t primitiveId, Intersection* isect, int* layerIndex) const
{
	return main->IntersectInternal(r, primitiveId, isect, layerIndex);
}

BSDF* LayeredGeometricPrimitivePatch::GetLayeredBSDF(int layerIndex,
	const DifferentialGeometry &dg,
	const Transform &ObjectToWorld, MemoryArena &arena) const
{
	return main->GetLayeredBSDF(layerIndex, dg, ObjectToWorld, arena);
}

BSSRDF* LayeredGeometricPrimitivePatch::GetLayeredBSSRDF(int layerIndex,
	const DifferentialGeometry &dg,
	const Transform &ObjectToWorld, MemoryArena &arena) const
{
	return main->GetLayeredBSSRDF(layerIndex, dg, ObjectToWorld, arena);
}


GeometricPrimitive* CreateGeometricPrimitive(
	const Reference<Shape>& s, const Reference<Material>& m,
	const AreaLight* a)
{
	// Create LayeredGeometricPrimitives for layered materials
	try {
		Reference<LayeredMaterial> lmtl = dynamic_cast<LayeredMaterial*>(m.GetPtr());
		Reference<ShrinkableShape> sshp = dynamic_cast<ShrinkableShape*>(s.GetPtr());
		if (!lmtl.GetPtr() || !sshp.GetPtr())
			throw std::bad_cast();
		return new LayeredGeometricPrimitiveMain(sshp, lmtl, a);
	} catch (const std::bad_cast&) {
		return new GeometricPrimitive(s, m, a);
	}
}
