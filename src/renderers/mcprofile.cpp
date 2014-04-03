
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

#include "mcprofile.h"
#include "paramset.h"
#include "parallel.h"
#include "progressreporter.h"
#include "reflection.h"
#include <fstream>
#include <sstream>
#include "multipole/MultipoleProfileCalculator/MultipoleProfileCalculator.h"

using namespace std;

// Use double precision in simulation
typedef VectorBase<double> DVector;
typedef PointBase<double> DPoint;
typedef NormalBase<double> DNormal;
typedef RayBase<double> DRay;

DVector UniformSampleSphereD(double u1, double u2) {
    double z = 1. - 2. * u1;
    double r = sqrt(max(0., 1. - z*z));
    double phi = 2. * M_PI * u2;
    double x = r * cos(phi);
    double y = r * sin(phi);
    return DVector(x, y, z);
}

DVector UniformSampleHemisphereD(double u1, double u2) {
    double z = u1;
    double r = sqrt(max(0., 1. - z*z));
    double phi = 2. * M_PI * u2;
    double x = r * cosf(phi);
    double y = r * sinf(phi);
    return DVector(x, y, z);
}

void ConcentricSampleDiskD(double u1, double u2, double *dx, double *dy) {
    double r, theta;
    // Map uniform random numbers to $[-1,1]^2$
    double sx = 2 * u1 - 1;
    double sy = 2 * u2 - 1;

    // Map square to $(r,\theta)$

    // Handle degeneracy at the origin
    if (sx == 0.0 && sy == 0.0) {
        *dx = 0.0;
        *dy = 0.0;
        return;
    }
    if (sx >= -sy) {
        if (sx > sy) {
            // Handle first region of disk
            r = sx;
            if (sy > 0.0) theta = sy/r;
            else          theta = 8.0 + sy/r;
        }
        else {
            // Handle second region of disk
            r = sy;
            theta = 2.0 - sx/r;
        }
    }
    else {
        if (sx <= sy) {
            // Handle third region of disk
            r = -sx;
            theta = 4.0 - sy/r;
        }
        else {
            // Handle fourth region of disk
            r = -sy;
            theta = 6.0 + sx/r;
        }
    }
    theta *= M_PI / 4.;
    *dx = r * cos(theta);
    *dy = r * sin(theta);
}

DVector CosineSampleHemisphereD(double u1, double u2) {
    DVector ret;
    ConcentricSampleDiskD(u1, u2, &ret.x, &ret.y);
    ret.z = sqrt(max(0., 1. - ret.x*ret.x - ret.y*ret.y));
    return ret;
}

// How layers and interfaces are labeled
//------------------------------------ Interface 0
//            Layer 0
//------------------------------------ Interface 1
//            Layer 1
//------------------------------------ Interface 2
//            Layer 2
//------------------------------------ Interface 3

struct MiniIsect {
	MiniIsect() {}
	MiniIsect(int iid, double invIOR, const DPoint& p) {
		interfaceId = iid;
		invRelativeIOR = invIOR;
		this->p = p;
	}

	int interfaceId;
	double invRelativeIOR;
	DPoint p;
};

// represents our mini scene
class MiniScene {
public:
	MiniScene(const Layer* pLayers, int nLayers) : layers(pLayers, pLayers + nLayers) {
		double depth = 0.;
		interfaceDepths.push_back(depth);
		for (const Layer& l : layers) {
			interfaceDepths.push_back(depth += l.thickness);
		}
	}

	bool Intersect(const DRay& ray, MiniIsect* isect, int currentLayer) const {
		double costheta = ray.d.z;
		if (costheta == 0.) return false;

		if (ray.d.z > 0.) {
			// Going inwards, check only the layer below currentLayer
			double nextDepth = interfaceDepths[currentLayer + 1];
			if (nextDepth - ray.o.z < costheta * ray.maxt) {
				// hit
				isect->interfaceId = currentLayer + 1;
				double t = (nextDepth - ray.o.z) / costheta;
				isect->p = ray(t);
				double nextIOR = (currentLayer + 1 == layers.size()) ? 1. : layers[currentLayer + 1].ior;
				double currentIOR = layers[currentLayer].ior;
				isect->invRelativeIOR = currentIOR / nextIOR;
				return true;
			}
		} else {
			// Going outwards, check only the layer above
			double nextDepth = interfaceDepths[currentLayer];
			if (nextDepth - ray.o.z > costheta * ray.maxt) {
				// hit
				isect->interfaceId = currentLayer;
				double t = (nextDepth - ray.o.z) / costheta;
				isect->p = ray(t);
				double nextIOR = (currentLayer == 0) ? 1. : layers[currentLayer - 1].ior;
				double currentIOR = layers[currentLayer].ior;
				isect->invRelativeIOR = currentIOR / nextIOR;
				return true;
			}
		}
		return false;
	}

	const Layer& getLayer(int layerId) const { return layers[layerId]; }
	int getNumLayers() const { return (int)layers.size(); }
private:
	vector<Layer> layers;
	vector<double> interfaceDepths;
};


class ProfileRendererTask : public Task {
public:
	ProfileRendererTask(const MiniScene* scene, MCProfile& profile,
		uint64_t photons, double extent, int segments, Mutex& mtx,
		ProgressReporter& pr, int taskId)
		: scene(scene), profile(profile), nPhotons(photons),
		  extent(extent), nSegments(segments), mutex(mtx),
		  reporter(pr), taskId(taskId)
	{
	}

	void Run() override;
private:
	const MiniScene* scene;
	MCProfile& profile;
	uint64_t nPhotons;
	double extent;
	int nSegments;
	Mutex& mutex;
	ProgressReporter& reporter;
	int taskId;

	void TraceSinglePhoton(MCProfile& profile, RNG& rng);
	void CombineProfile(const MCProfile& profile);
};


void ProfileRendererTask::Run() {
	RNG rng(89 * taskId);

	MCProfile localProfile(nSegments);

	for (uint64_t i = 0; i < nPhotons; i++) {
		TraceSinglePhoton(localProfile, rng);
	}

	CombineProfile(localProfile);
	reporter.Update();
}


void ProfileRendererTask::TraceSinglePhoton(MCProfile& profile, RNG& rng) {
	// The core part of tracing a photon through the mini scene

	// Setup
	DVector d0 = DVector(0, 0, 1);
	DRay ray(DPoint(0., 0., 0.), d0, 0.);
	int currentLayerId = 0;
	MiniIsect isect;
	double pathThroughput = 1.;
	double len = 0.;

	// Loop
	while (currentLayerId >= 0 && currentLayerId < scene->getNumLayers()) {
		const Layer& currentLayer = scene->getLayer(currentLayerId);
		int targetLayerId = currentLayerId;
		double scatterMfp = 1. / currentLayer.musp;
		len *= scatterMfp;

		do {
			// Sample mfp(s) according to the exponential distribution: musp * exp(-mfp * musp)
			if (len == 0.)
				len = min(-log(rng.RandomDouble()), 1e7) * scatterMfp;

			ray.maxt = len;
			// Test intersection
			if (scene->Intersect(ray, &isect, currentLayerId)) {
				// Attenuation due to absorption
				double distance = (isect.p - ray.o).Length();
				pathThroughput *= exp(-currentLayer.mua * distance);
				len = max(1e-7 * scatterMfp, len - distance);

				// Compute outgoing direction
				double cosi = Clamp(ray.d.z, -1., 1.);
				bool up = ray.d.z < 0.;

				// Compute _sint_ using Snell's law
				double sint2 = (1. - cosi * cosi) * isect.invRelativeIOR * isect.invRelativeIOR;
				double cost;
				double fresnel;
				bool reflect;
				if (sint2 >= 1.) {
					// Handle total internal reflection
					fresnel = 1.;
					reflect = true;
				}
				else {
					cost = sqrt(max(0., 1. - sint2));
					fresnel = FrDiel(abs(cosi), cost, isect.invRelativeIOR, 1.);
					reflect = rng.RandomDouble() < fresnel;
				}

				// Generate new ray
				DVector newd;
				if (reflect) {
					newd = DVector(ray.d.x, ray.d.y, -ray.d.z);
				} else {
					targetLayerId = up ? currentLayerId - 1 : currentLayerId + 1;
					newd = DVector(isect.invRelativeIOR * ray.d.x, isect.invRelativeIOR * ray.d.y, up ? -cost : cost);
				}
				ray = DRay(isect.p, newd, 0.);

				// Book-keeping
				if (currentLayerId != targetLayerId) {
					if (isect.interfaceId == 0) {
						double rd = DVector(isect.p.x, isect.p.y, 0.).Length();
						int segmentId = (int)(rd * nSegments / extent);
						if (rd < extent && segmentId < nSegments)
							profile[segmentId].reflectance += pathThroughput;
					} else if (isect.interfaceId == scene->getNumLayers()) {
						double rd = DVector(isect.p.x, isect.p.y, 0.).Length();
						int segmentId = (int)(rd * nSegments / extent);
						if (rd < extent && segmentId < nSegments)
							profile[segmentId].transmittance += pathThroughput;
					}
				}
			} else {
				DPoint p = ray(len);

				// Attenuation due to absorption
				pathThroughput *= exp(-currentLayer.mua * len);
				len = 0.;

				// Random scatter
				DVector newd = UniformSampleSphereD(rng.RandomDouble(), rng.RandomDouble());
				ray = DRay(p, newd, 0.);
			}
		} while (currentLayerId == targetLayerId);

		currentLayerId = targetLayerId;

		// Possibly terminate the path
		if (pathThroughput < 1e-5) {
			double continueProbability = pathThroughput * 1e5;
			if (rng.RandomDouble() > continueProbability)
				break;
			pathThroughput /= continueProbability;
		}

		// convert remaining length to go
		len *= currentLayer.musp;
	}
}


void ProfileRendererTask::CombineProfile(const MCProfile& profile) {
	// Using atomic memory ops probably won't win here, coz there
	// are so many operations (typically thousands), let alone the
	// false sharing
	MutexLock lock(mutex);
	for (int i = 0; i < nSegments; i++) {
		this->profile[i].reflectance += profile[i].reflectance;
		this->profile[i].transmittance += profile[i].transmittance;
	}
}


class MultipoleReferenceTask : public Task {
public:
	MultipoleReferenceTask(const Layer* pLayer, int nLayers, MCProfile& profile,
		double* ref, double* trans, double extent, int nSegments, bool lerpOnThinSlab,
		ProgressReporter& pr)
		: layers(pLayer, pLayer + nLayers), profile(profile), ref(ref), trans(trans),
		  extent(extent), nSegments(nSegments), lerpOnThinSlab(lerpOnThinSlab),
		  reporter(pr)
	{
		
	}

	void Run() override;
private:
	vector<Layer> layers;
	MCProfile& profile;
	double* ref;
	double* trans;
	double extent;
	int nSegments;
	bool lerpOnThinSlab;
	ProgressReporter& reporter;
};


void MultipoleReferenceTask::Run() {
	// Prepare _MPC_LayerSpec_s
	vector<MPC_LayerSpec> layerSpecs;
	for (size_t i = 0; i < layers.size(); i++) {
		const Layer& layer = layers[i];
		MPC_LayerSpec spec;
		spec.mua = layer.mua;
		spec.musp = layer.musp;
		spec.ior = layer.ior;
		spec.thickness = layer.thickness;
		spec.g_HG = 0.f;
		layerSpecs.push_back(spec);
	}
	// Prepare MPC_Options
	MPC_Options options;
	options.desiredLength = 1024; // able to afford this
	options.desiredStepSize = (float)(extent * 1.01 / options.desiredLength);
	options.lerpOnThinSlab = lerpOnThinSlab;

	// Compute
	MPC_Output* pOutput;
	MPC_ComputeDiffusionProfile((uint32)layerSpecs.size(), &layerSpecs[0], &options,
		&pOutput);
	*ref = pOutput->totalReflectance;
	*trans = pOutput->totalTransmittance;

	// Resample
	vector<float> distArray;
	distArray.reserve(nSegments);
	for (int i = 0; i < nSegments; i++) {
		distArray.push_back((float)((i + .5) * (extent / nSegments)));
	}
	MPC_ResampleDistribution(pOutput, nSegments, &distArray[0]);

	// Gather results
	profile.clear();
	profile.reserve(nSegments);
	for (int i = 0; i < nSegments; i++) {
		profile.push_back(MCProfileEntry(
			pOutput->pReflectance[i], pOutput->pTransmittance[i]));
	}

	MPC_FreeOutput(pOutput);
	reporter.Update();
}


static ostream& printProfile(ostream& os, const MCProfile& p, const double MCProfileEntry::*field,
	const vector<double>* rarray = NULL)
{
	if (!p.size()) return os;

	bool mulr = rarray != NULL;
	os << p[0].*field * (mulr ? rarray->at(0) : 1.);
	for (size_t i = 1; i < p.size(); i++) {
		os << "\t" << p[i].*field * (mulr ? rarray->at(i) : 1.);
	}

	return os;
}


void MonteCarloProfileRenderer::Render(const Scene* scene) {
	Info("Monte Carlo Profile Renderer");
	for (int i = 0; i < layers.size(); i++) {
		const Layer& l = layers[i];
		Info("Layer %d: mua=%f, musp=%f, ior=%f, thickness=%f",
			i, l.mua, l.musp, l.ior, l.thickness);
	}
	Info("Gather samples within %f mean-free-paths, divided into %d segments.",
		mfpRange, nSegments);
	Info("Output filename: %s", filename.c_str());

	// Compute extent of the profile
	double mfpTotal = 0.;
	for (const Layer& l : layers) {
		mfpTotal += 1. / (l.mua + l.musp);
	}
	double mfp = mfpTotal / (double)layers.size();

	double extent = mfpRange * mfp;

	// Create scene
	MiniScene miniScene(&layers[0], layers.size());
	MCProfile profile(nSegments);

	// Create _ProfileRendererTask_s
	int nTasks = (int)min<uint64_t>((nPhotons + 65535) / 65536, 1048576);
	Mutex* mutex = Mutex::Create();
	vector<Task*> tasks;
	tasks.reserve(nTasks);
	ProgressReporter reporter(nTasks, "Random Walk");
	for (int i = 0; i < nTasks; i++) {
		tasks.push_back(new ProfileRendererTask(&miniScene, profile,
			(uint64_t)(i + 1) * nPhotons / nTasks - (uint64_t)i * nPhotons / nTasks,
			extent, nSegments, *mutex, reporter, i));
	}
	EnqueueTasks(tasks);
	WaitForAllTasks();
	for (Task* task : tasks)
		delete task;
	Mutex::Destroy(mutex);
	reporter.Done();

	// Compute reference profiles genereated by multipole
	MCProfile noLerpProfile, lerpProfile;
	double totalNoLerpReflectance = 0.;
	double totalNoLerpTransmittance = 0.;
	double totalLerpReflectance = 0.;
	double totalLerpTransmittance = 0.;
	ProgressReporter multipoleReporter(2, "Profile");
	vector<Task*> multipoleTasks;
	multipoleTasks.push_back(new MultipoleReferenceTask(&layers[0], (int)layers.size(),
		noLerpProfile, &totalNoLerpReflectance, &totalNoLerpTransmittance,
		extent, nSegments, false, multipoleReporter));
	multipoleTasks.push_back(new MultipoleReferenceTask(&layers[0], (int)layers.size(),
		lerpProfile, &totalLerpReflectance, &totalLerpTransmittance,
		extent, nSegments, true, multipoleReporter));
	EnqueueTasks(multipoleTasks);
	WaitForAllTasks();
	for (Task* task : multipoleTasks)
		delete task;
	MPC_ClearCache();
	multipoleReporter.Done();

	// Normalize result
	Info("Gathering results...");
	double totalReflectance = 0.;
	double totalTransmittance = 0.;
	for (int i = 0; i < nSegments; i++) {
		MCProfileEntry& entry = profile[i];
		// normalize by total photons and area of ring
		double sum = (double)(2 * i + 1) * extent / nSegments;
		double h = extent / nSegments;
		double area = M_PI * sum * h;
		double factor = (double)nPhotons * area;
		totalReflectance += entry.reflectance;
		totalTransmittance += entry.transmittance;
		entry.reflectance /= factor;
		entry.transmittance /= factor;
	}
	totalReflectance /= (double)nPhotons;
	totalTransmittance /= (double)nPhotons;

	// Output result
	ofstream out(filename, ios::out | ios::trunc);
	if (!out) {
		Error("Failed to open output file: %s", filename.c_str());
		abort();
	}
	// header
	out << "Name\tTotal";
	vector<double> rarray;
	rarray.reserve(nSegments);
	for (int i = 0; i < nSegments; i++) {
		double distance = (i + .5) * extent / nSegments;
		out << "\t" << distance;
		rarray.push_back(distance);
	}
	out << endl;
	// content
	out << "Monte-Carlo Reflectance\t" << totalReflectance << "\t";
	printProfile(out, profile, &MCProfileEntry::reflectance) << endl;
	out << "Monte-Carlo Transmittance\t" << totalTransmittance << "\t";
	printProfile(out, profile, &MCProfileEntry::transmittance) << endl;
	out << "Multipole Reflectance\t" << totalNoLerpReflectance << "\t";
	printProfile(out, noLerpProfile, &MCProfileEntry::reflectance) << endl;
	out << "Multipole Transmittance\t" << totalNoLerpTransmittance << "\t";
	printProfile(out, noLerpProfile, &MCProfileEntry::transmittance) << endl;
	out << "Lerped Reflectance\t" << totalLerpReflectance << "\t";
	printProfile(out, lerpProfile, &MCProfileEntry::reflectance) << endl;
	out << "Lerped Transmittance\t" << totalLerpTransmittance << "\t";
	printProfile(out, lerpProfile, &MCProfileEntry::transmittance) << endl;
	// print again as r*Rd(r)
	out << "Monte-Carlo Reflectance\t" << totalReflectance << "\t";
	printProfile(out, profile, &MCProfileEntry::reflectance, &rarray) << endl;
	out << "Monte-Carlo Transmittance\t" << totalTransmittance << "\t";
	printProfile(out, profile, &MCProfileEntry::transmittance, &rarray) << endl;
	out << "Multipole Reflectance\t" << totalNoLerpReflectance << "\t";
	printProfile(out, noLerpProfile, &MCProfileEntry::reflectance, &rarray) << endl;
	out << "Multipole Transmittance\t" << totalNoLerpTransmittance << "\t";
	printProfile(out, noLerpProfile, &MCProfileEntry::transmittance, &rarray) << endl;
	out << "Lerped Reflectance\t" << totalLerpReflectance << "\t";
	printProfile(out, lerpProfile, &MCProfileEntry::reflectance, &rarray) << endl;
	out << "Lerped Transmittance\t" << totalLerpTransmittance << "\t";
	printProfile(out, lerpProfile, &MCProfileEntry::transmittance, &rarray) << endl;
	out.close();
}

Spectrum MonteCarloProfileRenderer::Li(const Scene* scene, const RayDifferential& ray,
	const Sample* sample, RNG& rng, MemoryArena& arena,
	Intersection* isect, Spectrum* T) const
{
	return Spectrum(0.f);
}

Spectrum MonteCarloProfileRenderer::Transmittance(const Scene* scene, const RayDifferential& ray,
	const Sample* sample, RNG& rng, MemoryArena& arena) const
{
	return Spectrum(0.f);
}

MonteCarloProfileRenderer* CreateMonteCarloProfileRenderer(const ParamSet& params) {
	int nLayers;
	const Layer* layers = params.FindLayer("layers", &nLayers);
	if (!layers || nLayers < 1) {
		Error("No layers param set for MCProfileRenderer.");
		abort();
	}
	float mfpRange = params.FindOneFloat("mfprange", 16.f);
	int segments = params.FindOneInt("segments", 1024);
	string strPhotons = params.FindOneString("photons", "100");
	uint64_t photons = _strtoui64(strPhotons.c_str(), NULL, 10);
	string filename = params.FindOneFilename("filename", "");

	return new MonteCarloProfileRenderer(layers, nLayers, mfpRange, segments,
		photons, filename);
}
