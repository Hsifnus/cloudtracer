#include "slab.h"
#include "ray.h"
#include "CGL/vector3D.h"
#include "CGL/matrix3x3.h"
#include "bsdf.h"

#include "static_scene/primitive.h"
#include "static_scene/sphere.h"
#include "static_scene/triangle.h"

#include <math.h>
#include <random>

namespace CGL {

// sample a ray given circle given collector's mean and variance
Vector3D Collector::sample() {

	std::normal_distribution<double> delta(0.0,variance);
	std::uniform_real_distribution<> angle(0.0, 2.0*M_PI);

	double alpha = delta(gen);
	double theta = angle(gen);

	Vector3D displacement = Vector3D(alpha * std::cos(theta), alpha * std::sin(theta), 0.0);
	return mean + origin2w * displacement;

}

bool Collector::same_mesh(const StaticScene::Primitive *other) {

	auto *tri_other = dynamic_cast<const StaticScene::Triangle *> (other);
	auto *sph_other = dynamic_cast<const StaticScene::Sphere *> (other);

	auto *tri = dynamic_cast<const StaticScene::Triangle *> (cloud);
	auto *sph = dynamic_cast<const StaticScene::Sphere *> (cloud);

	if (tri) {
		return tri_other ? tri->mesh == tri_other->mesh : false;
	} else if (sph) {
		return sph_other ? sph->object == sph_other->object : false;
	} else {
		return false;
	}

}

// project a ray down or up to the surface of the cloud the collector is attached to
bool Collector::project(const Ray &r, const StaticScene::Intersection& isect, double &delta) {
	if (same_mesh(isect.primitive)) {

		Vector3D new_mean = r.o + isect.t * r.d;
		delta = (isect.t * r.d).norm();
		mean = new_mean;
		make_coord_space(origin2w, isect.n);

		return true;
	} else {
		return false;
	}
}

Slab Collector::generate_slab(Vector3D origin) {
	Matrix3x3 w2o = origin2w.T();
	double thickness = 1.1 * std::abs((w2o * origin)[2]);
	return Slab(origin, thickness, origin2w, cloud, gen);
}

Collector Slab::basic_sample() {

	std::normal_distribution<double> delta(0.0, thickness);
	std::uniform_real_distribution<> angle(0.0, 2.0*M_PI);
	std::uniform_real_distribution<> var(0.1, 10.0);

	double alpha = delta(gen);
	double theta = angle(gen);
	double variance = var(gen);

	Vector3D displacement = Vector3D(alpha * std::cos(theta), alpha * std::sin(theta), 0.0);

	return Collector(origin + origin2w * displacement, variance, origin2w, cloud, gen);

}

// TODO: Implement light transport as described in:
// http://www-evasion.imag.fr/Membres/Antoine.Bouthors/research/phd/thesis/thesis.pdf
std::vector<Collector> Slab::transport_sample() {

	std::vector<Collector> result;
	return result;

}

}