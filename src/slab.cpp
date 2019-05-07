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

#define PATH_COUNT 100
#define EXT_COEFF 0.5

namespace CGL {

// sample a ray given circle given collector's mean and variance
Vector3D Collector::sample() {

	std::normal_distribution<double> deltaX(0.0, varianceX);
	std::normal_distribution<double> deltaY(0.0, varianceY);
	std::uniform_real_distribution<> angle(0.0, 2.0*M_PI);

	double alphaX = deltaX(gen);
	double alphaY = deltaY(gen);
	double theta = angle(gen);

	Vector3D displacement = Vector3D(alphaX * std::cos(theta), alphaY * std::sin(theta), 0.0);
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
	std::uniform_real_distribution<double> angle(0.0, 2.0*M_PI);
	std::uniform_real_distribution<double> var(0.1, 10.0);

	double alpha = delta(gen);
	double theta = angle(gen);
	double variance = var(gen);

	Vector3D displacement = Vector3D(alpha * std::cos(theta), alpha * std::sin(theta), 0.0);

	return Collector(origin + origin2w * displacement, variance, origin2w, cloud, gen);

}

// TODO: Implement light transport as described in:
// http://www-evasion.imag.fr/Membres/Antoine.Bouthors/research/phd/thesis/thesis.pdf
std::vector<Collector> Slab::transport_sample(const Ray &exit, Vector3D incident, std::vector<double> &L, BSDF *bsdf, int init_depth) {

	L.clear();
	// maximum scattering order
	int O = init_depth;
	// path count
	int N = PATH_COUNT;
	// assume we operate with 1 incident light direction for now

	// vector of collector means
	std::vector<double> x;	std::vector<double> y;
	// intermediate vectors used to compute variance
	std::vector<double> X_L2;	std::vector<double> Y_L2;

	std::uniform_real_distribution<double> u(0.0, 1.0);

	for (int i = 0; i < O; i++) {
		L.push_back(0);	x.push_back(0);	y.push_back(0); X_L2.push_back(0);	Y_L2.push_back(0);
	}

	Matrix3x3 w2o = origin2w.T();
	double phase = dot(origin2w[2], incident);

	for (int i = 0; i < N; i++) {
		Vector3D curr = w2o * (origin - exit.o);
		Vector3D v = w2o * exit.d.unit();
		Vector3D v_new = Vector3D();
		Vector3D c = Vector3D();
		float pdf;
		double l;
		double d;
		double ext_phase;
		for (int o = 0; curr[2] < 0 && o < O && curr[2] > -thickness; o++) {
			l = std::log(EXT_COEFF * u(gen));
			curr += l * v;
			// std::cout << o << ", " << curr << std::endl;

			d = -curr[2] / phase;
			ext_phase = std::exp(-EXT_COEFF * d) * hg_phase(v, incident);
			L[o] += ext_phase;
			c = curr + d * incident;
			x[o] += c[0] * ext_phase;
			y[o] += c[1] * ext_phase;
			X_L2[o] += c[0] * c[0] * ext_phase;
			Y_L2[o] += c[1] * c[1] * ext_phase;

			bsdf->sample_f(-v, &v_new, &pdf);
			v = v_new;
		}
	}

	std::vector<Collector> result;
	for (int o = 0; o < O; o++) {
		double sigma_x = std::sqrt((X_L2[o] - x[o] * x[o]) / L[o]);
		double sigma_y = std::sqrt((Y_L2[o] - y[o] * y[o]) / L[o]);
		x[o] /= L[o];	y[o] /= L[o];
		L[o] /= N;
		Vector3D mean = origin + origin2w * Vector3D(x[o], y[o], 0.0);
		Matrix3x3 o2w_new = Matrix3x3();
		Collector c_new = Collector(mean, sigma_x, sigma_y, origin2w, cloud, gen);
		c_new.order = o+1;
		result.push_back(c_new);
	}

	return result;

}

}