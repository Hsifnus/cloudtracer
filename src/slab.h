#ifndef CGL_SLAB_H
#define CGL_SLAB_H

#include "CGL/matrix3x3.h"
#include "CGL/vector3D.h"
#include "ray.h"
#include "static_scene/primitive.h"
#include <random>

namespace CGL {

class Collector;
class Slab;

class Collector {
 public:
  Collector(Vector3D m, double v, Matrix3x3 o2w, const StaticScene::Primitive *c, std::default_random_engine g) : mean(m), variance(v), origin2w(o2w), cloud(c), gen(g) { }

  Vector3D mean;
  double variance;
  Matrix3x3 origin2w;
  const StaticScene::Primitive* cloud;
  std::default_random_engine gen;

  Vector3D sample();
  bool same_mesh(const StaticScene::Primitive* other);
  bool project(const Ray &r, const StaticScene::Intersection& isect, double &delta);
  Slab generate_slab(Vector3D origin);

};

class Slab {
 public:
  Slab(Vector3D o, double t, Matrix3x3 o2w, const StaticScene::Primitive* c, std::default_random_engine g): origin(o), thickness(t), origin2w(o2w), cloud(c), gen(g) { }

  Vector3D origin;
  double thickness;
  Matrix3x3 origin2w;
  const StaticScene::Primitive* cloud;
  std::default_random_engine gen;

  Collector basic_sample();
  std::vector<Collector> transport_sample();

};

}

#endif  // CGL_SLAB_H