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
  Collector(Vector3D m, double v, Matrix3x3 o2w, const StaticScene::Primitive *c, std::default_random_engine g) : mean(m), varianceX(v), varianceY(v), origin2w(o2w), cloud(c), gen(g), order(0) { }
  Collector(Vector3D m, double vx, double vy, Matrix3x3 o2w, const StaticScene::Primitive *c, std::default_random_engine g) : mean(m), varianceX(vx), varianceY(vy), origin2w(o2w), cloud(c), gen(g), order(0) { }

  Vector3D mean;
  double varianceX, varianceY;
  Matrix3x3 origin2w;
  const StaticScene::Primitive* cloud;
  std::default_random_engine gen;
  int order;

  Vector3D sample();
  bool same_mesh(const StaticScene::Primitive* other);
  bool project(const Ray &r, const StaticScene::Intersection& isect, double &delta);
  Slab generate_slab(Vector3D origin);
};

inline bool operator==(const Collector &left, const Collector &right) {
  return left.mean == right.mean && left.varianceX == right.varianceX 
    && left.varianceY == right.varianceY && left.origin2w[0] == right.origin2w[0]
    && left.origin2w[1] == right.origin2w[1] && left.origin2w[2] == right.origin2w[2];
}

class Slab {
 public:
  Slab(Vector3D o, double t, Matrix3x3 o2w, const StaticScene::Primitive* c, std::default_random_engine g): origin(o), thickness(t), origin2w(o2w), cloud(c), gen(g) { }

  Vector3D origin;
  double thickness;
  Matrix3x3 origin2w;
  const StaticScene::Primitive* cloud;
  std::default_random_engine gen;

  Collector basic_sample();
  std::vector<Collector> transport_sample(const Ray &exit, const Vector3D incident, std::vector<double> &L, BSDF *bsdf, int init_depth);

};

}

#endif  // CGL_SLAB_H