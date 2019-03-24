#include "sphere.h"

#include <cmath>

#include  "../bsdf.h"
#include "../misc/sphere_drawing.h"

namespace CGL { namespace StaticScene {

bool Sphere::test(const Ray& r, double& t1, double& t2) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection test.
  // Return true if there are intersections and writing the
  // smaller of the two intersection times in t1 and the larger in t2.
  Vector3D roo = r.o-o;
  double a = dot(r.d, r.d); double b = dot(2*roo, r.d); double c = dot(roo, roo) - r2;
  double discrim = b*b - 4*a*c;
  
  if (discrim < 0) {
    return false;
  }

  t1 = (-b - sqrt(discrim)) / (2*a);
  t2 = (-b + sqrt(discrim)) / (2*a);

  bool intersects1 = t1 <= r.max_t && t1 >= r.min_t;  bool intersects2 = t2 <= r.max_t && t2 >= r.min_t;

  if (intersects1) {
    r.max_t = t1;
  } else if (intersects2) {
    r.max_t = t2;
  }

  return intersects1 || intersects2;

}

bool Sphere::intersect(const Ray& r) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection.
  // Note that you might want to use the the Sphere::test helper here.
  double t1, t2;
  return test(r, t1, t2);
}

bool Sphere::intersect(const Ray& r, Intersection *i) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection.
  // Note again that you might want to use the the Sphere::test helper here.
  // When an intersection takes place, the Intersection data should be updated
  // correspondingly.
  double t1, t2;
  bool intersects = test(r, t1, t2);
  
  if (!intersects) {
    return false;
  }

  i->t = t1;
  i->n = normal(r.o + t1 * r.d);
  i->primitive = this;
  i->bsdf = get_bsdf();

  return true;

}

void Sphere::draw(const Color& c, float alpha) const {
  Misc::draw_sphere_opengl(o, r, c);
}

void Sphere::drawOutline(const Color& c, float alpha) const {
    //Misc::draw_sphere_opengl(o, r, c);
}


} // namespace StaticScene
} // namespace CGL
