#include "bbox.h"

#include "GL/glew.h"

#include <algorithm>
#include <iostream>

namespace CGL {

bool BBox::intersect(const Ray& r, double& t0, double& t1) const {

  // TODO (Part 2.2):
  // Implement ray - bounding box intersection test
  // If the ray intersected the bounding box within the range given by
  // t0, t1, update t0 and t1 with the new intersection times.
  Vector3D tmin = Vector3D(); Vector3D tmax = Vector3D();
  for (int i = 0; i < 3; i++) {
    double oi = r.o[i];  double di = r.d[i];
    // Sometimes, the ray may be parallel to the bounding box, in which case we set arbitrarily high/low tmax/tmin values
    tmin[i] = di != 0 ? std::min((min[i] - oi) / di, (max[i] - oi) / di) : std::numeric_limits<double>::max();
    tmax[i] = di != 0 ? std::max((min[i] - oi) / di, (max[i] - oi) / di) : std::numeric_limits<double>::min();
  }

  double tminMax = tmin[0] >= tmin[1] && tmin[0] >= tmin[2] ? tmin[0]
          : tmin[1] >= tmin[2] && tmin[1] >= tmin[0] ? tmin[1] : tmin[2];
  double tmaxMin = tmax[0] <= tmax[1] && tmax[0] <= tmax[2] ? tmax[0]
          : tmax[1] <= tmax[2] && tmax[1] <= tmax[0] ? tmax[0] : tmax[2];

  if (tminMax > tmaxMin) { // no intersection found :(
    return false;
  }

  t0 = tminMax;
  t1 = tmaxMin;

  return true;

}

void BBox::draw(Color c, float alpha) const {

  glColor4f(c.r, c.g, c.b, alpha);

  // top
  glBegin(GL_LINE_STRIP);
  glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(max.x, max.y, max.z);
  glEnd();

  // bottom
  glBegin(GL_LINE_STRIP);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, min.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, min.y, min.z);
  glEnd();

  // side
  glBegin(GL_LINES);
  glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(min.x, min.y, max.z);
  glEnd();

}

std::ostream& operator<<(std::ostream& os, const BBox& b) {
  return os << "BBOX(" << b.min << ", " << b.max << ")";
}

} // namespace CGL
