#include "bvh.h"

#include "CGL/CGL.h"
#include "static_scene/triangle.h"

#include <limits>
#include <iostream>
#include <stack>

using namespace std;

namespace CGL { namespace StaticScene {

BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
                   size_t max_leaf_size) {

  root = construct_bvh(_primitives, max_leaf_size);

}

BVHAccel::~BVHAccel() {
  if (root) delete root;
}

BBox BVHAccel::get_bbox() const {
  return root->bb;
}

void BVHAccel::draw(BVHNode *node, const Color& c, float alpha) const {
  if (node->isLeaf()) {
    for (Primitive *p : *(node->prims))
      p->draw(c, alpha);
  } else {
    draw(node->l, c, alpha);
    draw(node->r, c, alpha);
  }
}

void BVHAccel::drawOutline(BVHNode *node, const Color& c, float alpha) const {
  if (node->isLeaf()) {
    for (Primitive *p : *(node->prims))
      p->drawOutline(c, alpha);
  } else {
    drawOutline(node->l, c, alpha);
    drawOutline(node->r, c, alpha);
  }
}

BVHNode *BVHAccel::construct_bvh(const std::vector<Primitive*>& prims, size_t max_leaf_size) {
  
  // TODO (Part 2.1):
  // Construct a BVH from the given vector of primitives and maximum leaf
  // size configuration. The starter code build a BVH aggregate with a
  // single leaf node (which is also the root) that encloses all the
  // primitives.

  bool naive = false;
  if (naive) {
    BBox centroid_box, bbox;

    for (Primitive *p : prims) {
        BBox bb = p->get_bbox();
        bbox.expand(bb);
        Vector3D c = bb.centroid();
        centroid_box.expand(c);
    }

    BVHNode *node = new BVHNode(bbox);


    node->prims = new vector<Primitive *>(prims);
    return node;
  } else {
    Vector3D cavg = Vector3D();
    BBox bbox;

    for (Primitive *p : prims) {
      BBox bb = p->get_bbox();
      bbox.expand(bb);
      cavg += bb.centroid() / prims.size();
    }
    BVHNode *node = new BVHNode(bbox);
    std::vector<Primitive*> left;
    std::vector<Primitive*> right;

    node->prims = new vector<Primitive *>(prims);
    // split into two child nodes if too many primitives in leaf
    if (prims.size() > max_leaf_size) {
      // find axis index corresponding to largest bbox dimension
      Vector3D shape = bbox.extent;
      int argmax = shape[0] >= shape[1] && shape[0] >= shape[2] 
        ? 0 : shape[1] >= shape[2] && shape[1] >= shape[0] 
        ? 1 : 2;
      // define split value
      double split = cavg[argmax];
      // place Primitives according to split
      for (Primitive *p : prims) {
        if (p->get_bbox().centroid()[argmax] <= split) {
          left.push_back(p);
        } else {
          right.push_back(p);
        }
      }
      // Set up child BVHNodes
      if (left.size() > 0) {
        node->l = construct_bvh(left, max_leaf_size);
      }
      if (right.size() > 0) {
        node->r = construct_bvh(right, max_leaf_size);
      }
    }
    return node;
  }

}


bool BVHAccel::intersect(const Ray& ray, BVHNode *node) const {

  // TODO (Part 2.3):
  // Fill in the intersect function.
  // Take note that this function has a short-circuit that the
  // Intersection version cannot, since it returns as soon as it finds
  // a hit, it doesn't actually have to find the closest hit.

  double t_min = ray.min_t; double t_max = ray.max_t;
  // ray doesn't intersect bbox within [max_t, min_t]
  if (!node->bb.intersect(ray, t_min, t_max) || std::max(ray.min_t, t_min) > std::min(ray.max_t, t_max)) {
    return false;
  }

  if (node->l == NULL && node->r == NULL) { // is a leaf
    for (Primitive *p : *(node->prims)) {
      if (p->intersect(ray))  {
        total_isects++;
        return true;
      }
    }
    return false;
  } else {
    return (node->l != NULL ? intersect(ray, node->l) : false) || (node->r != NULL ? intersect(ray, node->r) : false);
  }

}

bool BVHAccel::intersect(const Ray& ray, Intersection* i, BVHNode *node) const {

  // TODO (Part 2.3):
  // Fill in the intersect function.
  double t_min = ray.min_t; double t_max = ray.max_t;

  // ray doesn't intersect bbox within [max_t, min_t]
  if (!node->bb.intersect(ray, t_min, t_max) || std::max(ray.min_t, t_min) > std::min(ray.max_t, t_max)) {
    return false;
  }

  if (node->l == NULL && node->r == NULL) { // is a leaf
    bool found = false;
    for (Primitive *p : *(node->prims)) {
      if (p->intersect(ray, i)) {
        total_isects++;
        found = true;
      }
    }
    return found;
  } else {
    bool left = node->l != NULL ? intersect(ray, i, node->l) : false;
    bool right = node->r != NULL ? intersect(ray, i, node->r) : false;
    return left || right;
  }

}

}  // namespace StaticScene
}  // namespace CGL
