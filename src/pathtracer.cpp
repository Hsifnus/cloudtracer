#include "pathtracer.h"
#include "bsdf.h"
#include "ray.h"
#include "slab.h"
#include "polar_tex.h"

// #include "lenscamera.h"

#include <stack>
#include <random>
#include <algorithm>
#include <sstream>

#include "CGL/CGL.h"
#include "CGL/vector3D.h"
#include "CGL/matrix3x3.h"
#include "CGL/lodepng.h"

#include "GL/glew.h"

#include "static_scene/sphere.h"
#include "static_scene/triangle.h"
#include "static_scene/light.h"
#include "static_scene/object.h"

#define RR 0.7
#define COLLECTOR_NUM_SAMPLES 16
#define COLLECTOR_ITER_LIMIT 10
#define DELTA_THRESHOLD 1.0
#define MAX_SCATTER 20

using namespace CGL::StaticScene;

using std::min;
using std::max;

namespace CGL {

PathTracer::PathTracer(size_t ns_aa,
                       size_t max_ray_depth,
                       size_t ns_area_light,
                       size_t ns_diff,
                       size_t ns_glsy,
                       size_t ns_refr,
                       size_t num_threads,
                       size_t samples_per_batch,
                       float max_tolerance,
                       HDRImageBuffer* envmap,
                       bool direct_hemisphere_sample,
                       size_t mode,
                       string filename,
                       double lensRadius,
                       double focalDistance,
                       double deltaCeiling) : hypertex(PolarTex(42, 21, deltaCeiling, gen)) {
  state = INIT,
  this->ns_aa = ns_aa;
  this->max_ray_depth = max_ray_depth;
  this->ns_area_light = ns_area_light;
  this->ns_diff = ns_diff;
  this->ns_glsy = ns_diff;
  this->ns_refr = ns_refr;
  this->samplesPerBatch = samples_per_batch;
  this->maxTolerance = max_tolerance;
  this->lensRadius = lensRadius;
  this->focalDistance = focalDistance;
  this->direct_hemisphere_sample = direct_hemisphere_sample;
  this->mode = mode;
  this->filename = filename;
  this->deltaCeiling = deltaCeiling;

  if (envmap) {
    this->envLight = new EnvironmentLight(envmap);
  } else {
    this->envLight = NULL;
  }

  bvh = NULL;
  scene = NULL;
  camera = NULL;

  gridSampler = new UniformGridSampler2D();
  hemisphereSampler = new UniformHemisphereSampler3D();

  show_rays = true;

  imageTileSize = 32;
  numWorkerThreads = num_threads;
  workerThreads.resize(numWorkerThreads);

  tm_gamma = 2.2f;
  tm_level = 1.0f;
  tm_key = 0.18;
  tm_wht = 5.0f;

}

PathTracer::~PathTracer() {

  delete bvh;
  delete gridSampler;
  delete hemisphereSampler;

}

void PathTracer::set_scene(Scene *scene) {

  if (state != INIT) {
    return;
  }

  if (this->scene != nullptr) {
    delete scene;
    delete bvh;
    selectionHistory.pop();
  }

  if (this->envLight != nullptr) {
    scene->lights.push_back(this->envLight);
  }

  this->scene = scene;
  build_accel();

  if (has_valid_configuration()) {
    state = READY;
  }
}

void PathTracer::set_camera(Camera *camera) {

  if (state != INIT) {
    return;
  }

  this->camera = camera;

  this->camera->lensRadius = lensRadius;
  this->camera->focalDistance = focalDistance;
  
  if (has_valid_configuration()) {
    state = READY;
  }

}

void PathTracer::set_frame_size(size_t width, size_t height) {
  if (state != INIT && state != READY) {
    stop();
  }
  sampleBuffer.resize(width, height);
  frameBuffer.resize(width, height);
  cell_tl = Vector2D(0,0); 
  cell_br = Vector2D(width, height);
  render_cell = false;
  sampleCountBuffer.resize(width * height);
  if (has_valid_configuration()) {
    state = READY;
  }
}

bool PathTracer::has_valid_configuration() {
  return scene && camera && gridSampler && hemisphereSampler &&
         (!sampleBuffer.is_empty());
}

void PathTracer::update_screen() {
  switch (state) {
    case INIT:
    case READY:
      break;
    case VISUALIZE:
      visualize_accel();
      break;
    case RENDERING:
      glDrawPixels(frameBuffer.w, frameBuffer.h, GL_RGBA,
                   GL_UNSIGNED_BYTE, &frameBuffer.data[0]);
      if (render_cell)
        visualize_cell();
      break;
    case DONE:
        //sampleBuffer.tonemap(frameBuffer, tm_gamma, tm_level, tm_key, tm_wht);
      glDrawPixels(frameBuffer.w, frameBuffer.h, GL_RGBA,
                   GL_UNSIGNED_BYTE, &frameBuffer.data[0]);
      if (render_cell)
        visualize_cell();
      break;
  }
}

void PathTracer::stop() {
  switch (state) {
    case INIT:
    case READY:
      break;
    case VISUALIZE:
      while (selectionHistory.size() > 1) {
        selectionHistory.pop();
      }
      state = READY;
      break;
    case RENDERING:
      continueRaytracing = false;
    case DONE:
      for (int i=0; i<numWorkerThreads; i++) {
            workerThreads[i]->join();
            delete workerThreads[i];
        }
      state = READY;
      break;
  }
  render_silent = false;
}

void PathTracer::clear() {
  if (state != READY) return;
  delete bvh;
  bvh = NULL;
  scene = NULL;
  camera = NULL;
  selectionHistory.pop();
  sampleBuffer.resize(0, 0);
  frameBuffer.resize(0, 0);
  state = INIT;
  render_cell = false;
}

void PathTracer::start_visualizing() {
  if (state != READY) {
    return;
  }
  state = VISUALIZE;
}

void PathTracer::start_raytracing() {
  if (state != READY) return;

  // Intersection isect;
  // Ray r = camera->center_ray();
  // if (camera->lens_ind >= 0&& bvh->intersect(r, &isect)) {
  //   camera->focus_at(isect.t);
  // }

  rayLog.clear();
  workQueue.clear();

  state = RENDERING;
  continueRaytracing = true;
  workerDoneCount = 0;

  sampleBuffer.clear();
  if (!render_cell) {
    frameBuffer.clear();
    num_tiles_w = sampleBuffer.w / imageTileSize + 1;
    num_tiles_h = sampleBuffer.h / imageTileSize + 1;
    tilesTotal = num_tiles_w * num_tiles_h;
    tilesDone = 0;
    tile_samples.resize(num_tiles_w * num_tiles_h);
    memset(&tile_samples[0], 0, num_tiles_w * num_tiles_h * sizeof(int));

    // populate the tile work queue
    for (size_t y = 0; y < sampleBuffer.h; y += imageTileSize) {
        for (size_t x = 0; x < sampleBuffer.w; x += imageTileSize) {
            workQueue.put_work(WorkItem(x, y, imageTileSize, imageTileSize));
        }
    }
  } else {
    int w = (cell_br-cell_tl).x;
    int h = (cell_br-cell_tl).y;
    int imTS = imageTileSize / 4;
    num_tiles_w = w / imTS + 1;
    num_tiles_h = h / imTS + 1;
    tilesTotal = num_tiles_w * num_tiles_h;
    tilesDone = 0;
    tile_samples.resize(num_tiles_w * num_tiles_h);
    memset(&tile_samples[0], 0, num_tiles_w * num_tiles_h * sizeof(int));

    // populate the tile work queue
    for (size_t y = cell_tl.y; y < cell_br.y; y += imTS) {
      for (size_t x = cell_tl.x; x < cell_br.x; x += imTS) {
        workQueue.put_work(WorkItem(x, y, 
          min(imTS, (int)(cell_br.x-x)), min(imTS, (int)(cell_br.y-y)) ));
      }
    }
  }

  bvh->total_isects = 0; bvh->total_rays = 0;
  // launch threads
  fprintf(stdout, "[PathTracer] Rendering... "); fflush(stdout);
  for (int i=0; i<numWorkerThreads; i++) {
      workerThreads[i] = new std::thread(&PathTracer::worker_thread, this);
  }
}

void PathTracer::render_to_file(string filename, size_t x, size_t y, size_t dx, size_t dy) {
  if (x == -1) {
    unique_lock<std::mutex> lk(m_done);
    start_raytracing();
    cv_done.wait(lk, [this]{ return state == DONE; });
    lk.unlock();
    save_image(filename);
    fprintf(stdout, "[PathTracer] Job completed.\n");
  } else {
    render_cell = true;
    cell_tl = Vector2D(x,y);
    cell_br = Vector2D(x+dx,y+dy);
    ImageBuffer buffer;
    raytrace_cell(buffer);
    save_image(filename, &buffer);
    fprintf(stdout, "[PathTracer] Cell job completed.\n");
  }
}


void PathTracer::build_accel() {

  // collect primitives //
  fprintf(stdout, "[PathTracer] Collecting primitives... "); fflush(stdout);
  timer.start();
  vector<Primitive *> primitives;
  for (SceneObject *obj : scene->objects) {
    const vector<Primitive *> &obj_prims = obj->get_primitives();
    primitives.reserve(primitives.size() + obj_prims.size());
    primitives.insert(primitives.end(), obj_prims.begin(), obj_prims.end());
  }
  timer.stop();
  fprintf(stdout, "Done! (%.4f sec)\n", timer.duration());

  // build BVH //
  fprintf(stdout, "[PathTracer] Building BVH from %lu primitives... ", primitives.size()); 
  fflush(stdout);
  timer.start();
  bvh = new BVHAccel(primitives);
  timer.stop();
  fprintf(stdout, "Done! (%.4f sec)\n", timer.duration());

  // initial visualization //
  selectionHistory.push(bvh->get_root());
}

void PathTracer::visualize_accel() const {

  glPushAttrib(GL_ENABLE_BIT);
  glDisable(GL_LIGHTING);
  glLineWidth(1);
  glEnable(GL_DEPTH_TEST);

  // hardcoded color settings
  Color cnode = Color(.5, .5, .5); float cnode_alpha = 0.25f;
  Color cnode_hl = Color(1., .25, .0); float cnode_hl_alpha = 0.6f;
  Color cnode_hl_child = Color(1., 1., 1.); float cnode_hl_child_alpha = 0.6f;

  Color cprim_hl_left = Color(.6, .6, 1.); float cprim_hl_left_alpha = 1.f;
  Color cprim_hl_right = Color(.8, .8, 1.); float cprim_hl_right_alpha = 1.f;
  Color cprim_hl_edges = Color(0., 0., 0.); float cprim_hl_edges_alpha = 0.5f;

  BVHNode *selected = selectionHistory.top();

  // render solid geometry (with depth offset)
  glPolygonOffset(1.0, 1.0);
  glEnable(GL_POLYGON_OFFSET_FILL);

  if (selected->isLeaf()) {
    bvh->draw(selected, cprim_hl_left, cprim_hl_left_alpha);
  } else {
    bvh->draw(selected->l, cprim_hl_left, cprim_hl_left_alpha);
    bvh->draw(selected->r, cprim_hl_right, cprim_hl_right_alpha);
  }

  glDisable(GL_POLYGON_OFFSET_FILL);

  // draw geometry outline
  bvh->drawOutline(selected, cprim_hl_edges, cprim_hl_edges_alpha);

  // keep depth buffer check enabled so that mesh occluded bboxes, but
  // disable depth write so that bboxes don't occlude each other.
  glDepthMask(GL_FALSE);

  // create traversal stack
  stack<BVHNode *> tstack;

  // push initial traversal data
  tstack.push(bvh->get_root());

  // draw all BVH bboxes with non-highlighted color
  while (!tstack.empty()) {

    BVHNode *current = tstack.top();
    tstack.pop();

    current->bb.draw(cnode, cnode_alpha);
    if (current->l) tstack.push(current->l);
    if (current->r) tstack.push(current->r);
  }

  // draw selected node bbox and primitives
  if (selected->l) selected->l->bb.draw(cnode_hl_child, cnode_hl_child_alpha);
  if (selected->r) selected->r->bb.draw(cnode_hl_child, cnode_hl_child_alpha);

  glLineWidth(3.f);
  selected->bb.draw(cnode_hl, cnode_hl_alpha);

  // now perform visualization of the rays
  if (show_rays) {
      glLineWidth(1.f);
      glBegin(GL_LINES);

      for (size_t i=0; i<rayLog.size(); i+=500) {

          const static double VERY_LONG = 10e4;
          double ray_t = VERY_LONG;

          // color rays that are hits yellow
          // and rays this miss all geometry red
          if (rayLog[i].hit_t >= 0.0) {
              ray_t = rayLog[i].hit_t;
              glColor4f(1.f, 1.f, 0.f, 0.1f);
          } else {
              glColor4f(1.f, 0.f, 0.f, 0.1f);
          }

          Vector3D end = rayLog[i].o + ray_t * rayLog[i].d;

          glVertex3f(rayLog[i].o[0], rayLog[i].o[1], rayLog[i].o[2]);
          glVertex3f(end[0], end[1], end[2]);
      }
      glEnd();
  }

  glDepthMask(GL_TRUE);
  glPopAttrib();
}

void PathTracer::visualize_cell() const {
  glPushAttrib(GL_VIEWPORT_BIT);
  glViewport(0, 0, sampleBuffer.w, sampleBuffer.h);

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(0, sampleBuffer.w, sampleBuffer.h, 0, 0, 1);

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glTranslatef(0, 0, -1);

  glColor4f(1.0, 0.0, 0.0, 0.8);
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);

  // Draw the Red Rectangle.
  glBegin(GL_LINE_LOOP);
  glVertex2f(cell_tl.x, sampleBuffer.h-cell_br.y);
  glVertex2f(cell_br.x, sampleBuffer.h-cell_br.y);
  glVertex2f(cell_br.x, sampleBuffer.h-cell_tl.y);
  glVertex2f(cell_tl.x, sampleBuffer.h-cell_tl.y);
  glEnd();

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();

  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();

  glPopAttrib();

  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
}

void PathTracer::key_press(int key) {

  BVHNode *current = selectionHistory.top();
  switch (key) {
  case ']':
      ns_aa *=2;
      fprintf(stdout, "[PathTracer] Samples per pixel changed to %lu\n", ns_aa);
      //tm_key = clamp(tm_key + 0.02f, 0.0f, 1.0f);
      break;
  case '[':
      //tm_key = clamp(tm_key - 0.02f, 0.0f, 1.0f);
      ns_aa /=2;
      if (ns_aa < 1) ns_aa = 1;
      fprintf(stdout, "[PathTracer] Samples per pixel changed to %lu\n", ns_aa);
      break;
  case '=': case '+':
      ns_area_light *= 2;
      fprintf(stdout, "[PathTracer] Area light sample count increased to %zu.\n", ns_area_light);
      break;
  case '-': case '_':
      if (ns_area_light > 1) ns_area_light /= 2;
      fprintf(stdout, "[PathTracer] Area light sample count decreased to %zu.\n", ns_area_light);
      break;
  case '.': case '>':
      max_ray_depth++;
      fprintf(stdout, "[PathTracer] Max ray depth increased to %zu.\n", max_ray_depth);
      break;
  case ',': case '<':
      if (max_ray_depth) max_ray_depth--;
      fprintf(stdout, "[PathTracer] Max ray depth decreased to %zu.\n", max_ray_depth);
      break;
  case ';': case ':':
      focalDistance += .1;
      camera->focalDistance = focalDistance;
      fprintf(stdout, "[PathTracer] Focal distance increased to %f.\n", camera->focalDistance);
      break;
  case '\'': case '\"':
      focalDistance -= .1;
      camera->focalDistance = focalDistance;
      fprintf(stdout, "[PathTracer] Focal distance decreased to %f.\n", camera->focalDistance);
      break;
  case 'k': case 'K':
      if (lensRadius == 0)
        lensRadius = .03125f;
      else
        lensRadius *= sqrt(2.);
      camera->lensRadius = lensRadius;
      fprintf(stdout, "[PathTracer] Aperture increased to %f.\n", camera->lensRadius);
      break;
  case 'l': case 'L':
      if (lensRadius <= .03125f)
        lensRadius = 0.;
      else
        lensRadius /= sqrt(2.);
      camera->lensRadius = lensRadius;
      fprintf(stdout, "[PathTracer] Aperture decreased to %f.\n", camera->lensRadius);
      break;
  case KEYBOARD_UP:
      if (current != bvh->get_root()) {
          selectionHistory.pop();
      }
      break;
  case KEYBOARD_LEFT:
      if (current->l) {
          selectionHistory.push(current->l);
      }
      break;
  case KEYBOARD_RIGHT:
      if (current->l) {
          selectionHistory.push(current->r);
      }
      break;

  case 'C':
    render_cell = !render_cell;
    if (render_cell)
      fprintf(stdout, "[PathTracer] Now in cell render mode.\n");
    else
      fprintf(stdout, "[PathTracer] No longer in cell render mode.\n");
  break;

  default:
      return;
  }
}

Spectrum PathTracer::estimate_direct_lighting_hemisphere(const Ray& r, const Intersection& isect) {
  // Estimate the lighting from this intersection coming directly from a light.
  // For this function, sample uniformly in a hemisphere. 

  // make a coordinate system for a hit point
  // with N aligned with the Z direction.
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  // w_out points towards the source of the ray (e.g.,
  // toward the camera if this is a primary ray)
  const Vector3D& hit_p = r.o + r.d * isect.t;
  const Vector3D& w_out = w2o * (-r.d);

  // This is the same number of total samples as estimate_direct_lighting_importance (outside of delta lights). 
  // We keep the same number of samples for clarity of comparison.
  int num_samples = scene->lights.size() * ns_area_light;
  Spectrum L_out;

  // TODO (Part 3.2): 
  // Write your sampling loop here
  // COMMENT OUT `normal_shading` IN `est_radiance_global_illumination` BEFORE YOU BEGIN
  float pdf = 1.0f / (2.0f * PI);
  for (int i = 0; i < num_samples; i++) {
    const Vector3D w_in = hemisphereSampler->get_sample();
    const Vector3D sample = o2w * w_in;
    const Ray sampleRay = Ray(hit_p + EPS_D * sample, sample);
    Intersection intersection = isect;
    if (bvh->intersect(sampleRay, &intersection)) {
      float cosThetai = (w_in[2] / w_in.norm());
      L_out += isect.bsdf->f(w_in, w_out) * cosThetai * intersection.bsdf->get_emission() / (static_cast<float>(num_samples) * pdf);
    }
  }

  return L_out;

}

Spectrum PathTracer::estimate_direct_lighting_importance(const Ray& r, const Intersection& isect) {
  // Estimate the lighting from this intersection coming directly from a light.
  // To implement importance sampling, sample only from lights, not uniformly in a hemisphere. 

  // make a coordinate system for a hit point
  // with N aligned with the Z direction.
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  // w_out points towards the source of the ray (e.g.,
  // toward the camera if this is a primary ray)
  const Vector3D& hit_p = r.o + r.d * isect.t;
  const Vector3D& w_out = w2o * (-r.d);
  Spectrum L_out, rad_in;

  // TODO (Part 3.2): 
  // Here is where your code for looping over scene lights goes
  // COMMENT OUT `normal_shading` IN `est_radiance_global_illumination` BEFORE YOU BEGIN
  Vector3D wi = Vector3D(); Vector3D w_in = Vector3D();
  float distToLight, pdf;
  for (SceneLight *&light : scene->lights) {
    int sample_count = light->is_delta_light() ? 1 : ns_area_light;
    for (int i = 0; i < sample_count; i++) {
      pdf = direct_hemisphere_sample ? 0.0f : 1.0f;
      rad_in = light->sample_L(hit_p, &wi, &distToLight, &pdf);
      w_in = w2o * wi;
      if (w_in[2] >= 0) {
        const Ray sampleRay = Ray(hit_p + EPS_D * wi, wi, distToLight);
        Intersection intersection;
        if (!bvh->intersect(sampleRay, &intersection)) {
          float cosThetai = (fabs(w_in[2]) / w_in.norm());
          if (pdf != 0.0f) {
            L_out += isect.bsdf->f(w_in, w_out) * cosThetai * rad_in / (static_cast<float>(sample_count) * pdf);
          } else {
            i--;
          }
        }
      }
    }
  }

  return L_out;

}

Spectrum PathTracer::zero_bounce_radiance(const Ray&r, const Intersection& isect) {

  // TODO (Part 4.2):
  // Returns the light that results from no bounces of light

  return isect.bsdf->get_emission();

}

Spectrum PathTracer::one_bounce_radiance(const Ray&r, const Intersection& isect) {
  
  // TODO (Part 4.2):
  // Returns either the direct illumination by hemisphere or importance sampling
  // depending on `direct_hemisphere_sample`
  // (you implemented these functions in Part 3)

  return estimate_direct_lighting_importance(r, isect);
  
}

Spectrum PathTracer::at_least_one_bounce_radiance(const Ray&r, const Intersection& isect) {
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  // if (r.depth == max_ray_depth) {
  //   cout << isect.t * r.d.z << endl;
  // }

  Vector3D hit_p = r.o + r.d * isect.t;
  Vector3D w_out = w2o * (-r.d);
  Spectrum L_out = Spectrum();

  if (r.depth == 0) {
    return L_out;
  }

  if (!isect.bsdf->is_delta() && (this->mode <= 3 || max_ray_depth - r.depth + 3 == this->mode)) {
    L_out = one_bounce_radiance(r, isect);
  }

  if (isect.bsdf->is_cloud()) {
    L_out += canonical_cloud_slab_radiance(r, isect);
  } else {

    // TODO (Part 4.2): 
    // Here is where your code for sampling the BSDF,
    // performing Russian roulette step, and returning a recursively 
    // traced ray (when applicable) goes

    if (r.depth > 1) {
      Vector3D w_in = Vector3D(); float pdf;
      Spectrum sample_bsdf = isect.bsdf->sample_f(w_out, &w_in, &pdf);
      Vector3D wi = o2w * w_in;
      if (coin_flip(RR)) {
        Intersection next_isect;
        const Ray sampleRay = Ray(hit_p + EPS_D * wi, wi, INF_D, r.depth-1);
        if (bvh->intersect(sampleRay, &next_isect)) {
          Spectrum rad_in = at_least_one_bounce_radiance(sampleRay, next_isect);
          if (isect.bsdf->is_delta())
            rad_in += zero_bounce_radiance(sampleRay, next_isect);
          if (pdf != 0.0f) {
            float cosThetai = (fabs(w_in[2]) / w_in.norm());
            L_out += sample_bsdf * cosThetai * rad_in / (RR * pdf);
          }
        }
      }
    }
  }

  // cout << "Depth " << r.depth << " L_out: " << L_out << endl;

  return L_out;
}

Spectrum PathTracer::est_radiance_global_illumination(const Ray &r) {
  Intersection isect;
  Spectrum L_out;

  // You will extend this in assignment 3-2. 
  // If no intersection occurs, we simply return black.
  // This changes if you implement hemispherical lighting for extra credit.

  if (!bvh->intersect(r, &isect)) 
    return envLight ? envLight->sample_dir(r) : L_out;

  switch (this->mode) {
    // This line returns a color depending only on the normal vector 
    // to the surface at the intersection point.
    case 0:
      return normal_shading(isect.n);
    // TODO (Part 3): Return the direct illumination.
    case 1:
      return estimate_direct_lighting_hemisphere(r, isect);
    case 2:
      return estimate_direct_lighting_importance(r, isect);
    // TODO (Part 4): Accumulate the "direct" and "indirect" 
    // parts of global illumination into L_out rather than just direct
    case 3:
      return zero_bounce_radiance(r, isect) + at_least_one_bounce_radiance(r, isect);
    default:
      return zero_bounce_radiance(r, isect) + at_least_one_bounce_radiance(r, isect);
  }
  return normal_shading(isect.n);
}



Spectrum PathTracer::basic_cloud_slab_radiance(const Ray &r, const StaticScene::Intersection& isect) {

  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  Vector3D hit_p = r.o + r.d * isect.t;
  Vector3D w_out = w2o * (-r.d);
  Spectrum L_out = Spectrum();

  if (r.depth == 0 || coin_flip(1.0 - RR)) {
    return L_out;
  }

  Ray marched = Ray(hit_p + EPS_D * r.d, r.d);
  double distance;
  double alpha = ray_march(marched, 0.1, distance);
  if (coin_flip(alpha)) {
    hit_p = r.o + r.d * isect.t + r.d.unit() * distance;
    if (r.depth > 1) {
      float pdf = 0.0f;
      Vector3D wi = Vector3D();
      Vector3D w_in = Vector3D();
      float distToLight = 0.0f;
      scene->lights[0]->sample_L(hit_p, &wi, &distToLight, &pdf);

      Matrix3x3 o2w_wi;
      make_coord_space(o2w_wi, wi);

      Matrix3x3 o2w_inv;
      make_coord_space(o2w_inv, r.d.unit());

      Collector c = Collector(hit_p, 1.0, o2w_wi, isect.primitive, gen);
      Collector c_transmit = Collector(hit_p, 1.0, o2w_inv, isect.primitive, gen);
      Intersection next_isect;
      double delta = 0.0; double delta_transmit = 0.0;

      Ray next_ray = Ray(hit_p + EPS_D * wi, wi);
      bool intersects = bvh->intersect(next_ray, &next_isect);
      if (intersects && c.project(next_ray, next_isect, delta)) {

        Slab slab = c.generate_slab(c.mean);
        Collector c_next = slab.basic_sample();
        scene->lights[0]->sample_L(c_next.mean, &wi, &distToLight, &pdf);
        next_ray = Ray(c_next.mean + EPS_D * wi, wi);
        Ray next_inv_ray = Ray(c_next.mean - EPS_D * wi, -wi);

        int count = 0;
        while ((!bvh->intersect(next_ray, &next_isect) || !c.project(next_ray, next_isect, delta)) 
          && (!bvh->intersect(next_inv_ray, &next_isect) || !c.project(next_ray, next_isect, delta))) {

          if (count >= COLLECTOR_ITER_LIMIT) {
            c_next = c;
          } else {
            c_next = slab.basic_sample();
          }
          scene->lights[0]->sample_L(c_next.mean, &wi, &distToLight, &pdf);
          next_ray = Ray(c_next.mean + EPS_D * wi, wi);
          next_inv_ray = Ray(c_next.mean - EPS_D * wi, -wi);
          if (count >= COLLECTOR_ITER_LIMIT) {
            break;
          } else {
            count++;
          }
        }

        Spectrum L_direct = one_bounce_radiance(next_ray, next_isect);
        Spectrum L_collect = Spectrum();

        for (int i = 0; i < COLLECTOR_NUM_SAMPLES; i++) {
          
          Vector3D o_next = c_next.sample();
          Vector3D d_prev = (o_next - r.o).unit();
          Spectrum sample_bsdf = next_isect.bsdf->sample_f(-d_prev, &w_in, &pdf);
          Vector3D d_next = c_next.origin2w.T() * w_in;
          next_ray = Ray(o_next + EPS_D * d_next, d_next);

          if (bvh->intersect(next_ray, &next_isect)) {
            if (c_next.same_mesh(next_isect.primitive)) {
              next_ray.o += (next_isect.t + EPS_D) * next_ray.d;
              if (bvh->intersect(next_ray, &next_isect)) {
                next_ray.depth = r.depth-1;
                Spectrum rad_in = at_least_one_bounce_radiance(next_ray, next_isect);
                if (next_isect.bsdf->is_delta())
                  rad_in += zero_bounce_radiance(next_ray, next_isect);
                if (pdf != 0.0f) {
                  float cosThetai = (fabs(w_in[2]) / w_in.norm());
                  L_collect += hg_phase(w_out, d_prev) * sample_bsdf * cosThetai * rad_in / (pdf);
                  // std::cout << sample_bsdf << std::endl;
                }
              }
            } else {
              next_ray.depth = r.depth-1;
              Spectrum rad_in = at_least_one_bounce_radiance(next_ray, next_isect);
              if (next_isect.bsdf->is_delta())
                rad_in += zero_bounce_radiance(next_ray, next_isect);
              if (pdf != 0.0f) {
                float cosThetai = (fabs(w_in[2]) / w_in.norm());
                L_collect += hg_phase(w_out, d_prev) * sample_bsdf * cosThetai * rad_in / (pdf);
                // std::cout << sample_bsdf << std::endl;
              }
            }
          }
          
        }
        L_out += (L_collect / COLLECTOR_NUM_SAMPLES + L_direct) / RR;
      }
    }
  } else {      
    Spectrum L_transmit = Spectrum();
    Ray transmitted = Ray(hit_p + EPS_D * r.d, r.d);
    StaticScene::Intersection next_isect;
    if (bvh->intersect(transmitted, &next_isect)) {
      transmitted.o += r.d * (next_isect.t + EPS_D);
      transmitted.depth = r.depth;

      if (bvh->intersect(transmitted, &next_isect)) {
        L_transmit = at_least_one_bounce_radiance(transmitted, next_isect);
      } else {
        transmitted.o -= r.d * (next_isect.t + EPS_D);
        L_transmit = at_least_one_bounce_radiance(transmitted, next_isect);
      }
    }
    L_out += L_transmit / RR;
  }

  return L_out;
}

Spectrum PathTracer::canonical_cloud_slab_radiance(const Ray &r, const StaticScene::Intersection& isect) {

  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  Vector3D hit_p = r.o + r.d * isect.t;
  Vector3D w_out = w2o * (-r.d);
  Spectrum L_out = Spectrum();
  const Ray exit = Ray(hit_p, r.d);
  BSDF *cloud = isect.bsdf;

  if (r.depth == 0 || coin_flip(1.0 - RR)) {
    return L_out;
  }

  Ray marched = Ray(hit_p + EPS_D * r.d, r.d);
  double distance;
  double alpha = ray_march(marched, 0.1, distance);
  if (coin_flip(alpha)) {
    hit_p = r.o + r.d * isect.t + r.d.unit() * distance;
    if (r.depth > 1) {
      float pdf = 0.0f;
      Vector3D wi = Vector3D();
      Vector3D w_in = Vector3D();
      float distToLight = 0.0f;
      scene->lights[0]->sample_L(hit_p, &wi, &distToLight, &pdf);

      Matrix3x3 o2w_wi;
      make_coord_space(o2w_wi, wi);

      Matrix3x3 o2w_inv;
      make_coord_space(o2w_inv, r.d.unit());

      Collector c = Collector(hit_p, 1.0, o2w_wi, isect.primitive, gen);
      Collector c_transmit = Collector(hit_p, 1.0, o2w_inv, isect.primitive, gen);
      Intersection next_isect;
      double delta = 0.0; double delta_transmit = 0.0;

      Ray next_ray = Ray(hit_p + EPS_D * wi, wi);
      Ray next_inv_ray = Ray(hit_p - EPS_D * wi, -wi);
      bool intersects = bvh->intersect(next_ray, &next_isect);
      if (intersects && c.project(next_ray, next_isect, delta)) {

        std::vector<double> L;
        Collector argmin = c;
        std::vector<Collector> collectors;
        double deltaMin = 1000000.0;
        for (int i = 0; i < COLLECTOR_ITER_LIMIT && deltaMin <= DELTA_THRESHOLD; i++) {

          Slab slab = argmin.generate_slab(argmin.mean);
          collectors = slab.transport_sample(exit, wi, L, cloud, std::min(static_cast<int>(r.depth-1), MAX_SCATTER));

          for (Collector c_order : collectors) {
            // scene->lights[0]->sample_L(c_order.mean, &wi, &distToLight, &pdf);
            next_ray = Ray(c_order.mean + EPS_D * wi, wi);
            next_inv_ray = Ray(c_order.mean - EPS_D * wi, -wi);
            if ((bvh->intersect(next_ray, &next_isect) && c.project(next_ray, next_isect, delta)) ||
              (bvh->intersect(next_inv_ray, &next_isect) && c.project(next_inv_ray, next_isect, delta))) {
              if (argmin == c || deltaMin > delta) {
                deltaMin = delta;
                argmin = c_order;
              }
            }
          }
        }

        Spectrum L_direct = one_bounce_radiance(next_ray, next_isect);
        Spectrum L_collect = Spectrum();

        for (int o = 0; o < collectors.size(); o++) {
          for (int i = 0; i < COLLECTOR_NUM_SAMPLES; i++) {
            Vector3D o_next = collectors[o].sample();
            Vector3D d_prev = (o_next - r.o).unit();
            // cout << w_out << ", " << w_in << endl;
            Vector3D d_next = collectors[o].origin2w.T() * w_in;
            next_ray = Ray(o_next + EPS_D * d_next, d_next);

            if (bvh->intersect(next_ray, &next_isect)) {
              if (collectors[o].same_mesh(next_isect.primitive)) {
                next_ray.o += (next_isect.t + EPS_D) * next_ray.d;
                if (bvh->intersect(next_ray, &next_isect)) {
                  next_ray.depth = r.depth-collectors[o].order;
                  Spectrum rad_in = at_least_one_bounce_radiance(next_ray, next_isect);
                  if (next_isect.bsdf->is_delta())
                    rad_in += zero_bounce_radiance(next_ray, next_isect);
                  if (pdf != 0.0f) {
                    float cosThetai = (fabs(w_in[2]) / w_in.norm());
                    L_collect += L[o] * cosThetai * rad_in / (pdf);
                    // std::cout << sample_bsdf << std::endl;
                  }
                }
              } else {
                next_ray.depth = r.depth-1;
                Spectrum rad_in = at_least_one_bounce_radiance(next_ray, next_isect);
                if (next_isect.bsdf->is_delta())
                  rad_in += zero_bounce_radiance(next_ray, next_isect);
                if (pdf != 0.0f) {
                  float cosThetai = (fabs(w_in[2]) / w_in.norm());
                  L_collect += L[o] * cosThetai * rad_in / (pdf);
                  // std::cout << sample_bsdf << std::endl;
                }
              }
            }
          }
        }

        // std::cout << L_out << " + " << L_collect / COLLECTOR_NUM_SAMPLES << " + " << L_direct << std::endl;
        // cout << "COLLECT: " << alpha << ", " << distance << ", " << (L_collect / COLLECTOR_NUM_SAMPLES + L_direct) / RR << endl;
        // cout << "DIRECT: " << L_direct << endl;
        L_out += (L_collect / COLLECTOR_NUM_SAMPLES + L_direct) / RR;
      }
    }
  } else {      
    Spectrum L_transmit = Spectrum();
    Ray transmitted = Ray(hit_p + EPS_D * r.d, r.d);
    StaticScene::Intersection next_isect;
    if (bvh->intersect(transmitted, &next_isect)) {
      transmitted.o += r.d * (next_isect.t + EPS_D);
      transmitted.depth = r.depth;

      if (bvh->intersect(transmitted, &next_isect)) {
        L_transmit = at_least_one_bounce_radiance(transmitted, next_isect);
      } else {
        transmitted.o -= r.d * (next_isect.t + EPS_D);
        L_transmit = at_least_one_bounce_radiance(transmitted, next_isect);
      }
    }
    // std::cout << L_transmit << std::endl;
    L_out += L_transmit / RR;
    // cout << "TRANSMIT: " << alpha << ", " << distance << ", " << L_transmit / RR << endl; 
  }

  return L_out;
}

double PathTracer::ray_march(const Ray &r, const double step_size, double &distance) {
  StaticScene::Intersection isect;
  Ray transmitted = r;
  double d = 0.0; double alpha = 0.0; double c = 4.0;
  bool locked = false;
  Vector3D centroid = Vector3D();
  if (bvh->intersect(transmitted, &isect)) {
    auto *sph = dynamic_cast<const StaticScene::Sphere *>(isect.primitive);
    auto *tri = dynamic_cast<const StaticScene::Triangle *>(isect.primitive);
    if (sph) {
      centroid = sph->object->o;
    } else if (tri) {
      auto *mesh = tri->mesh;
      Vector3D sum = 0;
      for (int i = 0; i < mesh->vertexI; i++) {
        sum += mesh->positions[i];
      }
      centroid = sum / mesh->vertexI;
    } else {
      cout << "An error occurred!" << endl;
      throw 1;
    }
    while (d <= isect.t * r.d.norm() && alpha < 1.0) {

      Vector3D pos = transmitted.o - centroid;
      double phi = atan(pos.y / pos.x);
      if (pos.x < 0) { // quadrants II, III
        phi += M_PI;
      } else if (phi < 0) { // quadrant IV
        phi += 2 * M_PI;
      }
      double theta = atan(std::sqrt(pos.x * pos.x + pos.y * pos.y) / pos.z);
      double depth = min(isect.t * r.d.norm() - d, d);
      const Ray snapshot = transmitted;
      depth = estimate_depth(depth, snapshot);

      double density = hypertex.sample(depth, theta, phi);
      // cout << depth << ", " << density << endl;
      double alpha_k = 1.0 - pow(1.0 - density, c*step_size);
      alpha += alpha_k * (1 - alpha);
      d += step_size;
      transmitted.o += step_size * r.d.unit();
      if (!locked && coin_flip(alpha_k)) {
        distance = d;
        locked = true;
      }
    }
  }
  // cout << "Result: " << clamp(alpha, 0.0, 1.0) << ", " << distance << endl;
  return clamp(alpha, 0.0, 1.0);
}

double PathTracer::estimate_depth(double max_depth, const Ray &r) {
  Matrix3x3 o2w = Matrix3x3();
  make_coord_space(o2w, r.d);
  Vector3D x = o2w[0];
  Vector3D y = o2w[1];
  StaticScene::Intersection isect;
  double init_depth = max_depth;

  Ray x_ray = Ray(r.o, x);
  x_ray.max_t = max_depth;
  if (bvh->intersect(x_ray, &isect)) {
    max_depth = isect.t * x_ray.d.norm();
  }

  Ray y_ray = Ray(r.o, y);
  y_ray.max_t = max_depth;
  if (bvh->intersect(y_ray, &isect)) {
    max_depth = isect.t * y_ray.d.norm();
  }

  Ray neg_x_ray = Ray(r.o, -x);
  neg_x_ray.max_t = max_depth;
  if (bvh->intersect(neg_x_ray, &isect)) {
    max_depth = isect.t * neg_x_ray.d.norm();
  }

  Ray neg_y_ray = Ray(r.o, -y);
  neg_y_ray.max_t = max_depth;
  if (bvh->intersect(neg_y_ray, &isect)) {
    max_depth = isect.t * neg_y_ray.d.norm();
  }

  return max_depth;
}

Spectrum PathTracer::raytrace_pixel(size_t x, size_t y, bool useThinLens) {

  // TODO (Part 1.1):
  // Make a loop that generates num_samples camera rays and traces them 
  // through the scene. Return the average Spectrum. 
  // You should call est_radiance_global_illumination in this function.

  // TODO (Part 5):
  // Modify your implementation to include adaptive sampling.
  // Use the command line parameters "samplesPerBatch" and "maxTolerance"

  int num_samples = ns_aa;            // total samples to evaluate
  Vector2D origin = Vector2D(x,y);    // bottom left corner of the pixel
  Spectrum spectrum = Spectrum();
  double dbx = (double) x;  double dby = (double) y;
  double w = (double) sampleBuffer.w;  double h = (double) sampleBuffer.h;
  float sum = 0.0; float sq_sum = 0.0;
  float fsamplesPerBatch = static_cast<float>(samplesPerBatch);
  int totalSamplesCounted = 0;
  Vector2D samplesForLens = gridSampler->get_sample();

  if (num_samples == 1) {
    Ray currRay = camera->generate_ray_for_thin_lens((dbx+0.5) / w, (dby+0.5) / h, 
      samplesForLens.x, samplesForLens.y * 2.0 * PI);
    currRay.depth = max_ray_depth;
    spectrum = est_radiance_global_illumination(currRay);
    totalSamplesCounted = 1;
  } else if (num_samples > 1) {
    Vector2D offset;
    for (int i = 0; i < num_samples; i++) {
      samplesForLens = gridSampler->get_sample();
      offset = gridSampler->get_sample();
      Ray currRay = camera->generate_ray_for_thin_lens((dbx+offset.x) / w, (dby+offset.y) / h, 
        samplesForLens.x, samplesForLens.y * 2.0 * PI);
      currRay.depth = max_ray_depth;
      Spectrum sample = est_radiance_global_illumination(currRay);
      spectrum += sample;
      float illum = sample.illum();
      sum += illum;
      sq_sum += illum*illum;
      totalSamplesCounted++;
      if ((i+1) % samplesPerBatch == 0) {
        float mu = sum / totalSamplesCounted;
        float var = (1.0f/(totalSamplesCounted - 1.0f)) * fabs(sq_sum - sum*mu);
        if (1.96f*sqrt(var/totalSamplesCounted) <= maxTolerance * mu) {
          break;
        }
      }
    }
  }

  sampleCountBuffer[x + y * frameBuffer.w] = totalSamplesCounted;
  return spectrum / totalSamplesCounted;

}


void PathTracer::raytrace_tile(int tile_x, int tile_y,
                               int tile_w, int tile_h) {

  size_t w = sampleBuffer.w;
  size_t h = sampleBuffer.h;

  size_t tile_start_x = tile_x;
  size_t tile_start_y = tile_y;

  size_t tile_end_x = std::min(tile_start_x + tile_w, w);
  size_t tile_end_y = std::min(tile_start_y + tile_h, h);

  size_t tile_idx_x = tile_x / imageTileSize;
  size_t tile_idx_y = tile_y / imageTileSize;
  size_t num_samples_tile = tile_samples[tile_idx_x + tile_idx_y * num_tiles_w];

  for (size_t y = tile_start_y; y < tile_end_y; y++) {
    if (!continueRaytracing) return;
    for (size_t x = tile_start_x; x < tile_end_x; x++) {
      // TODO: 4.0
      // Change from false to true to enable thin lens
      Spectrum s = raytrace_pixel(x, y, false);
      sampleBuffer.update_pixel(s, x, y);
    }
  }

  tile_samples[tile_idx_x + tile_idx_y * num_tiles_w] += 1;
  sampleBuffer.toColor(frameBuffer, tile_start_x, tile_start_y, tile_end_x, tile_end_y);
}

void PathTracer::raytrace_cell(ImageBuffer& buffer) {
  size_t tile_start_x = cell_tl.x;
  size_t tile_start_y = cell_tl.y;

  size_t tile_end_x = cell_br.x;
  size_t tile_end_y = cell_br.y;

  size_t w = tile_end_x - tile_start_x;
  size_t h = tile_end_y - tile_start_y;
  HDRImageBuffer sb(w, h);
  buffer.resize(w,h);

  stop();
  render_cell = true;
  {
    unique_lock<std::mutex> lk(m_done);
    start_raytracing();
    cv_done.wait(lk, [this]{ return state == DONE; });
    lk.unlock();
  }

  for (size_t y = tile_start_y; y < tile_end_y; y++) {
    for (size_t x = tile_start_x; x < tile_end_x; x++) {
        buffer.data[w*(y-tile_start_y)+(x-tile_start_x)] = frameBuffer.data[x+y*sampleBuffer.w];
    }
  }
}

void PathTracer::worker_thread() {

  Timer timer;
  timer.start();

  WorkItem work;
  while (continueRaytracing && workQueue.try_get_work(&work)) {
    raytrace_tile(work.tile_x, work.tile_y, work.tile_w, work.tile_h);
    { 
      lock_guard<std::mutex> lk(m_done);
      ++tilesDone;
      if (!render_silent)  cout << "\r[PathTracer] Rendering... " << int((double)tilesDone/tilesTotal * 100) << '%';
      cout.flush();
    }
  }

  workerDoneCount++;
  if (!continueRaytracing && workerDoneCount == numWorkerThreads) {
    timer.stop();
    if (!render_silent)  fprintf(stdout, "\n[PathTracer] Rendering canceled!\n");
    state = READY;
  }

  if (continueRaytracing && workerDoneCount == numWorkerThreads) {
    timer.stop();
    if (!render_silent)  fprintf(stdout, "\r[PathTracer] Rendering... 100%%! (%.4fs)\n", timer.duration());
    if (!render_silent)  fprintf(stdout, "[PathTracer] BVH traced %llu rays.\n", bvh->total_rays);
    if (!render_silent)  fprintf(stdout, "[PathTracer] Averaged %f intersection tests per ray.\n", (((double)bvh->total_isects)/bvh->total_rays));

    lock_guard<std::mutex> lk(m_done);
    state = DONE;
    cv_done.notify_one();
  }
}

void PathTracer::save_image(string filename, ImageBuffer* buffer) {

  if (state != DONE) return;

  if (!buffer)
    buffer = &frameBuffer;

  if (filename == "") {
    time_t rawtime;
    time (&rawtime);

    time_t t = time(nullptr);
    tm *lt = localtime(&t);
    stringstream ss;
    ss << this->filename << "_screenshot_" << lt->tm_mon+1 << "-" << lt->tm_mday << "_" 
      << lt->tm_hour << "-" << lt->tm_min << "-" << lt->tm_sec << ".png";
    filename = ss.str();  
  }

  uint32_t* frame = &buffer->data[0];
  size_t w = buffer->w;
  size_t h = buffer->h;
  uint32_t* frame_out = new uint32_t[w * h];
  for(size_t i = 0; i < h; ++i) {
    memcpy(frame_out + i * w, frame + (h - i - 1) * w, 4 * w);
  }
  
  for (size_t i = 0; i < w * h; ++i) {
    frame_out[i] |= 0xFF000000;
  }

  fprintf(stderr, "[PathTracer] Saving to file: %s... ", filename.c_str());
  lodepng::encode(filename, (unsigned char*) frame_out, w, h);
  fprintf(stderr, "Done!\n");
  
  delete[] frame_out;

  save_sampling_rate_image(filename);
}

void PathTracer::save_sampling_rate_image(string filename) {
  size_t w = frameBuffer.w;
  size_t h = frameBuffer.h;
  ImageBuffer outputBuffer(w, h);

  for (int x = 0; x < w; x++) {
      for (int y = 0; y < h; y++) {
          float samplingRate = sampleCountBuffer[y * w + x] * 1.0f / ns_aa;

          Color c;
          if (samplingRate <= 0.5) {
              float r = (0.5 - samplingRate) / 0.5;
              c = Color(0.0f, 0.0f, 1.0f) * r + Color(0.0f, 1.0f, 0.0f) * (1.0 - r);
          } else {
              float r = (1.0 - samplingRate) / 0.5;
              c = Color(0.0f, 1.0f, 0.0f) * r + Color(1.0f, 0.0f, 0.0f) * (1.0 - r);
          }
          outputBuffer.update_pixel(c, x, h - 1 - y);
      }
  }
  uint32_t* frame_out = new uint32_t[w * h];
  
  for (size_t i = 0; i < w * h; ++i) {
    uint32_t out_color_hex = 0;
    frame_out[i] = outputBuffer.data.data()[i];
    frame_out[i] |= 0xFF000000;
  }

  lodepng::encode(filename.substr(0,filename.size()-4) + "_rate.png", (unsigned char*) frame_out, w, h);
  
  delete[] frame_out;
}

}  // namespace CGL
