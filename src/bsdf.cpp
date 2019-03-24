#include "bsdf.h"

#include <iostream>
#include <algorithm>
#include <utility>

#define EPS_B 1e-10;

using std::min;
using std::max;
using std::swap;

namespace CGL {

void make_coord_space(Matrix3x3& o2w, const Vector3D& n) {
  Vector3D z = Vector3D(n.x, n.y, n.z);
  Vector3D h = z;
  if (fabs(h.x) <= fabs(h.y) && fabs(h.x) <= fabs(h.z)) h.x = 1.0;
  else if (fabs(h.y) <= fabs(h.x) && fabs(h.y) <= fabs(h.z)) h.y = 1.0;
  else h.z = 1.0;

  z.normalize();
  Vector3D y = cross(h, z);
  y.normalize();
  Vector3D x = cross(z, y);
  x.normalize();

  o2w[0] = x;
  o2w[1] = y;
  o2w[2] = z;
}

// Diffuse BSDF //

Spectrum DiffuseBSDF::f(const Vector3D& wo, const Vector3D& wi) {

  // TODO (Part 3.1): 
  // This function takes in both wo and wi and returns the evaluation of
  // the BSDF for those two directions.

  return reflectance / PI;
}

Spectrum DiffuseBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO (Part 3.1): 
  // This function takes in only wo and provides pointers for wi and pdf,
  // which should be assigned by this function.
  // After sampling a value for wi, it returns the evaluation of the BSDF
  // at (wo, *wi).

  *wi = sampler.get_sample(pdf);
  return reflectance / PI;
}

// Mirror BSDF //

Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // set wi and pdf
  reflect(wo, wi);
  *pdf = 1.0f;
  return reflectance / abs_cos_theta(*wi);
}

// Microfacet BSDF //

double MicrofacetBSDF::G(const Vector3D& wo, const Vector3D& wi) {
    return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
}

double MicrofacetBSDF::D(const Vector3D& h) {
  // TODO: 2.2
  // Compute Beckmann normal distribution function (NDF) here.
  // You will need the roughness alpha.

  return std::pow(cos_theta(h), 100.0);;
}

Spectrum MicrofacetBSDF::F(const Vector3D& wi) {
  // TODO: 2.3
  // Compute Fresnel term for reflection on dielectric-conductor interface.
  // You will need both eta and etaK, both of which are Spectrum.

  return Spectrum();
}

Spectrum MicrofacetBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  // TODO: 2.1
  // Implement microfacet model here

  return Spectrum();
}

Spectrum MicrofacetBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO: 2.4
  // *Importance* sample Beckmann normal distribution function (NDF) here.
  // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
  //       and return the sampled BRDF value.

  *wi = cosineHemisphereSampler.get_sample(pdf); //placeholder
  return MicrofacetBSDF::f(wo, *wi);
}

// Refraction BSDF //

Spectrum RefractionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum RefractionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  return Spectrum();
}

// Glass BSDF //

Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO: 1.4
  // Compute Fresnel coefficient and either reflect or refract based on it.

  *pdf = 1.0;
  if (!refract(wo, wi, ior)) { // either possible refraction or total internal reflection
    reflect(wo, wi);
    return reflectance / abs_cos_theta(*wi);
  }
  double iord = wo[2] < 0 ? static_cast<float>(ior) : 1.0 / static_cast<float>(ior);
  double R0 = pow((iord-1.0)/(iord+1.0),2.0);
  double prob = R0 + (1.0-R0) * pow(1.0-abs_cos_theta(*wi),5.0);
  if (coin_flip(prob)) { // Schlick's approximation for reflection probability
    reflect(wo, wi);
    *pdf = prob;
    return prob * reflectance / abs_cos_theta(*wi);
  } else {
    *pdf = 1.0 - prob;
    Spectrum result = f(wo, *wi);
    // std::cout << result << ", " << *wi << " -> " << wo << std::endl;
    return (1.0 - prob) * transmittance / abs_cos_theta(*wi) / pow(iord, 2.0);
  }
}

void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {

  // TODO: 1.1
  // Implement reflection of wo about normal (0,0,1) and store result in wi.

  *wi = Vector3D(-wo[0],-wo[1],wo[2]);
}

bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) {

  // TODO: 1.3
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.

  double sthetao = sin_theta(wo);
  double iord = cos_theta(wo) < 0.0 ? static_cast<float>(ior) : 1.0 / static_cast<float>(ior);
  if (sthetao * iord > 1.0 || sthetao * iord < -1.0) { // internal reflection check
    return false;
  }
  double sthetai = sthetao * iord;
  *wi = Vector3D(-iord*wo[0],-iord*wo[1],sqrt(1-pow(iord, 2.0)*(1.0-pow(wo[2], 2.0)))*(wo[2]<0?1:-1));
  return true;
}

// Emission BSDF //

Spectrum EmissionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum EmissionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *pdf = 1.0 / PI;
  *wi  = sampler.get_sample(pdf);
  return Spectrum();
}

} // namespace CGL
