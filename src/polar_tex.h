#ifndef CGL_POLAR_TEX_H
#define CGL_POLAR_TEX_H

#include <algorithm>
#include <random>

namespace CGL {

double sigmoid(double mean, double weight);

double clamp(double val, double lo, double hi);

class PolarTex {
  public:
  	PolarTex(size_t w, size_t h, double d, std::default_random_engine gen);
  	double sample(double d, double theta, double phi);

  private:
  	double depth;
  	size_t width;
  	size_t height;
  	std::vector<std::vector<double>> polarMap;
  	std::vector<std::vector<double>> weightMap;
};

}

#endif  // CGL_POLAR_TEX_H