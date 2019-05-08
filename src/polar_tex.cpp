#include "polar_tex.h"

#include <random>
#include <algorithm>
#include <iostream>

namespace CGL {

double sigmoid(double mean, double curr, double weight) {
	double pow = weight * (mean - curr);
	return 1.0 / (1.0 + std::exp(pow));
}

double clamp(double val, double lo, double hi) {
	return std::max(std::min(val, hi), lo);
}

PolarTex::PolarTex(size_t w, size_t h, double d, std::default_random_engine gen) {
	depth = d;
	width = w;
	height = h;

	for (int i = 0; i < height; i++) {
		std::vector<double> row;
		std::vector<double> weightRow;
		double polarAngle = M_PI * (static_cast<double>(i) / static_cast<double>(height - 1) - 0.5);
		std::normal_distribution<double> mean(0.5, std::cos(polarAngle) / 5.0);
		for (int j = 0; j < width; j++) {
			if (j < width - 1) {
				std::normal_distribution<double> delta(clamp(mean(gen), 0.0, 1.0), std::cos(polarAngle));
				row.push_back(clamp(delta(gen), 0.0, 1.0));
				std::normal_distribution<double> weight(mean(gen), 2.0 * std::cos(polarAngle));
				weightRow.push_back(std::pow(weight(gen), 4.0));
			} else {
				row.push_back(row[0]);
				weightRow.push_back(weightRow[0]);
			}
		}
		polarMap.push_back(row);
		weightMap.push_back(weightRow);
	}
}

double PolarTex::sample(double d, double theta, double phi) {

	double curr = d / depth;
	// [0.0, 1.0] is the interesting interval for polar texture output
	// std::cout << theta << ", " << phi << std::endl;
	if (d < 0.0) {
		return 0.0;
	} else if (d > 1.0) {
		return 1.0;
	}

	// compute adjacent polar and azimuthal map indices
	double polarIndex = (theta + M_PI * 0.5) * static_cast<double>(height - 1) / M_PI;
	double azIndex = phi * static_cast<double>(width - 1) / (2.0 * M_PI);

	int polarFloor = (int) std::floor(polarIndex);	int polarCeil = (int) std::ceil(polarIndex);
	int azFloor = (int) std::floor(azIndex);	int azCeil = (int) std::ceil(azIndex);

	// Set up bilerp of sigmoids functions.
	double pAlpha00 = sigmoid(polarMap[polarFloor][azFloor], curr, weightMap[polarFloor][azFloor]);
	double pAlpha01 = sigmoid(polarMap[polarFloor][azCeil], curr, weightMap[polarFloor][azCeil]);
	double pAlpha10 = sigmoid(polarMap[polarCeil][azFloor], curr, weightMap[polarCeil][azFloor]);
	double pAlpha11 = sigmoid(polarMap[polarCeil][azCeil], curr, weightMap[polarCeil][azCeil]);
	// std::cout << pAlpha00 << ", " << pAlpha01 << ", " << pAlpha10 << ", " << pAlpha11 << std::endl;

	double lerp0 = (azCeil - azIndex) * pAlpha00 + (azIndex - azFloor) * pAlpha01;
	double lerp1 = (azCeil - azIndex) * pAlpha10 + (azIndex - azFloor) * pAlpha11;
	// std::cout << (polarCeil - polarIndex) * lerp0 + (polarIndex - polarFloor) * lerp1<< std::endl;
	return (polarCeil - polarIndex) * lerp0 + (polarIndex - polarFloor) * lerp1;
}

}