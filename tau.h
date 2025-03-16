#include <vector>
#include "input.h"
#include "geom.h"

double findSigma(const Problem& prob, int index, const Geometry& geom);

std::vector<std::vector<double>> createOpticalDepthMatrix(const Problem& prob, const Geometry& geom);