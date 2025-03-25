#include <vector>
#include "input.h"
#include "geom.h"
#include <Eigen/Dense>

double findSigma(const Problem& prob, int index, const Geometry& geom);

std::vector<std::vector<double>> createOpticalDepthMatrix(const Problem& prob,
    const Geometry& geom, Eigen::VectorXd Sigma);