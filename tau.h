#include <vector>
#include "input.h"
#include "geom.h"
#include <Eigen/Dense>


std::vector<std::vector<double>> createOpticalDepthMatrix(const Problem& prob,
    const Geometry& geom, Eigen::VectorXd Sigma);