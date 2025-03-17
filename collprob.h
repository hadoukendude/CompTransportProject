#pragma once
#include <vector>
#include "input.h"
#include "geom.h"
#include <Eigen/Dense>

Eigen::MatrixXd createCollProbMatrix(const Problem& prob, const std::vector<std::vector<double>>& tau,
    const Geometry& geom);