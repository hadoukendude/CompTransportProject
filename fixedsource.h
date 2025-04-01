#pragma once
#include <vector>
#include <Eigen/Dense>
#include "matrix.h"


void solveTransportFS(struct MatrixData& matData, Eigen::MatrixXd collProbMat);