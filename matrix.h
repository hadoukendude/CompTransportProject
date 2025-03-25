#pragma once
#include <vector>
#include <Eigen/Dense>


struct MatrixData
{
    Eigen::MatrixXd P; // collision probability matrix
    Eigen::VectorXd R; // collision rate
    Eigen::VectorXd D; // interval width vector
    Eigen::VectorXd G; // source vector

    Eigen::VectorXd X; // total xs vector
    Eigen::VectorXd S; // scattering xs vector
    Eigen::VectorXd F; // fission xs vector

    Eigen::MatrixXd H; // transport matrix
};

void initializeDelta(struct Problem& prob, struct Geometry& geom, struct MatrixData& matrix_data);
void initializeXS(struct Problem& prob, struct Geometry& geom, struct MatrixData& matrix_data);
void initializeTransportMatrix(struct MatrixData& matData, Eigen::MatrixXd collProbMat);