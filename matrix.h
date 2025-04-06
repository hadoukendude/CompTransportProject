#pragma once
#include <vector>
#include <Eigen/Dense>


struct MatrixData
{
    Eigen::MatrixXd A; // albedo matrix
    Eigen::MatrixXd P; // collision probability matrix (vol to vol)
    Eigen::MatrixXd P_sv; // collision probability matrix (vol to surf)
    Eigen::MatrixXd P_vs; // collision probability matrix (surf to vol)
    Eigen::MatrixXd P_ss; // collision probability matrix (surf to surf)

    Eigen::VectorXd R; // collision rate vector
    Eigen::VectorXd D; // interval width vector
    Eigen::VectorXd G; // source vector

    Eigen::VectorXd X; // total xs vector
    Eigen::VectorXd S; // scattering xs vector
    Eigen::VectorXd F; // fission xs vector

    Eigen::MatrixXd H; // transport matrix
    Eigen::VectorXd Flux; // final flux distribution

    double k{ 1.0 }; // eigenvalue
};

void initializeDelta(struct Problem& prob, struct Geometry& geom, struct MatrixData& matrix_data);
void initializeXS(struct Problem& prob, struct Geometry& geom, struct MatrixData& matrix_data);