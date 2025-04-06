#include "eigenvalue.h"
#include <Eigen/Eigenvalues>
#include <iostream>

void solveTransportEV(struct MatrixData& mat, Eigen::MatrixXd CPMat)
{
    
    
    Eigen::Index n{ mat.X.size() };
    mat.S.array() /= mat.X.array();
    mat.F.array() /= mat.X.array();
    double K0 = mat.k;

    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n,n);
    Eigen::MatrixXd S = CPMat * mat.S.asDiagonal();
    Eigen::MatrixXd F = CPMat * mat.F.asDiagonal();
    mat.H = (I-S);
    Eigen::MatrixXd A = (I - S).inverse() * F;

    //Eigen::EigenSolver<Eigen::MatrixXd> solver(A);
    //std::cout << "Largest Eigenvalue: " << solver.eigenvalues().real().maxCoeff() << std::endl; // Directly calculates eigenvalue, comment out for long runs

    Eigen::VectorXd R0{ mat.Flux.array() * mat.D.array() * mat.X.array() };
    Eigen::VectorXd R1{ A*R0 };

    double K1 = (R1.array() / mat.D.array() * mat.F.array()).sum() / (R0.array() / mat.D.array() * mat.F.array()).sum();
    std::cout << K1 << std::endl;
    R0 = R1;

    while (abs(K1 - K0) > 1e-6) // eigenvalue iteration to 0.1 pcm convergence
    {
        K0 = K1;
        R1 = A * R0;
        K1 = (R1.array() / mat.D.array() * mat.F.array()).sum() / (R0.array() / mat.D.array() * mat.F.array()).sum();
        std::cout << K1 << std::endl;
        R0 = R1;
    }

    mat.Flux = R1.array() / mat.X.array() / mat.D.array();
    mat.R = R1;
    mat.k = K1;
}