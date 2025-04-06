#include "fixedsource.h"

#include <iostream> // REMOVE LATER

// solves the system of equations for the fixed source problem
// H*r = s, r = inv(H)*s
void solveTransportFS(struct MatrixData& mat, Eigen::MatrixXd collProbMat)
{
    Eigen::VectorXd G{ mat.G.size() };
    G = mat.G.array() * mat.D.array();
    Eigen::VectorXd SF{ mat.S.size() };
    SF = mat.S + mat.F;

    Eigen::Index n{ collProbMat.rows() };
    mat.H = Eigen::MatrixXd::Zero(n, n);
    Eigen::MatrixXd I{ Eigen::MatrixXd::Identity(n, n) };

    mat.S.array() /= mat.X.array();
    mat.F.array() /= mat.X.array();
    mat.H = I - (collProbMat * (mat.S.asDiagonal() + mat.F.asDiagonal()));

    mat.G = mat.G.array() * mat.D.array(); // S_i = s_i * d_i
    mat.G = collProbMat * mat.G; // S_tilde = P_ii' * S_i
    mat.R = mat.H.inverse() * mat.G;
    mat.Flux = mat.R.array() / mat.X.array() / mat.D.array();

    /*
    {
        // lets see
        std::cout << mat.P << '\n' << '\n' << mat.P_ss << '\n' << '\n' << mat.P_vs << '\n' << '\n' << mat.P_sv << '\n' << '\n';
        //mat.A = Eigen::MatrixXd::Zero(2, 2);
        //mat.A(0, 0) = 1.;

        Eigen::MatrixXd Z{ Eigen::MatrixXd::Zero(n, n) };
        Z = collProbMat + (mat.P_vs * mat.A * (Eigen::MatrixXd::Identity(2, 2) - mat.P_ss * mat.A).inverse() * mat.P_sv);
        std::cout << "Z: " << '\n' << Z << std::endl;

        Eigen::MatrixXd Y{ Eigen::MatrixXd::Zero(n, n) };
        Y = I - (Z * (mat.S.asDiagonal() + mat.F.asDiagonal()));
        //Y = collProbMat * (mat.S.asDiagonal() + mat.F.asDiagonal()) 
        //    + (mat.P_vs * (Eigen::MatrixXd::Identity(2, 2) - mat.P_ss).inverse() * mat.P_sv) * SF.asDiagonal();
        Eigen::MatrixXd fux{ Eigen::MatrixXd::Zero(n, 1) };

        fux = Y.inverse() * Z * G;
        std::cout << fux << '\n';
        fux = fux.array() / mat.X.array() / mat.D.array();
        std::cout << fux << '\n';
    }
    */
}