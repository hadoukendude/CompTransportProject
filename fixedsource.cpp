#include "fixedsource.h"

// solves the system of equations for the fixed source problem
// H*r = s, r = inv(H)*s
void solveTransportFS(struct MatrixData& mat, Eigen::MatrixXd collProbMat)
{
    Eigen::Index n = collProbMat.rows();
    mat.H = Eigen::MatrixXd::Zero(n, n);
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n, n);

    mat.S.array() /= mat.X.array();
    mat.F.array() /= mat.X.array();
    mat.H = I - (collProbMat * (mat.S.asDiagonal() + mat.F.asDiagonal()));

    mat.G = mat.G.array() * mat.D.array(); // S_i = s_i * d_i
    mat.G = collProbMat * mat.G; // S_tilde = P_ii' * S_i
    mat.R = mat.H.inverse() * mat.G;
    mat.Flux = mat.R.array() / mat.X.array() / mat.D.array();
}