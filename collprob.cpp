#include <cmath>
#include <vector>
#include "input.h"
#include "tau.h"
#include "collprob.h"
#include <Eigen/Dense>
#include <boost/math/special_functions/expint.hpp>


using namespace boost::math;
double calculateCollProb(double sigma, double delta)
{
    return 0.0;
}

Eigen::MatrixXd createCollProbMatrix(const Problem& prob, const std::vector<std::vector<double>>& tau,
    const struct MatrixData& matData)
{
    size_t n{ tau.size() };
    Eigen::MatrixXd P{ Eigen::MatrixXd::Zero(n, n) };
    Eigen::MatrixXd P_diag{ Eigen::MatrixXd::Zero(n, n) };
    Eigen::MatrixXd P_upper{ Eigen::MatrixXd::Zero(n, n) };

    double sigma_i{};
    double sigma_j{};
    double delta_i{};
    double delta_j{};

    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            sigma_i = matData.X[i];
            sigma_j = matData.X[j];

            delta_i = matData.D[i];
            delta_j = matData.D[j];
            
            if (i == j)
            {
                P(i,j) = 1. - 1. / (2 * sigma_j * delta_j) * (1. - 2. * expint(3, sigma_j * delta_j));
            }
            else
                P(i,j) = 1. / (2. * sigma_j*delta_j) * (expint(3,tau[i][j]) - expint(3,tau[i][j] + sigma_i*delta_i)
                    - expint(3, tau[i][j] + sigma_j*delta_j) + expint(3, tau[i][j] + sigma_i*delta_i + sigma_j*delta_j));
        }
    }

    return P;
}