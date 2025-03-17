#include <cmath>
#include <vector>
#include "input.h"
#include "tau.h"
#include <Eigen/Dense>
#include <boost/math/special_functions/expint.hpp>


using namespace boost::math;
Eigen::MatrixXd createCollProbMatrix(const Problem& prob, const std::vector<std::vector<double>>& tau,
    const Geometry& geom)
{
    size_t n{ tau.size() };
    Eigen::MatrixXd P(n, n);
    double sigma_j{};
    double delta_j{};

    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            sigma_j = findSigma(prob, i, geom);
            sigma_j = (sigma_j == 0) ? 1e-8 : sigma_j;

            delta_j = geom.edges[j+1] - geom.edges[j];
            
            if (i == j)
            {
                std::cout << sigma_j << std::endl;
                std::cout << delta_j << std::endl;
                P(i, j) = 1.0 - 1.0 / (2 * sigma_j * delta_j) * (1.0 - 2.0 * expint(3, sigma_j * delta_j));
            }
            else if (i > j)
                P(i,j) = 0.0;
            else
                P(i,j) = 0.0;
        }
    }
    //std::cout << boost::math::expint(1.0);
    return P;
}