#include <cmath>
#include <vector>
#include "input.h"
#include "tau.h"
#include <Eigen/Dense>
//#include <boost/math/special_functions/expint.hpp>


double expInt(double x)
{
    return -std::expint(-x);
}

Eigen::MatrixXd createCollProbMatrix(const Problem& prob, const std::vector<std::vector<double>>& tau)
{
    size_t n{ tau.size() };
    Eigen::MatrixXd P(n, n);

    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            if (i == j)
                P(i,j) = 0.0;
            else if (i > j)
                P(i,j) = 0.0;
            else
                P(i,j) = 0.0;
        }
    }
    //std::cout << boost::math::expint(1.0);
    return P;
}