#include <iostream>
#include <vector>
#include "input.h"
#include "geom.h"
#include "tau.h"


double calculateOpticalDepth(const Problem& prob, int start, int end, const Geometry& geom, Eigen::VectorXd sigma_t)
{
    std::vector<int> mats{ prob.cellMaterial };
    double tau{ 0.0 };

    if (start == end)
        return tau; // No distance within the same cell
    else if (start+1 == end)
        return tau; // No distance between adjacent cells
    else if (start < end) // going forwards
    {
        double sigma{ 1.0 };
        for (size_t i = start+1; i < end; i++) // iterate 
        {
            tau += sigma_t[i]*(geom.edges[i+1] - geom.edges[i]);
        }
        return tau;
    }
    else // going backwards
        return -1;
}


std::vector<std::vector<double>> createOpticalDepthMatrix(const Problem& prob,
    const Geometry& geom, Eigen::VectorXd sigma_t)
{
    size_t n{ geom.indices.size() };
    std::vector<std::vector<double>> opticalDepthMatrix{
        std::vector<std::vector<double>>(n, std::vector<double>(n))};


    for (size_t i = 0; i < n; i++) // iterate over indices
        for (size_t j = 0; j < n; j++)
            opticalDepthMatrix[i][j] = calculateOpticalDepth(prob, i, j, geom, sigma_t);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < i; j++)
            opticalDepthMatrix[i][j] = opticalDepthMatrix[j][i];

    return opticalDepthMatrix;
}
