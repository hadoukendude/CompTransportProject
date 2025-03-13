#include <iostream>
#include <vector>
#include <cmath>
#include "input.h"

std::vector<std::vector<double>> createOpticalDepthMatrix(const Problem& prob, const std::vector<double>& intervals)
{
    size_t n{ intervals.size() };
    std::vector<std::vector<double>> opticalDepthMatrix{ std::vector<std::vector<double>>(n, std::vector<double>(n))};

    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < i; j++)
        {
            std::cout << j;
            double deltaX = intervals[i] - intervals[j];
            opticalDepthMatrix[i][j] = 1 * deltaX;
        }
    }
    return opticalDepthMatrix;
}