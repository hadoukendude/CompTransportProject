#include <iostream>
#include "input.h"
#include "geom.h"


std::vector<double> linspace(double start, double end, int num)
{
    std::vector<double> result{};
    double step = (end - start) / (num - 1);
    for (int i = 1; i < num; i++)
        result.push_back(start + i * step);
    return result;
}

std::vector<double> createIntervals(const Problem& prob)
{
    std::vector<double> mesh{ 0.0 };
    
    int nPoints{ 1 };
    //for (int i = 0; i < prob.dimensions.size(); i++)
    //{
    //    // UNCOMMENT AFTER TESTING
    //    //if (prob.sigma[i] == 0)
    //    //    nPoints = 2;
    //    //else
    //    //    nPoints = static_cast<int>(prob.sigma[i]*prob.dimensions[i]);
    //    nPoints = 1;
    //    //assume dimensions lists width of each region
    //    std::vector<double> newInterval = linspace( mesh.back(), mesh.back()+prob.dimensions[i], nPoints + 1);
    //    mesh.insert( mesh.end(), newInterval.begin(), newInterval.end() );
    //}
    for (auto e : prob.dimensions)
    {
        nPoints = 1;
        // assume dimensions lists cumulative distance
        std::vector<double> newInterval = linspace(mesh.back(), e, nPoints + 1);
        mesh.insert(mesh.end(), newInterval.begin(), newInterval.end());
    }

    return mesh;
}

std::vector<double> findMidpoints(const std::vector<double>& mesh)
{
    std::vector<double> midpoints{};
    for (int i = 0; i < mesh.size() - 1; i++)
        midpoints.push_back((mesh[i] + mesh[i + 1]) / 2.);
    return midpoints;
}