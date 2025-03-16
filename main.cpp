// CompTransportProject.cpp : Defines the entry point for the application.

#include <iostream>
#include <vector>
#include "main.h"
#include "input.h"
#include "geom.h"
#include "tau.h"
#include "collprob.h"

#include <Eigen/Dense>

template <typename V>
std::ostream& operator<<(std::ostream& out, const std::vector<V>& v)
{
    for (auto e : v)
        out << e << " ";
    return out << std::endl;
}
template <typename W>
std::ostream& operator<<(std::ostream& out, const std::vector<std::vector<W>>& v)
{
    for (auto e : v)
    { 
        for (auto f : e)
        {
            out << f << '\t';
        }
        out << std::endl;
    }
    return out << std::endl;
}

int main()
{
    Problem problem{ parseInput("problem1") };

    std::cout << "sigma_total: " << problem.sigma << std::endl;

    std::cout << "Materials: " << problem.materials << std::endl;
    std::cout << "Cells: " << problem.cellMaterial << std::endl;
    std::cout << "Dimensions: " << problem.dimensions << std::endl;

    std::vector<double> interfaces{ createIntervals(problem) };
    std::vector<double> midpoints{ findMidpoints(interfaces) };
    
    struct Geometry geometry { midpoints, interfaces };

    std::cout << "x_(i+1/2): " << std::endl << geometry.edges;
    std::cout << "x_i: " << std::endl << geometry.indices;

    std::vector<std::vector<double>> opticalDepthMatrix{ createOpticalDepthMatrix(problem, geometry) };

    std::cout << "t_(ii'): " << std::endl << opticalDepthMatrix;

    std::cout << "P_(ii'): " << std::endl << createCollProbMatrix(problem, opticalDepthMatrix);

    return 0;
}
