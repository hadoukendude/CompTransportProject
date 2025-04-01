// CompTransportProject.cpp : Defines the entry point for the application.

#include <iostream>
#include <vector>

#include "main.h"
#include "input.h"
#include "geom.h"
#include "tau.h"
#include "collprob.h"
#include "matrix.h"
#include "fixedsource.h"
#include "eigenvalue.h"
#include "plot.h"

#include <Eigen/Dense>

template <typename V>
std::ostream& operator<<(std::ostream& out, const std::vector<V>& v)
{
    for (auto e : v)
        out << e << " ";
    return out << '\n';
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
        out << '\n';
    }
    return out << '\n';
}

int main()
{
    Problem problem{ parseInput("problem1") };

    std::cout << "sigma_total: " << problem.sigma << '\n';

    for (size_t i = 0; i < problem.sigma.size(); i++) // absorption cross section perturbation
    {
        problem.sig_a[i] *= 1;
        problem.sigma[i] = problem.sig_a[i] + problem.sig_s[i];
    }

    std::cout << "Materials: " << problem.materials;
    std::cout << "Cells: " << problem.cellMaterial;
    std::cout << "Dimensions: " << problem.dimensions;

    std::vector<double> interfaces{ createIntervals(problem) };
    std::vector<double> midpoints{ findMidpoints(interfaces) };
    struct Geometry geometry { midpoints, interfaces };
    std::cout << "x_(i+1/2): " << '\n' << geometry.edges;
    std::cout << "x_i: " << '\n' << geometry.indices << '\n';

    struct MatrixData matrix_data;
    initializeDelta(problem, geometry, matrix_data);
    std::cout << "D_i: " << '\n' << matrix_data.D.adjoint() << '\n';

    initializeXS(problem, geometry, matrix_data);
    std::cout << "total: " << '\n' << matrix_data.X.adjoint() << '\n';
    std::cout << "scattering: " << '\n' << matrix_data.S.adjoint() << '\n';
    std::cout << "fission: " << '\n' << matrix_data.F.adjoint() << '\n';
    if (problem.problemType == ProblemType::fixedsource)
        std::cout << "source: " << '\n' << matrix_data.G.adjoint() << '\n';

    std::vector<std::vector<double>> opticalDepthMatrix{ createOpticalDepthMatrix(problem, geometry, matrix_data.X) };
    std::cout << '\n' << "t_(ii'): " << '\n' << opticalDepthMatrix;

    Eigen::MatrixXd collProbMat{ createCollProbMatrix(problem, opticalDepthMatrix, matrix_data) };
    std::cout << "P_(ii'): " << '\n' << collProbMat << '\n'; // ------- TODO add reflective BC ---------

    //int a, b;
    //std::cin >> a, b;
    //std::cout << collProbMat(b, a) * matrix_data.X(a) * matrix_data.D(a) << std::endl;  // reciprocity
    //std::cout << collProbMat(a, b) * matrix_data.X(b) * matrix_data.D(b) << std::endl;

    for (size_t i = 0; i < matrix_data.D.size() ; i++)
        std::cout << collProbMat.col(i).sum() << std::endl;  // herbert's conservation (=1, for specular, <1 for vacuum)

    if (problem.problemType == ProblemType::fixedsource)
        solveTransportFS(matrix_data, collProbMat);
    else
        solveTransportEV(matrix_data, collProbMat);

    //std::cout << '\n' << "H: " << '\n' << matrix_data.H;

    std::cout << matrix_data.Flux << '\n';
    //plotter(matrix_data);
    
    return 0;
}