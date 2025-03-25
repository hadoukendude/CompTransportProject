// CompTransportProject.cpp : Defines the entry point for the application.

#include <iostream>
#include <vector>

#include "main.h"
#include "input.h"
#include "geom.h"
#include "tau.h"
#include "collprob.h"
#include "matrix.h"

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

    std::cout << "Materials: " << problem.materials;
    std::cout << "Cells: " << problem.cellMaterial;
    std::cout << "Dimensions: " << problem.dimensions;

    std::vector<double> interfaces{ createIntervals(problem) };
    std::vector<double> midpoints{ findMidpoints(interfaces) };
    struct Geometry geometry { midpoints, interfaces };

    std::cout << "x_(i+1/2): " << std::endl << geometry.edges;
    std::cout << "x_i: " << std::endl << geometry.indices << std::endl;

    
    struct MatrixData matrix_data;
    initializeDelta(problem, geometry, matrix_data);
    std::cout << "D_i: " << std::endl << matrix_data.D.adjoint() << std::endl;
    initializeXS(problem, geometry, matrix_data);
    std::cout << "total: " << std::endl << matrix_data.X.adjoint() << std::endl;
    std::cout << "scattering: " << std::endl << matrix_data.S.adjoint() << std::endl;
    std::cout << "fission: " << std::endl << matrix_data.F.adjoint() << std::endl;
    std::cout << "source: " << std::endl << matrix_data.G.adjoint() << std::endl;

    
    std::vector<std::vector<double>> opticalDepthMatrix{ createOpticalDepthMatrix(problem, geometry, matrix_data.X) };
    std::cout << std::endl << "t_(ii'): " << std::endl << opticalDepthMatrix; // TODO replace geom with matdata and remove findSigma

    Eigen::MatrixXd collProbMat{ createCollProbMatrix(problem, opticalDepthMatrix, matrix_data) };
    std::cout << "P_(ii'): " << std::endl << collProbMat << std::endl; // likely needs more testing

    initializeTransportMatrix(matrix_data, collProbMat);
    std::cout << std::endl << "H: " << std::endl << matrix_data.H; // fix narrowing conversion

    // schizo 3 am bs
    matrix_data.R = matrix_data.H.inverse() * matrix_data.G;
    std::cout << std::endl << std::endl << matrix_data.R;
    std::cout << std::endl << std::endl << matrix_data.R.array()/matrix_data.X.array()/matrix_data.D.array();
    return 0;
}
