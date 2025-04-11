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
#include "timer.h"

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
    Timer timer;
    timer.start();

    std::cout << R"(
                   _  (`-')(`-')                   
         _         \-.(OO )(OO )_.->               
         \-,-----. _.'    \(_| \_)--. ,-.    ,-.   
          |  .--./(_...--''\  `.'  /,-| |-.,-| |-. 
         /_) (`-')|  |_.' | \    .')'-| |-''-| |-' 
         ||  |OO )|  .___.' .'    \   `-'    `-'   
        (_'  '--'\|  |     /  .'.  \               
           `-----'`--'    `--'   '--'                      
        )" << '\n';

    std::cout << "Reading input file..." << '\n';
    Problem problem{ parseInput("testproblem") };
    
    //problem.sig_a[2] *= 1.00; // absorption xs perturbation
    //problem.sigma[2] = problem.sig_a[2] + problem.sig_s[2];

    //std::cout << "Materials: " << problem.materials << "Cells: " 
    // << problem.cellMaterial << "Dimensions: " << problem.dimensions;

    std::cout << "Creating mesh..." << '\n';
    std::vector<double> interfaces{ createIntervals(problem) };
    std::vector<double> midpoints{ findMidpoints(interfaces) };
    struct Geometry geometry { midpoints, interfaces };
    //std::cout << "x_(i+1/2): " << '\n' << geometry.edges << "x_i: " << '\n' << geometry.indices << '\n';

    struct MatrixData matrix_data;
    initializeDelta(problem, geometry, matrix_data);
    //std::cout << "D_i: " << '\n' << matrix_data.D.adjoint() << '\n';

    initializeXS(problem, geometry, matrix_data);
    //std::cout << "total: " << '\n' << matrix_data.X.adjoint() << '\n' << "scattering: " << '\n' << matrix_data.S.adjoint() 
    //    << '\n' << "fission: " << '\n' << matrix_data.F.adjoint() << '\n';
    //if (problem.problemType == ProblemType::fixedsource)
        //std::cout << "source: " << '\n' << matrix_data.G.adjoint() << '\n';

    std::vector<std::vector<double>> opticalDepthMatrix{ createOpticalDepthMatrix(problem, geometry, matrix_data.X) };
    //std::cout << '\n' << "t_(ii'): " << '\n' << opticalDepthMatrix;

    std::cout << "Calculating collision probabilities..." << '\n';
    Eigen::MatrixXd collProbMat{ createCollProbMatrix(problem, opticalDepthMatrix, matrix_data) };
    //std::cout << "P_(ii'): " << '\n' << collProbMat << '\n'; // ------- TODO add reflective BC ---------

    //int a, b;
    //std::cin >> a, b;
    //std::cout << collProbMat(b, a) * matrix_data.X(a) * matrix_data.D(a) << std::endl;  // reciprocity
    //std::cout << collProbMat(a, b) * matrix_data.X(b) * matrix_data.D(b) << std::endl;
    //for (size_t i = 0; i < matrix_data.D.size(); i++)
    //    std::cout << collProbMat.col(i).sum() << '\t';  // herbert's conservation (=1, for specular, <1 for vacuum)
    { // from Stoke's theorem
        //std::cout << "Conservation 1: " << std::endl;
        //for (size_t i = 0; i < matrix_data.D.size(); i++)
        //{
        //    //std::cout << (matrix_data.D(i)*collProbMat.col(i)).sum() << '\t';  // herbert's conservation (=1, for specular, <1 for vacuum)
        //    //std::cout << 0.25 * (matrix_data.P_sv(0, i) + matrix_data.P_sv(1, i)) << std::endl;
        //    std::cout << ((matrix_data.D(i) * collProbMat.col(i)).sum() + 0.25 * (matrix_data.P_sv(0, i) + matrix_data.P_sv(1, i))) / matrix_data.D(i) << '\t';
        //}
        //std::cout << std::endl;
    }
    {
        //std::cout << "Conservation 2: " << std::endl;
        //std::cout << 0.25 * matrix_data.P_ss.col(0).sum() + (matrix_data.X.array() * matrix_data.D.array() * matrix_data.P_vs.col(0).array()).sum() << '\t'; // should be 1/4
        //std::cout << 0.25 * matrix_data.P_ss.col(1).sum() + (matrix_data.X.array() * matrix_data.D.array() * matrix_data.P_vs.col(1).array()).sum() << std::endl; // should also be 1/4 
    }


    if (problem.problemType == ProblemType::fixedsource)
        solveTransportFS(matrix_data, collProbMat);
    else
    {
        solveTransportEV(matrix_data, collProbMat);
    }

    /*
    // Eigenvalue perturbation
    //struct Problem pertprob = problem;
    //pertprob.sig_a[2] *= 1.00; // absorption xs perturbation
    //pertprob.sigma[2] = pertprob.sig_a[2] + pertprob.sig_s[2];

    //struct MatrixData pertdata {};
    //pertdata.D = matrix_data.D;
    //initializeXS(pertprob, geometry, pertdata);
    //std::vector<std::vector<double>> pertODM{ createOpticalDepthMatrix(pertprob, geometry, pertdata.X) };
    //Eigen::MatrixXd pertCP{ createCollProbMatrix(pertprob, pertODM, pertdata) };
    
    //pertdata.S.array() = pertdata.S.array() / pertdata.X.array();
    //pertdata.F.array() = pertdata.F.array() / pertdata.X.array();
    //pertdata.H = matrix_data.H.transpose();

    //Eigen::MatrixXd H = (collProbMat * matrix_data.S.asDiagonal()) - (pertCP * pertdata.S.asDiagonal());
    Eigen::MatrixXd F_0 = collProbMat * matrix_data.F.asDiagonal();
    //Eigen::MatrixXd F = (pertCP * pertdata.F.asDiagonal()) - F_0 ;
    Eigen::VectorXd R = matrix_data.R;

    // adjoint solver I guess
    std::cout << "Calculating adjoint flux..." << '\n';
    Eigen::MatrixXd Z{ matrix_data.H.transpose().inverse() * F_0.transpose() };
    //Eigen::VectorXd A0{ pertdata.Flux.array() * pertdata.D.array() * pertdata.X.array() };
    Eigen::VectorXd A0{ matrix_data.Flux.array() * matrix_data.D.array() * matrix_data.X.array() }; 
    A0.normalize();
    Eigen::VectorXd A1{ (Z - (matrix_data.k * Eigen::MatrixXd::Identity(Z.rows(), Z.cols()))).inverse() * A0 };

    //Eigen::EigenSolver<Eigen::MatrixXd> solver(Z);
    //std::cout << "Largest Eigenvalue: " << solver.eigenvalues().real().maxCoeff() << std::endl; // Directly calculates eigenvalue, comment out for long runs

    while ( abs(((Z * A0) - (matrix_data.k * A0)).maxCoeff()) > 1e-5)
    {
        A1 = (Z - (matrix_data.k * Eigen::MatrixXd::Identity(Z.rows(), Z.cols())) ).inverse() * A0; // rayleigh quotient iteration to find adjoint flux
        A0 = A1.normalized();
        //std::cout << ((Z * A0) - (matrix_data.k * A0)).maxCoeff() << std::endl;
    }
    A1 = A0.normalized();
    //std::cout << "adjie" << A1 << '\n';

    std::cout << "Adjoint relation: " << ((A1.transpose() * matrix_data.H * matrix_data.R) - (matrix_data.R.transpose() * matrix_data.H.transpose() * A1)) << std::endl;

    //auto lambda{ ((A1.transpose() * H * R) - (A1.transpose() * F * R / matrix_data.k)).array() / (A1.transpose() * F_0 * R).array()};
    auto lambda{ (problem.sig_a[2] * A1.transpose() * R).array() / (A1.transpose() * F_0 * R).array() };
    //std::cout << A1 << '\n';
    std::cout << "Eigenvalue perturbation: " << lambda << '\n';
    
    for (double i = 0; i <= 20; i++)
    {
        auto lambda{ ((i*0.02 - 0.2) * problem.sig_a[2] * A1.transpose() * R).array() / (A1.transpose() * F_0 * R).array() };
        std::cout << "Eigenvalue perturbation: " << lambda << '\n';
    }
    */

    timer.stop();
    std::cout << "Time elapsed: " << timer.elapsedSeconds() << " seconds" << '\n';

    plotter(matrix_data, problem.problemType);
    
    return 0;
}

// TODO ------------------------- create a universal n variable for the number of intervals so I don't
// keep saying Eigen::Index n{... blah blah blah