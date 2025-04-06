#include <cmath>
#include <vector>
#include "input.h"
#include "tau.h"
#include "collprob.h"
#include <Eigen/Dense>
#include <boost/math/special_functions/expint.hpp>


using namespace boost::math;

std::vector<std::vector<double>> calculateTauR(const struct MatrixData& matData)
{
    Eigen::Index n{ matData.D.size() };
    std::vector<std::vector<double>> tau_r{
        std::vector<std::vector<double>>(n, std::vector<double>(n)) };

    double tau{ 0.0 };
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            tau = 0.;
            for (size_t k = i + 1; k < n; k++)
            {
                tau += matData.X(k) * matData.D(k);
            }
            for (size_t k = j + 1; k < n; k++)
            {
                tau += matData.X(k) * matData.D(k);
            }
            tau_r[i][j] = tau;
        }
    }
    //std::cout << "t_R: " << '\n';
    //for (auto e : tau_r)
    //{
    //    for (auto f : e)
    //        std::cout << f << '\t';
    //    std::cout << '\n';
    //} // print tau_BC
    return tau_r;
}

std::vector<std::vector<double>> calculateTauBC(const struct MatrixData& matData)
{
    size_t n{static_cast<size_t> (matData.D.size()) };
    std::vector<std::vector<double>> tau_BC{
        std::vector<std::vector<double>>(n+1, std::vector<double>(2)) };

    tau_BC[0][0] = 0;
    for (size_t i = 1; i <= n; i++)
    {
        tau_BC[i][0] = matData.X(i-1)*matData.D(i-1) + tau_BC[i - 1][0];
    }

    double tau_cell{ (matData.D.array() * matData.X.array()).sum() }; // its over for taucells...
    tau_BC[0][1] = tau_cell;
    for (size_t i = 0; i < n; i++)
    {
        tau_BC[i+1][1] = tau_BC[i][1] - (matData.X(i)*matData.D(i));
    }

    //std::cout << "t_BC" << '\n';
    //for (auto e : tau_BC)
    //{
    //    for (auto f : e)
    //        std::cout << f << '\t';
    //    std::cout << '\n';
    //} // print tau_BC
    //std::cout << std::endl;
    return tau_BC; // TODO ------------ change name
}
void calcSurfToVol(struct MatrixData& matData, const std::vector<std::vector<double>>& tau)
{
    Eigen::Index n{ matData.D.size() };
    matData.P_vs = Eigen::MatrixXd::Zero(n, 2);

    Eigen::MatrixXd test = Eigen::MatrixXd::Zero(n, 1);

    for (size_t i = 0; i < n; i++)
    {
        if (matData.X(i) < 1e-20)
        {
            matData.P_vs(i, 0) = 1. / 2. * expint(3, tau[i][0]);
            matData.P_vs(i, 1) = 1. / 2. * expint(3, tau[i + 1][1]);
        }
        else
        {
            matData.P_vs(i, 0) = 1. / (2. * matData.X(i) * matData.D(i)) * (expint(3, tau[i][0]) - expint(3, tau[i + 1][0]));
            matData.P_vs(i, 1) = 1. / (2. * matData.X(i) * matData.D(i)) * (expint(3, tau[i + 1][1]) - expint(3, tau[i][1]));

            //matData.P_vs(i, 0) = 1. / (2. * matData.D(i)) * (expint(3, tau[i][0]) - expint(3, tau[i + 1][0]));
            //matData.P_vs(i, 1) = 1. / (2. * matData.D(i)) * (expint(3, tau[i + 1][1]) - expint(3, tau[i][1]));
        }
    }

    //std::cout << 1. / (2. * matData.X(0) * matData.D(0)) << '\n';
    //std::cout << 1. / (2. * matData.X(1) * matData.D(1)) << '\n';
    //std::cout << (expint(3, tau[0][0]) - expint(3, tau[1][0])) << '\n';
    //std::cout << 1. / (2. * matData.X(3) * matData.D(3)) << '\n';
    //std::cout << (expint(3, tau[2][0]) - expint(3, tau[3][0])) << '\n';
    //std::cout << test << '\n';
    //std::cout << matData.P_vs << '\n';
}
void calcVolToSurf(struct MatrixData& matData)
{
    Eigen::Index n{ matData.D.size() };
    matData.P_sv = Eigen::MatrixXd::Zero(2, n);
    for (size_t i = 0; i < n; i++)
    {
        matData.P_sv(0, i) = matData.P_vs(i, 0) * 4 * matData.D(i);
        matData.P_sv(1, i) = matData.P_vs(i, 1) * 4 * matData.D(i);
    }
    //std::cout << matData.P_sv << '\n';
}
void calcSurftoSurf(struct MatrixData& matData)
{
    double tau_LR{ (matData.D.array() * matData.X.array()).sum() }; // total optical depth of the problem
    
    matData.P_ss = Eigen::MatrixXd::Zero(2, 2);
    matData.P_ss(0, 1) = 2 * expint(3, tau_LR);
    matData.P_ss(1, 0) = 2 * expint(3, tau_LR);
    //std::cout << matData.P_ss; // test P_SS
}

Eigen::MatrixXd createCollProbMatrix(const Problem& prob, const std::vector<std::vector<double>>& tau,
    struct MatrixData& matData)
{
    size_t n{ tau.size() };
    Eigen::MatrixXd P{ Eigen::MatrixXd::Zero(n, n) };

    double sigma_i{};
    double sigma_j{};
    double delta_i{};
    double delta_j{};

    for (size_t i = 0; i < n; i++)
    {
            sigma_i = matData.X[i];     // row/destination
            delta_i = matData.D[i];

        for (size_t j = 0; j < n; j++)
        {
            sigma_j = matData.X[j];     // column/source
            delta_j = matData.D[j];
            
            if (i == j)
            {
                if (false)
                    P(i, j) = 1e10;
                else
                    P(i, j) = 1. - 1. / (2. * sigma_j * delta_j) * (1. - 2. * expint(3, sigma_j * delta_j));
            }
            else
            {
                P(i,j) = 1. / (2. * sigma_j*delta_j) * (expint(3,tau[i][j]) - expint(3,tau[i][j] + sigma_i*delta_i)
                    - expint(3, tau[i][j] + sigma_j*delta_j) + expint(3, tau[i][j] + sigma_i*delta_i + sigma_j*delta_j));
            }
        }
    }


    std::vector<std::vector<double>> tau_r{ calculateTauR(matData) };
    Eigen::MatrixXd P_r{ Eigen::MatrixXd::Zero(n, n) };
    for (size_t i = 0; i < n; i++)
    {
        sigma_i = matData.X[i];     // row/destination
        delta_i = matData.D[i];
        for (size_t j = 0; j < n; j++)
        {
            sigma_j = matData.X[j];     // column/source
            delta_j = matData.D[j];
            if (i == j)
            {
                P_r(i, j) = 1. / (2. * sigma_j * delta_j) * (expint(3, tau_r[i][j]) - expint(3, tau_r[i][j] + sigma_i * delta_i)
                    - expint(3, tau_r[i][j] + sigma_j * delta_j) + expint(3, tau_r[i][j] + sigma_i * delta_i + sigma_j * delta_j));
            }
            else
            {
                P_r(i, j) = 1. / (2. * sigma_j * delta_j) * (expint(3, tau_r[i][j]) - expint(3, tau_r[i][j] + sigma_i * delta_i)
                    - expint(3, tau_r[i][j] + sigma_j * delta_j) + expint(3, tau_r[i][j] + sigma_i * delta_i + sigma_j * delta_j));
            }
        }
    }
    //std::cout << P_r;


    //calcSurftoSurf(matData);
    //calcSurfToVol(matData, calculateTauBC(matData));
    //calcVolToSurf(matData);
    if (prob.rboundary == Rboundary::reflective)
        return P + P_r;
    else
        return P;
}