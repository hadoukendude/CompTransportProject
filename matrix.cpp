#include "input.h"
#include "geom.h"
#include "matrix.h"

void initializeDelta(struct Problem& prob, struct Geometry& geom, struct MatrixData& matrix_data)
{
    size_t n{ geom.edges.size() - 1 };
    matrix_data.D = Eigen::VectorXd::Zero(n);
    for (size_t i = 0; i < n; i++)
        matrix_data.D(i) = geom.edges[i + 1] - geom.edges[i];
}

void initializeXS(struct Problem& prob, struct Geometry& geom, struct MatrixData& matrix_data)
{
    size_t n{ geom.indices.size() };
    matrix_data.X = Eigen::VectorXd::Zero(n);
    matrix_data.S = Eigen::VectorXd::Zero(n);
    matrix_data.F = Eigen::VectorXd::Zero(n);

    matrix_data.G = Eigen::VectorXd::Zero(n);
    matrix_data.Flux = Eigen::VectorXd::Ones(n);

    int i{ 0 }; // cell number
    double x_i{ 0 };
    int matNo{ 0 };

    for (int index = 0; index < n; index++) // mesh point/subinterval number
    {
        i = 0;
        x_i = geom.indices[index];
        for (auto& e : prob.dimensions)
        {
            if (x_i < e)
                break;
            i++;
        }
        if (prob.problemType == ProblemType::fixedsource)
            matrix_data.G[index] = prob.source[i];
        matNo = prob.materials[prob.cellMaterial[i]];
        matrix_data.X[index] = prob.sigma[matNo] ? prob.sigma[matNo] : 1e-8;
        matrix_data.S[index] = prob.sig_s[matNo];
        matrix_data.F[index] = prob.nsig_f[matNo];
    }
}