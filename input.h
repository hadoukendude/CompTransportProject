#pragma once
#include <iostream>
#include <string>
#include <vector>

namespace ProblemType
{
    enum ProblemType
    {
        fixedsource,
        eigenvalue,
    };
}
namespace Lboundary
{
    enum Lboundary
    {
        vacuum,
        reflective,
    };
}
namespace Rboundary
{
    enum Rboundary
    {
        vacuum,
        reflective,
    };
}

struct Problem
{
    ProblemType::ProblemType problemType;
    Lboundary::Lboundary lboundary;
    Rboundary::Rboundary rboundary;

    std::vector<int> materials;
    std::vector<double> sig_a;
    std::vector<double> sig_s;
    std::vector<double> sig_f;
    std::vector<double> nsig_f;
    std::vector<double> sigma;

    std::vector<int> cellMaterial;
    std::vector<double> dimensions;
    std::vector<double> source;
};

struct Problem parseInput(std::string name);

std::ostream& operator<<(std::ostream& out, ProblemType::ProblemType probType);
std::ostream& operator<<(std::ostream& out, Lboundary::Lboundary lbound);
std::ostream& operator<<(std::ostream& out, Rboundary::Rboundary rbound);