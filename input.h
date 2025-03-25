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

    std::vector<int> materials; // list of materials
    std::vector<double> sig_a; // xs corresponding to the above material
    std::vector<double> sig_s;
    std::vector<double> sig_f;
    std::vector<double> nsig_f;
    std::vector<double> sigma;

    std::vector<int> cellMaterial; // list of materials for each cell in order
    std::vector<double> dimensions; // position of each material interface
    std::vector<double> source; // source term for each cell
};

struct Problem parseInput(std::string name);

std::ostream& operator<<(std::ostream& out, ProblemType::ProblemType probType);
std::ostream& operator<<(std::ostream& out, Lboundary::Lboundary lbound);
std::ostream& operator<<(std::ostream& out, Rboundary::Rboundary rbound);