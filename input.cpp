#include <iostream>
#include <string>
#include <fstream>
#include "input.h"

constexpr std::string_view getProblemType(ProblemType::ProblemType probType)
{
    switch (probType)
    {
    case ProblemType::eigenvalue: return "eigenvalue";
    case ProblemType::fixedsource: return "fixedsource";
    default:    return "invalid";
    }
}

constexpr std::string_view getLboundary(Lboundary::Lboundary lbound)
{
    switch (lbound)
    {
    case Lboundary::reflective: return "reflective";
    case Lboundary::vacuum: return "vacuum";
    default:    return "invalid";
    }
}

constexpr std::string_view getRboundary(Rboundary::Rboundary rbound)
{
    switch (rbound)
    {
    case Rboundary::reflective: return "reflective";
    case Rboundary::vacuum: return "vacuum";
    default:    return "invalid";
    }
}

std::ostream& operator<<(std::ostream& out, ProblemType::ProblemType probType)
{
    return out << getProblemType(probType);
}

std::ostream& operator<<(std::ostream& out, Lboundary::Lboundary lbound)
{
    return out << getLboundary(lbound);
}

std::ostream& operator<<(std::ostream& out, Rboundary::Rboundary rbound)
{
    return out << getRboundary(rbound);
}

struct Problem parseInput(std::string name)
{
    Problem problem;
    
    std::ifstream file(name);
    if (!file)
    {
        std::cerr << "Error: file not found" << std::endl;
        std::exit(1);
    }

    std::string line;

    // assign problem type
    std::getline(file, line);
    line = std::strtok(&line[0], "=");
    problem.problemType = (line == "fixedsource" ? ProblemType::fixedsource : ProblemType::eigenvalue);

    // assign boundary conditions
    std::getline(file, line);
    line = std::strtok(&line[0], "=");
    if (line == "specular")
        problem.lboundary = Lboundary::reflective;
    else if (line == "periodic")
        std::exit(1);
    else
        problem.lboundary = Lboundary::vacuum;

    std::getline(file, line);
    line = std::strtok(&line[0], "=");
    if (line == "specular")
        problem.rboundary = Rboundary::reflective;
    else if (line == "periodic")
        std::exit(1);
    else
        problem.rboundary = Rboundary::vacuum;

    // assign materials
    std::getline(file, line); // skip line
    while (std::getline(file, line))
    {
        if (line[0] == '#')
            break;

        char* token = std::strtok(&line[0], " ");
        problem.materials.push_back(std::stoi(token));

        problem.sig_a.push_back(std::stod(std::strtok(NULL, " ")));
        problem.sig_s.push_back(std::stod(std::strtok(NULL, " ")));
        problem.sig_f.push_back(std::stod(std::strtok(NULL, " ")));
        problem.nsig_f.push_back(std::stod(std::strtok(NULL, " ")));
        if (std::strtok(NULL, " ") != NULL)
            std::cerr << "Error: too many cross-sections" << std::endl;
    }

    // build geometry
    while (std::getline(file, line))
    {
        char* token = std::strtok(&line[0], " ");
        problem.cellMaterial.push_back(std::stoi(token));
        problem.dimensions.push_back(std::stod(std::strtok(NULL, " ")));
        if (problem.problemType == ProblemType::fixedsource)
            problem.source.push_back(std::stod(std::strtok(NULL, " ")));
        if (std::strtok(NULL, " ") != NULL)
            std::cerr << "Error: too many arguments" << std::endl;
    }

    for (size_t i = 0; i < problem.sig_a.size(); i++)
        problem.sigma.push_back(problem.sig_a[i] + problem.sig_s[i]);

    return problem;
}