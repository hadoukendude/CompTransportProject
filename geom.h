#pragma once
#include <vector>
#include "input.h"

struct Geometry
{
    std::vector<double> indices;
    std::vector<double> edges;
};

std::vector<double> createIntervals(Problem prob);

std::vector<double> findMidpoints(std::vector<double> mesh);