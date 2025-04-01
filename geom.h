#pragma once
#include <vector>
#include "input.h"

struct Geometry
{
    std::vector<double> indices;
    std::vector<double> edges;
};

std::vector<double> createIntervals(const Problem& prob);

std::vector<double> findMidpoints(const std::vector<double>& mesh);