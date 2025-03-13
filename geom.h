#pragma once
#include <vector>
#include "input.h"


std::vector<double> createIntervals(Problem prob);

std::vector<double> findMidpoints(std::vector<double> mesh);