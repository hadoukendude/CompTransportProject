// CompTransportProject.cpp : Defines the entry point for the application.
//

#include <iostream>
#include "main.h"
#include "input.h"
#include "geom.h"
#include "transport.h"

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
    
    //std::cout << findMidpoints(createIntervals(problem));

    std::cout << createOpticalDepthMatrix(problem, findMidpoints(createIntervals(problem)));

    return 0;
}
