#include <vector>
#include <sciplot/sciplot.hpp>
#include <Eigen/Dense>
#include "input.h"
#include "matrix.h"
#include <cmath>

using namespace sciplot;
void plotter(const MatrixData& mat, ProblemType::ProblemType pt)
{
    Eigen::Index n{ mat.D.size() };
    std::vector<double> x(n);
    for (size_t i = 1; i < n; i++)
        x[i] = mat.D(i) + x[i - 1];

    Plot2D plot;
    plot.legend();

    plot.xlabel("z [cm]");
    if (pt == ProblemType::fixedsource)
    {
        plot.ylabel("Flux [n/cm**2/s]");
        plot.drawCurve(x, mat.Flux.array());
    }
    else if (pt == ProblemType::eigenvalue)
    {
        plot.ylabel("Normalized Flux");
        plot.drawCurve(x, mat.Flux.array()/mat.Flux.maxCoeff()).label("CPM");
    }

    //std::vector <double> y(n);
    //for (size_t i=0;i<n;i++)
    //    y[i] = std::sin(0.197484176581 * (x[i]+0.4));
    //plot.drawCurve(x, y).label("Diffusion Equation");


    // Create figure to hold plot
    Figure fig = { {plot} };
    // Create canvas to hold figure
    Canvas canvas = { {fig} };

    canvas.show();
    canvas.save("flux.png");
}