#include <sciplot/sciplot.hpp>
#include <Eigen/Dense>
#include "matrix.h"

using namespace sciplot;
void plotter(const MatrixData& mat)
{
    Vec x = linspace(0.0, mat.D.sum(), (mat.D.size()-1));

    Plot2D plot;

    plot.xlabel("z [cm]");
    plot.ylabel("Flux [n/cm**2/s]");

    plot.drawCurve(x, mat.Flux.array());
    plot.legend().hide();

    // Create figure to hold plot
    Figure fig = { {plot} };
    // Create canvas to hold figure
    Canvas canvas = { {fig} };

    canvas.show();
}