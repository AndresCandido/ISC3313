#include "SciProgLib.h"
#include <iostream>
#include <cmath>
#include <math.h>
#include <iomanip>      // std::setprecision

using namespace std;

double f(double m) {

    double cd = 0.25, v = 36, t = 4, g = 9.81, Fm;

    Fm = sqrt(g * m / cd) * tanh(sqrt(g * cd / m) * t) - v;

    return Fm;
}

double df(double m) {

    double cd = 0.25, v = 36, t = 4, g = 9.81, FPm;

    FPm = (0.5) * sqrt(g / (m * cd)) * tanh(sqrt(g * cd / m) * t) - (g * t / (2 * m)) * (1 - pow(tanh(sqrt(g * cd / m) * t), 2));

    return FPm;
}

double f2(double x) {

    return pow(x, 2) / 10.0 - 2.0 * sin(x);
}

int main() {

    double root = SciProgLib::RootFinding::Bisect(140, 150, 0.1, f);

    cout << "The root is " << setprecision(9) << root << endl;

    double root2 = SciProgLib::RootFinding::NewtonRaphson(140, 0.1, 100, f, df);

    cout << "The root is " << setprecision(9) << root2 << endl;

    double optima = SciProgLib::Optimization::GoldenSectionSearch(4, 0, 0.00001, f2);

    cout << "The optima is " << setprecision(9) << optima << endl;

    double optima2 = SciProgLib::Optimization::Parabolic_Interpolation(0, 1, 4, 0.00001, f2);

    cout << "The optima at x is " << setprecision(9) << optima2 << endl;

    return 0;  
}