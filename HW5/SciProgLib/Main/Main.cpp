#include "SciProgLib.h"
#include <iostream>
#include <cmath>
#include <math.h>
#include <iomanip>      // std::setprecision
#include <fstream>
#include <vector>

using namespace std;

double f(double t) {

    double cd = 0.25, g = 9.81, m = 68.1, Fm;

    Fm = sqrt(g * m / cd) * tanh(sqrt(g * cd / m) * t);

    return Fm;
}

double df(double m) {

    double cd = 0.25, v = 36, t = 3, g = 9.81, FPm;

    FPm = (0.5) * sqrt(g / (m * cd)) * tanh(sqrt(g * cd / m) * t) - (g * t / (2 * m)) * (1 - pow(tanh(sqrt(g * cd / m) * t), 2));

    return FPm;
}

double f2(double x) {

    return pow(x, 2) / 10.0 - 2.0 * sin(x);
}

double f3(double m) {

    double Fm;

    Fm = 0.2 + (25.0 * m) - (200.0 * pow(m, 2)) + (675.0 * pow(m, 3)) - (900.0 * pow(m, 4)) + (400.0 * pow(m, 5));

    return Fm;
}

double f4(double m) {

    return -0.1 * pow(m, 4) - 0.15 * pow(m, 3) - 0.5 * pow(m, 2) - 0.25 * m + 1.2;
}

double f5(double m) {

    return pow(m, 4) + 2 * pow(m, 3) + 8 * pow(m, 2) + 5 * m;
}

int main() {

    //double root = SciProgLib::RootFinding::Bisect(140, 150, 0.1, f);
    //cout << "The root is " << setprecision(9) << root << endl;

    //double root2 = SciProgLib::RootFinding::NewtonRaphson(140, 0.1, 100, f, df);
    //cout << "The root is " << setprecision(9) << root2 << endl;

    double optima = SciProgLib::Optimization::GoldenSectionSearch(-2, 1, 0.001, f5); 
    cout << "The optima is " << setprecision(9) << optima << endl;

    double optima2 = SciProgLib::Optimization::ParabolicInterpolation(-2, -1, 1, 0.00001, f5); 
    cout << "The optima at x is " << setprecision(9) << optima2 << endl;

    double integrate = SciProgLib::Integration::Comp_Trapezoid_rule(0, 3, 2, -0.8668, 0.00001, f);
    cout << "The area under the curve is " << setprecision(9) << integrate << endl;

    double integrate2 = SciProgLib::Integration::Comp_Simpson_rule(0, 3, 2, 0.2001, 0.00001, f);
    cout << "The area under the curve is " << setprecision(9) << integrate2 << endl;

    //double data[] = { 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7 };  // data points must be equispaced for data simpson to work properly
    //double f[] = { 1.543, 1.669, 1.811, 1.971, 2.151, 2.352, 2.577, 2.828 };
    //unsigned int n = sizeof(data) / sizeof(data[0]) - 1;   // length(x) -1
    //double integrate3 = SciProgLib::Integration::DataSimpsonsRules(data, f, n);
    //cout << "The area under the curve is " << setprecision(9) << integrate3 << endl;

    ifstream xdata("xdata.dat", ifstream::binary); // turn file into binary object,  what is a binary object?
    ifstream fdata("fdata.dat", ifstream::binary);
    vector<double> X, Fx;
    double x, fx;
    while (xdata >> x) {  // files must be in the main folder of the library so they can be read 
        X.push_back(x);
    }
    xdata.close();
    while (fdata >> fx) {
        Fx.push_back(fx);
    }
    fdata.close();   // why do we need to put parenthesis after this type of command?
    int n2 = X.size() - 1;


    double integrate4 = SciProgLib::Integration::Data_trapezoid(X, Fx, n2);
    cout << "The area under the curve is (using Data Trapezoid) " << setprecision(9) << integrate4 << endl;

    //double derivative = SciProgLib::Differentiation::BackwardDifference(0.5, 0.5, f4);
    //cout << "The derivative at x=0.5 is " << derivative << endl;

    //double derivative2 = SciProgLib::Differentiation::ForwardDifference(0.5, 0.5, f4);
    //cout << "The derivative at x=0.5 is " << derivative2 << endl;

    //double derivative3 = SciProgLib::Differentiation::CenteredDifference(0.5, 0.5, f4);
    //cout << "The derivative at x=0.5 is " << derivative3 << endl;

    return 0;  
}