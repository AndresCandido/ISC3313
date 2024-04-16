#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

using namespace std;

namespace SciProgLib 
{
	class RootFinding
	{
	public:
		static double Bisect(double lower, double upper, double tol, double(*f)(double m));

		static double NewtonRaphson(double x, double tol, double Maxit, double(*f)(double m), double(*df)(double m));
	};

	class Optimization 
	{
	public:
		static double GoldenSectionSearch(double xlower, double xupper, double tol, double(*f)(double x));

		static double ParabolicInterpolation(double x1, double x2, double x3, double tol, double(*f)(double x));

	};

	class Integration
	{
	public:
		static double Comp_Trapezoid_rule(double xLower, double xUpper, double n, double avgVal, double tol, double(*f)(double x));

		static double Comp_Simpson_rule(double xLower, double xUpper, double n, double avgVal, double tol, double(*f)(double x));

		static double DataSimpsonsRules(double x[], double fx[], int n);

		static double Data_trapezoid(vector<double> X, vector<double> Fx, int n);

	};

	class Differentiation
	{
	public:
		static double ForwardDifference(double x, double h, double(*f)(double x));

		static double BackwardDifference(double x, double h, double(*f)(double x));

		static double CenteredDifference(double x, double h, double(*f)(double x));

	};
}

