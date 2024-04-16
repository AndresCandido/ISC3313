#include <iostream>
#include <cmath>

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
		static double GoldenSectionSearch(double a, double b, double tol, double(*f)(double x));

		static double Parabolic_Interpolation(double x1, double x2, double x3, double tol, double(*f)(double x));
	};
}

