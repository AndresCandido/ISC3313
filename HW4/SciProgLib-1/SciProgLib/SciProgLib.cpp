#include "SciProgLib.h"

namespace SciProgLib
{
	double RootFinding::Bisect(double lower, double upper, double tol, double(*f)(double m))       // a = upper , b = lower 
	{
		if (f(upper) * f(lower) >= 0)
		{
			cout << "No roots.";

			return 0;
		}

		double error = 1, c_old = 0, c_new = 0;

		while (error > tol)
		{
			c_new = (upper + lower) / 2;

			if (f(upper) * f(c_new) < 0)
			{
				upper = c_new;
			}
			else
			{
				lower = c_new;
			}

			error = abs((c_new - c_old) / c_new);
			c_old = c_new;
		}

		return c_new; 
	}


	double RootFinding::NewtonRaphson(double x, double tol, double Maxit, double(*f)(double m), double(*df)(double m))      // Newton Raphson fucntion
	{
		double newx, iter = 0, RelError = 1;

		while (RelError > tol && iter <= Maxit && abs(f(x)) > 0) {   

			newx = x - (f(x) / df(x));

			RelError = abs(newx - x);    
			iter++;

			x = newx;
		}

		return x;
	}


	double Optimization::GoldenSectionSearch(double upperx, double lowerx, double tol, double(*f)(double x)) 
	{
		//specify a bracket [a,b] that you think contains the minimum aka the domain of the function
		double goldenR, x1, x2, Xopt = 0, Error = 1, phi = 1.61803398874989;

		goldenR = (phi - 1) * (upperx - lowerx); // calculate golden ratio

		x1 = lowerx + goldenR;
		x2 = upperx - goldenR;

		while (Error > tol){

			if (f(x1) < f(x2))
			{
				lowerx = x2;
				x2 = x1;
				goldenR = (phi - 1) * (upperx - lowerx);
				x1 = lowerx + goldenR;
				Xopt = x1;
			}
			else
			{
				upperx = x1;
				x1 = x2;
				goldenR = (phi - 1) * (upperx - lowerx);
				x2 = upperx - goldenR;
				Xopt = x2;
			}

			Error = (2 - phi) * abs((upperx - lowerx) / Xopt);
		}    

		return Xopt; // Xopt is the minimum of the function
	}


	double Optimization::Parabolic_Interpolation(double x1, double x2, double x3, double tol, double(*f)(double x)) {

		double x4, x4old = 0, Error = 1;

		x4 = x2 - 0.5 * ( (pow(x2 - x1, 2) * (f(x2) - f(x3)) - pow(x2 - x3, 2) * (f(x2) - f(x1))) /
			((x2 - x1) * (f(x2) - f(x3)) - (x2 - x3) * (f(x2) - f(x1)) ) );

		while (Error > tol) {

			x4old = x4;

			if (x2 > x4) {
				x3 = x2;
				x2 = x4;
			}
			else {
				x1 = x2;
				x2 = x4;
			}

			x4 = x2 - 0.5 * ((pow(x2 - x1, 2) * (f(x2) - f(x3)) - pow(x2 - x3, 2) * (f(x2) - f(x1))) /
				((x2 - x1) * (f(x2) - f(x3)) - (x2 - x3) * (f(x2) - f(x1))));

			Error = abs((x4 - x4old) / x4);
		}

		return x4;
	}

}

