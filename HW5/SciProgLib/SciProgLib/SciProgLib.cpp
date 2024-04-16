#include "SciProgLib.h"
#include <fstream>
#include <vector>

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

	double Optimization::GoldenSectionSearch(double lowerx, double upperx, double tol, double(*f)(double x)) 
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

	double Optimization::ParabolicInterpolation(double x1, double x2, double x3, double tol, double(*f)(double x)) {

		double x4, x4old = 0, Error = 1, iter = 0;

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

			iter++;           // this couple of lines are to fit the criteria given in HW5
			if (iter == 5) {  // criteria = run this code for 5 iterations and then return the minimum
				cout << "Iterations = " << iter << endl;
				return x4;
			}

		}

		return x4;
	}

	double Integration::Comp_Trapezoid_rule(double xLower, double xUpper, double n, double avgVal, double tol, double(*f)(double x)) {

		double h, x1, sum, Error;

		Error = -pow(xUpper - xLower, 3) / (12 * pow(n, 2)) * avgVal;
		if (Error > tol) {

			n = 2 * n;
			return  Comp_Trapezoid_rule(xLower, xUpper, n, avgVal, tol, f);
		}

		h = abs(xUpper - xLower) / n;
		sum = h / 2.0 * (f(xLower) + f(xUpper));
		x1 = xLower + h;

		for (int i = 1; i < n; i++) {

			sum = sum + (h * f(x1));
			x1 = x1 + h;
		}

		cout << "The number of segments needed for CT = " << n << endl;
		return sum;
	}

	double Integration::Comp_Simpson_rule(double xLower, double xUpper, double n, double avgVal, double tol, double(*f)(double x)) {

		double h, x0, x1, sum, Error;

		Error = -pow(xUpper - xLower, 5) / (180 * pow(n, 4)) * avgVal;
		if (abs(Error) > tol) {

			n = 2 * n;
			return  Comp_Simpson_rule(xLower, xUpper, n, avgVal, tol, f);
		}

		h = abs(xUpper - xLower) / n;
		sum = h / 3.0 * (f(xLower) + f(xUpper));
		x0 = 0;
		x1 = xLower + h;


		for (int i = 1; i < n; i++) {


			if (i % 2 == 0) { // evens

				sum = sum + (h / 3) * (2 * f(x1));
			}
			else {            // odds

				sum = sum + (h / 3) * (4 * f(x1));
			}

			x1 = x1 + h;
		}

		cout << "The number of segments needed for CS = " << n << endl;
		return sum;
	}

	double Integration::DataSimpsonsRules(double x[], double fx[], int n) {

		double h = (x[n] - x[0]) / n;
		double sum = 0.0;

		if (n % 2 == 0) { // first case procedure

			sum = (h / 3.0) * (fx[0] + fx[n]);

			for (int i = 1; i < n; i++) {

				if (i % 2 == 0) {

					sum = sum + (2.0 * (h / 3.0) * fx[i]);
				}
				else {

					sum = sum + (4.0 * (h / 3.0) * fx[i]);
				}
			}
		}

		else if (n % 3 == 0) {    // second case procedure

			sum = (h / 3.0) * (fx[0] + fx[n]);

			for (int i = 1; i < n; i++) {

				sum = sum + (3.0 * h / 8.0) * (3.0 * fx[i]);

			}
		}

		else {  // third case procedure

			sum = (h / 3) * (fx[0] + fx[n - 3]);

			for (int i = 1; i < n - 3; i++) {

				if (i % 2 == 0) {

					sum = sum + (2.0 * h / 3.0 * fx[i]);
				}

				else {

					sum = sum + (4.0 * h / 3.0 * fx[i]);
				}
			}

			sum = sum + (3.0 * h / 8.0) * (fx[n - 3] + fx[n]);
			for (int i = n - 2; i < n; i++) {

				sum = sum + ((3.0 * h / 8.0) * 3.0 * fx[i]);;
			}
		}

		return sum;
	}

	double Integration::Data_trapezoid(vector<double> X, vector<double> Fx, int n) {
		double sum = 0, h = 0;

		for (int i = 1; i <= n; i++) {

			h = X[i] - X[i - 1];
			sum = sum + ((h / 2) * (Fx[i - 1] + Fx[i]));
		}

		return sum;
	}

	double Differentiation::ForwardDifference(double x, double h, double(*f)(double x)) {
	
		return (f(x + h) - f(x)) / h;
	}

	double Differentiation::BackwardDifference(double x, double h, double(*f)(double x)) {
	
		return (f(x) - f(x -h)) / h;
	}

	double Differentiation::CenteredDifference(double x, double h, double(*f)(double x)) {
	
		return (f(x + h) - f(x - h)) / (2 * h);
	}
}

