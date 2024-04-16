#include "SciProgLib.h"
#include <fstream>
#include <vector>

namespace SciProgLib
{
	double RootFinding::Bisect(double lower, double upper, double tol, double(*f)(double m))       
	{
		double iter = 0;

		if (f(upper) * f(lower) >= 0) // check if points are greater or equal to zero 
		{
			cout << "No roots.";

			return 0;
		}

		double error = 1, c_old = 0, c_new = 0;

		while (error > tol) // runs loop while the error is greater than the stablished tolerance
		{
			c_new = (upper + lower) / 2;

			if (f(upper) * f(c_new) < 0) //bisection occurs in the if / else statement
			{
				lower = c_new;
			}
			else
			{
				upper = c_new;
			}

			error = abs((c_new - c_old) / c_new); //calculate new error
			c_old = c_new;                        //update value of c_new

			iter++;
		}

		cout << "Iterations to find root with Bisect method: " << iter << endl;

		return c_new; 
	}

	double RootFinding::NewtonRaphson(double x, double tol, double Maxit, double(*f)(double m), double(*df)(double m))      // Newton Raphson fucntion
	{
		double newx, iter = 0, RelError = 1;

		while (RelError > tol && iter <= Maxit && abs(f(x)) > 0) {   // runs loop while the error is greater than the stablished tolerance, iterations are less than maxit, and abs(f(x)) > 0.

			newx = x - (f(x) / df(x)); //find new value of x

			RelError = abs(newx - x);   //update error 
			iter++;

			x = newx; //update x
		}

		cout << "Iterations to find root with Newton Raphson method: " << iter << endl;

		return x;
	}

	double Optimization::GoldenSectionSearch(double lowerx, double upperx, double tol, double(*f)(double x), bool min) 
	{
		//specify a bracket [a,b] that you think contains the minimum/maximum aka the domain of the function
		double goldenR, x1, x2, Xopt = 0, Error = 1, phi = 1.61803398874989, coef = 1;

		goldenR = (phi - 1) * (upperx - lowerx); // calculate golden ratio

		if (min == false) { //here we are able to change whether we find the minimum or maximum
			coef = -1;
		}
		
		x1 = lowerx + goldenR;
		x2 = upperx - goldenR;

		while (Error > tol){   // runs loop while the error is greater than the stablished tolerance

			if (coef * (f(x1)) < coef * (f(x2)))  //GoldenSectionSearch occurs in if / else statement
			{
				lowerx = x2;
				x2 = x1;
				goldenR = (phi - 1) * (upperx - lowerx); //update golden ratio
				x1 = lowerx + goldenR;
				Xopt = x1; //update xopt
			}
			else
			{
				upperx = x1;
				x1 = x2;
				goldenR = (phi - 1) * (upperx - lowerx); //update golden ratio
				x2 = upperx - goldenR;
				Xopt = x2; //update xopt
			}

			Error = (2 - phi) * abs((upperx - lowerx) / Xopt); //update error
		}    

		return Xopt; // Xopt is the minimum of the function
	}

	double Optimization::ParabolicInterpolation(double x1, double x2, double x3, double tol, double(*f)(double x)) {

		double x4, x4old = 0, Error = 1;

		double den = ((x2 - x1) * (f(x2) - f(x3)) - (x2 - x3) * (f(x2) - f(x1))); //get denominator to check if it is equal to 0

		if (den == 0) { //check if den = 0

			den = numeric_limits<double>::min(); // if it is then den becomes a very small number instead of zero
		}

		x4 = x2 - 0.5 * ( (pow(x2 - x1, 2) * (f(x2) - f(x3)) - pow(x2 - x3, 2) * (f(x2) - f(x1))) /  //get value of x4
			den );



		while (Error > tol) {   // runs loop while the error is greater than the stablished tolerance

			x4old = x4;

			if (x2 > x4) { //interval gets smaller ( closer to the optima )
				x3 = x2;
				x2 = x4;
			}
			else {
				x1 = x2;
				x2 = x4;
			}
			
			den = ((x2 - x1) * (f(x2) - f(x3)) - (x2 - x3) * (f(x2) - f(x1))); //update den

			if (den == 0) { //check again

				den = numeric_limits<double>::min();
			}

			x4 = x2 - 0.5 * ((pow(x2 - x1, 2) * (f(x2) - f(x3)) - pow(x2 - x3, 2) * (f(x2) - f(x1))) /    //update x4
				den);


			Error = abs((x4 - x4old) / x4); //update error


		}

		return x4;
	}

	double Integration::Comp_Trapezoid_rule(double xLower, double xUpper, double n, double avgVal, double tol, double(*f)(double x)) {

		double h, x1, sum, Error;

		Error = -pow(xUpper - xLower, 3) / (12 * pow(n, 2)) * avgVal; //compute error
		if (Error > tol) {   //if the error is greater than the stablished tolerance runs function again with new n

			n = 2 * n;
			return  Comp_Trapezoid_rule(xLower, xUpper, n, avgVal, tol, f);  
		}

		h = abs(xUpper - xLower) / n; //compute h for getting values of the sum
		sum = h / 2.0 * (f(xLower) + f(xUpper)); //get current value of sum
		x1 = xLower + h; //update x1

		for (int i = 1; i < n; i++) { //for loop adds up the sums, also updates x1

			sum = sum + (h * f(x1));
			x1 = x1 + h;
		}

		cout << "The number of segments needed for CT = " << n << endl;
		return sum;
	}

	double Integration::Comp_Simpson_rule(double xLower, double xUpper, double n, double avgVal, double tol, double(*f)(double x)) {

		double h, x0, x1, sum, Error;

		Error = -pow(xUpper - xLower, 5) / (180 * pow(n, 4)) * avgVal; //comput error
		if (abs(Error) > tol) {  // if the absolute error is greater than the stablished tolerance runs function again with new n

			n = 2 * n;
			return  Comp_Simpson_rule(xLower, xUpper, n, avgVal, tol, f);
		}

		h = abs(xUpper - xLower) / n;                      //similar to comp_trapezoid
		sum = h / 3.0 * (f(xLower) + f(xUpper));
		x0 = 0;
		x1 = xLower + h;


		for (int i = 1; i < n; i++) {  //loop with diffferent cases based on if i is an odd or even number


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

		for (int i = 1; i <= n; i++) { //for loop adds up the sums, also updates h

			h = X[i] - X[i - 1];
			sum = sum + ((h / 2) * (Fx[i - 1] + Fx[i]));
		}

		return sum;
	}

	double Differentiation::ForwardDifference(double x, double h, double(*f)(double x)) {

		return (f(x + h) - f(x)) / h;
	}

	double Differentiation::BackwardDifference(double x, double h, double(*f)(double x)) {

		return (f(x) - f(x - h)) / h;
	}

	double Differentiation::CenteredDifference(double x, double h, double(*f)(double x)) {

		return (f(x + h) - f(x - h)) / (2 * h);
	}

	double RootFinding::Secant(double x0, double x1, double tol, int Maxit, double(*f)(double x)) {

		double newx, iter = 0, RelError = 1;

		while (RelError > tol && iter <= Maxit && abs(f(x1)) > 0) {

			newx = x1 - (f(x1) / Differentiation::BackwardDifference(x1, x1 - x0, f));  //use backwards difference based on the nature of the secant method

			RelError = abs(newx - x1); //update error
			iter++;

			x0 = x1;
			x1 = newx;
		}

		cout << "Iterations to find root with Secant method: " << iter << endl;

		return x1;

	}

	double Integration::GaussLegendre(int n, double(*f)(double(x))){
		double I = 0, Error;
		double c[6], x[6];

		if (n == 1) {   // diffierent weights (c) and arguments (x) arrays for each case of n (1 through 6)
			c[0] = 2;

			x[0] = 0.0;
		}
		else if(n == 2) {
			c[0] = 1;
			c[1] = 1;

			x[0] = -1 / sqrt(3);
			x[1] = 1 / sqrt(3);
		}
		else if (n == 3) { //It would appear the equation f(x) used returns an I equal to 0 when n=3 
			c[0] = 5/9;
			c[1] = 8/9;
			c[2] = 5/9;

			x[0] = -sqrt(3/5);
			x[1] = 0.0;
			x[2] = sqrt(3/5);
		}
		else if (n == 4) {
			c[0] = (18 - sqrt(30)) / 36;
			c[1] = (18 + sqrt(30)) / 36;
			c[2] = (18 + sqrt(30)) / 36;
			c[3] = (18 - sqrt(30)) / 36;

			x[0] = -sqrt(525 + 70 * sqrt(30)) / 35;
			x[1] = -sqrt(525 - 70 * sqrt(30)) / 35;
			x[2] = sqrt(525 - 70 * sqrt(30)) / 35;
			x[3] = sqrt(525 + 70 * sqrt(30)) / 35;
		}
		else if (n == 5) {
			c[0] = (322 - 13 * sqrt(70)) / 900;
			c[1] = (322 + 13 * sqrt(70)) / 900;
			c[2] = 128 / 225;
			c[3] = (322 + 13 * sqrt(70)) / 900;
			c[4] = (322 - 13 * sqrt(70)) / 900;

			x[0] = -sqrt(245 + 14 * sqrt(70)) / 21;
			x[1] = -sqrt(245 - 14 * sqrt(70)) / 21;
			x[2] = 0.0;
			x[3] = sqrt(245 - 14 * sqrt(70)) / 21;
			x[4] = sqrt(245 + 14 * sqrt(70)) / 21;
		}
		else if (n == 6) {
			c[0] = 0.171324492379170;
			c[1] = 0.360761573048139;
			c[2] = 0.467913934572691;
			c[3] = 0.467913934572691;
			c[4] = 0.360761573048139;
			c[5] = 0.171324492379170;

			x[0] = -0.932469514203152;
			x[1] = -0.661209386466265;
			x[2] = -0.238619186083197;
			x[3] = 0.238619186083197;
			x[4] = 0.661209386466265;
			x[5] = 0.932469514203152;
		}
		else {
			cout << "n must be a number from 1 through 6" << endl;
			return 0;
		}

		for (int i = 0; i < n; i++) { //Find the sum to get the vaule of I by adding up all Is

			I += c[i] * f(x[i]);
		}

		Error = (abs(I - 1.640533333333341) / 1.640533333333341) * 100; //Calculate percent relative error

		cout << "The percent relative error is " << Error << endl;

		return I;
}

}

