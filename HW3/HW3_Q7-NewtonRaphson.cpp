#include <iostream>
#include <cmath>
#include <math.h>
#include <iomanip>      // std::setprecision

using namespace std;



double f(double cd) {

    double m = 95, v = 46, t = 9, g = 9.81, Fm;

    Fm = sqrt(g * m / cd) * tanh(sqrt(g * cd / m) * t) - v;

    return Fm;
}

double fprime(double cd) {

    double v = 46, t = 9, g = 9.81, m = 95, FPm;

    FPm = -0.5 * sqrt((g * m) / pow(cd, 3)) * tanh(sqrt((g * cd) / m) * t) + (g * t) / (2 * m) * (1 - (pow(tanh(sqrt((g * cd) / m) * t), 2)));  

    return FPm;
}

void NewtonRaphson(double x, double RelError, double MaxIter, double Iter) {

    double newx;

    while (RelError > 0.00001 && Iter <= MaxIter && abs(f(x)) > 0) {    //start while loop that evaluates the newton-raphson method here, use percent relative error as condition for the loop

         //use the f and fprime functions above to find the values of the function and its prime at certain values of x
        newx = x - (f(x) / fprime(x));

        RelError = abs(newx - x);    //stopping criteria being if the function value is close to zero
        Iter++;

        x = newx;

    }

    cout << "The root is " << setprecision(9) << x << ", the number of iterations was: " << Iter;

    return;
}

int main() {

    double x;

    cout << "Input a value for cd: ";
    cin >> x;

    double RelError = 1, MaxIter = 3, Iter = 0;

    NewtonRaphson(x, RelError, MaxIter, Iter);

    return 0;
}