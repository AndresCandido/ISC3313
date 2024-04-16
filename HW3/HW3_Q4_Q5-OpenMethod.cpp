#include <iostream>
#include <cmath>
#include <math.h>
#include <iomanip>      // std::setprecision

using namespace std;



double f(double m) {

    double Fm;

    Fm = sin(sqrt(m)) - m;

    return Fm;
}

double fprime(double m) {

    double FPm;

    FPm = ((cos(pow(m,0.5)))/(2*pow(m,0.5)))-1;

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

        cout << "The root is " << setprecision(9) << x << ", the number of iterations is: " << Iter << endl;

    }

    cout << "The final root is " << setprecision(9) << x;

    return;
}

int main() {

    double x;

    cout << "Input a value for x: ";
    cin >> x;

    double RelError = 1, MaxIter = 100, Iter = 0;

    NewtonRaphson(x, RelError, MaxIter, Iter);

    return 0;
}