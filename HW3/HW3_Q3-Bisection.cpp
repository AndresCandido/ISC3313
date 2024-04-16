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

int main() {

    //specify a bracket [a,b] that you think contains the root
    double a, b, c;

    cout << "Enter upper and lower brackets:" << endl << "Upper bracket:";
    cin >> a;
    cout << "Lower bracket:";
    cin >> b;

    double RelError = 1, Iter = 1;

    while (RelError > 0.01) {        //start while loop here, use percent relative error as condition for the loop

         //evaluate f(a) and f(b)

        double Fa, Fb; //plug a and b into given function

        Fa = f(a);
        Fb = f(b);

        //if f(a)*f(b)<0 then you have at least one root

        if (Fa * Fb < 0) {

            //compute the midpoint within the interval,c=1/2(a+b)
            c = (a + b) / 2;

            //check one of the new half - intervals for the root, f(a)f(c) < 0.0
            double Fc; //plug c into function

            Fc = f(c);

            if (Fa * Fc < 0) {

                //if found, then set b=c and start at step 2 - otherwise set a=c and start at step 2
                b = c;
            }
            else {
                a = c;
            }
        }

        else {
            cout << "No roots";

            return 0;
        }
        RelError = abs((b - a) / b);//check/update percent relative error condition here

        Iter++;

    }   // end while loop here

    cout << "The root is " << setprecision(9) << c << ", the number of iterations was: " << Iter;

    return 0;  //exit the algorithm if the percent relative error is less than your tolerance
}