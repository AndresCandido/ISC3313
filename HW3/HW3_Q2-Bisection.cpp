#include <iostream>
#include <cmath>
#include <math.h>
#include <iomanip>      // std::setprecision

using namespace std;

double f(double m) {
    double Fm;

    Fm = (2*m)-0.5;

    return Fm;
}

int main() {

    //specify a bracket [a,b] that you think contains the root
    double a, b, c;

    cout << "Enter upper and lower brackets:" << endl << "Upper bracket:";
    cin >> a;
    cout << "Lower bracket:";
    cin >> b;

    double RelError = 1, MaxIter = 5, Iter = 0;;

    while (RelError > 0.00001 && Iter <= MaxIter) {        //start while loop here, use percent relative error as condition for the loop

        Iter++;
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

        else if (Fa * Fb == 0) { 
        
            cout << "The final root is " << setprecision(9) << c;

            return 0;
        }

        else {
            cout << "No roots";

            return 0;
        }

        if (b != 0) {     //check/update percent relative error condition here

            RelError = abs((b - a) / b);

        }

        else {
            RelError = abs(a - b) / a;
        }

        cout << "The root for iteration " << Iter <<" is: " << setprecision(9) << c << " and the rel. error is: " << RelError << endl;
  
    }   // end while loop here

    return 0; 
}