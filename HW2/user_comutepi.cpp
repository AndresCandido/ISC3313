#include <iostream>  //libraries included in the original computepi.cpp
#include <random>
#include <ctime>
#include <cmath>

using namespace std;

//Use Define to define pi
#define MYPI 3.141592653589793238
int main()
{
    double e; // e is used to calculate relative absolute error and it is determined by user input (10^-2, 10^-3, or 10^-4)
    cout << "Enter value of e: ";
    cin >> e; 

    double Ntot_Sum = 0;

    // initialize random seed
    srand((unsigned)time(0));

    for (int i = 0; i < 5; i++) {      // loop the program 5 times to get the average Ntot

        int Nhit = 0, Ntot = 0;        // these are needed for the approximation of pi and the final result
        double E = 1;                  // E, Nhit, and Ntot reset to 0 every time the for loop starts
                                       // E starts at 1 so that the first while loop can start
        while (E > e)
        {
            Ntot++;           // every time the while loop runs Ntot increases by 1           

            double x = (double)rand() / RAND_MAX;
            double y = (double)rand() / RAND_MAX;

            // check if the random numbers lie in the circle using the formula in the slides, 
            if (sqrt(pow(x - 0.5, 2) + pow(y - 0.5, 2)) < 0.5)
            {
                // if the dart landed in the circle increase the number of darts that made it in the 
                // circle by 1

                Nhit++;

            }

            double pi = 4.0 * Nhit / Ntot;  // pi is determined with the given formula

            E = (abs(MYPI - pi)) / MYPI;    // new E is determined using the given formula
        }


        // print the Ntot of each for loop run, useful for checking if average is correct
        cout << Ntot << endl;

        Ntot_Sum = Ntot_Sum + Ntot; // Add all Ntots

    }

    double average_Ntot = Ntot_Sum / 5; // get Ntot average 

    cout << "The average number of iterations is: " << average_Ntot << endl;    // print average

    return 0;
       

}
