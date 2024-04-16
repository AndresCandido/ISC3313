//This program takes two floats from the user, multiplies them, and shows their product

#include <iostream>  //iostream is needed for cin & cout
#include <iomanip>   // iomanip is needed for fixed and setprecision, which show the desired decimal point
using namespace std; 

int main()
{
	double x, y, prod;        // Here we declare all variables needed in the program
	

	cout << "Enter the value of x:" << endl;    // ask user for x
	cin >> x;

	cout << "Enter the value of y:" << endl;   // ask user for y
	cin >> y;

	prod = x * y;                              // multiply x by y

	cout << "The product of x and y is: " << fixed << setprecision(6) << prod << endl; //display the product of x and y. setprecision is used here to show 6 decimal space

	return 0;
}