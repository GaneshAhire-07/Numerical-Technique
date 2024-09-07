#include "BisectionMethod.hpp"
#include <bits/stdc++.h>
using namespace std;
// Constructor initializes the function, tolerance, and maximum iterations
BisectionMethod::BisectionMethod(function<double(double)> f, double tol, int maxIter) {
    func = f;
    tolerance = tol;
    maxIterations = maxIter;
}

// Method to find the root using the Bisection Method
void BisectionMethod::findRoot(double a, double b) {
    // Open output file to log iterations
    ofstream outFile("output.txt");
    if (!outFile.is_open()) {
        cout << "Failed to open output file." << endl;
        return;
    }

    // Check if initial interval is valid
    if (func(a) * func(b) >= 0) {
        cout << "Invalid interval: The function must have different signs at the endpoints a and b." << endl;
        return;
    }

    double c; // Midpoint
    for (int iter = 1; iter <= maxIterations; ++iter) {
        c = (a + b) / 2; // Calculate midpoint
        double fc = func(c); // Evaluate function at midpoint

        // Log current iteration to file
        outFile << "Iteration : " << iter << " +ve = " << fixed  << c << ", -ve = " << fc << endl;

        // Check if the current midpoint is a root
        if (fabs(fc) < tolerance) { //check the absolute value tolerance level
            cout << "Root found: " << c << " after " << iter << " iterations." << endl;
            outFile.close();
            return; //close the function after found the root
        }

        // Determine the subinterval that contains the root
        if (func(a) * fc < 0) { //therefore is a is -ve value to swap the a and b
            b = c; // Root is in the left
        } else {
            a = c; // Root is in the right subinterval
        }
    }

    // If the loop exits without finding the root within tolerance
    cout << "Max iterations reached without finding the root to the desired tolerance." << endl;
    outFile.close();
}
