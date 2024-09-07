#include <bits/stdc++.h>
#include "NewtonRaphson.hpp"
using namespace std;

// Define the function f(x) for which we want to find the root
double f(double X)
{
    return (X * X - 4 * X + 2); // Quadratic function
}

// Define the derivative of the function f'(x)
double fdash(double X)
{
    return (2 * X - 4); // Derivative of the quadratic function
}

// Function to perform Newton-Raphson method to find the root of the equation
double NewtonRaphson::Newton_Raphson()
{
    // Repeat the process until the difference between consecutive values is within the tolerance
    do
    {
        prev = X0; // Store the previous value of X0
        X0 = X0 - f(X0) / fdash(X0); // Update X0 using the Newton-Raphson formula
    } while (fabs(X0 - prev) > Tol); // Continue until the difference is less than the tolerance

    root = X0; // The root of the equation is the final value of X0
    cout << "Roots of the given equation: " << root; // Print the root
    return root; // Return the root
}
