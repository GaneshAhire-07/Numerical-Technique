#include "BisectionMethod.hpp"
#include <bits/stdc++.h>

// Define the function for which we want to find the root
double sampleFunction(double x)
{
    return x * x * x - x - 1; // Example function: f(x) = x^2 - 4(0,5),f(x) = x^3 - x - 1(1,2)
}

int main()
{
    // Create a BisectionMethod object with the function, tolerance, and max iterations
    BisectionMethod bisection(sampleFunction, 1e-5, 100); // tolerance = 10^-5, max iterations = 100
    bisection.findRoot(1, 2);// Initial interval [a, b]

    return 0;
}