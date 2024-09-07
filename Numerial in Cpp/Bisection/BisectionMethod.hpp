#include <bits/stdc++.h>
using namespace std;

class BisectionMethod {
private:
    function<double(double)> func; // Function to find root of
    double tolerance; // Desired tolerance for root finding
    int maxIterations; // Maximum number of iterations

public:
    BisectionMethod(function<double(double)> f, double tol, int maxIter); // Constructor
    void findRoot(double a, double b); // Method to find the root in interval [a, b]
};