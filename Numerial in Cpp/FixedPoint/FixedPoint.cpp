#include"FixedPoint.hpp"
#include <bits/stdc++.h>
using namespace std;
double FixedPoint::Function(double x)
{
    return ((x*x*x+2)/4.0);
}

double FixedPoint::RootFinding()
{
   do
    {
        prev = X0;
        X0 =  Function(X0);
    } while (fabs(X0 - prev) < tol);

    root = X0;
    cout << "roots of given equation: " << root;
}