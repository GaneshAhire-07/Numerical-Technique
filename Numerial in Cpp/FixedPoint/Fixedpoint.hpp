#include <bits/stdc++.h>
using namespace std;
class FixedPoint
{
public:
    double X0 = 0;
    // cout << "enter the x0 value: " << Fx.X0 << endl;
    double prev;
    double root;
    double tol = 1e-6;
    int maxIter = 100;

    double RootFinding();
    double Function(double);
};