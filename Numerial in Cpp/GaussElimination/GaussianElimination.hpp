#include <bits/stdc++.h>
using namespace std;
class GaussianElimination {
    int row;
    int col;
    vector<vector<double>> mat;

public:
    GaussianElimination(vector<vector<double>>& matrix);
    void printMatrix();
    void rowReduction();
    void backSubstitution();
    void gaussianElimination();
};
