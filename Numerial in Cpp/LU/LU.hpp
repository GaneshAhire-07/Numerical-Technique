#include <bits/stdc++.h>
using namespace std;

class Matrix
{
public:
    int Rows, Cols;
    vector<vector<double>> Mat;

    Matrix();
    Matrix(int, int);
    void init();
    void ReadMatrixFromFile(string);
    void PrintMatrix();
    void LUDecomposition();
    double GetLowerTriangular(Matrix, Matrix, int, int);
    double GetUpperTriangular(Matrix, Matrix, int, int);
    double CalculateDeterminant();
    double LUDecompositionAndGetDeterminant();
};
