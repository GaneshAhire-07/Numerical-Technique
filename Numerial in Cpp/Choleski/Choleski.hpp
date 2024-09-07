#include <bits/stdc++.h>
using namespace std;

class Matrix
{
public:
    int Rows, Cols;
    vector<vector<double>> Mat;

    Matrix();
    Matrix(int, int);

    void ReadMatrixFromFile(string);
    void PrintMatrix();
    bool isSymmetric();
    double getDiagonalElement(Matrix, int);
    double getLowerTriangularMatrix(Matrix, int, int);
    void CholeskyDecomposition();
    Matrix LowerMatrixTranspose(Matrix);

    // void getLowerTriangularMatrix();
};
