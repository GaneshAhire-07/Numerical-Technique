#include <bits/stdc++.h>

class Matrix
{
public:
    int rows;
    int cols;
    int index;
    double **mat;
    void printMatrix();

    Matrix();
    Matrix(int r, int c);
    ~Matrix(); // Destructor Declaration

    void ReadFiles();
    void Jacobi();
    bool makeDiagonallyDominant();
    bool isDiagonallyDominant();
    void SwapRows(int r1, int r2);
    int getDDRowAt(int r);
};