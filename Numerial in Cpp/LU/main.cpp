#include <iostream>
#include "LU.hpp"
using namespace std;

int main()
{
    Matrix matrix;

    string MatrixFile = "input.txt";
    matrix.ReadMatrixFromFile(MatrixFile);
    matrix.PrintMatrix();
    matrix.LUDecomposition();

    matrix.LUDecompositionAndGetDeterminant();

    return 0;
}