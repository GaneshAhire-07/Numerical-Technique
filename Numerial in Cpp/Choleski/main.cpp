#include <iostream>
#include "Choleski.hpp"
using namespace std;

int main()
{
    Matrix matrix;

    string MatrixFile = "input.txt";
    matrix.ReadMatrixFromFile(MatrixFile);
    matrix.PrintMatrix();
    matrix.CholeskyDecomposition();

    // matrix.PrintMatrix();

    return 0;
}