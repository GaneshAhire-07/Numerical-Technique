#include "jacobi.hpp"
#include <iostream>

using namespace std;


int main() {
    Matrix matrix;
    matrix.ReadFiles();

    if (matrix.makeDiagonallyDominant()) {
        matrix.Jacobi();
    } else {
        cerr << "The matrix is not suitable for the Jacobi method." << endl;
    }

    return 0;
}
