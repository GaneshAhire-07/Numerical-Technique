#include "GaussianElimination.hpp"
#include<bits/stdc++.h>
using namespace std;

int main() {
    ifstream inFile("input.txt");
    if (!inFile.is_open()) {
        cout << "Failed to open input file." << endl;
        return 1;
    }
    int rows, cols;
    inFile >> rows >> cols;

    vector<vector<double>> matrix(rows, vector<double>(cols));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            inFile >> matrix[i][j];
        }
    }
    inFile.close();

    GaussianElimination sol(matrix);
    sol.gaussianElimination();

    return 0;
}
