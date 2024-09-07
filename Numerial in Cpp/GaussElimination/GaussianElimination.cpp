#include <bits/stdc++.h>
using namespace std;
class GaussianElimination {
private:
    vector<vector<double>> mat;
    int row;
    int col;

    void rowReduction();
    void backSubstitution();

public:
    GaussianElimination(vector<vector<double>> &matrix);
    void printMatrix();
    void gaussianElimination();
};

// Constructor implementation
GaussianElimination::GaussianElimination(vector<vector<double>> &matrix)
{
    mat = matrix;
    row = matrix.size();
    col = matrix[0].size();
}

// Function to print the matrix
void GaussianElimination::printMatrix() {
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            cout << " " << mat[i][j] << " ";
        }
        cout << "\n";
    }
}

// Function for row reduction
void GaussianElimination::rowReduction() {
    for (int r = 0; r < row; r++) {
        double pivot = mat[r][r];
        if (pivot == 0) {
            for (int i = r + 1; i < row; i++) {
                if (mat[i][r] != 0) {     
                    swap(mat[r], mat[i]);
                    pivot = mat[r][r];
                    break;
                }
            }
        }
        if (pivot == 0) continue;
        for (int c = r; c < col; c++) {
            mat[r][c] /= pivot;
        }
        for (int i = 0; i < row; i++) {
            if (i != r) {
                double factor = mat[i][r];
                for (int c = r; c < col; c++) {
                    mat[i][c] -= factor * mat[r][c];
                }
            }
        }
    }
}
//function for column reduction
// Function for back substitution
void GaussianElimination::backSubstitution() {
    vector<double> ans(row);
    for (int r = row - 1; r >= 0; r--) {
        double sum = mat[r][col - 1];
        for (int c = r + 1; c < col - 1; c++) {
            sum -= mat[r][c] * ans[c];
        }
        ans[r] = sum;
    }
    ofstream outFile("output.txt");
    if (outFile.is_open()) {
        for (int i = 0; i < row; i++) {
            outFile << fixed << setprecision(5) << ans[i] << endl;
        }
        outFile.close();
    } else {
        cout << "Failed to open output file." << endl;
    }
}

// Function to perform Gaussian elimination
void GaussianElimination::gaussianElimination() {
    rowReduction();
    backSubstitution();
}

