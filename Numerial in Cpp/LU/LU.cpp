#include <bits/stdc++.h>
#include "LU.hpp"
using namespace std;

// Default constructor for the Matrix class
Matrix::Matrix()
{
    cout << "L U Decomposition " << endl;
}

// Parameterized constructor for the Matrix class
Matrix::Matrix(int r, int c)
{
    Rows = r; // Set number of rows
    Cols = c; // Set number of columns
    Mat.resize(Rows, vector<double>(Cols, 0)); // Initialize matrix with zeros
    init(); // Initialize the matrix elements
}

// Initialize the matrix elements to zero
void Matrix::init()
{
    for (size_t i = 0; i < Rows; i++)
    {
        for (size_t j = 0; j < Cols; j++)
        {
            Mat[i][j] = 0; // Set each element to zero
        }
    }
}

// Read matrix from a file
void Matrix::ReadMatrixFromFile(string MatrixFile)
{
    ifstream Fin; // Input file stream
    Fin.open(MatrixFile); // Open the file
    if (!Fin) // Check if the file was successfully opened
    {
        cout << "Error: Unable to open matrix file!" << endl;
        exit(1); // Exit the program with an error
    }

    Fin >> Rows >> Cols; // Read the number of rows and columns
    Mat = vector<vector<double>>(Rows, vector<double>(Cols)); // Resize the matrix

    cout << "\nrows : " << Rows << "\tcols : " << Cols << endl;

    for (int i = 0; i < Rows; i++) // Read the matrix elements from the file
    {
        for (int j = 0; j < Cols; j++)
        {
            Fin >> Mat[i][j];
        }
    }
    Fin.close(); // Close the file
}

// Print the matrix to the console
void Matrix::PrintMatrix()
{
    for (int i = 0; i < Mat.size(); i++) // Loop through rows
    {
        for (int j = 0; j < Mat[i].size(); j++) // Loop through columns
        {
            cout << Mat[i][j] << " "; // Print each element
        }
        cout << endl; // Newline after each row
    }
}

// Calculate the element of the upper triangular matrix
double Matrix::GetUpperTriangular(Matrix L, Matrix U, int rowIndex, int colIndex)
{
    if (rowIndex == colIndex) // Diagonal elements of U
    {
        return 1;
    }
    double ans = Mat[rowIndex][colIndex]; // Start with the original matrix element
    double temp = 0;
    for (int k = 0; k < rowIndex; k++) // Sum the products of corresponding elements in L and U
    {
        temp += L.Mat[rowIndex][k] * U.Mat[k][colIndex];
    }
    ans = (ans - temp) / L.Mat[rowIndex][rowIndex]; // Adjust the element for the current row
    return ans;
}

// Calculate the element of the lower triangular matrix
double Matrix::GetLowerTriangular(Matrix L, Matrix U, int rowIndex, int colIndex)
{
    double ans = Mat[rowIndex][colIndex]; // Start with the original matrix element
    double temp = 0.0;
    for (int k = 0; k < colIndex; k++) // Sum the products of corresponding elements in L and U
    {
        temp += L.Mat[rowIndex][k] * U.Mat[k][colIndex];
    }
    ans = ans - temp; // Adjust the element for the current row
    return ans;
}

// Perform LU Decomposition of the matrix
void Matrix::LUDecomposition()
{
    Matrix L(this->Rows, this->Cols); // Create L matrix
    Matrix U(this->Rows, this->Cols); // Create U matrix

    for (int rowIndex = 0; rowIndex < Rows; rowIndex++) // Loop through rows
    {
        for (int colIndex = 0; colIndex < Cols; colIndex++) // Loop through columns
        {
            if (colIndex >= rowIndex)
                L.Mat[colIndex][rowIndex] = GetLowerTriangular(L, U, colIndex, rowIndex); // Fill L
        }

        for (int colIndex = 0; colIndex < Cols; colIndex++) // Loop through columns
        {
            if (colIndex >= rowIndex)
                U.Mat[rowIndex][colIndex] = GetUpperTriangular(L, U, rowIndex, colIndex); // Fill U
        }
    }

    cout << "Lower triangular matrix:" << endl;
    L.PrintMatrix(); // Print L matrix
    cout << "\nUpper triangular matrix:" << endl;
    U.PrintMatrix(); // Print U matrix
}

// Perform LU Decomposition and calculate the determinant of the matrix
double Matrix::LUDecompositionAndGetDeterminant()
{
    Matrix L(Rows, Cols); // Create L matrix
    Matrix U(Rows, Cols); // Create U matrix

    for (int rowIndex = 0; rowIndex < Rows; rowIndex++) // Loop through rows
    {
        for (int colIndex = 0; colIndex < Cols; colIndex++) // Loop through columns
        {
            if (colIndex >= rowIndex)
                L.Mat[colIndex][rowIndex] = GetLowerTriangular(L, U, colIndex, rowIndex); // Fill L
        }

        for (int colIndex = 0; colIndex < Cols; colIndex++) // Loop through columns
        {
            if (colIndex >= rowIndex)
                U.Mat[rowIndex][colIndex] = GetUpperTriangular(L, U, rowIndex, colIndex); // Fill U
        }
    }

    double detL = 1.0; // Initialize determinant of L
    double detU = 1.0; // Initialize determinant of U

    for (int i = 0; i < Rows; i++)
    {
        detL *= U.Mat[i][i]; // Multiply the diagonal elements of U
        detU *= L.Mat[i][i]; // Multiply the diagonal elements of L
    }

    cout << "Determinant of Lower : " << detL << endl;
    cout << "Determinant Of Upper : " << detU << endl;

    double detA = detL * detU; // Determinant of A is the product of determinants of L and U
    cout << "Determinant of A : " << detA << endl;
}
