#include "Choleski.hpp" // Include the header file for the Choleski class definitions
#include <bits/stdc++.h> // Include all standard libraries (useful for competitive programming)
using namespace std; // Use the standard namespace

// Default constructor for the Matrix class
Matrix::Matrix()
{
    cout << "Given The Matrix" << endl; // Print a message when the matrix is created
}

// Parameterized constructor for the Matrix class
Matrix::Matrix(int r, int c)
{
    Rows = r; // Set the number of rows
    Cols = c; // Set the number of columns
    Mat.resize(Rows, vector<double>(Cols, 0)); // Resize the matrix and initialize all elements to 0
}

// Method to read a matrix from a file
void Matrix::ReadMatrixFromFile(string MatrixFile)
{
    ifstream Fin; // Create an input file stream object
    Fin.open(MatrixFile); // Open the matrix file
    if (!Fin) // Check if the file opened successfully
    {
        cout << "Unable to open matrix file!" << endl; // Print an error message if the file can't be opened
       // exit(1); // Exit the program with an error code
    }
    Fin >> Rows >> Cols; // Read the number of rows and columns from the file
    this->Mat.resize(Rows, vector<double>(Cols, 0)); // Resize the matrix and initialize all elements to 0

    cout << "rows : " << Rows << "\tcols : " << Cols << endl; // Print the number of rows and columns

    for (int i = 0; i < this->Rows; i++) // Loop through each row
        for (int j = 0; j < this->Cols; j++) // Loop through each column
            Fin >> Mat[i][j]; // Read each element from the file into the matrix

    Fin.close(); // Close the file
}

// Method to print the matrix to the console
void Matrix::PrintMatrix()
{
    for (int i = 0; i < Mat.size(); i++) // Loop through each row
    {
        for (int j = 0; j < Mat[i].size(); j++) // Loop through each column
        {
            cout << Mat[i][j] << " "; // Print the element followed by a space
        }
        cout << endl; // Print a newline after each row
    }
}

// Method to check if the matrix is symmetric
bool Matrix::isSymmetric()
{
    for (int i = 0; i < Rows; i++) // Loop through each row
    {
        for (int j = 0; j < Cols; j++) // Loop through each column
        {
            if (i > j) // Only check elements below the diagonal
            {
                if (Mat[i][j] != Mat[j][i]) // Check if the element is equal to its transpose
                    return false; // Return false if any element is not equal to its transpose
            }
        }
    }
    return true; // Return true if all elements are symmetric
}

// Method to get the diagonal element for Cholesky decomposition
double Matrix::getDiagonalElement(Matrix LowerTriangularMatrix, int RowIndex)
{
    double Ans1 = Mat[RowIndex][RowIndex]; // Get the diagonal element from the original matrix
    double temp = 0.0; // Initialize a temporary variable
    for (int k = 0; k < RowIndex; k++) // Loop through each element in the current row up to the diagonal
    {
        temp += LowerTriangularMatrix.Mat[RowIndex][k] * LowerTriangularMatrix.Mat[RowIndex][k]; // Sum the squares of the elements in the lower triangular matrix
    }
    double Ans2 = sqrt(Ans1 - temp); // Calculate the diagonal element for the lower triangular matrix
    return Ans2; // Return the calculated diagonal element
}

// Method to get the off-diagonal element for Cholesky decomposition
double Matrix::getLowerTriangularMatrix(Matrix LowerTriangularMatrix, int RowIndex, int ColIndex)
{
    double Ans1 = Mat[RowIndex][ColIndex]; // Get the off-diagonal element from the original matrix
    double temp = 0.0; // Initialize a temporary variable
    double Ans2 = LowerTriangularMatrix.Mat[ColIndex][ColIndex]; // Get the diagonal element from the lower triangular matrix
    for (int k = 0; k < ColIndex; k++) // Loop through each element in the current row up to the column index
    {
        temp += LowerTriangularMatrix.Mat[RowIndex][k] * LowerTriangularMatrix.Mat[ColIndex][k]; // Sum the products of corresponding elements in the lower triangular matrix
    }
    double Ans = (Ans1 - temp) / Ans2; // Calculate the off-diagonal element for the lower triangular matrix
    return Ans; // Return the calculated off-diagonal element
}

// Method to get the transpose of the lower triangular matrix
Matrix Matrix::LowerMatrixTranspose(Matrix LowerTriangularMatrix)
{
    Matrix LowerTriangularMatrixTranspose(this->Rows, this->Cols); // Create a new matrix for the transpose
    for (int RowIndex = 0; RowIndex < Rows; RowIndex++) // Loop through each row
    {
        for (int ColIndex = 0; ColIndex < Rows; ColIndex++) // Loop through each column
        {
            LowerTriangularMatrixTranspose.Mat[RowIndex][ColIndex] = LowerTriangularMatrix.Mat[ColIndex][RowIndex]; // Set the element in the transpose matrix
        }
    }
    return LowerTriangularMatrixTranspose; // Return the transpose matrix
}

// Method to perform Cholesky decomposition
void Matrix::CholeskyDecomposition()
{
    if (isSymmetric()) // Check if the matrix is symmetric
    {
        cout << "The matrix is symmetric." << endl; // Print a message if the matrix is symmetric
    }
    else
    {
        cout << "The matrix is not symmetric." << endl; // Print a message if the matrix is not symmetric
        return; // Return without performing the decomposition
    }

    Matrix LowerTriangularMatrix(Rows, Cols); // Create a matrix for the lower triangular matrix
    int RowIndex, ColIndex; // Declare row and column indices

    for (RowIndex = 0; RowIndex < Rows; RowIndex++) // Loop through each row
    {
        for (ColIndex = 0; ColIndex <= RowIndex; ColIndex++) // Loop through each column up to the current row index
        {
            if (RowIndex == ColIndex) // If it's a diagonal element
            {
                LowerTriangularMatrix.Mat[RowIndex][RowIndex] = getDiagonalElement(LowerTriangularMatrix, RowIndex); // Calculate and set the diagonal element
            }
            else if (ColIndex < RowIndex) // If it's an off-diagonal element
            {
                LowerTriangularMatrix.Mat[RowIndex][ColIndex] = getLowerTriangularMatrix(LowerTriangularMatrix, RowIndex, ColIndex); // Calculate and set the off-diagonal element
            }
        }
    }

    cout << "Lower Triangular Matrix:" << endl; // Print a message
    LowerTriangularMatrix.PrintMatrix(); // Print the lower triangular matrix

    Matrix LowerTriangularMatrixTranspose = LowerMatrixTranspose(LowerTriangularMatrix); // Get the transpose of the lower triangular matrix
    cout << "Transpose of Lower Triangular Matrix:" << endl; // Print a message
    LowerTriangularMatrixTranspose.PrintMatrix(); // Print the transpose of the lower triangular matrix
}
