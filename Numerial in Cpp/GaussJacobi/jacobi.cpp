#include "jacobi.hpp"
#include <bits/stdc++.h>
using namespace std;

// Parameterized constructor for the Matrix class
Matrix::Matrix(int r, int c) : rows(r), cols(c)
{
    mat = new double *[rows]; // Allocate memory for rows
    for (int i = 0; i < rows; i++)
    {
        mat[i] = new double[cols](); // Allocate memory for columns and initialize to zero
    }
}

// Destructor for the Matrix class
Matrix::~Matrix()
{
    for (int i = 0; i < rows; i++)
    {
        delete[] mat[i]; // Deallocate memory for each row
    }
    delete[] mat; // Deallocate memory for the array of pointers
}

// Function to read matrix data from files
void Matrix::ReadFiles()
{
    ifstream fin; // Input file stream

    // Open the input file and check for errors
    fin.open("input.txt");
    if (!fin.is_open())
    {
        cerr << "Error: Unable to open input file." << endl;
        exit(1);
    }

    // Read the number of rows and columns
    fin >> rows >> cols;

    // Allocate memory for the matrix with the new dimensions
    mat = new double *[rows];
    for (int i = 0; i < rows; i++)
    {
        mat[i] = new double[cols]();
    }

    // Read the matrix elements except the last column
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols - 1; j++)
        {
            fin >> mat[i][j];
        }
    }
    fin.close(); // Close the input file

    // Open the temp file and check for errors
    fin.open("temp.txt");
    if (!fin.is_open())
    {
        cerr << "Error: Unable to open temp file." << endl;
        exit(1);
    }

    // Read the last column elements from the temp file
    for (int i = 0; i < rows; i++)
    {
        fin >> mat[i][cols - 1];
    }
    fin.close(); // Close the temp file

    cout << "Combined matrix:" << endl;
    printMatrix(); // Print the combined matrix
}

// Function to print the matrix
void Matrix::printMatrix()
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            cout << mat[i][j] << "   ";
        }
        cout << endl;
    }
}

// Function to check if the matrix is diagonally dominant
bool Matrix::isDiagonallyDominant()
{
    for (int r = 0; r < rows; r++)
    {
        double sum = 0;
        for (int c = 0; c < cols - 1; c++)
        {
            if (r != c)
            {
                sum += fabs(mat[r][c]); // Sum the absolute values of non-diagonal elements
            }
        }
        if (fabs(mat[r][r]) < sum) // Check if diagonal element is greater than the sum
        {
            return false;
        }
    }
    return true;
}

// Function to find a row that can be swapped to make the matrix diagonally dominant
int Matrix::getDDRowAt(int r)
{
    for (int i = r + 1; i < rows; i++)
    {
        double sum = 0.0;
        for (int c = 0; c < cols - 1; c++)
        {
            if (r != c)
            {
                sum += fabs(mat[i][c]); // Sum the absolute values of non-diagonal elements
            }
        }
        if (fabs(mat[r][r]) >= sum) // Check if diagonal element is greater than the sum
        {
            return i;
        }
    }
    return -1; // Return -1 if no suitable row is found
}

// Function to swap two rows of the matrix
void Matrix::SwapRows(int r1, int r2)
{
    for (int i = 0; i < cols; i++)
    {
        double tmp = mat[r1][i];
        mat[r1][i] = mat[r2][i];
        mat[r2][i] = tmp;
    }
}

// Function to make the matrix diagonally dominant if possible
bool Matrix::makeDiagonallyDominant()
{
    if (isDiagonallyDominant())
    {
        cout << "Matrix is already diagonally dominant." << endl;
        return true;
    }

    for (int r = 0; r < rows; r++)
    {
        double sum = 0;
        for (int c = 0; c < cols - 1; c++)
        {
            if (r != c)
            {
                sum += fabs(mat[r][c]);
            }
        }
        if (sum > fabs(mat[r][r]))
        {
            int index = getDDRowAt(r);
            if (index != -1)
            {
                SwapRows(r, index);
                cout << "Swapped rows " << r + 1 << " and " << index + 1 << endl;
            }
            else
            {
                cout << "Not possible to make the matrix diagonally dominant." << endl;
                return false;
            }
        }
    }
    cout << "Matrix after making it diagonally dominant:" << endl;
    printMatrix();
    return true;
}

// Function to perform the Jacobi iterative method
void Matrix::Jacobi()
{
    double *ans = new double[cols - 1](); // Solution vector
    double *prev = new double[cols - 1](); // Previous iteration values
    int maxIter = 100; // Maximum number of iterations
    double tol = 0.0001; // Tolerance for convergence

    ofstream fout("output.txt"); // Output file stream
    if (!fout.is_open())
    {
        cerr << "Error: Unable to open output file." << endl;
        exit(1);
    }

    int iter = 0; // Iteration counter
    while (true)
    {
        iter++;
        bool done = true;

        for (int r = 0; r < rows; r++) // Loop through rows
        {
            double sum = 0.0;
            for (int c = 0; c < cols - 1; c++) // Calculate the sum of the product of non-diagonal elements and their previous iteration values
            {
                if (r != c)
                {
                    sum += mat[r][c] * prev[c];
                }
            }
            ans[r] = (mat[r][cols - 1] - sum) / mat[r][r]; // Calculate the current iteration value
            if (fabs(ans[r] - prev[r]) > tol) // Check if the difference is greater than the tolerance
            {
                done = false;
            }
        }

        for (int i = 0; i < cols - 1; i++)
        {
            prev[i] = ans[i]; // Update previous iteration values
        }

        if (done || iter >= maxIter) // Break the loop if done or maximum iterations reached
        {
            break;
        }
    }

    cout << "Solution after " << iter << " iterations:" << endl;
    for (int i = 0; i < cols - 1; i++)
    {
        cout << "x" << i + 1 << " = " << ans[i] << endl; // Print the solution
        fout << ans[i] << endl; // Write the solution to the output file
    }

    fout.close(); // Close the output file
    delete[] ans; // Deallocate memory for the solution vector
    delete[] prev; // Deallocate memory for the previous iteration values
}
