#include "seidel.hpp"
#include <bits/stdc++.h>
using namespace std;

double tol = 0.0001;

string fileName_L = "output.txt";
string fileName_R = "input.txt";

// Default constructor implementation
Matrix::Matrix()
{
    // Initialize rows and cols to zero
    rows = 0;
    cols = 0;
    mat = nullptr; // Ensure mat is initialized to nullptr
}

// Constructor implementation
Matrix::Matrix(int r, int c) : rows(r), cols(c)
{
    // Dynamically allocate memory for the array
    mat = new double *[rows];
    for (int i = 0; i < rows; i++)
    {
        mat[i] = new double[cols]();
    }
}

// Destructor implementation
Matrix::~Matrix()
{
    // Deallocate memory for the matrix
    for (int i = 0; i < rows; i++)
    {
        delete[] mat[i];
    }
    delete[] mat;
}

// Function to read matrix data from files
void Matrix::ReadFiles()
{
    ifstream fin;
    fin.open("input.txt");

    if (!fin.is_open())
    {
        cerr << "Error: Unable to open input file." << endl;
        exit(1);
    }

    fin >> rows >> cols;

    // Dynamically allocate memory for the array
    mat = new double *[rows];
    for (int i = 0; i < rows; i++)
    {
        mat[i] = new double[cols]();
    }
    cout << "\n Rows: " << rows << "\t cols: " << cols << endl;

    // Read the matrix from the first file
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < (cols - 1); j++)
        {
            fin >> mat[i][j];
        }
    }
    fin.close();

    // Open the second file
    fin.open("temp.txt");

    if (!fin.is_open())
    {
        cerr << "Error: Unable to open temp file." << endl;
        exit(1);
    }

    // Read the matrix elements from the second file
    for (int i = 0; i < rows; i++)
    {
        fin >> mat[i][cols - 1];
    }
    fin.close();

    cout << "The final matrix after combining elements from both files:" << endl;
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            cout << mat[i][j] << "   ";
        }
        cout << endl;
    }
}

// Function to get the index of the diagonally dominant row
int Matrix::getDDRowAt(int r)
{
    for (int i = r + 1; i < rows; i++)
    {
        double sum = 0.0;
        for (int c = 0; c < cols - 1; c++)
        {
            if (r != c)
            {
                sum += fabs(mat[i][c]);
            }
        }
        if (fabs(mat[r][r]) >= sum)
        {
            return i;
        }
    }
    cout << "No diagonally dominant row found after row " << r << endl;
    return -1;
}

// Function to swap two rows in the matrix
void Matrix::SwapRows(int r1, int r2)
{
    for (int i = 0; i < cols; i++)
    {
        double tmp = mat[r1][i];
        mat[r1][i] = mat[r2][i];
        mat[r2][i] = tmp;
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
                sum += fabs(mat[r][c]);
            }
        }
        if (fabs(mat[r][r]) < sum)
        {
            return false;
        }
    }
    return true;
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
            cout << "Row " << r + 1 << " can be made diagonally dominant." << endl;
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

    cout << "Diagonally dominant matrix:" << endl;
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            cout << mat[i][j] << "   ";
        }
        cout << endl;
    }

    return true;
}

// Function to perform Gauss-Seidel iteration
void Matrix::GaussSeidel()
{
    double *ans = new double[cols - 1]();
    double *prev = new double[cols - 1]();

    ofstream fout;
    fout.open("output.txt");

    if (!fout.is_open())
    {
        cerr << "Error: Unable to open output file." << endl;
        exit(1);
    }

    while (true)
    {
        bool done = true;

        for (int r = 0; r < rows; r++)
        {
            double sum = 0.0;

            for (int c = 0; c < cols - 1; c++)
            {
                if (r != c)
                {
                    sum += mat[r][c] * ans[c];
                }
            }

            ans[r] = (mat[r][cols - 1] - sum) / mat[r][r];

            if (fabs(ans[r] - prev[r]) > tol)
            {
                done = false;
            }

            prev[r] = ans[r];
        }

        if (done)
        {
            break;
        }
    }

    cout << "Solution:" << endl;
    for (int i = 0; i < cols - 1; i++)
    {
        cout << "x" << i + 1 << " = " << ans[i] << endl;
        fout << ans[i] << endl;
    }

    fout.close();
    delete[] ans;
    delete[] prev;
}
