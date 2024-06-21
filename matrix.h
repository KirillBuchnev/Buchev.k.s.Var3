#ifndef MATRIX
#define MATRIX
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <chrono>
using namespace std;
struct Matrix {
    int rows;
    int cols;
    vector<vector<double>> data;
};
Matrix createMatrix(int rows, int cols);
void fillMatrix(Matrix& matrix);
double calculateDeterminant(const Matrix& matrix);
Matrix getInverseMatrix(const Matrix& matrix);
double getValidInput(const string& prompt);

vector<double> solveCramer(const Matrix& A, const vector<double>& b);

vector<double> solveInverseMatrix(const Matrix& A, const vector<double>& b);

void showMenu();

#endif