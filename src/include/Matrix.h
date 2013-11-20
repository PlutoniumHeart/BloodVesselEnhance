#ifndef MATRIX_H
#define MATRIX_H


#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <cassert>


class Matrix
{
public:
    Matrix();
    Matrix(int row, int col);
    Matrix(const Matrix& other);
    virtual ~Matrix();

    bool Init(int row, int col);
    bool JacobiEigenv(double* dblEigenValue, Matrix& mtxEigenVector, int iMaxIt = 60, double eps = 0.000001);
    int GetColumns() const;
    int GetRows() const;
    void SetData(double* value);
    double* GetData() const;
    bool SetElement(int row, int col, double value);
    double GetElement(int row, int col) const;

    Matrix& operator=(const Matrix& other);
    bool operator==(const Matrix& other) const;
    bool operator!=(const Matrix& other) const;
    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;
    Matrix operator*(double value) const;
    Matrix operator*(const Matrix& other) const;
protected:
    int m_iColumns;
    int m_iRows;
    double* m_pdData;
};


#endif // MATRIX_H
