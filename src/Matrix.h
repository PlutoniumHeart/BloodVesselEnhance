#ifndef MATRIX_H
#define MATRIX_H


#include <cmath>
#include <iostream>
#include <cassert>


class Matrix
{
public:
    Matrix();
    Matrix(int row, int col);
    virtual ~Matrix();

    bool Init(int row, int col);
    bool JacobiEigenv(double* dblEigenValue, Matrix& mtxEigenVector, int iMaxIt = 60, double eps = 0.000001);
    double* GetData() const;
    bool SetElement(int row, int col, double value);
    double GetElement(int row, int col) const;
protected:
    int m_iColumns;
    int m_iRows;
    double* m_pdData;
};


#endif // !MATRIX_H
