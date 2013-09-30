#include "Matrix.h"


/* ===========================================================================================
// Basic constructor
============================================================================================*/
Matrix::Matrix()
    : m_iColumns(1)
    , m_iRows(1)
    , m_pdData(NULL)
{
    bool bSuccess = Init(m_iRows, m_iColumns);
    assert(bSuccess);
}


/* ===========================================================================================
// Matrix constructor
============================================================================================*/
Matrix::Matrix(int row, int col)
    : m_iColumns(col)
    , m_iRows(row)
    , m_pdData(NULL)
{
    bool bSuccess = Init(m_iRows, m_iColumns);
    assert(bSuccess);
}


/* ===========================================================================================
// Matrix copy constructor
============================================================================================*/
Matrix::Matrix(const Matrix& other)
    : m_iColumns(other.GetColumns())
    , m_iRows(other.GetRows())
    , m_pdData(NULL)
{
    bool bSuccess = Init(m_iRows, m_iColumns);
    assert(bSuccess);

    memcpy(m_pdData, other.m_pdData, sizeof(double)*m_iColumns*m_iRows);
}


/* ===========================================================================================
// Matrix destructor
============================================================================================*/
Matrix::~Matrix()
{
    if(m_pdData != NULL)
    {
        delete [] m_pdData;
        m_pdData = NULL;
    }
}


/* ===========================================================================================
// Initialize the matrix
============================================================================================*/
bool Matrix::Init(int row, int col)
{
    if(m_pdData!=NULL)
    {
        delete [] m_pdData;
        m_pdData = NULL;
    }

    m_iRows = row;
    m_iColumns = col;
    int size = row*col;

    if(size<0)
    {
        return false;
    }

    m_pdData = new double[size];

    if(m_pdData == NULL)
    {
        return false;
    }

    memset(m_pdData, 0, sizeof(double)*size);

    return true;
}


/* ===========================================================================================
// Jocobi method for calculating eigen values of a real symmetrix matrix
============================================================================================*/
bool Matrix::JacobiEigenv(double* dblEigenValue, Matrix& mtxEigenVector, int nMaxIt, double eps)
{
    int i = 0, j = 0, p = 0, q = 0, u = 0, w = 0, t = 0, s = 0, l = 0;
    double fm = 0.0, cn = 0.0, sn = 0.0, omega = 0.0, x = 0.0, y = 0.0, d = 0.0;

    if(!mtxEigenVector.Init(m_iColumns, m_iRows))
    {
        return false;
    }

    l = 1;

    for(i=0;i<=m_iColumns-1;i++)
    {
        mtxEigenVector.m_pdData[i*m_iColumns]=1.0;
        for(j=0;j<=m_iColumns-1;j++)
        {
            if(i!=j)
            {
                mtxEigenVector.m_pdData[i*m_iColumns] = 0.0;
            }
        }
    }

    while(true)
    {
        fm = 0.0;
        for(i=0;i<=m_iColumns-1;i++)
        {
            for(j=0;j<=i-1;j++)
            {
                d=fabs(m_pdData[i*m_iColumns+j]);
                if((i!=j) && (d>fm))
                {
                    fm=d;
                    p=i;
                    q=j;
                }
            }
        }

        // Error smaller than the eps, use values on the diagnoal as Eigen value and return true
        if(fm<eps)
        {
            for(i=0;i<m_iColumns;++i)
            {
                dblEigenValue[i] = GetElement(i, i);
            }
            return true;
        }

        if(l>nMaxIt)
        {
            return false;
        }

        l++;
        u = p*m_iColumns+q;
        w = p*m_iColumns+p;
        t = q*m_iColumns+p;
        s = q*m_iColumns+q;
        x = -m_pdData[u];
        y = (m_pdData[s]-m_pdData[w])/2.0;
        omega = x/sqrt(x*x + y*y);

        if(y<0.0)
        {
            omega = -omega;
        }

        sn = 1.0 + sqrt(1.0-omega*omega);
        sn = omega/sqrt(2.0*sn);
        cn = sqrt(1.0 - sn*sn);
        fm = m_pdData[w];
        m_pdData[w] = fm*cn*cn + m_pdData[s]*sn*sn + m_pdData[u]*omega;
        m_pdData[s] = fm*sn*sn + m_pdData[s]*cn*cn - m_pdData[u]*omega;
        m_pdData[u] = 0.0;
        m_pdData[t] = 0.0;

        for(j=0;j<=m_iColumns-1;j++)
        {
            if((j!=p) && (j!=q))
            {
                u = p*m_iColumns+j;
                w = q*m_iColumns+j;
                fm = m_pdData[u];
                m_pdData[u] = fm*cn + m_pdData[w]*sn;
                m_pdData[w] = fm*sn + m_pdData[w]*cn;
            }
        }

        for(i=0;i<=m_iColumns-1;i++)
        {
            if((i!=p) && (i!=q))
            {
                u = i*m_iColumns+p;
                w = i*m_iColumns+q;
                fm = m_pdData[u];
                m_pdData[u] = fm*cn + m_pdData[w]*sn;
                m_pdData[w] = -fm*sn + m_pdData[w]*cn;
            }
        }

        for(i=0;i<m_iColumns-1;i++)
        {
            u = i*m_iColumns+p;
            w = i*m_iColumns+q;
            fm = mtxEigenVector.m_pdData[u];
            mtxEigenVector.m_pdData[u] = fm*cn + mtxEigenVector.m_pdData[w]*sn;
            mtxEigenVector.m_pdData[w] = -fm*sn + mtxEigenVector.m_pdData[w]*cn;
        }
    }

    for(i=0;i<m_iColumns;++i)
    {
        dblEigenValue[i] = GetElement(i,i);
    }

    return true;
}


/* ===========================================================================================
// Operator '=' overload
============================================================================================*/
Matrix& Matrix::operator=(const Matrix& other)
{
    if(&other!=this)
    {
        bool bSuccess = Init(other.GetRows(), other.GetColumns());
        assert(bSuccess);

        memcpy(m_pdData, other.m_pdData, sizeof(double)*m_iColumns*m_iRows);
    }

    return *this;
}


/* ===========================================================================================
// Operator '==' overload
============================================================================================*/
bool Matrix::operator==(const Matrix& other) const
{
    int i = 0, j = 0;
    if(m_iColumns!=other.GetColumns() ||
       m_iRows!=other.GetRows())
    {
        return false;
    }

    for(i=0;i<m_iRows;++i)
    {
        for(j=0;j<m_iColumns;++j)
        {
            if(GetElement(i, j) != other.GetElement(i, j))
            {
                return false;
            }
        }
    }

    return true;
}


/* ===========================================================================================
// Operator '!=' overload
============================================================================================*/
bool Matrix::operator!=(const Matrix& other) const
{
    return !(*this == other);
}


/* ===========================================================================================
// Operator '+' overload
============================================================================================*/
Matrix Matrix::operator+(const Matrix& other) const
{
    int i = 0, j = 0;
    assert(m_iColumns == other.GetColumns() &&
           m_iRows == other.GetRows());

    Matrix result(*this);

    for(i=0;i<m_iRows;++i)
    {
        for(j=0;j<m_iColumns;++j)
        {
            result.SetElement(i, j, result.GetElement(i,j) + other.GetElement(i,j));
        }
    }

    return result;
}


/* ===========================================================================================
// Operator '-' overload
============================================================================================*/
Matrix Matrix::operator-(const Matrix& other) const
{
    int i = 0, j = 0;
    assert(m_iColumns == other.GetColumns() &&
           m_iRows == other.GetRows());

    Matrix result(*this);

    for(i=0;i<m_iRows;++i)
    {
        for(j=0;j<m_iColumns;++j)
        {
            result.SetElement(i, j, result.GetElement(i, j) - other.GetElement(i, j));
        }
    }

    return result;
}


/* ===========================================================================================
// Operator '*' overload (Scalar multiplication)
============================================================================================*/
Matrix Matrix::operator*(double value) const
{
    int i = 0, j = 0;
    Matrix result(*this);

    for(i=0;i<m_iRows;++i)
    {
        for(j=0;j<m_iColumns;++j)
        {
            result.SetElement(i, j, result.GetElement(i, j)*value);
        }
    }

    return result;
}


/* ===========================================================================================
// Operator '=' overload (Matrix multiplication)
============================================================================================*/
Matrix Matrix::operator*(const Matrix& other) const
{
    assert(m_iColumns == other.GetRows());

    int i = 0, j = 0, k = 0;
    Matrix result(m_iRows, other.GetColumns());

    double value;
    for(i=0;i<result.GetRows();++i)
    {
        for(j=0;j<other.GetColumns();++j)
        {
            value = 0.0;
            for(k=0;k<m_iColumns;++k)
            {
                value += GetElement(i, k)*other.GetElement(k, j);
            }
            result.SetElement(i, j, value);
        }
    }

    return result;
}


/* ===========================================================================================
// Get data pointer of the matrix
============================================================================================*/
void Matrix::SetData(double* value)
{
    memset(m_pdData, 0, sizeof(double)*m_iColumns*m_iRows);
    memcpy(m_pdData, value, sizeof(double)*m_iColumns*m_iRows);
}


/* ===========================================================================================
// Get data pointer of the matrix
============================================================================================*/
double* Matrix::GetData() const
{
    return m_pdData;
}


/* ===========================================================================================
// Set matrix element
============================================================================================*/
bool Matrix::SetElement(int row, int col, double value)
{
    if(col<0 || col>=m_iColumns || row<0 || row>=m_iRows)
    {
        return false;
    }

    if(m_pdData == NULL)
    {
        return false;
    }

    m_pdData[col+row*m_iColumns] = value;

    return true;
}


/* ===========================================================================================
// Get matrix element
============================================================================================*/
double Matrix::GetElement(int row, int col) const
{
    if(col<0 || col>=m_iColumns || row<0 || row>=m_iRows)
    {
        return false;
    }

    if(m_pdData == NULL)
    {
        return false;
    }

    return m_pdData[col+row*m_iColumns];
}


/* ===========================================================================================
// Get matrix column number
============================================================================================*/
int Matrix::GetColumns() const
{
    return m_iColumns;
}


/* ===========================================================================================
// Get matrix row number
============================================================================================*/
int Matrix::GetRows() const
{
    return m_iRows;
}
