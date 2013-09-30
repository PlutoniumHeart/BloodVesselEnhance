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
// Get data pointer of the matrix
============================================================================================*/
double* Matrix::GetData() const
{
    return m_pdData;
}


double Matrix::GetElement(int row, int col) const
{
    assert(col>=0 && col<m_iColumns && row>=0 && row<m_iRows);
    assert(m_pdData);
    return m_pdData[col+row*m_iColumns];
}