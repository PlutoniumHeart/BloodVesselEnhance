#include "StructureClassify.h"


StructureClassify::StructureClassify()
    : m_iWidth(0)
    , m_iHeight(0)
    , m_iDepth(0)
    , m_pInput(NULL)
    , m_ipHist(NULL)
    , m_pResult(NULL)
    , m_pTemp(NULL)
    , m_pOutput(NULL)
    , m_Filter(NULL)
{
    m_Filter = new CFilter<unsigned char, short>;
}


StructureClassify::~StructureClassify()
{
    if(m_pInput != NULL)
    {
        delete [] m_pInput;
        m_pInput = NULL;
    }

    if(m_ipHist != NULL)
    {
        delete [] m_ipHist;
        m_ipHist = NULL;
    }

    if(m_pResult != NULL)
    {
        delete [] m_pResult;
        m_pResult = NULL;
    }

    if(m_pTemp != NULL)
    {
        delete [] m_pTemp;
        m_pTemp = NULL;
    }

    if(m_pOutput != NULL)
    {
        delete [] m_pOutput;
        m_pOutput = NULL;
    }

    if(m_Filter != NULL)
    {
        delete m_Filter;
        m_Filter = NULL;
    }
}


void StructureClassify::SetInput(BaseImageSeries* input)
{
    m_iWidth = input->GetImageWidth();
    m_iHeight = input->GetImageHeight();
    m_iDepth = input->GetImageDepth();
    m_pInput = new unsigned char[m_iWidth*m_iHeight*m_iDepth];
    for(int i=0;i<m_iWidth*m_iHeight*m_iDepth;i++)
    {
        m_pInput[i] = input->GetUnsignedCharPixelData()[i];
    }
    Histogram(input);
    m_Filter->SetDim(m_iWidth, m_iHeight, m_iDepth);
}


void StructureClassify::Histogram(BaseImageSeries* input)
{
    int i = 0;
    int w = input->GetImageWidth(), h = input->GetImageHeight(), d = input->GetImageDepth();

    if(input->GetShortPixelData() != NULL)
    {
        m_ipHist = new int[0xFFF]; // Only care about DICOM with 12 bits

        for(i=0;i<SHRT_MAX+abs(SHRT_MIN);i++)
        {
            m_ipHist[i] = 0;
        }

        for(i=0;i<w*h;i++)
        {
            m_ipHist[input->GetShortPixelData()[i] + abs(SHRT_MIN)]++;
        }
    }
    else if(input->GetUnsignedCharPixelData() != NULL)
    {
        m_ipHist = new int[UCHAR_MAX+1];

        for(i=0;i<=UCHAR_MAX;i++)
        {
            m_ipHist[i] = 0;
        }

        for(i=0;i<w*h*d;i++)
        {
            m_ipHist[input->GetUnsignedCharPixelData()[i]]++;
        }
    }
    else
    {
        std::cerr<<"Input image is an empty image. "<<std::endl;
        return;
    }
}


int* StructureClassify::GetHistogram()
{
    return m_ipHist;
}


int StructureClassify::ClassifyBrightTube(int iterations, float rMax, float rMin)
{
    int i = 0, j = 0;
    int secondPeak = 1;

    m_pResult = new unsigned char[m_iWidth*m_iHeight*m_iDepth];
    m_pOutput = new unsigned char[m_iWidth*m_iHeight*m_iDepth];

    for(i=1;i<255;i++)
    {
        if(m_ipHist[i]>m_ipHist[secondPeak])
        {
            secondPeak = i;
        }
    }

#ifdef _DEBUG
    std::cout<<"Max grey value: "<<secondPeak<<std::endl;
#endif

    double r1 = rMin, r2 = rMax;
    float sigma = 0.0;

    short* xx = new short[m_iWidth*m_iHeight*m_iDepth];
    short* yy = new short[m_iWidth*m_iHeight*m_iDepth];
    short* zz = new short[m_iWidth*m_iHeight*m_iDepth];
    short* xy = new short[m_iWidth*m_iHeight*m_iDepth];
    short* xz = new short[m_iWidth*m_iHeight*m_iDepth];
    short* yz = new short[m_iWidth*m_iHeight*m_iDepth];

    Matrix i_Matrix(3, 3);
    Matrix i_MatrixEigenVector(3, 3);

    double lamada[3], temp = 0.0;
    double lineValue = 0.0, tempValue = 0.0;

    for(j=0;j<iterations;j++)
    {
        sigma = pow(r2/r1, (float)j/(iterations-1))*r1/4.0;

#ifdef _DEBUG
        std::cout<<"Sigma: "<<sigma<<std::endl;
#endif

        m_Filter->RecursiveGaussianFilter(2,0,0,m_pInput,xx,sigma);
        m_Filter->RecursiveGaussianFilter(0,0,2,m_pInput,zz,sigma);
        m_Filter->RecursiveGaussianFilter(1,0,1,m_pInput,xz,sigma);
        m_Filter->RecursiveGaussianFilter(0,2,0,m_pInput,yy,sigma);
        m_Filter->RecursiveGaussianFilter(1,1,0,m_pInput,xy,sigma);
        m_Filter->RecursiveGaussianFilter(0,1,1,m_pInput,yz,sigma);

        for(i=0;i<m_iWidth*m_iHeight*m_iDepth;i++)
        {
            if(m_pInput[i]>0 && m_pInput[i]>secondPeak)
            {
                *(i_Matrix.GetData()  ) = xx[i];
                *(i_Matrix.GetData()+1) = xy[i];
                *(i_Matrix.GetData()+2) = xz[i];
                *(i_Matrix.GetData()+3) = xy[i];
                *(i_Matrix.GetData()+4) = yy[i];
                *(i_Matrix.GetData()+5) = yz[i];
                *(i_Matrix.GetData()+6) = xz[i];
                *(i_Matrix.GetData()+7) = yz[i];
                *(i_Matrix.GetData()+8) = zz[i];

                i_Matrix.JacobiEigenv(lamada, i_MatrixEigenVector);
                if(lamada[0]>lamada[1])
                {
                    temp = lamada[0];
                    lamada[0] = lamada[1];
                    lamada[1] = temp;
                }
                if(lamada[0]>lamada[2])
                {
                    temp = lamada[0];
                    lamada[0] = lamada[2];
                    lamada[2] = temp;
                }
                if(lamada[1]>lamada[2])
                {
                    temp = lamada[1];
                    lamada[1] = lamada[2];
                    lamada[2] = temp;
                }

                // Improved version
                if(lamada[0]<0 && lamada[1]<0)
                    lineValue = -lamada[1]*(1-exp(-(m_pInput[i])/255.0))*sigma*sigma;
                else
                    lineValue = 0.0;

                if(j == 0)
                {
                    if(lineValue>255.0)
                        m_pResult[i] = 255;
                    else if(lineValue<0.0)
                        m_pResult[i] = 0;
                    else
                        m_pResult[i] = lineValue;
                }
                else
                {
                    tempValue = m_pResult[i];
                    if(tempValue<lineValue)
                    {
                        if(lineValue>255.0)
                            m_pResult[i] = 255;
                        else if(lineValue<0.0)
                            m_pResult[i] = 0;
                        else
                            m_pResult[i] = lineValue;
                    }
                }
            }
            else
            {
                m_pResult[i] = 0.0;
            }
        }
    }

    delete [] xx;
    delete [] yy;
    delete [] zz;
    delete [] xy;
    delete [] xz;
    delete [] yz;

    return 0;
}


unsigned char* StructureClassify::GetResult()
{
    return m_pResult;
}
