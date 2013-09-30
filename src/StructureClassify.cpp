#include "StructureClassify.h"


StructureClassify::StructureClassify()
    : m_ipHist(NULL)
{
}


StructureClassify::~StructureClassify()
{
}


void StructureClassify::Histogram(BaseImage* input)
{
    int i = 0;
    int w = input->GetImageWidth(), h = input->GetImageHeight();

    if(input->GetShortPixelData() != NULL)
    {
        m_ipHist = new int[SHRT_MAX+abs(SHRT_MIN)];

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
        m_ipHist = new int[UCHAR_MAX];

        for(i=0;i<UCHAR_MAX;i++)
        {
            m_ipHist[i] = 0;
        }

        for(i=0;i<w*h;i++)
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


void StructureClassify::Histogram(BaseImageSeries* input)
{
    int i = 0;
    int w = input->GetImageWidth(), h = input->GetImageHeight(), d = input->GetImageDepth();

    if(input->GetShortPixelData() != NULL)
    {
        m_ipHist = new int[SHRT_MAX+abs(SHRT_MIN)];

        for(i=0;i<SHRT_MAX+abs(SHRT_MIN);i++)
        {
            m_ipHist[i] = 0;
        }

        for(i=0;i<w*h*d;i++)
        {
            m_ipHist[input->GetShortPixelData()[i] + abs(SHRT_MIN)]++;
        }
    }
    else if(input->GetUnsignedCharPixelData() != NULL)
    {
        m_ipHist = new int[UCHAR_MAX];

        for(i=0;i<UCHAR_MAX;i++)
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