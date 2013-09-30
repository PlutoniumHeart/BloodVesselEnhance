#include "BaseImageSeries.h"


BaseImageSeries::BaseImageSeries()
    : m_lDepth(0)
{
}


BaseImageSeries::~BaseImageSeries()
{
}


unsigned long BaseImageSeries::GetImageDepth()
{
    return m_lDepth;
}


int BaseImageSeries::CastShortToUnsignedChar(BaseImageSeries* output, int window_pos, int window_half_size)
{
    int i = 0;

    if(output->GetUnsignedCharPixelData() != NULL)
    {
        for(i=0;i<m_lWidth*m_lHeight*m_lDepth;i++)
        {
            if(m_sPixelData[i]>=window_pos+window_half_size)
            {
                output->GetUnsignedCharPixelData()[i] = 255;
            }
            else if(m_sPixelData[i]<=window_pos-window_half_size)
            {
                output->GetUnsignedCharPixelData()[i] = 0;
            }
            else
            {
                output->GetUnsignedCharPixelData()[i] = 255*(m_sPixelData[i]-window_pos-window_half_size-1)/(2.0*window_half_size);
            }
        }
    }
    else
    {
        std::cerr<<"Input is likely not a unsigned char image type. "<<std::endl;
        return -1;
    }
    output->Update();
    return 0;
}


int BaseImageSeries::CastUnsignedCharToShort(BaseImageSeries* output)
{
    int i = 0;

    if(output->GetShortPixelData() != NULL)
    {
        for(i=0;i<m_lWidth*m_lHeight*m_lDepth;i++)
        {
            output->GetShortPixelData()[i] = m_cPixelData[i];
        }
    }
    else
    {
        std::cerr<<"Input is likely not a short image type. "<<std::endl;
        return -1;
    }
    output->Update();
    return 0;
}
