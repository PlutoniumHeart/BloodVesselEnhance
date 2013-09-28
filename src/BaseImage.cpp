#include "BaseImage.h"


BaseImage::BaseImage()
    : m_lWidth(0)
    , m_lHeight(0)
    , m_sPixelData(NULL)
    , m_cPixelData(NULL)
{
}


BaseImage::~BaseImage()
{
    if(m_sPixelData!=NULL)
    {
        delete [] m_sPixelData;
    }
    if(m_cPixelData!=NULL)
    {
        delete [] m_cPixelData;
    }
}


short* BaseImage::GetShortPixelData()
{
    return m_sPixelData;
}


unsigned char* BaseImage::GetUnsignedCharPixelData()
{
    return m_cPixelData;
}


unsigned long BaseImage::GetImageWidth()
{
    return m_lWidth;
}


unsigned long BaseImage::GetImageHeight()
{
    return m_lHeight;
}


int BaseImage::CastShortToUnsignedChar(BaseImage* output, int window_pos, int window_half_size)
{
    int i = 0;

    if(output->GetUnsignedCharPixelData() != NULL)
    {
        for(i=0;i<m_lWidth*m_lHeight;i++)
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


int BaseImage::CastUnsignedCharToShort(BaseImage* output)
{
    int i = 0;

    if(output->GetShortPixelData() != NULL)
    {
        for(i=0;i<m_lWidth*m_lHeight;i++)
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