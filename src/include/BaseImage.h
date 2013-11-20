#ifndef BASEIMAGE_H
#define BASEIMAGE_H


#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>


const unsigned char D2 = 2;

typedef itk::ImageRegion<D2>            ImageRegion2D;


class BaseImage
{
public:
    BaseImage();
    virtual ~BaseImage();

    virtual short* GetShortPixelData();
    virtual unsigned char* GetUnsignedCharPixelData();
    unsigned long GetImageWidth();
    unsigned long GetImageHeight();
    unsigned long GetImageDepth();
    virtual void Update() = 0;
    virtual int CastShortToUnsignedChar(BaseImage* output, int window_pos, int window_half_size);
    virtual int CastUnsignedCharToShort(BaseImage* output);
protected:
    template<typename T1, typename T2, typename T3> 
    int Reshape221(T1 input, T2* destination);

    template<typename T1, typename T2, typename T3> 
    int Reshape122(T2* input, T1 destination);
protected:
    unsigned long m_lWidth;
    unsigned long m_lHeight;
    unsigned long m_lDepth;
    short* m_sPixelData;
    unsigned char* m_cPixelData;
};

template<typename T1, typename T2, typename T3>
int BaseImage::Reshape221(T1 input, T2* destination)
{
    int i = 0;
    T3 in(input, input->GetLargestPossibleRegion());
    in.GoToBegin();
    while(!in.IsAtEnd())
    {
        destination[i++] = in.Get();
        ++in;
    }
    return 0;
}

template<typename T1, typename T2, typename T3>
int BaseImage::Reshape122(T2* input, T1 destination)
{
    int i = 0;
    T3 out(destination, destination->GetLargestPossibleRegion());
    out.GoToBegin();
    while(!out.IsAtEnd())
    {
        out.Set(*(input+i++));
        ++out;
    }
    return 0;
}


#endif // BASEIMAGEIO_H
