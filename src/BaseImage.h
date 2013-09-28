#ifndef IMAGEIO_H
#define IMAGEIO_H


#include <itkImage.h>
#include <itkGDCMImageIO.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMetaDataDictionary.h>
#include <itkMetaDataObject.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkPNGImageIO.h>
#include <itkJPEGImageIO.h>
//#include <itkCastImageFilter.h>


const unsigned char D2 = 2;

typedef itk::ImageRegion<D2>                                ImageRegion2D;

typedef itk::Image<short, D2>                               DICOMImageType2D;
typedef itk::GDCMImageIO                                    DICOMIOType;
typedef itk::ImageFileReader<DICOMImageType2D>              DICOM2DReadType;
typedef itk::ImageFileWriter<DICOMImageType2D>              DICOM2DWriteType;
typedef itk::ImageRegionConstIterator<DICOMImageType2D>     DICOMConstIterator2DType;
typedef itk::ImageRegionIterator<DICOMImageType2D>          DICOMIterator2DType;
typedef itk::MetaDataDictionary                             DICOMDictionaryType;

typedef itk::Image<unsigned char, D2>                       PNGImageType;
typedef itk::PNGImageIO                                     PNGIOType;
typedef itk::ImageFileReader<PNGImageType>                  PNGReadType;
typedef itk::ImageFileWriter<PNGImageType>                  PNGWriteType;
typedef itk::ImageRegionConstIterator<PNGImageType>         PNGConstIteratorType;
typedef itk::ImageRegionIterator<PNGImageType>              PNGIteratorType;

typedef itk::Image<unsigned char, D2>                       JPEGImageType;
typedef itk::JPEGImageIO                                    JPEGIOType;
typedef itk::ImageFileReader<JPEGImageType>                 JPEGReadType;
typedef itk::ImageFileWriter<JPEGImageType>                 JPEGWriteType;
typedef itk::ImageRegionConstIterator<JPEGImageType>        JPEGConstIteratorType;
typedef itk::ImageRegionIterator<JPEGImageType>             JPEGIteratorType;


class BaseImage
{
public:
    BaseImage();
    ~BaseImage();

    short* GetShortPixelData();
    unsigned char* GetUnsignedCharPixelData();
    unsigned long GetImageWidth();
    unsigned long GetImageHeight();
    virtual void Update() = 0;
    int CastShortToUnsignedChar(BaseImage* output, int window_pos, int window_half_size);
    int CastUnsignedCharToShort(BaseImage* output);
protected:
    template<typename T1, typename T2, typename T3> 
    int Reshape221(T1 input, T2* destination);

    template<typename T1, typename T2, typename T3> 
    int Reshape122(T2* input, T1 destination);
protected:
    unsigned long m_lWidth;
    unsigned long m_lHeight;
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


#endif // !BASEIMAGEIO_H
