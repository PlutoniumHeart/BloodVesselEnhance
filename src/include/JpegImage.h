#ifndef JPEGIMAGE_H
#define JPEGIMAGE_H


#include <itkJPEGImageIO.h>
#include "BaseImage.h"


typedef itk::Image<unsigned char, D2>                       JPEGImageType;
typedef itk::JPEGImageIO                                    JPEGIOType;
typedef itk::ImageFileReader<JPEGImageType>                 JPEGReadType;
typedef itk::ImageFileWriter<JPEGImageType>                 JPEGWriteType;
typedef itk::ImageRegionConstIterator<JPEGImageType>        JPEGConstIteratorType;
typedef itk::ImageRegionIterator<JPEGImageType>             JPEGIteratorType;


class JpegImage : public BaseImage
{
public:
    JpegImage(int width, int height);
    JpegImage(std::string fileName);
    virtual ~JpegImage();

    int ReadJpegFile(std::string inputFile);
    int WriteJpegFile(std::string outputFile);
    JPEGImageType::Pointer GetImageObject();
    virtual void Update();
protected:
    JPEGIOType::Pointer m_pJpegIO;
    JPEGImageType::Pointer m_ImageObject;
};


#endif // JPEGIMAGE_H
