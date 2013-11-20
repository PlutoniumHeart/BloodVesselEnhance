#ifndef PNGIMAGE_H
#define PNGIMAGE_H


#include <itkPNGImageIO.h>
#include "BaseImage.h"


typedef itk::Image<unsigned char, D2>                       PNGImageType;
typedef itk::PNGImageIO                                     PNGIOType;
typedef itk::ImageFileReader<PNGImageType>                  PNGReadType;
typedef itk::ImageFileWriter<PNGImageType>                  PNGWriteType;
typedef itk::ImageRegionConstIterator<PNGImageType>         PNGConstIteratorType;
typedef itk::ImageRegionIterator<PNGImageType>              PNGIteratorType;


class PngImage : public BaseImage
{
public:
    PngImage(int width, int height);
    PngImage(std::string fileName);
    virtual ~PngImage();

    int ReadPngFile(std::string inputFile);
    int WritePngFile(std::string outputFile);
    PNGImageType::Pointer GetImageObject();
    virtual void Update();
protected:
    PNGIOType::Pointer m_pPngIO;
    PNGImageType::Pointer m_ImageObject;
};


#endif // PNGIMAGE_H
