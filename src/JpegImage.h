#ifndef JPEGIMAGE_H
#define JPEGIMAGE_H


#include "BaseImage.h"


class JpegImage : public BaseImage
{
public:
    JpegImage(int width, int height);
    JpegImage(std::string fileName);
    ~JpegImage();

    int ReadJpegFile(std::string inputFile);
    int WriteJpegFile(std::string outputFile);
    JPEGImageType::Pointer GetImageObject();
    virtual void Update();
protected:
    JPEGIOType::Pointer m_pJpegIO;
    JPEGImageType::Pointer m_ImageObject;
};


#endif // !JPEGIMAGE_H
