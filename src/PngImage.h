#ifndef PNGIMAGE_H
#define PNGIMAGE_H


#include "BaseImage.h"


class PngImage : public BaseImage
{
public:
    PngImage(int width, int height);
    PngImage(std::string fileName);
    ~PngImage();

    int ReadPngFile(std::string inputFile);
    int WritePngFile(std::string outputFile);
    PNGImageType::Pointer GetImageObject();
    virtual void Update();
protected:
    PNGIOType::Pointer m_pPngIO;
    PNGImageType::Pointer m_ImageObject;
};


#endif // !PNGIMAGE_H
