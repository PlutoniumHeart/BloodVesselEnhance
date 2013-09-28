#ifndef DICOM2DIO_H
#define DICOM2DIO_H


#include "BaseImage.h"


class DicomSlice : public BaseImage
{
public:
    DicomSlice(int width, int height);
    DicomSlice(std::string fileName);
    ~DicomSlice();

    int ReadDicomFile(std::string inputFile);
    int WriteDicomFile(std::string outputFile);
    DICOMImageType2D::Pointer GetImageObject();
    virtual void Update();
protected:
    DICOMDictionaryType m_HeaderDictionary;
    DICOMIOType::Pointer m_pDicomIO;
    DICOMImageType2D::Pointer m_ImageObject;
};


#endif // !DICOM2DIO_H
