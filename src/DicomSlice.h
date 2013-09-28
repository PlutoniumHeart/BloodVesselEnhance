#ifndef DICOM2DIO_H
#define DICOM2DIO_H


#include <itkGDCMImageIO.h>
#include <itkMetaDataDictionary.h>
#include <itkMetaDataObject.h>
#include "BaseImage.h"


typedef itk::Image<short, D2>                               DICOMImageType2D;
typedef itk::GDCMImageIO                                    DICOMIOType;
typedef itk::ImageFileReader<DICOMImageType2D>              DICOM2DReadType;
typedef itk::ImageFileWriter<DICOMImageType2D>              DICOM2DWriteType;
typedef itk::ImageRegionConstIterator<DICOMImageType2D>     DICOMConstIterator2DType;
typedef itk::ImageRegionIterator<DICOMImageType2D>          DICOMIterator2DType;
typedef itk::MetaDataDictionary                             DICOMDictionaryType;


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
