#ifndef DICOMSERIES_H
#define DICOMSERIES_H


#include <itkGDCMSeriesFileNames.h>

#include "BaseImageSeries.h"
#include "DicomSlice.h"


typedef itk::Image<short, D3>                                           DICOMSeriesType;
typedef itk::ImageSeriesReader<DICOMSeriesType>                         DICOMSeriesReadType;
typedef itk::ImageSeriesWriter<DICOMSeriesType, DICOMImageType2D>       DICOMSeriesWriteType;
typedef itk::ImageRegionConstIterator<DICOMSeriesType>                  DICOMSeriesConstIteratorType;
typedef itk::ImageRegionIterator<DICOMSeriesType>                       DICOMSeriesIteratorType;
typedef itk::GDCMSeriesFileNames                                        DICOMSeriesNameGeneratorType;


class DicomSeries : public BaseImageSeries
{
public:
    DicomSeries(int width, int height, int depth);
    DicomSeries(std::string folderName, int seriesNumber);
    virtual ~DicomSeries();

    virtual void Update();
    int ReadDicomSeries(std::string folderName, std::string seriesUID);
    int WriteDicomSeries(std::string folderName, std::string seriesFormat, int begin, int end);
protected:
    DICOMIOType::Pointer m_pDicomIO;
    DICOMSeriesType::Pointer m_ImageObject;
private:
    int m_iNumberOfSeries;
    std::string *m_pSerieName;
    std::vector<std::string> m_vFileNameContainer;
private:
    DICOMSeriesNameGeneratorType::Pointer nameGenerator;
    int SetDataDirectory(std::string folderName);
    std::string GetSerieUID(int number);
};


#endif // DICOMSERIES_H
