#ifndef JPEGSERIES_H
#define JPEGSERIES_H


#include "BaseImageSeries.h"
#include "JpegImage.h"


typedef itk::Image<unsigned char, D3>                           JPEGSeriesType;
typedef itk::ImageSeriesReader<JPEGSeriesType>                  JPEGSeriesReadType;
typedef itk::ImageSeriesWriter<JPEGSeriesType, JPEGImageType>   JPEGSeriesWriteType;
typedef itk::ImageRegionConstIterator<JPEGSeriesType>           JPEGSeriesConstIteratorType;
typedef itk::ImageRegionIterator<JPEGSeriesType>                JPEGSeriesIteratorType;


class JpegSeries : public BaseImageSeries
{
public:
    JpegSeries(int width, int height, int depth);
    JpegSeries(std::string folderName, std::string seriesFormat, int begin, int end);
    virtual ~JpegSeries();

    virtual void Update();
    int ReadJpegSeries(std::string folderName, std::string seriesFormat, int begin, int end);
    int WriteJpegSeries(std::string folderName, std::string seriesFormat, int begin, int end);
    JPEGSeriesType::Pointer GetImageObject();
protected:
    JPEGIOType::Pointer m_pJpegIO;
    JPEGSeriesType::Pointer m_ImageObject;
};


#endif // JPEGSERIES_H
