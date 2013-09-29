#ifndef PNGSERIES_H
#define PNGSERIES_H


#include "BaseImageSeries.h"
#include "PngImage.h"


typedef itk::Image<unsigned char, D3>                           PNGSeriesType;
typedef itk::ImageSeriesReader<PNGSeriesType>                   PNGSeriesReadType;
typedef itk::ImageSeriesWriter<PNGSeriesType, PNGImageType>     PNGSeriesWriteType;
typedef itk::ImageRegionConstIterator<PNGSeriesType>            PNGSeriesConstIteratorType;
typedef itk::ImageRegionIterator<PNGSeriesType>                 PNGSeriesIteratorType;


class PngSeries : public BaseImageSeries
{
public:
    PngSeries(int width, int height, int depth);
    PngSeries(std::string folderName, std::string seriesFormat, int begin, int end);
    ~PngSeries();

    virtual void Update();
    int ReadPngSeries(std::string folderName, std::string seriesFormat, int begin, int end);
    int WritePngSeries(std::string folderName, std::string seriesFormat, int begin, int end);
    PNGSeriesType::Pointer GetImageObject();
protected:
    PNGIOType::Pointer m_pPngIO;
    PNGSeriesType::Pointer m_ImageObject;
};



#endif // !PNGSERIES_H
