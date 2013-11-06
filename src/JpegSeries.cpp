#include "JpegSeries.h"


JpegSeries::JpegSeries(int width, int height, int depth)
{
    m_lWidth = width;
    m_lHeight = height;
    m_lDepth = depth;

    m_pJpegIO = JPEGIOType::New();

    m_ImageObject = JPEGSeriesType::New();
    JPEGSeriesType::IndexType index;
    JPEGSeriesType::SizeType size;
    JPEGSeriesType::RegionType region;
    index[0] = 0;
    index[1] = 0;
    index[2] = 0;
    size[0] = m_lWidth;
    size[1] = m_lHeight;
    size[2] = m_lDepth;
    region.SetSize(size);
    region.SetIndex(index);
    m_ImageObject->SetRegions(region);
    m_ImageObject->Allocate();
    m_ImageObject->FillBuffer(0);

    m_cPixelData = new unsigned char[m_lWidth*m_lHeight*m_lDepth];
    Reshape321<JPEGSeriesType::Pointer, unsigned char, JPEGSeriesConstIteratorType>(m_ImageObject, m_cPixelData);
}


JpegSeries::JpegSeries(std::string folderName, std::string seriesFormat, int begin, int end)
{
    m_pJpegIO = JPEGIOType::New();
    m_ImageObject = JPEGSeriesType::New();
    ReadJpegSeries(folderName, seriesFormat, begin, end);
}


JpegSeries::~JpegSeries()
{
}


int JpegSeries::ReadJpegSeries(std::string folderName, std::string seriesFormat, int begin, int end)
{
    JPEGSeriesReadType::Pointer reader = JPEGSeriesReadType::New();
    SeriesNameGeneratorType::Pointer nameGenerator = SeriesNameGeneratorType::New();

    nameGenerator->SetSeriesFormat(folderName+seriesFormat);
    nameGenerator->SetStartIndex(begin);
    nameGenerator->SetEndIndex(end);
    nameGenerator->SetIncrementIndex(1);

    reader->SetImageIO(m_pJpegIO);
    reader->SetFileNames(nameGenerator->GetFileNames());
    try
    {
        reader->Update();
    }
    catch(itk::ExceptionObject &e)
    {
        std::cerr<<"Exception caught in reading file series: "<<std::endl;
        std::cerr<<e<<std::endl;
    }
    m_ImageObject = reader->GetOutput();
    ImageRegion3D::SizeType tmp = m_ImageObject->GetLargestPossibleRegion().GetSize();
    m_lWidth = tmp.GetElement(0);
    m_lHeight = tmp.GetElement(1);
    m_lDepth = tmp.GetElement(2);
    m_cPixelData = new unsigned char[m_lWidth*m_lHeight*m_lDepth];
    Reshape321<JPEGSeriesType::Pointer, unsigned char, JPEGSeriesConstIteratorType>(m_ImageObject, m_cPixelData);
    return 0;
}


int JpegSeries::WriteJpegSeries(std::string folderName, std::string seriesFormat, int begin, int end)
{
    JPEGSeriesWriteType::Pointer writer = JPEGSeriesWriteType::New();
    SeriesNameGeneratorType::Pointer nameGenerator = SeriesNameGeneratorType::New();

    nameGenerator->SetSeriesFormat(folderName+seriesFormat);
    nameGenerator->SetStartIndex(begin);
    nameGenerator->SetEndIndex(end);
    nameGenerator->SetIncrementIndex(1);

    writer->SetInput(m_ImageObject);
    writer->SetFileNames(nameGenerator->GetFileNames());
    try
    {
        writer->Update();
    }
    catch(itk::ExceptionObject &e)
    {
        std::cerr<<"Execption caught in writing files: "<<std::endl;
        std::cerr<<e<<std::endl;
    }
    return 0;
}


void JpegSeries::Update()
{
    Reshape123<JPEGSeriesType::Pointer, unsigned char, JPEGSeriesIteratorType>(m_cPixelData, m_ImageObject);
}


JPEGSeriesType::Pointer JpegSeries::GetImageObject()
{
    return m_ImageObject;
}
