#include "PngSeries.h"


PngSeries::PngSeries(int width, int height, int depth)
{
    m_lWidth = width;
    m_lHeight = height;
    m_lDepth = depth;

    m_pPngIO = PNGIOType::New();

    m_ImageObject = PNGSeriesType::New();
    PNGSeriesType::IndexType index;
    PNGSeriesType::SizeType size;
    PNGSeriesType::RegionType region;
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
    Reshape321<PNGSeriesType::Pointer, unsigned char, PNGSeriesConstIteratorType>(m_ImageObject, m_cPixelData);
}


PngSeries::PngSeries(std::string folderName, std::string seriesFormat, int begin, int end)
{
    m_pPngIO = PNGIOType::New();
    m_ImageObject = PNGSeriesType::New();
    ReadPngSeries(folderName, seriesFormat, begin, end);
}


PngSeries::~PngSeries()
{
}


int PngSeries::ReadPngSeries(std::string folderName, std::string seriesFormat, int begin, int end)
{
    PNGSeriesReadType::Pointer reader = PNGSeriesReadType::New();
    SeriesNameGeneratorType::Pointer nameGenerator = SeriesNameGeneratorType::New();

    nameGenerator->SetSeriesFormat(folderName+seriesFormat);
    nameGenerator->SetStartIndex(begin);
    nameGenerator->SetEndIndex(end);
    nameGenerator->SetIncrementIndex(1);

    reader->SetImageIO(m_pPngIO);
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
    Reshape321<PNGSeriesType::Pointer, unsigned char, PNGSeriesConstIteratorType>(m_ImageObject, m_cPixelData);
    return 0;
}


int PngSeries::WritePngSeries(std::string folderName, std::string seriesFormat, int begin, int end)
{
    PNGSeriesWriteType::Pointer writer = PNGSeriesWriteType::New();
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


void PngSeries::Update()
{
    Reshape123<PNGSeriesType::Pointer, unsigned char, PNGSeriesIteratorType>(m_cPixelData, m_ImageObject);
}


PNGSeriesType::Pointer PngSeries::GetImageObject()
{
    return m_ImageObject;
}
