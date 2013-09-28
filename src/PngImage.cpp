#include "PngImage.h"


PngImage::PngImage(int width, int height)
{
    m_lWidth = width;
    m_lHeight = height;

    m_pPngIO = PNGIOType::New();

    m_ImageObject = PNGImageType::New();
    PNGImageType::IndexType index;
    PNGImageType::SizeType size;
    PNGImageType::RegionType region;
    index[0] = 0;
    index[1] = 0;
    size[0] = m_lWidth;
    size[1] = m_lHeight;
    region.SetSize(size);
    region.SetIndex(index);
    m_ImageObject->SetRegions(region);
    m_ImageObject->Allocate();
    m_ImageObject->FillBuffer(0);

    m_cPixelData = new unsigned char[m_lWidth*m_lHeight];
    Reshape221<PNGImageType::Pointer, unsigned char, PNGConstIteratorType>(m_ImageObject, m_cPixelData);
}


PngImage::PngImage(std::string fileName)
{
    m_pPngIO = PNGIOType::New();
    m_ImageObject = PNGImageType::New();
    ReadPngFile(fileName);
}


PngImage::~PngImage()
{
}


int PngImage::ReadPngFile(std::string inputFile)
{
    PNGReadType::Pointer reader = PNGReadType::New();
    reader->SetImageIO(m_pPngIO);
    reader->SetFileName(inputFile);
    try
    {
        reader->Update();
    }
    catch(itk::ExceptionObject &e)
    {
        std::cerr<<"Exception in reading PNG file: "<<std::endl;
        std::cerr<<e<<std::endl;
        return -1;
    }
    m_ImageObject = reader->GetOutput();
    ImageRegion2D::SizeType tmp = m_ImageObject->GetLargestPossibleRegion().GetSize();
    m_lWidth = tmp.GetElement(0);
    m_lHeight = tmp.GetElement(1);
    m_cPixelData = new unsigned char[m_lWidth*m_lHeight];
    Reshape221<PNGImageType::Pointer, unsigned char, PNGConstIteratorType>(m_ImageObject, m_cPixelData);
    return 0;
}


int PngImage::WritePngFile(std::string outputFile)
{
    PNGWriteType::Pointer writer = PNGWriteType::New();

    writer->SetImageIO(m_pPngIO);
    writer->SetFileName(outputFile);
    writer->SetInput(m_ImageObject);
    try
    {
        writer->Update();
    }
    catch(itk::ExceptionObject &e)
    {
        std::cerr<<"Exception in file writing: "<<std::endl;
        std::cerr<<e<<std::endl;
        return -1;
    }
    return 0;
}


PNGImageType::Pointer PngImage::GetImageObject()
{
    return m_ImageObject;
}


void PngImage::Update()
{
    Reshape122<PNGImageType::Pointer, unsigned char, PNGIteratorType>(m_cPixelData, m_ImageObject);
}