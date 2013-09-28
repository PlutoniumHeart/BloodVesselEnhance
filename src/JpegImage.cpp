#include "JpegImage.h"


JpegImage::JpegImage(int width, int height)
{
    m_lWidth = width;
    m_lHeight = height;

    m_pJpegIO = JPEGIOType::New();

    m_ImageObject = JPEGImageType::New();
    JPEGImageType::IndexType index;
    JPEGImageType::SizeType size;
    JPEGImageType::RegionType region;
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
    Reshape221<JPEGImageType::Pointer, unsigned char, JPEGConstIteratorType>(m_ImageObject, m_cPixelData);
}


JpegImage::JpegImage(std::string fileName)
{
    m_pJpegIO = JPEGIOType::New();
    m_ImageObject = JPEGImageType::New();
    ReadJpegFile(fileName);
}


JpegImage::~JpegImage()
{
}


int JpegImage::ReadJpegFile(std::string inputFile)
{
    JPEGReadType::Pointer reader = JPEGReadType::New();
    reader->SetImageIO(m_pJpegIO);
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
    Reshape221<JPEGImageType::Pointer, unsigned char, JPEGConstIteratorType>(m_ImageObject, m_cPixelData);
    return 0;
}


int JpegImage::WriteJpegFile(std::string outputFile)
{
    JPEGWriteType::Pointer writer = JPEGWriteType::New();

    writer->SetImageIO(m_pJpegIO);
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


JPEGImageType::Pointer JpegImage::GetImageObject()
{
    return m_ImageObject;
}


void JpegImage::Update()
{
    Reshape122<JPEGImageType::Pointer, unsigned char, JPEGIteratorType>(m_cPixelData, m_ImageObject);
}