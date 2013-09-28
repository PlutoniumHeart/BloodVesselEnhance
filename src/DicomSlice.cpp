#include "DicomSlice.h"


DicomSlice::DicomSlice(int width, int height)
{
    m_lWidth = width;
    m_lHeight = height;

    m_pDicomIO = DICOMIOType::New();

    m_ImageObject = DICOMImageType2D::New();
    DICOMImageType2D::IndexType index;
    DICOMImageType2D::SizeType size;
    DICOMImageType2D::RegionType region;
    index[0] = 0;
    index[1] = 0;
    size[0] = m_lWidth;
    size[1] = m_lHeight;
    region.SetSize(size);
    region.SetIndex(index);
    m_ImageObject->SetRegions(region);
    m_ImageObject->Allocate();
    m_ImageObject->FillBuffer(0);

    m_sPixelData = new short[m_lWidth*m_lHeight];
    Reshape221<DICOMImageType2D::Pointer, short, DICOMConstIterator2DType>(m_ImageObject, m_sPixelData);
}


DicomSlice::DicomSlice(std::string fileName)
{
    m_pDicomIO = DICOMIOType::New();
    ReadDicomFile(fileName);
}


DicomSlice::~DicomSlice()
{
}


int DicomSlice::ReadDicomFile(std::string inputFile)
{
    DICOM2DReadType::Pointer reader = DICOM2DReadType::New();
    reader->SetImageIO(m_pDicomIO);
    reader->SetFileName(inputFile);
    try
    {
        reader->Update();
    }
    catch(itk::ExceptionObject &e)
    {
        std::cerr<<"Exception in reading DICOM file: "<<std::endl;
        std::cerr<<e<<std::endl;
        return -1;
    }
    m_ImageObject = reader->GetOutput();
    m_HeaderDictionary = m_ImageObject->GetMetaDataDictionary();
    ImageRegion2D::SizeType tmp = m_ImageObject->GetLargestPossibleRegion().GetSize();
    m_lWidth = tmp.GetElement(0);
    m_lHeight = tmp.GetElement(1);
    m_sPixelData = new short[m_lWidth*m_lHeight];
    Reshape221<DICOMImageType2D::Pointer, short, DICOMConstIterator2DType>(m_ImageObject, m_sPixelData);
    return 0;
}


int DicomSlice::WriteDicomFile(std::string outputFile)
{
    DICOM2DWriteType::Pointer writer = DICOM2DWriteType::New();

    //Reshape122<DICOMImageType2D::Pointer, short, DICOMIterator2DType>(m_sPixelData, m_ImageObject);
    writer->SetMetaDataDictionary(m_HeaderDictionary);
    writer->UseInputMetaDataDictionaryOff();
    writer->SetImageIO(m_pDicomIO);
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


DICOMImageType2D::Pointer DicomSlice::GetImageObject()
{
    return m_ImageObject;
}


void DicomSlice::Update()
{
    Reshape122<DICOMImageType2D::Pointer, short, DICOMIterator2DType>(m_sPixelData, m_ImageObject);
}