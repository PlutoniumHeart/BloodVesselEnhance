#include "DicomSeries.h"


DicomSeries::DicomSeries(int width, int height, int depth)
    : m_iNumberOfSeries(0)
    , m_pSerieName(NULL)
{
    m_lWidth = width;
    m_lHeight = height;
    m_lDepth = depth;

    m_pDicomIO = DICOMIOType::New();

    m_ImageObject = DICOMSeriesType::New();
    DICOMSeriesType::IndexType index;
    DICOMSeriesType::SizeType size;
    DICOMSeriesType::RegionType region;
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

    m_sPixelData = new short[m_lWidth*m_lHeight*m_lDepth];
    Reshape321<DICOMSeriesType::Pointer, short, DICOMSeriesConstIteratorType>(m_ImageObject, m_sPixelData);
}


DicomSeries::DicomSeries(std::string folderName, int seriesNumber)
    : m_iNumberOfSeries(0)
    , m_pSerieName(NULL)
{
    m_pDicomIO = DICOMIOType::New();
    m_ImageObject = DICOMSeriesType::New();
    SetDataDirectory(folderName);
    if(seriesNumber>=m_iNumberOfSeries)
    {
        std::cerr<<"seriesNumber is greater than what's in the folder."<<std::endl;
        return;
    }
    ReadDicomSeries(folderName, GetSerieUID(seriesNumber));
}


DicomSeries::~DicomSeries()
{
}


int DicomSeries::SetDataDirectory(std::string folderName)
{
    short i = 0;
    DICOMSeriesNameGeneratorType::Pointer nameGenerator = DICOMSeriesNameGeneratorType::New();
    nameGenerator->SetUseSeriesDetails(true);
    nameGenerator->AddSeriesRestriction("0020|0012"); // Determine number of series
    nameGenerator->SetDirectory(folderName);

    const std::vector<std::string> &seriesUID = nameGenerator->GetSeriesUIDs();
    std::vector<std::string>::const_iterator seriesItr;
    m_iNumberOfSeries = seriesUID.end() - seriesUID.begin();
    m_pSerieName = new std::string[m_iNumberOfSeries];
    for(seriesItr=seriesUID.begin();seriesItr!=seriesUID.end();seriesItr++)
    {
        m_pSerieName[i++] = seriesItr->c_str();
    }

    return 0;
}


int DicomSeries::ReadDicomSeries(std::string folderName, std::string seriesName)
{
    DICOMSeriesType::RegionType region;
    DICOMSeriesType::SizeType size;

    if(folderName == "")
    {
        std::cerr<<"Error! No directory specified. "<<std::endl;
        return -1;
    }
    std::string seriesIdentifier;

    DICOMSeriesReadType::Pointer seriesReader = DICOMSeriesReadType::New();
    seriesReader->SetImageIO(m_pDicomIO);
    
    nameGenerator = DICOMSeriesNameGeneratorType::New();
    nameGenerator->SetUseSeriesDetails(true);
    nameGenerator->AddSeriesRestriction("0020|0012"); // Determine number of series
    nameGenerator->SetDirectory(folderName);

    const std::vector<std::string> &seriesUID = nameGenerator->GetSeriesUIDs();
    std::vector<std::string>::const_iterator seriesItr = seriesUID.begin();

    if(!seriesName.empty())
    {
        seriesIdentifier = seriesName;
    }
    else
    {
        seriesIdentifier = seriesUID.begin()->c_str();
    }

    m_vFileNameContainer = nameGenerator->GetFileNames(seriesIdentifier);
    seriesReader->SetFileNames(m_vFileNameContainer);

    try
    {
        seriesReader->Update();
    }
    catch(itk::ExceptionObject &e)
    {
        std::cerr<<"Exceptions caught in reading file series: "<<std::endl;
        std::cerr<<e<<std::endl;
    }
    region = seriesReader->GetOutput()->GetLargestPossibleRegion();
    size = region.GetSize();
    m_ImageObject = seriesReader->GetOutput();
    ImageRegion3D::SizeType tmp = m_ImageObject->GetLargestPossibleRegion().GetSize();
    m_lWidth = tmp.GetElement(0);
    m_lHeight = tmp.GetElement(1);
    m_lDepth = tmp.GetElement(2);
    m_sPixelData = new short[m_lWidth*m_lHeight*m_lDepth];
    Reshape321<DICOMSeriesType::Pointer, short, DICOMSeriesConstIteratorType>(m_ImageObject, m_sPixelData);
    return 0;
}


int DicomSeries::WriteDicomSeries(std::string folderName)
{
    DICOMSeriesWriteType::Pointer seriesWriter = DICOMSeriesWriteType::New();

    seriesWriter->SetInput(m_ImageObject);
    seriesWriter->SetImageIO(m_pDicomIO);

    nameGenerator->SetOutputDirectory(folderName);
    seriesWriter->SetFileNames(nameGenerator->GetOutputFileNames());

    try
    {
        seriesWriter->Update();
    }
    catch(itk::ExceptionObject &e)
    {
        std::cerr<<"Exeption caught in writing image series: "<<std::endl;
        std::cerr<<e<<std::endl;
    }

    return 0;
}


std::string DicomSeries::GetSerieUID(int number)
{
    return m_pSerieName[number];
}


void DicomSeries::Update()
{
    Reshape123<DICOMSeriesType::Pointer, short, DICOMSeriesIteratorType>(m_sPixelData, m_ImageObject);
}