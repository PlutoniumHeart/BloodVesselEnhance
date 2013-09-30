#include <chrono>
#include <fstream>

#include "DicomSlice.h"
#include "PngImage.h"
#include "JpegImage.h"
#include "JpegSeries.h"
#include "PngSeries.h"
#include "DicomSeries.h"
#include "StructureClassify.h"


int main()
{
    /*DicomSlice* pDicomSlice = new DicomSlice("../data/2d_DICOM/LI.MR.RESEARCH_PHANTOM.0007.0047.2013.08.06.19.22.15.906250.27825948.IMA");
    PngImage* pPng = new PngImage("../data/2d_DICOM/278.png");
    JpegImage* pJpeg = new JpegImage("../data/2d_DICOM/215.jpg");
    DicomSlice* pEmptyDicom = new DicomSlice(512, 512);
    PngImage* pEmptyPng = new PngImage(320, 320);
    JpegImage* pEmptyJpeg = new JpegImage(320, 320);

    pDicomSlice->WriteDicomFile("../data/output/2DDICOMTest.IMA");
    pPng->WritePngFile("../data/output/278.png");
    pJpeg->WriteJpegFile("../data/output/215.jpeg");

    pPng->CastUnsignedCharToShort(pEmptyDicom);
    pDicomSlice->CastShortToUnsignedChar(pEmptyPng,153,193);
    pDicomSlice->CastShortToUnsignedChar(pEmptyJpeg,153,193);

    pEmptyDicom->WriteDicomFile("../data/output/Empty.IMA");
    pEmptyPng->WritePngFile("../data/output/Empty.png");
    pEmptyJpeg->WriteJpegFile("../data/output/Empty.jpeg");


    delete pDicomSlice;
    delete pPng;
    delete pJpeg;
    delete pEmptyDicom;
    delete pEmptyPng;
    delete pEmptyJpeg;


    JpegSeries* pJpegSeries = new JpegSeries("../data/SegmentedLiver/", "%03d.jpg", 85, 325);
    PngSeries* pPngSeries = new PngSeries("../data/3d_png/", "%03d.png", 1, 10);
    DicomSeries* pDicomSeries = new DicomSeries("../data/AX_RADIAL_PRE_0007/", 0);
    JpegSeries* pJpegEmpty = new JpegSeries(320, 320, 96);

    pJpegSeries->WriteJpegSeries("../data/output/3Djpeg/", "%03d.jpg", 85, 325);
    pPngSeries->WritePngSeries("../data/output/3Dpng/", "%03d.png", 1, 10);
    pDicomSeries->WriteDicomSeries("../data/output/3Ddicom/");

    auto begin = std::chrono::high_resolution_clock::now();
    pDicomSeries->CastShortToUnsignedChar(pJpegEmpty, 72, 102);
    pJpegEmpty->WriteJpegSeries("../data/output/3DJpeg/", "%03d.jpg", 1, 96);
    auto end = std::chrono::high_resolution_clock::now();

    auto ticks = std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count();
    std::cout<<"Cast operation took "<<ticks<<" micro seconds."<<std::endl;

    delete pJpegSeries;
    delete pPngSeries;
    delete pDicomSeries;
    delete pJpegEmpty;*/

    std::fstream file;
    JpegSeries* pJpegSeries = new JpegSeries("../data/SegmentedLiver/", "%03d.jpg", 85, 325);
    StructureClassify* pStructureClassify = new StructureClassify();

    pStructureClassify->Histogram(pJpegSeries);
    /*file.open("../data/output/histogram.txt", std::ios::out);
    if(file.is_open())
    {
        for(int i=0;i<UCHAR_MAX;i++)
        {
            file<<pStructureClassify->GetHistogram()[i]<<", ";
        }
    }
    else
    {
        std::cerr<<"Failed to open file."<<std::endl;
        return -1;
    }
    file.close();*/

    delete pJpegSeries;
    delete pStructureClassify;
    return 0;
}
