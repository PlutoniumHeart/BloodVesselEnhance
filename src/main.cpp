#include "DicomSlice.h"
#include "PngImage.h"
#include "JpegImage.h"
#include "JpegSeries.h"


int main()
{
    DicomSlice* pDicomSlice = new DicomSlice("../data/2d_DICOM/LI.MR.RESEARCH_PHANTOM.0007.0047.2013.08.06.19.22.15.906250.27825948.IMA");
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
    pJpegSeries->WriteJpegSeries("../data/output/3Djpeg/", "%03d.jpg", 85, 325);

    delete pJpegSeries;


    return 0;
}
