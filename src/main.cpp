#include "DicomSlice.h"
#include "PngImage.h"
#include "JpegImage.h"
#include "JpegSeries.h"
#include "PngSeries.h"
#include "DicomSeries.h"


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
    delete pEmptyJpeg;*/


    JpegSeries* pJpegSeries = new JpegSeries("../data/SegmentedLiver/", "%03d.jpg", 85, 325);
    PngSeries* pPngSeries = new PngSeries("../data/3d_png/", "%03d.png", 1, 10);
    DicomSeries* pDicomSeries = new DicomSeries("../data/AX_RADIAL_PRE_0007/", 0);

    pJpegSeries->WriteJpegSeries("../data/output/3Djpeg/", "%03d.jpg", 85, 325);
    pPngSeries->WritePngSeries("../data/output/3Dpng/", "%03d.png", 1, 10);
    pDicomSeries->WriteDicomSeries("../data/output/3Ddicom/");

    delete pJpegSeries;
    delete pPngSeries;
    delete pDicomSeries;

    return 0;
}
