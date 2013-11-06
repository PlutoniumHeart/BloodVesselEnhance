//#include <chrono>
#include <fstream>

#include "DicomSlice.h"
#include "PngImage.h"
#include "JpegImage.h"
#include "JpegSeries.h"
#include "PngSeries.h"
#include "DicomSeries.h"
#include "StructureClassify.h"
#include "RecursiveGaussian.h"


int main()
{
    std::fstream file;
    JpegSeries* pJpegSeries = new JpegSeries("../data/SegmentedLiver/", "%03d.jpg", 85, 325);
    JpegSeries* pResult = new JpegSeries(pJpegSeries->GetImageWidth(), pJpegSeries->GetImageHeight(), pJpegSeries->GetImageDepth());

    StructureClassify* pStructureClassify = new StructureClassify();
    pStructureClassify->SetInput(pJpegSeries);

    //auto begin = std::chrono::high_resolution_clock::now();
    pStructureClassify->ClassifyBrightTube(5, 80, 2);
    //auto end = std::chrono::high_resolution_clock::now();

    unsigned char *p = pStructureClassify->GetResult();
    unsigned char *pp = pResult->GetUnsignedCharPixelData();
    for(int i=0;i<pJpegSeries->GetImageWidth()*pJpegSeries->GetImageHeight()*pJpegSeries->GetImageDepth();i++)
    {
        pp[i] = p[i];
    }

    /*int max = 0;
    for(int i=0;i<pJpegSeries->GetImageWidth()*pJpegSeries->GetImageHeight()*pJpegSeries->GetImageDepth();i++)
    {
        if(pp[i]>max)
        {
            max = pp[i];
        }
    }*/

    pResult->Update();

    //auto ticks = std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count();
    //std::cout<<"Cast operation took "<<ticks<<" micro seconds."<<std::endl;
    pResult->WriteJpegSeries("../data/output/3Djpeg/", "%03d.jpg", 85, 325);

    delete pJpegSeries;
    delete pStructureClassify;
    return 0;
}
