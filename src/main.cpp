#include <fstream>

#include "DicomSlice.h"
#include "PngImage.h"
#include "JpegImage.h"
#include "JpegSeries.h"
#include "PngSeries.h"
#include "DicomSeries.h"
#include "StructureClassify.h"
#include "RecursiveGaussian.h"
#include "Util.h"


int main()
{
    std::fstream file;
    timespec time1, time2;
    
    JpegSeries* pJpegSeries = new JpegSeries("../data/SegmentedLiver/", "%03d.jpg", 85, 325);
    DicomSeries* pDicomResult = new DicomSeries(pJpegSeries->GetImageWidth(), pJpegSeries->GetImageHeight(), pJpegSeries->GetImageDepth());

    StructureClassify* pStructureClassify = new StructureClassify();
    pStructureClassify->SetInput(pJpegSeries);

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    pStructureClassify->ClassifyBrightTube(5, 80, 2);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);

    unsigned char *p = pStructureClassify->GetResult();
    short* pp = pDicomResult->GetShortPixelData();
    for(int i=0;i<pJpegSeries->GetImageWidth()*pJpegSeries->GetImageHeight()*pJpegSeries->GetImageDepth();i++)
    {
        pp[i] = p[i];
    }

    pDicomResult->Update();
#ifdef _MSC_VER
    std::cout<<"Operation took: "<<diff(time1, time2).tv_sec<<"."<<diff(time1, time2).tv_usec<<"s"<<std::endl;
#else
    std::cout<<"Operation took: "<<diff(time1, time2).tv_sec<<"."<<diff(time1, time2).tv_nsec<<"s"<<std::endl;
#endif
    pDicomResult->WriteDicomSeries("../data/output/Dicom_out/", "%03d.IMG", 85, 325);

    delete pJpegSeries;
    delete pStructureClassify;
    return 0;
}
