#ifndef STRUCTURECLASSIFY_H
#define STRUCTURECLASSIFY_H


#include <climits>

#include "BaseImage.h"
#include "BaseImageSeries.h"
#include "RecursiveGaussian.h"
#include "Matrix.h"


class StructureClassify
{
public:
    StructureClassify();
    ~StructureClassify();

    void SetInput(BaseImageSeries* input);
    int* GetHistogram();
    int ClassifyBrightTube(int iterations, float rMax, float rMin);
    unsigned char* GetResult();
protected:
    void Histogram(BaseImageSeries* input);
private:
    int m_iWidth;
    int m_iHeight;
    int m_iDepth;
    unsigned char* m_pInput;
    int* m_ipHist;
    short* m_pTemp;
    unsigned char* m_pResult;
    unsigned char* m_pOutput;

    CFilter<unsigned char, short>* m_Filter;
};


#endif // STRUCTURECLASSIFY_H
