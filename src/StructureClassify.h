#ifndef STRUCTURECLASSIFY_H
#define STRUCTURECLASSIFY_H


#include <climits>

#include "BaseImage.h"
#include "BaseImageSeries.h"


class StructureClassify
{
public:
    StructureClassify();
    ~StructureClassify();

    void Histogram(BaseImage* input);
    void Histogram(BaseImageSeries* input);
    int* GetHistogram();
private:
    int* m_ipHist;
};


#endif // !STRUCTURECLASSIFY_H
