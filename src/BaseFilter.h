#ifndef BASEFILTER_H
#define BASEFILTER_H


#include "BaseImage.h"
#include "BaseImageSeries.h"


enum Direction
{
    XDirection=1,
    YDirection=2,
    ZDirection=4,

    Undefined,
};


enum class DataType
{
    Short,
    UChar,
    
    Undefined,
};


class BaseFilter
{
public:
    BaseFilter();
    ~BaseFilter();

    virtual void Filter(int order, float sigma, int dir) = 0;
protected:
    unsigned long m_lWidth;
    unsigned long m_lHeight;
    unsigned long m_lDepth;
};


#endif // !BASEFILTER_H
