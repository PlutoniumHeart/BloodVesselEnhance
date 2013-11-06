#ifndef BASEIMAGESERIES_H
#define BASEIMAGESERIES_H


#include <itkImageSeriesReader.h>
#include <itkImageSeriesWriter.h>
#include <itkNumericSeriesFileNames.h>

#include "BaseImage.h"


const unsigned char D3 = 3;

typedef itk::ImageRegion<D3>                    ImageRegion3D;
typedef itk::NumericSeriesFileNames             SeriesNameGeneratorType;


class BaseImageSeries : public BaseImage
{
public:
    BaseImageSeries();
    virtual ~BaseImageSeries();
protected:
    template<typename T1, typename T2, typename T3> 
    int Reshape321(T1 input, T2* destination);

    template<typename T1, typename T2, typename T3> 
    int Reshape123(T2* input, T1 destination);
};

template<typename T1, typename T2, typename T3>
int BaseImageSeries::Reshape321(T1 input, T2* destination)
{
    int i = 0;
    T3 in(input, input->GetLargestPossibleRegion());
    in.GoToBegin();
    while(!in.IsAtEnd())
    {
        destination[i++] = in.Get();
        ++in;
    }
    return 0;
}

template<typename T1, typename T2, typename T3>
int BaseImageSeries::Reshape123(T2* input, T1 destination)
{
    int i = 0;
    T3 out(destination, destination->GetLargestPossibleRegion());
    out.GoToBegin();
    while(!out.IsAtEnd())
    {
        out.Set(*(input+i++));
        ++out;
    }
    return 0;
}


#endif // !BASEIMAGESERIES_H
