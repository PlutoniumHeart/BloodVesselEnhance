#ifndef RECURSIVEGAUSSIAN_H
#define RECURSIVEGAUSSIAN_H


#include "BaseFilter.h"


class RecursiveGaussian : public BaseFilter
{
public:
    RecursiveGaussian();
    ~RecursiveGaussian();

    void SetInput(BaseImage* image);
    void WriteToImage(BaseImage* output);
    virtual void Filter(int order, float sigma, int dir);
protected:
    int RecursiveGaussianX(float* pData, float * outData, int order, float sigma);
    int RecursiveGaussianY(float* pData, float * outData, int order, float sigma);
    int RecursiveGaussianZ(float* pData, float * outData, int order, float sigma);
protected:
    DataType m_eDataType;
    BaseImage* m_pInput;
    float* m_fInputData;
    float* m_fPixelData;
    float m_Max;
    float m_Min;
    bool m_bFlip;
};

#endif // !RECURSIVEGAUSSIAN_H
