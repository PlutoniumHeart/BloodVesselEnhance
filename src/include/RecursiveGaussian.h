#ifndef RECURSIVEGAUSSIAN_H
#define RECURSIVEGAUSSIAN_H


#define BLOCK_DIM 16


#include "GPUMemOps.cuh"


/*************************************************************************
*                                                                        *
*       External links to CUDA kernel calling functions.                 *
*                                                                        *
**************************************************************************/
extern "C" void 
Call_d_recursiveGaussianY( /*Pointers point to input and output data*/
    int *d_src, int *d_dest, 
    /*Image Dimension*/
    int width, int height, int depth, 
    /*Recursive Gaussian filter parameters*/
    float b0, float b1, float b2, float b3, float B, 
    /*CUDA kernel calling parameters*/
    int order, int n, int nthread, 
    /*Boundary condition parameters*/
    float M11, float M12, float M13, 
    float M21, float M22, float M23, 
    float M31, float M32, float M33);

extern "C" void 
Call_d_recursiveGaussianX( /*Pointers point to input and output data*/
    int *d_src, int *d_dest,
    /*Image Dimension*/
    int width, int height, int depth,
    /*Recursive Gaussian filter parameters*/
    float b0, float b1, float b2, float b3, float B, 
    /*CUDA kernel calling parameters*/
    int order, int n, int nthread, 
    /*Boundary condition parameters*/
    float M11, float M12, float M13, 
    float M21, float M22, float M23, 
    float M31, float M32, float M33);
extern "C" void 
Call_d_recursiveGaussianZ(/*Pointers point to input and output data*/
    int *d_src, int *d_dest, 
    /*Image Dimension*/
    int width, int height, int depth, 
    /*Recursive Gaussian filter parameters*/
    float b0, float b1, float b2, float b3, float B, 
    /*CUDA kernel calling parameters*/
    int order, int n, int nthread, 
    /*Boundary condition parameters*/
    float M11, float M12, float M13, 
    float M21, float M22, float M23, 
    float M31, float M32, float M33);
extern "C" void
Call_d_transpose(/*Pointers point to input and output data*/
    int *odata, int *idata,
    /*Image Dimension*/
    int width, int height, int depth, 
    int x, int y);
extern "C" void
Call_d_SobelFilter(/*Pointers point to input and output data*/
    int* d_src, int * d_dest,
    /*Image Dimension*/
    int width, int height, int depth,
    /*CUDA kernel calling parameters*/
    int n, int nthread);

//extern "C" void
//Call_d_SobelNewX(/*Pointers point to input and output data*/
//  int* d_src, int * d_dest,
//  /*Image Dimension*/
//  int width, int height, int depth);
//
//extern "C" void
//Call_d_SobelNewY(/*Pointers point to input and output data*/
//  int* d_src, int * d_dest,
//  /*Image Dimension*/
//  int width, int height, int depth);

extern "C" void setupTexture(int* h_volume, int width, int height, int depth);
extern "C" void deleteTexture();
extern "C" void SobelNew(int* odata, unsigned int pitch, int width, int height, int depth, float fScale);
/************************************************************************
 *                                                                      *
*************************************************************************/




/************************************************************************
 *                                                                      *
 *        CFilter Class template                                        *
 *                                                                      *
 ************************************************************************/
template <class T1, class T2>
class CFilter : public CGPUMemOps
{
public:
    CFilter(void);
    ~CFilter(void);
protected:
    int width;
    int height;
    int depth;
    
    int* h_input;
    int* h_output;
    
    bool flag;
public:
/*Recursive Gaussian filter with CPU*/
    int RecursiveGaussian(int xDirection, int yDirection, int zDirection, T1* d_src, T2* d_dest, /*int order,*/ float sigma);
    int RecursiveGaussianX(int* d_src, int* d_dest, int order, float sigma);
    int RecursiveGaussianY(int* d_src, int* d_dest, int order, float sigma);
    int RecursiveGaussianZ(int* d_src, int* d_dest, int order, float sigma);
    
    int iDivUp(int a, int b);
    
/*****Recursive Gaussian filter with CUDA****************************************/
    
    void CudaRecursiveGaussian(	/*Choose which direction need to be filtered,
                                        mark "1" otherwise mark "0"*/
        int xDirection, int yDirection, int zDirection, 
        /*Pointers point to the input and output data*/
        T1* d_src, T2* d_dest, 
        /*Gaussian filter parameters, order could be 1,
                2, 3; best acuracy could be reached with a 
                sigma no smaller than 1.*/
        /*int order,*/ float sigma);
/********************************************************************************/
    int CudaRecursiveGaussianX(int* d_src, int* d_dest, int order, float sigma);
    int CudaRecursiveGaussianY(int* d_src, int* d_dest, int order, float sigma);
    int CudaRecursiveGaussianZ(int* d_src, int* d_dest, int order, float sigma);
    
/********************************************************************************/
    void ToHost(T2* d_dest);
    
    /*Sobel filter with CPU*/
    void SobelFilter(T1* in, T2* out);
    
    /*Sobel filter with CUDA*/
    void CudaSobelFilter(T1* d_src, T2* d_dest);
    void CudaSobelNew(T1* d_src, T2* d_dest);
    
    //void CudaSobelNew(short* d_src, short* d_dest);
    
    void CudaClear();
    void ResetFlag();
    //void ToDisplay(int c, int w);
public:
    void SetDim(int w, int h, int d);
    int GetWidth();
    int GetHeight();
    int GetDepth();
    
    void SetInput(T1* input);
    int* GetInput();
    int* GetOutput();
    
    void RecursiveGaussianFilter(int xDirection, int yDirection, int zDirection,
                                 T1* d_src, T2* d_dest, /*int order,*/ float sigma);
    void SobelEdgeDetect(T1* d_src, T2* d_dest);
};


#include "../RecursiveGaussian.cpp"


#endif // RECURSIVEGAUSSIAN_H
