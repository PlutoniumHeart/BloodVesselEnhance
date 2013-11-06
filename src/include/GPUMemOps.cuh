#pragma once

#include <stdio.h>
#include <iostream>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_functions.h>

/************************************************************************
*																		*
*	CGPUMemOps Class													*
*																		*
*************************************************************************/

class CGPUMemOps
{
public:
	CGPUMemOps();
	~CGPUMemOps();
protected:
	int* d_image;//Pointer point to input image in device
	int* d_destination;//Pointer point to destination image in device
	short* d_shortImage;
	short* d_shortDest;
	cudaArray* cuArray;
public:
	void PassToDevice(int* h_data, int size);
	void PassToHost(int* h_data, int size, int s);
	void PassToDevice(short* h_data, int size);
	void PassToHost(short* h_data, int size, int s);
	//void PassToTexture(int* h_data, int width, int height, int depth, cudaChannelFormatDesc channelDesc);

	int* GetInputImage();
	int* GetOutputImage();
	short* GetShortInput();
	short* GetShortOutput();
	void ClearDevice();
};