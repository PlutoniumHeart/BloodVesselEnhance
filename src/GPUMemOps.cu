#include "GPUMemOps.cuh"

CGPUMemOps::CGPUMemOps()
{
	d_image = NULL;
	d_destination = NULL;
}

CGPUMemOps::~CGPUMemOps()
{
}

void CGPUMemOps::PassToDevice(int* h_data, int size)
{
	size_t memSize = size * sizeof(int);
	checkCudaErrors(cudaMalloc((void**)&d_image, memSize));
	checkCudaErrors(cudaMalloc((void**)&d_destination, memSize));
	checkCudaErrors(cudaMemcpy(d_image, h_data, memSize, cudaMemcpyHostToDevice));
}

void CGPUMemOps::PassToHost(int* h_data, int size, int s)
{
	size_t memSize = size * sizeof(int);
	if(s == 0)
	{
		checkCudaErrors(cudaMemcpy(h_data, d_destination, memSize, cudaMemcpyDeviceToHost));
	}
	else
	{
		checkCudaErrors(cudaMemcpy(h_data, d_image, memSize, cudaMemcpyDeviceToHost));
	}
}

void CGPUMemOps::PassToDevice(short* h_data, int size)
{
	size_t memSize = size * sizeof(short);
	checkCudaErrors(cudaMalloc((void**)&d_image, memSize));
	checkCudaErrors(cudaMalloc((void**)&d_destination, memSize));
	checkCudaErrors(cudaMemcpy(d_image, h_data, memSize, cudaMemcpyHostToDevice));
}

void CGPUMemOps::PassToHost(short* h_data, int size, int s)
{
	size_t memSize = size * sizeof(short);
	if(s == 0)
	{
		checkCudaErrors(cudaMemcpy(h_data, d_destination, memSize, cudaMemcpyDeviceToHost));
	}
	else
	{
		checkCudaErrors(cudaMemcpy(h_data, d_image, memSize, cudaMemcpyDeviceToHost));
	}
}

//void CGPUMemOps::PassToTexture(int* h_data, int width, int height, int depth, cudaChannelFormatDesc channelDesc)
//{
//	size_t memSize = width * height * depth * sizeof(int);
//	cudaExtent volumeSize;
//	volumeSize.width = width;
//	volumeSize.height = height;
//	volumeSize.depth = depth;
//	cudaMalloc3DArray(&cuArray, &channelDesc, width, volumeSize);
//	cudaMemcpy3D();
//	cudaMemcpyToArray(cuArray, 0, 0, h_data, memSize, cudaMemcpyHostToDevice);
//}

int* CGPUMemOps::GetInputImage()
{
	return d_image;
}

int* CGPUMemOps::GetOutputImage()
{
	return d_destination;
}

short* CGPUMemOps::GetShortInput()
{
	return d_shortImage;
}

short* CGPUMemOps::GetShortOutput()
{
	return d_shortDest;
}

void CGPUMemOps::ClearDevice()
{
	cudaFree(d_destination);
	cudaFree(d_image);
}