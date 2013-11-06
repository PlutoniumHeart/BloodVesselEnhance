#ifndef CUDARECURSIVEGAUSSIAN_KERNEL_CU
#define CUDARECURSIVEGAUSSIAN_KERNEL_CU

#include <helper_math.h>
#include <stdio.h>
#include <iostream>

#define BLOCK_DIM 16

__global__ void d_transpose(int *odata, int *idata, int width, int height, int depth)
{
    __shared__ int block[BLOCK_DIM][BLOCK_DIM+1];
    //int blockIdx_y, blockIdx_x;
    // do diagonal reordering
	/*if (width == height)
	{
		blockIdx_y = blockIdx.x;
		blockIdx_x = (blockIdx.x+blockIdx.y)%gridDim.x;
	}
	else
	{
		int bid = blockIdx.x + gridDim.x*blockIdx.y;
		blockIdx_y = bid%gridDim.y;
		blockIdx_x = ((bid/gridDim.y)+blockIdx_y)%gridDim.x;
	}*/

    int xIndex, yIndex, zIndex;
    
    xIndex = blockIdx.x * BLOCK_DIM + threadIdx.x;
    yIndex = (blockIdx.y * BLOCK_DIM + threadIdx.y) % height;
    zIndex = (blockIdx.y * BLOCK_DIM + threadIdx.y) / height;
    
    int slice = height * width;
    
    if((xIndex < width) && (yIndex < height) && (zIndex < depth))
    {
		block[threadIdx.y][threadIdx.x] = idata[zIndex * slice + yIndex * width + xIndex];
    }
    __syncthreads();
    
    xIndex = (blockIdx.y * BLOCK_DIM + threadIdx.y) % height;
    yIndex = (blockIdx.y * BLOCK_DIM + threadIdx.y) / height;
    zIndex = blockIdx.x * BLOCK_DIM + threadIdx.x;
    
    
    slice = height * depth;
    
    if ((xIndex < height) && (yIndex < depth) && (zIndex < width))
    {
		odata[zIndex * slice + yIndex * height + xIndex] = block[threadIdx.y][threadIdx.x];
    }
}


extern "C" void Call_d_transpose(int *odata, int *idata, int width, int height, int depth, int x, int y)
{
    dim3 grid(x, y, 1);
    dim3 threads(BLOCK_DIM, BLOCK_DIM, 1);
    d_transpose<<< grid, threads >>>(odata, idata, width, height, depth);
}
































__global__ void d_recursiveGaussianY(int *d_src, int *d_dest, int depth, int height, int width, float b0, float b1, float b2, float b3, float B, int order, float M11, float M12, float M13, float M21, float M22, float M23, float M31, float M32, float M33)
{
	float wP1 = 0.f, wP2 = 0.f, wP3 = 0.f;
	int y = 0;
	float outF1 = 0.f, outF2 = 0.f, outF3 = 0.f;
	unsigned int x = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int yy = blockIdx.y*width*height;
	if(x > width)
		return;
	
	d_src += x + yy;
	d_dest += x + yy;
	
	wP1 = (float)*d_src/sqrt(B); wP2 = wP1; wP3 = wP1;
	
	switch (order) 
	{
		case 0: 
		{
			for(y=0;y<height;y++)
			{
				float xC = (float)*d_src;
				float wC = (float)(xC - b1*wP1 - b2*wP2 - b3*wP3)/b0;
				*d_dest = (int)wC;
				d_src += width; d_dest += width;
				wP3 = wP2; wP2 = wP1; wP1 = wC;
			}
			
			d_src -= width;
			d_dest -= width;
			
			
			float up = (float)*d_src/(1.0+b1+b2+b3);
			float vp = (float)up/(1.0+b1+b2+b3);
			
			float out = 0.f;
			out = (float)M11*(*d_dest-up) + M12*(*(d_dest-width)-up) + M13*(*(d_dest-2*width)-up)+vp;
			outF1 = (float)M21*(*d_dest-up) + M22*(*(d_dest-width)-up) + M23*(*(d_dest-2*width)-up)+vp;
			outF2 = (float)M31*(*d_dest-up) + M32*(*(d_dest-width)-up) + M33*(*(d_dest-2*width)-up)+vp;
			out *= B; outF1 *= B; outF2 *= B;
			outF3 = outF2; outF2 = outF1; outF1 = out;
			
			*d_dest = (int)out;
			d_src -= width;
			d_dest -= width;
						
			for(y=height-1-1;y>=0;y--)
			{
				float wC = (float)*d_dest;
				out = (float)(B*wC - b1*outF1 - b2*outF2 - b3*outF3)/b0;
				*d_dest = (int)out;
				d_src -= width; d_dest -= width;
				outF3 = outF2; outF2 = outF1; outF1 = out;
			}
		} break;

		case 1:
		{
			float xP1 = (float)*(d_src);
			float xF1 = (float)*(d_src + 1*width);
			
			wP1 = (float)(xF1 - xP1)/2.0*(1.0+b1+b2+b3); wP3 = wP2 = wP1;
			
			float wC = ((xF1- xP1)/2.0 - b1*wP1 - b2*wP1 - b3*wP1)/b0;
			*d_dest = (int)wC;
			d_src += width; d_dest += width;
			wP3 = wP2; wP2 = wP1; wP1 = wC;
			
			for(y=0;y<height-1-1;y++)
			{
				xP1 = (float)*(d_src - width);
				xF1 = (float)*(d_src + width);
				wC = (float)((xF1- xP1)/2.0 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
				*d_dest = (int)wC;
				d_src += width; d_dest += width;
				wP3 = wP2; wP2 = wP1; wP1 = wC;
			}
			
			xP1 = (float)*(d_src - width);
			xF1 = (float)*d_src;
			wC = (float)((xF1- xP1)/2.0 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
			*d_dest = (int)wC;
			wP3 = wP2; wP2 = wP1; wP1 = wC;
			
			float up = (*d_src - *(d_src-width))/2.0*(1.0+b1+b2+b3);
			float vp = up/(1.0+b1+b2+b3);
			
			float out = 0.f;
			out = (float)M11*(*d_dest-up) + M12*(*(d_dest-width)-up) + M13*(*(d_dest-2*width)-up)+vp;
			outF1 = (float)M21*(*d_dest-up) + M22*(*(d_dest-width)-up) + M23*(*(d_dest-2*width)-up)+vp;
			outF2 = (float)M31*(*d_dest-up) + M32*(*(d_dest-width)-up) + M33*(*(d_dest-2*width)-up)+vp;
			out *= B; outF1 *= B; outF2 *= B;
			outF3 = outF2; outF2 = outF1; outF1 = out;
			
			*d_dest = (int)out;
			
			d_src -= width;
			d_dest -= width;
						
			for(y=height-2;y>=0;y--)
			{
				float wC = (float)*d_dest;
				out = (float)(B*wC - b1*outF1 - b2*outF2 - b3*outF3)/b0;
				*d_dest = (int)out;
				d_src -= width; d_dest -= width;
				outF3 = outF2; outF2 = outF1; outF1 = out;
			}		
		} break;

		case 2: 
		{
			float xP1 = (float)*d_src;
			float xC = (float)*d_src;
			float xF1 = (float)*(d_src+width);
			wP1 = 0.0/(1.0+b1+b2+b3); wP3 = wP2 = wP1;
			
			float wC = (float)((xF1 - 2*xC + xP1) - b1*wP1 - b2*wP2 - b3*wP3)/b0;
			*d_dest = (int)wC;
			d_src += width; d_dest += width;
			wP3 = wP2; wP2 = wP1; wP1 = wC;
			
			for(y=0;y<height-2;y++)
			{
				xC = (float)*d_src;
				xP1 = (float)*(d_src-width);
				xF1 = (float)*(d_src+width);
				wC = (float)(xF1 - 2*xC + xP1 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
				*d_dest = (int)wC;
				d_src += width; d_dest += width;
				wP3 = wP2; wP2 = wP1; wP1 = wC;
			}
			
			xC = (float)*d_src;
			xP1 = (float)*(d_src-width);
			xF1 = (float)*d_src;
			wC = (float)(xF1 - 2*xC + xP1 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
			*d_dest = (int)wC;
			wP3 = wP2; wP2 = wP1; wP1 = wC;
			
			float up = 0;
			float vp = up/(1.0+b1+b2+b3);
			
			float out = 0.f;
			out = (float)M11*(*d_dest-up) + M12*(*(d_dest-width)-up) + M13*(*(d_dest-2*width)-up)+vp;
			outF1 = (float)M21*(*d_dest-up) + M22*(*(d_dest-width)-up) + M23*(*(d_dest-2*width)-up)+vp;
			outF2 = (float)M31*(*d_dest-up) + M32*(*(d_dest-width)-up) + M33*(*(d_dest-2*width)-up)+vp;
			out *= B; outF1 *= B; outF2 *= B;
			outF3 = outF2; outF2 = outF1; outF1 = out;
			
			*d_dest = (int)out;
			
			d_src -= width;
			d_dest -= width;
			
			for(y=height-2;y>=0;y--)
			{
				float wC = (float)*d_dest;
				out = (float)(B*wC - b1*outF1 - b2*outF2 - b3*outF3)/b0;
				*d_dest = (int)out;
				d_src -= width; d_dest -= width;
				outF3 = outF2; outF2 = outF1; outF1 = out;
			}
		} break;
	}
}

extern "C" void Call_d_recursiveGaussianY(int *d_src, int *d_dest, int width, int height, int depth, float b0, float b1, float b2, float b3, float B, int order, int n, int nthread, float M11, float M12, float M13, float M21, float M22, float M23, float M31, float M32, float M33)
{
	d_recursiveGaussianY<<<dim3(n, depth), nthread>>>(d_src, d_dest, depth, height, width, b0, b1, b2, b3, B, order, M11, M12, M13, M21, M22, M23, M31, M32, M33);
}

__global__ void d_recursiveGaussianX(int *d_src, int *d_dest, int depth, int height, int width, float b0, float b1, float b2, float b3, float B, int order, float M11, float M12, float M13, float M21, float M22, float M23, float M31, float M32, float M33)
{
	float wP1 = 0.f, wP2 = 0.f, wP3 = 0.f;
	int y = 0;
	float outF1 = 0.f, outF2 = 0.f, outF3 = 0.f;
	unsigned int x = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int yy = blockIdx.y*width*height;
	if(x > width)
		return;
		
	d_src += x*width + yy;
	d_dest += x*width + yy;

	wP1 = (float)*d_src/sqrt(B); wP2 = wP1; wP3 = wP1;
	
	switch (order) 
	{
		case 0: 
		{
			for(y=0;y<width;y++)
			{
				float xC = (float)*d_src;
				float wC = (float)(xC - b1*wP1 - b2*wP2 - b3*wP3)/b0;
				*d_dest = (int)wC;
				d_src += 1; d_dest += 1;
				wP3 = wP2; wP2 = wP1; wP1 = wC;
			}
			
			d_src -= 1;
			d_dest -= 1;
			
			
			float up = (float)*d_src/(1.0+b1+b2+b3);
			float vp = (float)up/(1.0+b1+b2+b3);
			
			float out = 0.f;
			out = (float)M11*(*d_dest-up) + M12*(*(d_dest-1)-up) + M13*(*(d_dest-2*1)-up)+vp;
			outF1 = (float)M21*(*d_dest-up) + M22*(*(d_dest-1)-up) + M23*(*(d_dest-2*1)-up)+vp;
			outF2 = (float)M31*(*d_dest-up) + M32*(*(d_dest-1)-up) + M33*(*(d_dest-2*1)-up)+vp;
			out *= B; outF1 *= B; outF2 *= B;
			outF3 = outF2; outF2 = outF1; outF1 = out;
			
			*d_dest = (int)out;
			d_src -= 1;
			d_dest -= 1;
						
			for(y=width-1-1;y>=0;y--)
			{
				float wC = (float)*d_dest;
				out = (float)(B*wC - b1*outF1 - b2*outF2 - b3*outF3)/b0;
				*d_dest = (int)out;
				d_src -= 1; d_dest -= 1;
				outF3 = outF2; outF2 = outF1; outF1 = out;
			}
		} break;

		case 1:
		{
			float xP1 = (float)*(d_src);
			float xF1 = (float)*(d_src + 1*1);
			
			wP1 = (float)(xF1 - xP1)/2.0*(1.0+b1+b2+b3); wP3 = wP2 = wP1;
			
			float wC = ((xF1- xP1)/2.0 - b1*wP1 - b2*wP1 - b3*wP1)/b0;
			*d_dest = (int)wC;
			d_src += 1; d_dest += 1;
			wP3 = wP2; wP2 = wP1; wP1 = wC;
			
			for(y=0;y<width-2;y++)
			{
				xP1 = (float)*(d_src - 1);
				xF1 = (float)*(d_src + 1);
				wC = (float)((xF1- xP1)/2.0 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
				*d_dest = (int)wC;
				d_src += 1; d_dest += 1;
				wP3 = wP2; wP2 = wP1; wP1 = wC;
			}
			
			xP1 = (float)*(d_src - 1);
			xF1 = (float)*d_src;
			wC = (float)((xF1- xP1)/2.0 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
			*d_dest = (int)wC;
			wP3 = wP2; wP2 = wP1; wP1 = wC;
			
			float up = (*d_src - *(d_src-1))/2.0*(1.0+b1+b2+b3);
			//float up = 0;
			float vp = up/(1.0+b1+b2+b3);
			
			float out = 0.f;
			out = (float)M11*(*d_dest-up) + M12*(*(d_dest-1)-up) + M13*(*(d_dest-2*1)-up)+vp;
			outF1 = (float)M21*(*d_dest-up) + M22*(*(d_dest-1)-up) + M23*(*(d_dest-2*1)-up)+vp;
			outF2 = (float)M31*(*d_dest-up) + M32*(*(d_dest-1)-up) + M33*(*(d_dest-2*1)-up)+vp;
			out *= B; outF1 *= B; outF2 *= B;
			outF3 = outF2; outF2 = outF1; outF1 = out;
			
			*d_dest = (int)out;
			
			d_src -= 1;
			d_dest -= 1;
						
			for(y=width-2;y>=0;y--)
			{
				float wC = (float)*d_dest;
				out = (float)(B*wC - b1*outF1 - b2*outF2 - b3*outF3)/b0;
				*d_dest = (int)out;
				d_src -= 1; d_dest -= 1;
				outF3 = outF2; outF2 = outF1; outF1 = out;
			}		
		} break;

		case 2: 
		{
			float xP1 = (float)*d_src;
			float xC = (float)*d_src;
			float xF1 = (float)*(d_src+1);
			wP1 = 0.0/(1.0+b1+b2+b3); wP3 = wP2 = wP1;
			
			float wC = (int)((xF1 - 2*xC + xP1) - b1*wP1 - b2*wP2 - b3*wP3)/b0;
			*d_dest = (int)wC;
			d_src += 1; d_dest += 1;
			wP3 = wP2; wP2 = wP1; wP1 = wC;
			
			for(y=0;y<width-2;y++)
			{
				xC = (float)*d_src;
				xP1 = (float)*(d_src-1);
				xF1 = (float)*(d_src+1);
				wC = (float)(xF1 - 2*xC + xP1 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
				*d_dest = (int)wC;
				d_src += 1; d_dest += 1;
				wP3 = wP2; wP2 = wP1; wP1 = wC;
			}
			
			xC = (float)*d_src;
			xP1 = (float)*(d_src-1);
			xF1 = (float)*d_src;
			wC = (float)(xF1 - 2*xC + xP1 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
			*d_dest = (int)wC;
			wP3 = wP2; wP2 = wP1; wP1 = wC;
			
			float up = 0;
			float vp = up/(1.0+b1+b2+b3);
			
			float out = 0.f;
			out = (float)M11*(*d_dest-up) + M12*(*(d_dest-1)-up) + M13*(*(d_dest-2*1)-up)+vp;
			outF1 = (float)M21*(*d_dest-up) + M22*(*(d_dest-1)-up) + M23*(*(d_dest-2*1)-up)+vp;
			outF2 = (float)M31*(*d_dest-up) + M32*(*(d_dest-1)-up) + M33*(*(d_dest-2*1)-up)+vp;
			out *= B; outF1 *= B; outF2 *= B;
			outF3 = outF2; outF2 = outF1; outF1 = out;
			
			*d_dest = (int)out;
			
			d_src -= 1;
			d_dest -= 1;
			
			for(y=width-2;y>=0;y--)
			{
				float wC = (float)*d_dest;
				out = (float)(B*wC - b1*outF1 - b2*outF2 - b3*outF3)/b0;
				*d_dest = (int)out;
				d_src -= 1; d_dest -= 1;
				outF3 = outF2; outF2 = outF1; outF1 = out;
			}
		} break;
	}
}

extern "C" void Call_d_recursiveGaussianX(int *d_src, int *d_dest, int width, int height, int depth, float b0, float b1, float b2, float b3, float B, int order, int n, int nthread, float M11, float M12, float M13, float M21, float M22, float M23, float M31, float M32, float M33)
{
	d_recursiveGaussianX<<<dim3(n, depth), nthread>>>(d_src, d_dest, depth, height, width, b0, b1, b2, b3, B, order, M11, M12, M13, M21, M22, M23, M31, M32, M33);
}

__global__ void d_recursiveGaussianZ(int* d_src, int* d_dest, int depth, int height, int width, float b0, float b1, float b2, float b3, float B, int order, float M11, float M12, float M13, float M21, float M22, float M23, float M31, float M32, float M33)
{
	float wP1 = 0.f, wP2 = 0.f, wP3 = 0.f;
	int y = 0;
	float outF1 = 0.f, outF2 = 0.f, outF3 = 0.f;
	
	unsigned int x = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int yy = blockIdx.y*width;
	if(x > width)
		return;
	
	d_src += x + yy;
	d_dest += x + yy;
	
	wP1 = (float)*d_src/sqrt(B); wP2 = wP1; wP3 = wP1;
	
	switch (order) 
	{
		case 0: 
		{
			for(y=0;y<depth;y++)
			{
				float xC = (float)*d_src;
				float wC = (float)(xC - b1*wP1 - b2*wP2 - b3*wP3)/b0;
				*d_dest = (int)wC;
				d_src += width*height; d_dest += width*height;
				wP3 = wP2; wP2 = wP1; wP1 = wC;
			}
			
			d_src -= width*height;
			d_dest -= width*height;
			
			
			float up = (float)*d_src/(1.0+b1+b2+b3);
			float vp = (float)up/(1.0+b1+b2+b3);
			
			float out = 0.f;
			out = (float)M11*(*d_dest-up) + M12*(*(d_dest-width*height)-up) + M13*(*(d_dest-2*width*height)-up)+vp;
			outF1 = (float)M21*(*d_dest-up) + M22*(*(d_dest-width*height)-up) + M23*(*(d_dest-2*width*height)-up)+vp;
			outF2 = (float)M31*(*d_dest-up) + M32*(*(d_dest-width*height)-up) + M33*(*(d_dest-2*width*height)-up)+vp;
			out *= B; outF1 *= B; outF2 *= B;
			outF3 = outF2; outF2 = outF1; outF1 = out;
			
			*d_dest = (int)out;
			d_src -= width*height;
			d_dest -= width*height;
						
			for(y=depth-1-1;y>=0;y--)
			{
				float wC = (float)*d_dest;
				out = (float)(B*wC - b1*outF1 - b2*outF2 - b3*outF3)/b0;
				*d_dest = (int)out;
				d_src -= width*height; d_dest -= width*height;
				outF3 = outF2; outF2 = outF1; outF1 = out;
			}
		} break;

		case 1:
		{
			float xP1 = (float)*(d_src);
			float xF1 = (float)*(d_src + 1*width*height);
			
			wP1 = (float)(xF1 - xP1)/2.0*(1.0+b1+b2+b3); wP3 = wP2 = wP1;
			
			float wC = ((xF1- xP1)/2.0 - b1*wP1 - b2*wP1 - b3*wP1)/b0;
			*d_dest = (int)wC;
			d_src += width*height; d_dest += width*height;
			wP3 = wP2; wP2 = wP1; wP1 = wC;
			
			for(y=0;y<depth-2;y++)
			{
				xP1 = (float)*(d_src - width*height);
				xF1 = (float)*(d_src + width*height);
				wC = (float)((xF1- xP1)/2.0 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
				*d_dest = (int)wC;
				d_src += width*height; d_dest += width*height;
				wP3 = wP2; wP2 = wP1; wP1 = wC;
			}
			
			xP1 = (float)*(d_src - width*height);
			xF1 = (float)*d_src;
			wC = (float)((xF1- xP1)/2.0 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
			*d_dest = (int)wC;
			wP3 = wP2; wP2 = wP1; wP1 = wC;
						
			float up = (*d_src - *(d_src-width*height))/2.0*(1.0+b1+b2+b3);
			//float up = 0;
			float vp = up/(1.0+b1+b2+b3);
			
			float out = 0.f;
			out = (float)M11*(*d_dest-up) + M12*(*(d_dest-width*height)-up) + M13*(*(d_dest-2*width*height)-up)+vp;
			outF1 = (float)M21*(*d_dest-up) + M22*(*(d_dest-width*height)-up) + M23*(*(d_dest-2*width*height)-up)+vp;
			outF2 = (float)M31*(*d_dest-up) + M32*(*(d_dest-width*height)-up) + M33*(*(d_dest-2*width*height)-up)+vp;
			out *= B; outF1 *= B; outF2 *= B;
			outF3 = outF2; outF2 = outF1; outF1 = out;
			
			*d_dest = (int)out;
			
			d_src -= width*height;
			d_dest -= width*height;
						
			for(y=depth-2;y>=0;y--)
			{
				float wC = (float)*d_dest;
				out = (float)(B*wC - b1*outF1 - b2*outF2 - b3*outF3)/b0;
				*d_dest = (int)out;
				d_src -= width*height; d_dest -= width*height;
				outF3 = outF2; outF2 = outF1; outF1 = out;
			}		
		} break;

		case 2: 
		{
			float xP1 = (float)*d_src;
			float xC = (float)*d_src;
			float xF1 = (float)*(d_src+width*height);
			wP1 = 0.0/(1.0+b1+b2+b3); wP3 = wP2 = wP1;
			
			float wC = (float)((xF1 - 2*xC + xP1) - b1*wP1 - b2*wP2 - b3*wP3)/b0;
			*d_dest = (int)wC;
			d_src += width*height; d_dest += width*height;
			wP3 = wP2; wP2 = wP1; wP1 = wC;
			
			for(y=0;y<depth-1-1;y++)
			{
				xC = (float)*d_src;
				xP1 = (float)*(d_src-width*height);
				xF1 = (float)*(d_src+width*height);
				wC = (float)(xF1 - 2*xC + xP1 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
				*d_dest = (int)wC;
				d_src += width*height; d_dest += width*height;
				wP3 = wP2; wP2 = wP1; wP1 = wC;
			}
			
			xC = (float)*d_src;
			xP1 = (float)*(d_src-width*height);
			xF1 = (float)*d_src;
			wC = (float)(xF1 - 2*xC + xP1 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
			*d_dest = (int)wC;
			wP3 = wP2; wP2 = wP1; wP1 = wC;
					
			float up = 0;
			float vp = up/(1.0+b1+b2+b3);
			
			float out = 0.f;
			out = (float)M11*(*d_dest-up) + M12*(*(d_dest-width*height)-up) + M13*(*(d_dest-2*width*height)-up)+vp;
			outF1 = (float)M21*(*d_dest-up) + M22*(*(d_dest-width*height)-up) + M23*(*(d_dest-2*width*height)-up)+vp;
			outF2 = (float)M31*(*d_dest-up) + M32*(*(d_dest-width*height)-up) + M33*(*(d_dest-2*width*height)-up)+vp;
			out *= B; outF1 *= B; outF2 *= B;
			outF3 = outF2; outF2 = outF1; outF1 = out;
			
			*d_dest = (int)out;
			
			d_src -= width*height;
			d_dest -= width*height;
			
			for(y=depth-2;y>=0;y--)
			{
				float wC = (float)*d_dest;
				out = (float)(B*wC - b1*outF1 - b2*outF2 - b3*outF3)/b0;
				*d_dest = (int)out;
				d_src -= width*height; d_dest -= width*height;
				outF3 = outF2; outF2 = outF1; outF1 = out;
			}
		} break;
	}
}

extern "C" void Call_d_recursiveGaussianZ(int *d_src, int *d_dest, int width, int height, int depth, float b0, float b1, float b2, float b3, float B, int order, int n, int nthread, float M11, float M12, float M13, float M21, float M22, float M23, float M31, float M32, float M33)
{
	d_recursiveGaussianZ<<<dim3(n, height), nthread>>>(d_src, d_dest, depth, height, width, b0, b1, b2, b3, B, order, M11, M12, M13, M21, M22, M23, M31, M32, M33);
}


#endif //CUDARECURSIVEGAUSSIAN_KERNEL_CU

