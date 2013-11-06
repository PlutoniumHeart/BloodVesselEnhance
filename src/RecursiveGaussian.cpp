#ifndef Filter_CPP
#define Filter_CPP
#include "RecursiveGaussian.h"


template <class T1, class T2> CFilter<T1,T2>::CFilter(void)
{
	width = 0;
	height = 0;
	depth = 0;
	h_output = NULL;
	h_input = NULL;
	flag  = true;
}

template <class T1,class T2> CFilter<T1,T2>::~CFilter(void)
{
}

template <class T1,class T2> void CFilter<T1,T2>::SetDim(int w, int h, int d)
{
	width = w;
	height = h;
	depth = d;
}

template <class T1,class T2> int CFilter<T1,T2>::GetWidth()
{
	return width;
}

template <class T1,class T2> int CFilter<T1,T2>::GetHeight()
{
	return height;
}

template <class T1,class T2> int CFilter<T1,T2>::GetDepth()
{
	return depth;
}

template <class T1, class T2> void CFilter<T1,T2>::SetInput(T1* input)
{
	int i = 0;
	delete [] h_input;
	h_input = new int[width*height*depth];
	h_output = new int[width*height*depth];
	for(i=0;i<width*height*depth;i++)
	{
		h_input[i] = input[i];
	}
}

template <class T1, class T2> int* CFilter<T1,T2>::GetInput()
{
	return h_input;
}

template <class T1, class T2> int* CFilter<T1, T2>::GetOutput()
{
	return h_output;
}

template <class T1, class T2> int CFilter<T1, T2>::RecursiveGaussian(int xDirection, int yDirection, int zDirection, T1* d_src, T2* d_dest,/* int order,*/ float sigma)
{
	SetInput(d_src);
	if(xDirection >= 0)
	{
		if(flag)
			RecursiveGaussianX(h_input, h_output, xDirection, sigma);
		else
			RecursiveGaussianX(h_output, h_input, xDirection, sigma);
		flag = !flag;
	}
	if(yDirection >= 0)
	{
		if(flag)
			RecursiveGaussianY(h_input, h_output, yDirection, sigma);
		else
			RecursiveGaussianY(h_output, h_input, yDirection, sigma);
		flag = !flag;
	}
	if(zDirection >= 0)
	{
		if(flag)
			RecursiveGaussianZ(h_input, h_output, zDirection, sigma);
		else
			RecursiveGaussianZ(h_output, h_input, zDirection, sigma);
		flag = !flag;
	}
	int i = 0;
	if(flag)
	{
		for(i=0;i<width*height*depth;i++)
		{
			d_dest[i] = h_input[i];
		}
	}
	else
	{
		for(i=0;i<width*height*depth;i++)
		{
			d_dest[i] = h_output[i];
		}
	}
	delete [] h_output;
	ResetFlag();
	return 0;
}

template <class T1, class T2> int CFilter<T1, T2>::RecursiveGaussianX(int* d_src, int* d_dest, int order, float sigma)
{
	//const float nsigma = sigma < 0.5f ? 0.5f : sigma;
	float q = 0, m0 = 1.16680, m1 = 1.10783, m2 = 1.40586;
	float m1Sqr = m1*m1, m2Sqr = m2*m2;
	if(sigma < 1)
	{
		//std::cout<<"sigma value too small, result may be inaccurate! considering sigma as allowed minmun 0.5!"<<std::endl;
		q = 1.31564*(sqrt(1+0.490811*sigma*sigma)-1);
	}
	else
	{
		q = 1.31564*(sqrt(1+0.490811*sigma*sigma)-1);
	}
	
	const float
	qSqr = q*q,
	scale = (m0+q)*(m1Sqr+m2Sqr+2*m1*q+qSqr),
	b0 = 1.0,
    b1 = -q*((2*m0*m1+m1Sqr+m2Sqr+(2*m0+4*m1)*q+3*qSqr)/scale),
    b2 = qSqr*(m0+2*m1+3*q)/scale,
	b3 = -qSqr*q/scale,

	B = (m0*(m1Sqr+m2Sqr)/scale)*(m0*(m1Sqr+m2Sqr)/scale),
	a1= -b1, a2 = -b2, a3 = -b3,
	scale1 = 1.0/((1.0+a1-a2+a3)*(1.0-a1-a2-a3)*(1.0+a2+(a1-a3)*a3)),

	M11 = scale1*(-a3*a1+1.0-a3*a3-a2),
	M12 = scale1*(a3+a1)*(a2+a3*a1),
	M13 = scale1*a3*(a1+a3*a2),
	M21 = scale1*(a1+a3*a2),
	M22 = -scale1*(a2-1.0)*(a2+a3*a1),
	M23 = -scale1*a3*(a3*a1+a3*a3+a2-1.0),
	M31 = scale1*(a3*a1+a2+a1*a1-a2*a2),
	M32 = scale1*(a1*a2+a3*a2*a2-a1*a3*a3-a3*a3*a3-a3*a2+a3),
	M33 = scale1*a3*(a1+a3*a2);

	int i = 0, j = 0, k = 0; bool flag = true;
	int size = width*height*depth;

	float wP1 = 0.f, wP2 = 0.f, wP3 = 0.f;
	float outF1 = 0.f, outF2 = 0.f, outF3 = 0.f, out = 0.f;
	float xP1 = 0.f;
	float xF1 = 0.f;

	switch (order) 
	{
		case 0: 
		{
			for(k=0;k<depth;k++)
			{
				for(j=0;j<height;j++)
				{
					for(i=0;i<width;i++)
					{
						if(flag)
						{
							wP1 = d_src[height*width*k+height*j+i]/sqrt(B);
							wP2 = wP1;
							wP3 = wP1;
							flag = false;
						}

						float xC = (float)d_src[height*width*k+height*j+i];
						float wC = (float)(xC - b1*wP1 - b2*wP2 - b3*wP3)/b0;
						d_dest[height*width*k+height*j+i] = (int)wC;
						wP3 = wP2; wP2 = wP1; wP1 = wC;

						if(i == width-1)
						{
							flag = true;

							float up = (float)d_src[height*width*k+height*j+i]/(1.0+b1+b2+b3);
							float vp = (float)up/(1.0+b1+b2+b3);
							out = (float)M11*(d_dest[height*width*k+height*j+i]-up) + M12*(d_dest[height*width*k+height*j+i-1]-up) + M13*(d_dest[height*width*k+height*j+i-2]-up)+vp;
							outF1 = (float)M21*(d_dest[height*width*k+height*j+i]-up) + M22*(d_dest[height*width*k+height*j+i-1]-up) + M23*(d_dest[height*width*k+height*j+i-2]-up)+vp;
							outF2 = (float)M31*(d_dest[height*width*k+height*j+i]-up) + M32*(d_dest[height*width*k+height*j+i-1]-up) + M33*(d_dest[height*width*k+height*j+i-2]-up)+vp;
							out *= B; outF1 *= B; outF2 *= B;
							outF3 = outF2; outF2 = outF1; outF1 = out;

							d_dest[height*width*k+height*j+i] = (short)out;
						}
					}
					for(i=width-1-1;i>=0;i--)
					{
						float wC = (float)d_dest[height*width*k+height*j+i];
						out = (float)(B*wC - b1*outF1 - b2*outF2 - b3*outF3)/b0;
						d_dest[height*width*k+height*j+i] = (int)out;
						outF3 = outF2; outF2 = outF1; outF1 = out;
					}
				}
			}
		} break;
		case 1:
		{
			float wC = 0.f;
			float xC = 0.f;

			for(k=0;k<depth;k++)
			{
				for(j=0;j<height;j++)
				{
					for(i=0;i<width;i++)
					{
						if(flag)
						{
							xP1 = (float)d_src[height*width*k+height*j+i];
							xF1 = (float)d_src[height*width*k+height*j+i+1];
							wP1 = (float)(xF1 - xP1)/2.0*(1.0+b1+b2+b3);
							wP2 = wP1;
							wP3 = wP1;

							wC = ((xF1- xP1)/2.0 - b1*wP1 - b2*wP1 - b3*wP1)/b0;
							d_dest[height*width*k+height*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;

							flag = false;

							continue;
						}

						if(i == width-1)
						{
							flag = true;

							xP1 = (float)d_src[height*width*k+height*j+i-1];
							xF1 = (float)d_src[height*width*k+height*j+i];
							wC = (float)((xF1- xP1)/2.0 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
							d_dest[height*width*k+height*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;

							float up = (d_src[height*width*k+height*j+i] - d_src[height*width*k+height*j+i-1])/2.0*(1.0+b1+b2+b3);////
							float vp = (float)up/(1.0+b1+b2+b3);
							out = (float)M11*(d_dest[height*width*k+height*j+i]-up) + M12*(d_dest[height*width*k+height*j+i-1]-up) + M13*(d_dest[height*width*k+height*j+i-2]-up)+vp;
							outF1 = (float)M21*(d_dest[height*width*k+height*j+i]-up) + M22*(d_dest[height*width*k+height*j+i-1]-up) + M23*(d_dest[height*width*k+height*j+i-2]-up)+vp;
							outF2 = (float)M31*(d_dest[height*width*k+height*j+i]-up) + M32*(d_dest[height*width*k+height*j+i-1]-up) + M33*(d_dest[height*width*k+height*j+i-2]-up)+vp;
							out *= B; outF1 *= B; outF2 *= B;
							outF3 = outF2; outF2 = outF1; outF1 = out;

							d_dest[height*width*k+height*j+i] = (int)out;
						}
						else
						{
							xP1 = (float)d_src[height*width*k+height*j+i-1];
							xF1 = (float)d_src[height*width*k+height*j+i+1];
							wC = (float)((xF1- xP1)/2.0 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
							d_dest[height*width*k+height*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;
						}
					}
					for(i=width-1-1;i>=0;i--)
					{
						wC = (float)d_dest[height*width*k+height*j+i];
						out = (float)(B*wC - b1*outF1 - b2*outF2 - b3*outF3)/b0;
						d_dest[height*width*k+height*j+i] = (int)out;
						outF3 = outF2; outF2 = outF1; outF1 = out;
					}
				}
			}
		} break;
		case 2:
		{
			float xC = 0.f;
			float wC = 0.f;
			for(k=0;k<depth;k++)
			{
				for(j=0;j<height;j++)
				{
					for(i=0;i<width;i++)
					{
						if(flag)
						{
							xC = (float)d_src[height*width*k+height*j+i];
							xP1 = (float)d_src[height*width*k+height*j+i];
							float xF1 = (float)d_src[height*width*k+height*j+i+1];
							wP1 = 0;
							wP2 = wP1;
							wP3 = wP1;

							float wC = (float)((xF1 - 2*xC + xP1) - b1*wP1 - b2*wP2 - b3*wP3)/b0;
							d_dest[height*width*k+height*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;

							flag = false;

							continue;
						}

						if(i == width-1)
						{
							flag = true;

							xC = (float)d_src[height*width*k+height*j+i];
							xP1 = (float)d_src[height*width*k+height*j+i-1];
							xF1 = (float)d_src[height*width*k+height*j+i];
							wC = (float)(xF1 - 2*xC + xP1 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
							d_dest[height*width*k+height*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;

							float up = 0;
							float vp = (float)up/(1.0+b1+b2+b3);
							out = (float)M11*(d_dest[height*width*k+height*j+i]-up) + M12*(d_dest[height*width*k+height*j+i-1]-up) + M13*(d_dest[height*width*k+height*j+i-2]-up)+vp;
							outF1 = (float)M21*(d_dest[height*width*k+height*j+i]-up) + M22*(d_dest[height*width*k+height*j+i-1]-up) + M23*(d_dest[height*width*k+height*j+i-2]-up)+vp;
							outF2 = (float)M31*(d_dest[height*width*k+height*j+i]-up) + M32*(d_dest[height*width*k+height*j+i-1]-up) + M33*(d_dest[height*width*k+height*j+i-2]-up)+vp;
							out *= B; outF1 *= B; outF2 *= B;
							outF3 = outF2; outF2 = outF1; outF1 = out;

							d_dest[height*width*k+height*j+i] = (int)out;
						}
						else
						{
							xC = (float)d_src[height*width*k+height*j+i];
							xP1 = (float)d_src[height*width*k+height*j+i-1];
							xF1 = (float)d_src[height*width*k+height*j+i+1];
							wC = (float)(xF1 - 2*xC + xP1 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
							d_dest[height*width*k+height*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;
						}
					}
					for(i=width-1-1;i>=0;i--)
					{
						float wC = (float)d_dest[height*width*k+height*j+i];
						out = (float)(B*wC - b1*outF1 - b2*outF2 - b3*outF3)/b0;
						d_dest[height*width*k+height*j+i] = (int)out;
						outF3 = outF2; outF2 = outF1; outF1 = out;
					}
				}
			}
		} break;
	}
	return 0;
}

template <class T1, class T2> int CFilter<T1, T2>::RecursiveGaussianY(int* d_src, int* d_dest, int order, float sigma)
{
	//const float nsigma = sigma < 0.5f ? 0.5f : sigma;
	float q = 0, m0 = 1.16680, m1 = 1.10783, m2 = 1.40586;
	float m1Sqr = m1*m1, m2Sqr = m2*m2;
	if(sigma < 1)
	{
		//std::cout<<"sigma value too small, result may be inaccurate! considering sigma as allowed minmun 0.5!"<<std::endl;
		q = 1.31564*(sqrt(1+0.490811*sigma*sigma)-1);
	}
	else
	{
		q = 1.31564*(sqrt(1+0.490811*sigma*sigma)-1);
	}
	
	const float
	qSqr = q*q,
	scale = (m0+q)*(m1Sqr+m2Sqr+2*m1*q+qSqr),
	b0 = 1.0,
    b1 = -q*((2*m0*m1+m1Sqr+m2Sqr+(2*m0+4*m1)*q+3*qSqr)/scale),
    b2 = qSqr*(m0+2*m1+3*q)/scale,
	b3 = -qSqr*q/scale,

	B = (m0*(m1Sqr+m2Sqr)/scale)*(m0*(m1Sqr+m2Sqr)/scale),
	a1= -b1, a2 = -b2, a3 = -b3,
	scale1 = 1.0/((1.0+a1-a2+a3)*(1.0-a1-a2-a3)*(1.0+a2+(a1-a3)*a3)),

	M11 = scale1*(-a3*a1+1.0-a3*a3-a2),
	M12 = scale1*(a3+a1)*(a2+a3*a1),
	M13 = scale1*a3*(a1+a3*a2),
	M21 = scale1*(a1+a3*a2),
	M22 = -scale1*(a2-1.0)*(a2+a3*a1),
	M23 = -scale1*a3*(a3*a1+a3*a3+a2-1.0),
	M31 = scale1*(a3*a1+a2+a1*a1-a2*a2),
	M32 = scale1*(a1*a2+a3*a2*a2-a1*a3*a3-a3*a3*a3-a3*a2+a3),
	M33 = scale1*a3*(a1+a3*a2);

	int i = 0, j = 0, k = 0; bool flag = true;
	int size = width*height*depth;

	float wP1 = 0.f, wP2 = 0.f, wP3 = 0.f;
	float outF1 = 0.f, outF2 = 0.f, outF3 = 0.f, out = 0.f;
	float xP1 = 0.f;
	float xF1 = 0.f;

	switch (order) 
	{
		case 0: 
		{
			for(k=0;k<depth;k++)
			{
				for(i=0;i<width;i++)
				{
					for(j=0;j<height;j++)
					{
						if(flag)
						{
							wP1 = d_src[height*width*k+height*j+i]/sqrt(B);
							wP2 = wP1;
							wP3 = wP1;
							flag = false;
						}

						float xC = (float)d_src[height*width*k+height*j+i];
						float wC = (float)(xC - b1*wP1 - b2*wP2 - b3*wP3)/b0;
						d_dest[height*width*k+height*j+i] = (int)wC;
						wP3 = wP2; wP2 = wP1; wP1 = wC;

						if(j == height-1)
						{
							flag = true;

							float up = (float)d_src[height*width*k+height*j+i]/(1.0+b1+b2+b3);
							float vp = (float)up/(1.0+b1+b2+b3);
							out = (float)M11*(d_dest[height*width*k+height*j+i]-up) + M12*(d_dest[height*width*k+height*(j-1)+i]-up) + M13*(d_dest[height*width*k+height*(j-2)+i]-up)+vp;
							outF1 = (float)M21*(d_dest[height*width*k+height*j+i]-up) + M22*(d_dest[height*width*k+height*(j-1)+i]-up) + M23*(d_dest[height*width*k+height*(j-2)+i]-up)+vp;
							outF2 = (float)M31*(d_dest[height*width*k+height*j+i]-up) + M32*(d_dest[height*width*k+height*(j-1)+i]-up) + M33*(d_dest[height*width*k+height*(j-2)+i]-up)+vp;
							out *= B; outF1 *= B; outF2 *= B;
							outF3 = outF2; outF2 = outF1; outF1 = out;

							d_dest[height*width*k+height*j+i] = (int)out;
						}
					}
					for(j=height-1-1;j>=0;j--)
					{
						float wC = (float)d_dest[height*width*k+height*j+i];
						out = (float)(B*wC - b1*outF1 - b2*outF2 - b3*outF3)/b0;
						d_dest[height*width*k+height*j+i] = (int)out;
						outF3 = outF2; outF2 = outF1; outF1 = out;
					}
				}
			}
		} break;
		case 1:
		{
			float wC = 0.f;
			float xC = 0.f;

			for(k=0;k<depth;k++)
			{
				for(i=0;i<width;i++)
				{
					for(j=0;j<height;j++)
					{
						if(flag)
						{
							xP1 = (float)d_src[height*width*k+height*j+i];
							xF1 = (float)d_src[height*width*k+height*(j+1)+i];
							wP1 = (float)(xF1 - xP1)/2.0*(1.0+b1+b2+b3);
							wP2 = wP1;
							wP3 = wP1;

							wC = ((xF1- xP1)/2.0 - b1*wP1 - b2*wP1 - b3*wP1)/b0;
							d_dest[height*width*k+height*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;

							flag = false;

							continue;
						}

						if(j == height-1)
						{
							flag = true;

							xP1 = (float)d_src[height*width*k+height*(j-1)+i];
							xF1 = (float)d_src[height*width*k+height*j+i];
							wC = (float)((xF1- xP1)/2.0 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
							d_dest[height*width*k+height*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;

							float up = (d_src[height*width*k+height*j+i] - d_src[height*width*k+height*(j-1)+i])/2.0*(1.0+b1+b2+b3);;
							float vp = (float)up/(1.0+b1+b2+b3);
							out = (float)M11*(d_dest[height*width*k+height*j+i]-up) + M12*(d_dest[height*width*k+height*(j-1)+i]-up) + M13*(d_dest[height*width*k+height*(j-2)+i]-up)+vp;
							outF1 = (float)M21*(d_dest[height*width*k+height*j+i]-up) + M22*(d_dest[height*width*k+height*(j-1)+i]-up) + M23*(d_dest[height*width*k+height*(j-2)+i]-up)+vp;
							outF2 = (float)M31*(d_dest[height*width*k+height*j+i]-up) + M32*(d_dest[height*width*k+height*(j-1)+i]-up) + M33*(d_dest[height*width*k+height*(j-2)+i]-up)+vp;
							out *= B; outF1 *= B; outF2 *= B;
							outF3 = outF2; outF2 = outF1; outF1 = out;

							d_dest[height*width*k+height*j+i] = (int)out;
						}
						else
						{
							xP1 = (float)d_src[height*width*k+height*(j-1)+i];
							xF1 = (float)d_src[height*width*k+height*(j+1)+i];
							wC = (float)((xF1- xP1)/2.0 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
							d_dest[height*width*k+height*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;
						}
					}
					for(j=height-1-1;j>=0;j--)
					{
						wC = (float)d_dest[height*width*k+height*j+i];
						out = (float)(B*wC - b1*outF1 - b2*outF2 - b3*outF3)/b0;
						d_dest[height*width*k+height*j+i] = (int)out;
						outF3 = outF2; outF2 = outF1; outF1 = out;
					}
				}
			}
		} break;
		case 2:
		{
			float xC = 0.f;
			float wC = 0.f;
			for(k=0;k<depth;k++)
			{
				for(i=0;i<width;i++)
				{
					for(j=0;j<height;j++)
					{
						if(flag)
						{
							xC = (float)d_src[height*width*k+height*j+i];
							xP1 = (float)d_src[height*width*k+height*j+i];
							float xF1 = (float)d_src[height*width*k+height*(j+1)+i];
							wP1 = 0;
							wP2 = wP1;
							wP3 = wP1;

							float wC = (float)((xF1 - 2*xC + xP1) - b1*wP1 - b2*wP2 - b3*wP3)/b0;
							d_dest[height*width*k+height*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;

							flag = false;

							continue;
						}

						if(j == height-1)
						{
							flag = true;

							xC = (float)d_src[height*width*k+height*j+i];
							xP1 = (float)d_src[height*width*k+height*(j-1)+i];
							xF1 = (float)d_src[height*width*k+height*j+i];
							wC = (float)(xF1 - 2*xC + xP1 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
							d_dest[height*width*k+height*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;

							float up = 0;
							float vp = (float)up/(1.0+b1+b2+b3);
							out = (float)M11*(d_dest[height*width*k+height*j+i]-up) + M12*(d_dest[height*width*k+height*(j-1)+i]-up) + M13*(d_dest[height*width*k+height*(j-2)+i]-up)+vp;
							outF1 = (float)M21*(d_dest[height*width*k+height*j+i]-up) + M22*(d_dest[height*width*k+height*(j-1)+i]-up) + M23*(d_dest[height*width*k+height*(j-2)+i]-up)+vp;
							outF2 = (float)M31*(d_dest[height*width*k+height*j+i]-up) + M32*(d_dest[height*width*k+height*(j-1)+i]-up) + M33*(d_dest[height*width*k+height*(j-2)+i]-up)+vp;
							out *= B; outF1 *= B; outF2 *= B;
							outF3 = outF2; outF2 = outF1; outF1 = out;

							d_dest[height*width*k+height*j+i] = (int)out;
						}
						else
						{
							xC = (float)d_src[height*width*k+height*j+i];
							xP1 = (float)d_src[height*width*k+height*(j-1)+i];
							xF1 = (float)d_src[height*width*k+height*(j+1)+i];
							wC = (float)(xF1 - 2*xC + xP1 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
							d_dest[height*width*k+height*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;
						}
					}
					for(j=height-1-1;j>=0;j--)
					{
						float wC = (float)d_dest[height*width*k+height*j+i];
						out = (float)(B*wC - b1*outF1 - b2*outF2 - b3*outF3)/b0;
						d_dest[height*width*k+height*j+i] = (int)out;
						outF3 = outF2; outF2 = outF1; outF1 = out;
					}
				}
			}
		} break;
	}
	return 0;
}

template <class T1, class T2> int CFilter<T1, T2>::RecursiveGaussianZ(int* d_src, int* d_dest, int order, float sigma)
{
	//const float nsigma = sigma < 0.5f ? 0.5f : sigma;
	float q = 0, m0 = 1.16680, m1 = 1.10783, m2 = 1.40586;
	float m1Sqr = m1*m1, m2Sqr = m2*m2;
	if(sigma < 1)
	{
		//std::cout<<"sigma value too small, result may be inaccurate! considering sigma as allowed minmun 0.5!"<<std::endl;
		q = 1.31564*(sqrt(1+0.490811*sigma*sigma)-1);
	}
	else
	{
		q = 1.31564*(sqrt(1+0.490811*sigma*sigma)-1);
	}
	
	const float
	qSqr = q*q,
	scale = (m0+q)*(m1Sqr+m2Sqr+2*m1*q+qSqr),
	b0 = 1.0,
    b1 = -q*((2*m0*m1+m1Sqr+m2Sqr+(2*m0+4*m1)*q+3*qSqr)/scale),
    b2 = qSqr*(m0+2*m1+3*q)/scale,
	b3 = -qSqr*q/scale,

	B = (m0*(m1Sqr+m2Sqr)/scale)*(m0*(m1Sqr+m2Sqr)/scale),
	a1= -b1, a2 = -b2, a3 = -b3,
	scale1 = 1.0/((1.0+a1-a2+a3)*(1.0-a1-a2-a3)*(1.0+a2+(a1-a3)*a3)),

	M11 = scale1*(-a3*a1+1.0-a3*a3-a2),
	M12 = scale1*(a3+a1)*(a2+a3*a1),
	M13 = scale1*a3*(a1+a3*a2),
	M21 = scale1*(a1+a3*a2),
	M22 = -scale1*(a2-1.0)*(a2+a3*a1),
	M23 = -scale1*a3*(a3*a1+a3*a3+a2-1.0),
	M31 = scale1*(a3*a1+a2+a1*a1-a2*a2),
	M32 = scale1*(a1*a2+a3*a2*a2-a1*a3*a3-a3*a3*a3-a3*a2+a3),
	M33 = scale1*a3*(a1+a3*a2);

	int i = 0, j = 0, k = 0; bool flag = true;
	int size = width*height*depth;

	float wP1 = 0.f, wP2 = 0.f, wP3 = 0.f;
	float outF1 = 0.f, outF2 = 0.f, outF3 = 0.f, out = 0.f;
	float xP1 = 0.f;
	float xF1 = 0.f;

	switch (order) 
	{
		case 0: 
		{
			for(i=0;i<width;i++)
			{
				for(j=0;j<height;j++)
				{
					for(k=0;k<depth;k++)
					{
						if(flag)
						{
							wP1 = d_src[height*width*k+height*j+i]/sqrt(B);
							wP2 = wP1;
							wP3 = wP1;
							flag = false;
						}

						float xC = (float)d_src[height*width*k+height*j+i];
						float wC = (float)(xC - b1*wP1 - b2*wP2 - b3*wP3)/b0;
						d_dest[height*width*k+height*j+i] = (int)wC;
						wP3 = wP2; wP2 = wP1; wP1 = wC;

						if(k == depth-1)
						{
							flag = true;

							float up = (float)d_src[height*width*k+height*j+i]/(1.0+b1+b2+b3);
							float vp = (float)up/(1.0+b1+b2+b3);
							out = (float)M11*(d_dest[height*width*k+height*j+i]-up) + M12*(d_dest[height*width*(k-1)+height*j+i]-up) + M13*(d_dest[height*width*(k-2)+height*j+i]-up)+vp;
							outF1 = (float)M21*(d_dest[height*width*k+height*j+i]-up) + M22*(d_dest[height*width*(k-1)+height*j+i]-up) + M23*(d_dest[height*width*(k-2)+height*j+i]-up)+vp;
							outF2 = (float)M31*(d_dest[height*width*k+height*j+i]-up) + M32*(d_dest[height*width*(k-1)+height*j+i]-up) + M33*(d_dest[height*width*(k-2)+height*j+i]-up)+vp;
							out *= B; outF1 *= B; outF2 *= B;
							outF3 = outF2; outF2 = outF1; outF1 = out;

							d_dest[height*width*k+height*j+i] = (int)out;
						}
					}
					for(k=depth-1-1;k>=0;k--)
					{
						float wC = (float)d_dest[height*width*k+height*j+i];
						out = (float)(B*wC - b1*outF1 - b2*outF2 - b3*outF3)/b0;
						d_dest[height*width*k+height*j+i] = (int)out;
						outF3 = outF2; outF2 = outF1; outF1 = out;
					}
				}
			}
		} break;
		case 1:
		{
			float wC = 0.f;
			float xC = 0.f;

			for(i=0;i<width;i++)
			{
				for(j=0;j<height;j++)
				{
					for(k=0;k<depth;k++)
					{
						if(flag)
						{
							xP1 = (float)d_src[height*width*k+height*j+i];
							xF1 = (float)d_src[height*width*(k+1)+height*j+i];
							wP1 = (float)(xF1 - xP1)/2.0*(1.0+b1+b2+b3);
							wP2 = wP1;
							wP3 = wP1;

							wC = ((xF1- xP1)/2.0 - b1*wP1 - b2*wP1 - b3*wP1)/b0;
							d_dest[height*width*k+height*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;

							flag = false;

							continue;
						}

						if(k == depth-1)
						{
							flag = true;

							xP1 = (float)d_src[height*width*(k-1)+height*j+i];
							xF1 = (float)d_src[height*width*k+height*j+i];
							wC = (float)((xF1- xP1)/2.0 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
							d_dest[height*width*k+height*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;

							float up = (d_src[height*width*k+height*j+i] - d_src[height*width*(k-1)+height*j+i])/2.0*(1.0+b1+b2+b3);;
							float vp = (float)up/(1.0+b1+b2+b3);
							out = (float)M11*(d_dest[height*width*k+height*j+i]-up) + M12*(d_dest[height*width*(k-1)+height*j+i]-up) + M13*(d_dest[height*width*(k-2)+height*j+i]-up)+vp;
							outF1 = (float)M21*(d_dest[height*width*k+height*j+i]-up) + M22*(d_dest[height*width*(k-1)+height*j+i]-up) + M23*(d_dest[height*width*(k-2)+height*j+i]-up)+vp;
							outF2 = (float)M31*(d_dest[height*width*k+height*j+i]-up) + M32*(d_dest[height*width*(k-1)+height*j+i]-up) + M33*(d_dest[height*width*(k-2)+height*j+i]-up)+vp;
							out *= B; outF1 *= B; outF2 *= B;
							outF3 = outF2; outF2 = outF1; outF1 = out;

							d_dest[height*width*k+height*j+i] = (int)out;
						}
						else
						{
							xP1 = (float)d_src[height*width*(k-1)+height*j+i];
							xF1 = (float)d_src[height*width*(k+1)+height*j+i];
							wC = (float)((xF1- xP1)/2.0 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
							d_dest[height*width*k+height*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;
						}
					}
					for(k=depth-1-1;k>=0;k--)
					{
						wC = (float)d_dest[height*width*k+height*j+i];
						out = (float)(B*wC - b1*outF1 - b2*outF2 - b3*outF3)/b0;
						d_dest[height*width*k+height*j+i] = (int)out;
						outF3 = outF2; outF2 = outF1; outF1 = out;
					}
				}
			}
		} break;
		case 2:
		{
			float xC = 0.f;
			float wC = 0.f;
			for(i=0;i<width;i++)
			{
				for(j=0;j<height;j++)
				{
					for(k=0;k<depth;k++)
					{
						if(flag)
						{
							xC = (float)d_src[height*width*k+height*j+i];
							xP1 = (float)d_src[height*width*k+height*j+i];
							float xF1 = (float)d_src[height*width*(k+1)+height*j+i];
							wP1 = 0;
							wP2 = wP1;
							wP3 = wP1;

							float wC = (float)((xF1 - 2*xC + xP1) - b1*wP1 - b2*wP2 - b3*wP3)/b0;
							d_dest[height*width*k+height*j+i] = (short)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;

							flag = false;

							continue;
						}

						if(k == depth-1)
						{
							flag = true;

							xC = (float)d_src[height*width*k+height*j+i];
							xP1 = (float)d_src[height*width*(k-1)+height*j+i];
							xF1 = (float)d_src[height*width*k+height*j+i];
							wC = (float)(xF1 - 2*xC + xP1 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
							d_dest[height*width*k+height*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;

							float up = 0;
							float vp = (float)up/(1.0+b1+b2+b3);
							out = (float)M11*(d_dest[height*width*k+height*j+i]-up) + M12*(d_dest[height*width*(k-1)+height*j+i]-up) + M13*(d_dest[height*width*(k-2)+height*j+i]-up)+vp;
							outF1 = (float)M21*(d_dest[height*width*k+height*j+i]-up) + M22*(d_dest[height*width*(k-1)+height*j+i]-up) + M23*(d_dest[height*width*(k-2)+height*j+i]-up)+vp;
							outF2 = (float)M31*(d_dest[height*width*k+height*j+i]-up) + M32*(d_dest[height*width*(k-1)+height*j+i]-up) + M33*(d_dest[height*width*(k-2)+height*j+i]-up)+vp;
							out *= B; outF1 *= B; outF2 *= B;
							outF3 = outF2; outF2 = outF1; outF1 = out;

							d_dest[height*width*k+height*j+i] = (int)out;
						}
						else
						{
							xC = (float)d_src[height*width*k+height*j+i];
							xP1 = (float)d_src[height*width*(k-1)+height*j+i];
							xF1 = (float)d_src[height*width*(k+1)+height*j+i];
							wC = (float)(xF1 - 2*xC + xP1 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
							d_dest[height*width*k+height*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;
						}
					}
					for(k=depth-1-1;k>=0;k--)
					{
						float wC = (float)d_dest[height*width*k+height*j+i];
						out = (float)(B*wC - b1*outF1 - b2*outF2 - b3*outF3)/b0;
						d_dest[height*width*k+height*j+i] = (int)out;
						outF3 = outF2; outF2 = outF1; outF1 = out;
					}
				}
			}
		} break;
	}
	return 0;
}

template <class T1, class T2> int CFilter<T1, T2>::iDivUp(int a, int b)
{
	return (a % b != 0) ? (a / b + 1) : (a / b);
}

template <class T1, class T2> void CFilter<T1, T2>::CudaRecursiveGaussian(int xDirection, int yDirection, int zDirection, T1* d_src, T2* d_dest,/* int order, */float sigma)
{
	//int controle = 0;
	int size = width*height*depth;
	SetInput(d_src);//T* d_src
	PassToDevice(h_input, size);//int* h_input
	if(xDirection >= 0)
	{
		if(flag)
			CudaRecursiveGaussianX(GetInputImage(), GetOutputImage(), xDirection, sigma);
		else
			CudaRecursiveGaussianX(GetOutputImage(), GetInputImage(), xDirection, sigma);
		flag = !flag;
	}
	if(yDirection >= 0)
	{
		if(flag)
			CudaRecursiveGaussianY(GetInputImage(), GetOutputImage(), yDirection, sigma);
		else
			CudaRecursiveGaussianY(GetOutputImage(), GetInputImage(), yDirection, sigma);
		flag = !flag;
	}
	if(zDirection >= 0)
	{
		if(flag)
			CudaRecursiveGaussianZ(GetInputImage(), GetOutputImage(), zDirection, sigma);
		else
			CudaRecursiveGaussianZ(GetOutputImage(), GetInputImage(), zDirection, sigma);
		flag = !flag;
	}
	ToHost(d_dest);
	delete [] h_output;
}

template <class T1, class T2> int CFilter<T1, T2>::CudaRecursiveGaussianX(int* d_src, int* d_dest, int order, float sigma)
{
	int size = width*height*depth;

	int nthread = 2;
	int n = iDivUp(height, nthread);

    //const float nsigma = sigma < 0.5f ? 0.5f : sigma;
	float q = 0, m0 = 1.16680, m1 = 1.10783, m2 = 1.40586;
	float m1Sqr = m1*m1, m2Sqr = m2*m2;
	if(sigma < 1)
	{
		//std::cout<<"sigma value too small, result may be inaccurate! considering sigma as allowed minmun 0.5!"<<std::endl;
		q = 1.31564*(sqrt(1+0.490811*sigma*sigma)-1);
	}
	else
	{
		q = 1.31564*(sqrt(1+0.490811*sigma*sigma)-1);
	}
	
	const float
	qSqr = q*q,
	scale = (m0+q)*(m1Sqr+m2Sqr+2*m1*q+qSqr),
	b0 = 1.0,
    b1 = -q*((2*m0*m1+m1Sqr+m2Sqr+(2*m0+4*m1)*q+3*qSqr)/scale),
    b2 = qSqr*(m0+2*m1+3*q)/scale,
	b3 = -qSqr*q/scale,

	B = (m0*(m1Sqr+m2Sqr)/scale)*(m0*(m1Sqr+m2Sqr)/scale),
	a1= -b1, a2 = -b2, a3 = -b3,
	scale1 = 1.0/((1.0+a1-a2+a3)*(1.0-a1-a2-a3)*(1.0+a2+(a1-a3)*a3)),

	M11 = scale1*(-a3*a1+1.0-a3*a3-a2),
	M12 = scale1*(a3+a1)*(a2+a3*a1),
	M13 = scale1*a3*(a1+a3*a2),
	M21 = scale1*(a1+a3*a2),
	M22 = -scale1*(a2-1.0)*(a2+a3*a1),
	M23 = -scale1*a3*(a3*a1+a3*a3+a2-1.0),
	M31 = scale1*(a3*a1+a2+a1*a1-a2*a2),
	M32 = scale1*(a1*a2+a3*a2*a2-a1*a3*a3-a3*a3*a3-a3*a2+a3),
	M33 = scale1*a3*(a1+a3*a2);

	Call_d_recursiveGaussianX(d_src, d_dest, width, height, depth, b0, b1, b2, b3, B, order, n, nthread, M11, M12, M13, M21, M22, M23, M31, M32, M33);
	return 0;
}

template <class T1, class T2> int CFilter<T1, T2>::CudaRecursiveGaussianY(int* d_src, int* d_dest, int order, float sigma)
{
	int size = width*height*depth;

	int nthread = 256;
	int n = iDivUp(width, nthread);

    //const float nsigma = sigma < 0.5f ? 0.5f : sigma;
	float q = 0, m0 = 1.16680, m1 = 1.10783, m2 = 1.40586;
	float m1Sqr = m1*m1, m2Sqr = m2*m2;
	if(sigma < 1)
	{
		//std::cout<<"sigma value too small, result may be inaccurate! considering sigma as allowed minmun 0.5!"<<std::endl;
		q = 1.31564*(sqrt(1+0.490811*sigma*sigma)-1);
	}
	else
	{
		q = 1.31564*(sqrt(1+0.490811*sigma*sigma)-1);
	}
	
	const float
	qSqr = q*q,
	scale = (m0+q)*(m1Sqr+m2Sqr+2*m1*q+qSqr),
	b0 = 1.0,
    b1 = -q*((2*m0*m1+m1Sqr+m2Sqr+(2*m0+4*m1)*q+3*qSqr)/scale),
    b2 = qSqr*(m0+2*m1+3*q)/scale,
	b3 = -qSqr*q/scale,

	B = (m0*(m1Sqr+m2Sqr)/scale)*(m0*(m1Sqr+m2Sqr)/scale),
	a1= -b1, a2 = -b2, a3 = -b3,
	scale1 = 1.0/((1.0+a1-a2+a3)*(1.0-a1-a2-a3)*(1.0+a2+(a1-a3)*a3)),

	M11 = scale1*(-a3*a1+1.0-a3*a3-a2),
	M12 = scale1*(a3+a1)*(a2+a3*a1),
	M13 = scale1*a3*(a1+a3*a2),
	M21 = scale1*(a1+a3*a2),
	M22 = -scale1*(a2-1.0)*(a2+a3*a1),
	M23 = -scale1*a3*(a3*a1+a3*a3+a2-1.0),
	M31 = scale1*(a3*a1+a2+a1*a1-a2*a2),
	M32 = scale1*(a1*a2+a3*a2*a2-a1*a3*a3-a3*a3*a3-a3*a2+a3),
	M33 = scale1*a3*(a1+a3*a2);

	Call_d_recursiveGaussianY(d_src, d_dest, width, height, depth, b0, b1, b2, b3, B, order, n, nthread, M11, M12, M13, M21, M22, M23, M31, M32, M33);
	return 0;
}

template <class T1, class T2> int CFilter<T1, T2>::CudaRecursiveGaussianZ(int* d_src, int* d_dest, int order, float sigma)
{
	int size = width*height*depth;

	int nthread = 256;
	int n = iDivUp(width, nthread);

    //const float nsigma = sigma < 0.5f ? 0.5f : sigma;
	float q = 0, m0 = 1.16680, m1 = 1.10783, m2 = 1.40586;
	float m1Sqr = m1*m1, m2Sqr = m2*m2;
	if(sigma < 1)
	{
		//std::cout<<"sigma value too small, result may be inaccurate! considering sigma as allowed minmun 0.5!"<<std::endl;
		q = 1.31564*(sqrt(1+0.490811*sigma*sigma)-1);
	}
	else
	{
		q = 1.31564*(sqrt(1+0.490811*sigma*sigma)-1);
	}
	
	const float
	qSqr = q*q,
	scale = (m0+q)*(m1Sqr+m2Sqr+2*m1*q+qSqr),
	b0 = 1.0,
    b1 = -q*((2*m0*m1+m1Sqr+m2Sqr+(2*m0+4*m1)*q+3*qSqr)/scale),
    b2 = qSqr*(m0+2*m1+3*q)/scale,
	b3 = -qSqr*q/scale,

	B = (m0*(m1Sqr+m2Sqr)/scale)*(m0*(m1Sqr+m2Sqr)/scale),
	a1= -b1, a2 = -b2, a3 = -b3,
	scale1 = 1.0/((1.0+a1-a2+a3)*(1.0-a1-a2-a3)*(1.0+a2+(a1-a3)*a3)),

	M11 = scale1*(-a3*a1+1.0-a3*a3-a2),
	M12 = scale1*(a3+a1)*(a2+a3*a1),
	M13 = scale1*a3*(a1+a3*a2),
	M21 = scale1*(a1+a3*a2),
	M22 = -scale1*(a2-1.0)*(a2+a3*a1),
	M23 = -scale1*a3*(a3*a1+a3*a3+a2-1.0),
	M31 = scale1*(a3*a1+a2+a1*a1-a2*a2),
	M32 = scale1*(a1*a2+a3*a2*a2-a1*a3*a3-a3*a3*a3-a3*a2+a3),
	M33 = scale1*a3*(a1+a3*a2);

	Call_d_recursiveGaussianZ(d_src, d_dest, width, height, depth, b0, b1, b2, b3, B, order, n, nthread, M11, M12, M13, M21, M22, M23, M31, M32, M33);
	return 0;
}

template <class T1, class T2> void CFilter<T1, T2>::ToHost(T2 *d_dest)
{
	int s = 1;
	if(!flag)
		s = 0;
	int size = width*height*depth;
	PassToHost(h_output, size, s);
	int i = 0;
	for(i=0;i<width*height*depth;i++)
	{
		d_dest[i] = h_output[i];
	}
}

template <class T1, class T2> void CFilter<T1, T2>::SobelFilter(T1* in, T2* out)
{
	SetInput(in);
	int i = 0, j = 0, k = 0;
	int ul,um,ur,
		  ml,mm,mr,
		  ll,lm,lr;
	int Horz, Vert, Sum;

	for(k=0;k<depth;k++)
	{
		for(j=0;j<height;j++)
		{
			for(i=0;i<width;i++)
			{
				if(i == 0 && j != 0)
				{
					ul = 0;
					um = h_input[height*width*k+height*(j-1)+i];
					ur = h_input[height*width*k+height*(j-1)+i+1];
					ml = 0;
					mm = h_input[height*width*k+height*j+i];
					mr = h_input[height*width*k+height*j+i+1];
					ll = 0;
					lm = h_input[height*width*k+height*(j+1)+i];
					lr = h_input[height*width*k+height*(j+1)+i+1];
				}
				else if(i != 0 && j == 0)
				{
					ul = 0;
					um = 0;
					ur = 0;
					ml = h_input[height*width*k+height*j+i-1];
					mm = h_input[height*width*k+height*j+i];
					mr = h_input[height*width*k+height*j+i+1];
					ll = h_input[height*width*k+height*(j+1)+i-1];
					lm = h_input[height*width*k+height*(j+1)+i];
					lr = h_input[height*width*k+height*(j+1)+i+1];
				}
				else if(i == 0 && j == 0)
				{
					ul = 0;
					um = 0;
					ur = 0;
					ml = 0;
					mm = h_input[height*width*k+height*j+i];
					mr = h_input[height*width*k+height*j+i+1];
					ll = 0;
					lm = h_input[height*width*k+height*(j+1)+i];
					lr = h_input[height*width*k+height*(j+1)+i+1];
				}
				else
				{
					ul = h_input[height*width*k+height*(j-1)+i-1];
					um = h_input[height*width*k+height*(j-1)+i];
					ur = h_input[height*width*k+height*(j-1)+i+1];
					ml = h_input[height*width*k+height*j+i-1];
					mm = h_input[height*width*k+height*j+i];
					mr = h_input[height*width*k+height*j+i+1];
					ll = h_input[height*width*k+height*(j+1)+i-1];
					lm = h_input[height*width*k+height*(j+1)+i];
					lr = h_input[height*width*k+height*(j+1)+i+1];
				}
				Horz = ur + 2*mr + lr - ul - 2*ml - ll;
				Vert = ul + 2*um + ur - ll - 2*lm - lr;
				Sum = (int)(abs(Horz)+abs(Vert));

				//if ( Sum < 0 )
				//	h_output[height*width*k+height*j+i] = 0;
				//else if ( Sum > 0xffff ) 
				//	h_output[height*width*k+height*j+i] = 0xffff;
				h_output[height*width*k+height*j+i] = Sum;
			}
		}
	}
	int n = 0;
	for(n=0;n<width*height*depth;n++)
	{
		out[n] = h_output[n];
	}
	delete [] h_output;
}

template <class T1, class T2> void CFilter<T1, T2>::CudaSobelFilter(T1* d_src, T2* d_dest)
{
	int size = width*height*depth;
	SetInput(d_src);
	PassToDevice(h_input, size);

	int nthread = 256;
	int n = iDivUp(width, nthread);
	if(flag)
		Call_d_SobelFilter(GetInputImage(), GetOutputImage(), width, height, depth, n, nthread);
	else
		Call_d_SobelFilter(GetOutputImage(), GetInputImage(), width, height, depth, n, nthread);
	flag = !flag;
	ToHost(d_dest);
	delete [] h_output;
}

template <class T1, class T2> void CFilter<T1, T2>::CudaSobelNew(T1* d_src, T2* d_dest)
{
	//int size = width*height*depth;
	int* data = NULL;
	SetInput(d_src);
	setupTexture(h_input, width, height, depth);
	checkCudaErrors(cudaMalloc((void**)&data, width*height*depth*sizeof(int)));
	SobelNew(data, width, width, height, depth, 1.0);
	checkCudaErrors(cudaMemcpy(h_output, data, width*height*depth*sizeof(int), cudaMemcpyDeviceToHost));
	deleteTexture();
	for(int i=0;i<width*height*depth;i++)
	{
		d_dest[i] = h_output[i];
	}
	delete [] h_output;
}

template <class T1, class T2> void CFilter<T1, T2>::CudaClear()
{
	ClearDevice();
	ResetFlag();
}

template <class T1, class T2> void CFilter<T1, T2>::ResetFlag()
{
	flag = true;
}

//void CFilter::ToDisplay(int c, int w)
//{
//}
template <class T1, class T2> void CFilter<T1, T2>::RecursiveGaussianFilter(int xDirection, int yDirection, int zDirection, T1* d_src, T2* d_dest, /*int order, */float sigma)
{
	bool cudaFlag = false;
	int deviceCount;
	cudaGetDeviceCount(&deviceCount);
	int dev;
	for(dev = 0; dev < deviceCount; ++dev)
	{
		cudaDeviceProp deviceProp;
		cudaGetDeviceProperties(&deviceProp, dev);
		if(dev == 0)
		{
			if (deviceProp.major == 9999 && deviceProp.minor == 9999)
				cudaFlag = false;
			else
				cudaFlag = true;
		}
	}
	if(!cudaFlag)
		RecursiveGaussian(xDirection, yDirection, zDirection, d_src, d_dest, sigma);
	else
	{
		CudaRecursiveGaussian(xDirection, yDirection, zDirection, d_src, d_dest, sigma);
		CudaClear();
	}
}

template <class T1, class T2> void CFilter<T1, T2>::SobelEdgeDetect(T1* d_src, T2* d_dest)
{
	bool cudaFlag = false;
	int deviceCount;
	cudaGetDeviceCount(&deviceCount);
	int dev;
	for(dev = 0; dev < deviceCount; ++dev)
	{
		cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, dev);
		if(dev == 0)
		{
			if (deviceProp.major == 9999 && deviceProp.minor == 9999)
				cudaFlag = false;
			else
				cudaFlag = true;
		}
	}
	if(!cudaFlag)
		SobelFilter(d_src, d_dest);
	else
	{
		CudaSobelFilter(d_src, d_dest);
		CudaClear();
	}
}

#endif
