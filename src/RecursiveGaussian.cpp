#include "RecursiveGaussian.h"


RecursiveGaussian::RecursiveGaussian()
    : m_eDataType(DataType::Undefined)
    , m_fInputData(NULL)
    , m_fPixelData(NULL)
    , m_Max(0)
    , m_Min(10000000)
    , m_bFlip(false)
{
}


RecursiveGaussian::~RecursiveGaussian()
{
    if(m_fInputData != NULL)
    {
        delete [] m_fInputData;
        m_fInputData = NULL;
    }

    if(m_fPixelData != NULL)
    {
        delete [] m_fPixelData;
        m_fPixelData = NULL;
    }
}


void RecursiveGaussian::SetInput(BaseImage* input)
{
    m_pInput = input;
}


void RecursiveGaussian::WriteToImage(BaseImage* output)
{
    int i = 0, j = 0;
    m_bFlip = !m_bFlip;

    if(m_bFlip)
    {
        for(i=0;i<m_lWidth*m_lHeight*m_lDepth;i++)
        {
        
            if(m_fInputData[i]>m_Max)
            {
                m_Max = m_fInputData[i];
            }
            if(m_fInputData[i]<m_Min)
            {
                m_Min = m_fInputData[i];
            }
        }
    }
    else
    {
        for(i=0;i<m_lWidth*m_lHeight*m_lDepth;i++)
        {
            if(m_fPixelData[i]>m_Max)
            {
                m_Max = m_fPixelData[i];
            }
            if(m_fPixelData[i]<m_Min)
            {
                m_Min = m_fPixelData[i];
            }
        }
    }

    if(m_bFlip)
    {
        for(i=0;i<m_lWidth*m_lHeight*m_lDepth;i++)
        {
            output->GetUnsignedCharPixelData()[i] = 255*(m_fInputData[i]/(m_Max-m_Min));
        }
    }
    else
    {
        for(i=0;i<m_lWidth*m_lHeight*m_lDepth;i++)
        {
            output->GetUnsignedCharPixelData()[i] = 255*(m_fPixelData[i]/(m_Max-m_Min));
        }
    }
    output->Update();
}


void RecursiveGaussian::Filter(int order, float sigma, int dir)
{
    int i = 0;
    m_lWidth = m_pInput->GetImageWidth();
    m_lHeight = m_pInput->GetImageHeight();
    m_lDepth = m_pInput->GetImageDepth();

    m_fInputData = new float[m_lWidth*m_lHeight*m_lDepth];

    if(m_pInput->GetUnsignedCharPixelData() != NULL)
    {
        m_eDataType = DataType::UChar;
        for(i=0;i<m_lWidth*m_lHeight*m_lDepth;i++)
        {
            m_fInputData[i] = m_pInput->GetUnsignedCharPixelData()[i];
        }
    }
    else if(m_pInput->GetShortPixelData() != NULL)
    {
        m_eDataType = DataType::Short;
        for(i=0;i<m_lWidth*m_lHeight*m_lDepth;i++)
        {
            m_fInputData[i] = m_pInput->GetShortPixelData()[i];
        }
    }
    else
    {
        std::cerr<<"Input Image is NULL! "<<std::endl;
        return;
    }

    m_fPixelData = new float[m_lWidth*m_lHeight*m_lDepth];

    if( dir & (Direction::XDirection) )
    {
        RecursiveGaussianX(m_fInputData, m_fPixelData, order, sigma);
        m_bFlip = !m_bFlip;
    }
    if( dir & (Direction::YDirection) )
    {
        if(m_bFlip)
        {
            RecursiveGaussianY(m_fPixelData, m_fInputData, order, sigma);
        }
        else
        {
            RecursiveGaussianY(m_fInputData, m_fPixelData, order, sigma);
        }
        m_bFlip = !m_bFlip;
    }
    if( dir & (Direction::ZDirection) )
    {
        if(m_bFlip)
        {
            RecursiveGaussianZ(m_fPixelData, m_fInputData, order, sigma);
        }
        else
        {
            RecursiveGaussianZ(m_fInputData, m_fPixelData, order, sigma);
        }
        m_bFlip = !m_bFlip;
    }
}


int RecursiveGaussian::RecursiveGaussianX(float* pData, float * outData, int order, float sigma)
{
    const float nsigma = sigma < 0.5f ? 0.5f : sigma;
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
	int size = m_lWidth*m_lHeight*m_lDepth;

	float wP1 = 0.f, wP2 = 0.f, wP3 = 0.f;
	float outF1 = 0.f, outF2 = 0.f, outF3 = 0.f, out = 0.f;
	float xP1 = 0.f;
	float xF1 = 0.f;

	switch (order) 
	{
		case 0: 
		{
			for(k=0;k<m_lDepth;k++)
			{
				for(j=0;j<m_lHeight;j++)
				{
					for(i=0;i<m_lWidth;i++)
					{
						if(flag)
						{
							wP1 = pData[m_lHeight*m_lWidth*k+m_lHeight*j+i]/sqrt(B);
							wP2 = wP1;
							wP3 = wP1;
							flag = false;
						}

						float xC = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
						float wC = (float)(xC - b1*wP1 - b2*wP2 - b3*wP3)/b0;
						outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)wC;
						wP3 = wP2; wP2 = wP1; wP1 = wC;

						if(i == m_lWidth-1)
						{
							flag = true;

							float up = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i]/(1.0+b1+b2+b3);
							float vp = (float)up/(1.0+b1+b2+b3);
							out = (float)M11*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i]-up) + 
                                         M12*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i-1]-up) + 
                                         M13*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i-2]-up)+vp;

							outF1 = (float)M21*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i]-up) + 
                                           M22*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i-1]-up) + 
                                           M23*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i-2]-up)+vp;

							outF2 = (float)M31*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i]-up) + 
                                           M32*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i-1]-up) + 
                                           M33*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i-2]-up)+vp;

							out *= B; outF1 *= B; outF2 *= B;
							outF3 = outF2; outF2 = outF1; outF1 = out;

							outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (short)out;
						}
					}
					for(i=m_lWidth-1-1;i>=0;i--)
					{
						float wC = (float)outData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
						out = (float)(B*wC - b1*outF1 - b2*outF2 - b3*outF3)/b0;
						outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)out;
						outF3 = outF2; outF2 = outF1; outF1 = out;
					}
				}
			}
		} break;
		case 1:
		{
			float wC = 0.f;
			float xC = 0.f;

			for(k=0;k<m_lDepth;k++)
			{
				for(j=0;j<m_lHeight;j++)
				{
					for(i=0;i<m_lWidth;i++)
					{
						if(flag)
						{
							xP1 = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
							xF1 = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i+1];
							wP1 = (float)(xF1 - xP1)/2.0*(1.0+b1+b2+b3);
							wP2 = wP1;
							wP3 = wP1;

							wC = ((xF1- xP1)/2.0 - b1*wP1 - b2*wP1 - b3*wP1)/b0;
							outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;

							flag = false;

							continue;
						}

						if(i == m_lWidth-1)
						{
							flag = true;

							xP1 = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i-1];
							xF1 = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
							wC = (float)((xF1- xP1)/2.0 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
							outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;

							float up = (pData[m_lHeight*m_lWidth*k+m_lHeight*j+i] - pData[m_lHeight*m_lWidth*k+m_lHeight*j+i-1])/2.0*(1.0+b1+b2+b3);////
							float vp = (float)up/(1.0+b1+b2+b3);
							out = (float)M11*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i]-up) + 
                                         M12*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i-1]-up) + 
                                         M13*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i-2]-up)+vp;

							outF1 = (float)M21*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i]-up) + 
                                           M22*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i-1]-up) + 
                                           M23*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i-2]-up)+vp;

							outF2 = (float)M31*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i]-up) + 
                                           M32*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i-1]-up) + 
                                           M33*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i-2]-up)+vp;

							out *= B; outF1 *= B; outF2 *= B;
							outF3 = outF2; outF2 = outF1; outF1 = out;

							outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)out;
						}
						else
						{
							xP1 = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i-1];
							xF1 = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i+1];
							wC = (float)((xF1- xP1)/2.0 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
							outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;
						}
					}
					for(i=m_lWidth-1-1;i>=0;i--)
					{
						wC = (float)outData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
						out = (float)(B*wC - b1*outF1 - b2*outF2 - b3*outF3)/b0;
						outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)out;
						outF3 = outF2; outF2 = outF1; outF1 = out;
					}
				}
			}
		} break;
		case 2:
		{
			float xC = 0.f;
			float wC = 0.f;
			for(k=0;k<m_lDepth;k++)
			{
				for(j=0;j<m_lHeight;j++)
				{
					for(i=0;i<m_lWidth;i++)
					{
						if(flag)
						{
							xC = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
							xP1 = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
							float xF1 = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i+1];
							wP1 = 0;
							wP2 = wP1;
							wP3 = wP1;

							float wC = (float)((xF1 - 2*xC + xP1) - b1*wP1 - b2*wP2 - b3*wP3)/b0;
							outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;

							flag = false;

							continue;
						}

						if(i == m_lWidth-1)
						{
							flag = true;

							xC = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
							xP1 = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i-1];
							xF1 = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
							wC = (float)(xF1 - 2*xC + xP1 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
							outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;

							float up = 0;
							float vp = (float)up/(1.0+b1+b2+b3);
							out = (float)M11*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i]-up) + 
                                         M12*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i-1]-up) + 
                                         M13*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i-2]-up)+vp;

							outF1 = (float)M21*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i]-up) + 
                                           M22*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i-1]-up) + 
                                           M23*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i-2]-up)+vp;

							outF2 = (float)M31*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i]-up) + 
                                           M32*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i-1]-up) + 
                                           M33*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i-2]-up)+vp;

							out *= B; outF1 *= B; outF2 *= B;
							outF3 = outF2; outF2 = outF1; outF1 = out;

							outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)out;
						}
						else
						{
							xC = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
							xP1 = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i-1];
							xF1 = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i+1];
							wC = (float)(xF1 - 2*xC + xP1 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
							outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;
						}
					}
					for(i=m_lWidth-1-1;i>=0;i--)
					{
						float wC = (float)outData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
						out = (float)(B*wC - b1*outF1 - b2*outF2 - b3*outF3)/b0;
						outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)out;
						outF3 = outF2; outF2 = outF1; outF1 = out;
					}
				}
			}
		} break;
	}
	return 0;
}


int RecursiveGaussian::RecursiveGaussianY(float* pData, float * outData, int order, float sigma)
{
    const float nsigma = sigma < 0.5f ? 0.5f : sigma;
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
	int size = m_lWidth*m_lHeight*m_lDepth;

	float wP1 = 0.f, wP2 = 0.f, wP3 = 0.f;
	float outF1 = 0.f, outF2 = 0.f, outF3 = 0.f, out = 0.f;
	float xP1 = 0.f;
	float xF1 = 0.f;

	switch (order) 
	{
		case 0: 
		{
			for(k=0;k<m_lDepth;k++)
			{
				for(i=0;i<m_lWidth;i++)
				{
					for(j=0;j<m_lHeight;j++)
					{
						if(flag)
						{
							wP1 = pData[m_lHeight*m_lWidth*k+m_lHeight*j+i]/sqrt(B);
							wP2 = wP1;
							wP3 = wP1;
							flag = false;
						}

						float xC = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
						float wC = (float)(xC - b1*wP1 - b2*wP2 - b3*wP3)/b0;
						outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)wC;
						wP3 = wP2; wP2 = wP1; wP1 = wC;

						if(j == m_lHeight-1)
						{
							flag = true;

							float up = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i]/(1.0+b1+b2+b3);
							float vp = (float)up/(1.0+b1+b2+b3);
							out = (float)M11*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i]-up) + M12*(outData[m_lHeight*m_lWidth*k+m_lHeight*(j-1)+i]-up) + M13*(outData[m_lHeight*m_lWidth*k+m_lHeight*(j-2)+i]-up)+vp;
							outF1 = (float)M21*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i]-up) + M22*(outData[m_lHeight*m_lWidth*k+m_lHeight*(j-1)+i]-up) + M23*(outData[m_lHeight*m_lWidth*k+m_lHeight*(j-2)+i]-up)+vp;
							outF2 = (float)M31*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i]-up) + M32*(outData[m_lHeight*m_lWidth*k+m_lHeight*(j-1)+i]-up) + M33*(outData[m_lHeight*m_lWidth*k+m_lHeight*(j-2)+i]-up)+vp;
							out *= B; outF1 *= B; outF2 *= B;
							outF3 = outF2; outF2 = outF1; outF1 = out;

							outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)out;
						}
					}
					for(j=m_lHeight-1-1;j>=0;j--)
					{
						float wC = (float)outData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
						out = (float)(B*wC - b1*outF1 - b2*outF2 - b3*outF3)/b0;
						outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)out;
						outF3 = outF2; outF2 = outF1; outF1 = out;
					}
				}
			}
		} break;
		case 1:
		{
			float wC = 0.f;
			float xC = 0.f;

			for(k=0;k<m_lDepth;k++)
			{
				for(i=0;i<m_lWidth;i++)
				{
					for(j=0;j<m_lHeight;j++)
					{
						if(flag)
						{
							xP1 = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
							xF1 = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*(j+1)+i];
							wP1 = (float)(xF1 - xP1)/2.0*(1.0+b1+b2+b3);
							wP2 = wP1;
							wP3 = wP1;

							wC = ((xF1- xP1)/2.0 - b1*wP1 - b2*wP1 - b3*wP1)/b0;
							outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;

							flag = false;

							continue;
						}

						if(j == m_lHeight-1)
						{
							flag = true;

							xP1 = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*(j-1)+i];
							xF1 = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
							wC = (float)((xF1- xP1)/2.0 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
							outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;

							float up = (pData[m_lHeight*m_lWidth*k+m_lHeight*j+i] - pData[m_lHeight*m_lWidth*k+m_lHeight*(j-1)+i])/2.0*(1.0+b1+b2+b3);;
							float vp = (float)up/(1.0+b1+b2+b3);
							out = (float)M11*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i]-up) + M12*(outData[m_lHeight*m_lWidth*k+m_lHeight*(j-1)+i]-up) + M13*(outData[m_lHeight*m_lWidth*k+m_lHeight*(j-2)+i]-up)+vp;
							outF1 = (float)M21*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i]-up) + M22*(outData[m_lHeight*m_lWidth*k+m_lHeight*(j-1)+i]-up) + M23*(outData[m_lHeight*m_lWidth*k+m_lHeight*(j-2)+i]-up)+vp;
							outF2 = (float)M31*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i]-up) + M32*(outData[m_lHeight*m_lWidth*k+m_lHeight*(j-1)+i]-up) + M33*(outData[m_lHeight*m_lWidth*k+m_lHeight*(j-2)+i]-up)+vp;
							out *= B; outF1 *= B; outF2 *= B;
							outF3 = outF2; outF2 = outF1; outF1 = out;

							outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)out;
						}
						else
						{
							xP1 = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*(j-1)+i];
							xF1 = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*(j+1)+i];
							wC = (float)((xF1- xP1)/2.0 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
							outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;
						}
					}
					for(j=m_lHeight-1-1;j>=0;j--)
					{
						wC = (float)outData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
						out = (float)(B*wC - b1*outF1 - b2*outF2 - b3*outF3)/b0;
						outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)out;
						outF3 = outF2; outF2 = outF1; outF1 = out;
					}
				}
			}
		} break;
		case 2:
		{
			float xC = 0.f;
			float wC = 0.f;
			for(k=0;k<m_lDepth;k++)
			{
				for(i=0;i<m_lWidth;i++)
				{
					for(j=0;j<m_lHeight;j++)
					{
						if(flag)
						{
							xC = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
							xP1 = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
							float xF1 = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*(j+1)+i];
							wP1 = 0;
							wP2 = wP1;
							wP3 = wP1;

							float wC = (float)((xF1 - 2*xC + xP1) - b1*wP1 - b2*wP2 - b3*wP3)/b0;
							outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;

							flag = false;

							continue;
						}

						if(j == m_lHeight-1)
						{
							flag = true;

							xC = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
							xP1 = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*(j-1)+i];
							xF1 = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
							wC = (float)(xF1 - 2*xC + xP1 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
							outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;

							float up = 0;
							float vp = (float)up/(1.0+b1+b2+b3);
							out = (float)M11*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i]-up) + M12*(outData[m_lHeight*m_lWidth*k+m_lHeight*(j-1)+i]-up) + M13*(outData[m_lHeight*m_lWidth*k+m_lHeight*(j-2)+i]-up)+vp;
							outF1 = (float)M21*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i]-up) + M22*(outData[m_lHeight*m_lWidth*k+m_lHeight*(j-1)+i]-up) + M23*(outData[m_lHeight*m_lWidth*k+m_lHeight*(j-2)+i]-up)+vp;
							outF2 = (float)M31*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i]-up) + M32*(outData[m_lHeight*m_lWidth*k+m_lHeight*(j-1)+i]-up) + M33*(outData[m_lHeight*m_lWidth*k+m_lHeight*(j-2)+i]-up)+vp;
							out *= B; outF1 *= B; outF2 *= B;
							outF3 = outF2; outF2 = outF1; outF1 = out;

							outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)out;
						}
						else
						{
							xC = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
							xP1 = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*(j-1)+i];
							xF1 = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*(j+1)+i];
							wC = (float)(xF1 - 2*xC + xP1 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
							outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;
						}
					}
					for(j=m_lHeight-1-1;j>=0;j--)
					{
						float wC = (float)outData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
						out = (float)(B*wC - b1*outF1 - b2*outF2 - b3*outF3)/b0;
						outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)out;
						outF3 = outF2; outF2 = outF1; outF1 = out;
					}
				}
			}
		} break;
	}
	return 0;
}


int RecursiveGaussian::RecursiveGaussianZ(float* pData, float * outData, int order, float sigma)
{
    const float nsigma = sigma < 0.5f ? 0.5f : sigma;
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
	int size = m_lWidth*m_lHeight*m_lDepth;

	float wP1 = 0.f, wP2 = 0.f, wP3 = 0.f;
	float outF1 = 0.f, outF2 = 0.f, outF3 = 0.f, out = 0.f;
	float xP1 = 0.f;
	float xF1 = 0.f;

	switch (order) 
	{
		case 0: 
		{
			for(i=0;i<m_lWidth;i++)
			{
				for(j=0;j<m_lHeight;j++)
				{
					for(k=0;k<m_lDepth;k++)
					{
						if(flag)
						{
							wP1 = pData[m_lHeight*m_lWidth*k+m_lHeight*j+i]/sqrt(B);
							wP2 = wP1;
							wP3 = wP1;
							flag = false;
						}

						float xC = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
						float wC = (float)(xC - b1*wP1 - b2*wP2 - b3*wP3)/b0;
						outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)wC;
						wP3 = wP2; wP2 = wP1; wP1 = wC;

						if(k == m_lDepth-1)
						{
							flag = true;

							float up = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i]/(1.0+b1+b2+b3);
							float vp = (float)up/(1.0+b1+b2+b3);
							out = (float)M11*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i]-up) + M12*(outData[m_lHeight*m_lWidth*(k-1)+m_lHeight*j+i]-up) + M13*(outData[m_lHeight*m_lWidth*(k-2)+m_lHeight*j+i]-up)+vp;
							outF1 = (float)M21*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i]-up) + M22*(outData[m_lHeight*m_lWidth*(k-1)+m_lHeight*j+i]-up) + M23*(outData[m_lHeight*m_lWidth*(k-2)+m_lHeight*j+i]-up)+vp;
							outF2 = (float)M31*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i]-up) + M32*(outData[m_lHeight*m_lWidth*(k-1)+m_lHeight*j+i]-up) + M33*(outData[m_lHeight*m_lWidth*(k-2)+m_lHeight*j+i]-up)+vp;
							out *= B; outF1 *= B; outF2 *= B;
							outF3 = outF2; outF2 = outF1; outF1 = out;

							outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)out;
						}
					}
					for(k=m_lDepth-1-1;k>=0;k--)
					{
						float wC = (float)outData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
						out = (float)(B*wC - b1*outF1 - b2*outF2 - b3*outF3)/b0;
						outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)out;
						outF3 = outF2; outF2 = outF1; outF1 = out;
					}
				}
			}
		} break;
		case 1:
		{
			float wC = 0.f;
			float xC = 0.f;

			for(i=0;i<m_lWidth;i++)
			{
				for(j=0;j<m_lHeight;j++)
				{
					for(k=0;k<m_lDepth;k++)
					{
						if(flag)
						{
							xP1 = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
							xF1 = (float)pData[m_lHeight*m_lWidth*(k+1)+m_lHeight*j+i];
							wP1 = (float)(xF1 - xP1)/2.0*(1.0+b1+b2+b3);
							wP2 = wP1;
							wP3 = wP1;

							wC = ((xF1- xP1)/2.0 - b1*wP1 - b2*wP1 - b3*wP1)/b0;
							outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;

							flag = false;

							continue;
						}

						if(k == m_lDepth-1)
						{
							flag = true;

							xP1 = (float)pData[m_lHeight*m_lWidth*(k-1)+m_lHeight*j+i];
							xF1 = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
							wC = (float)((xF1- xP1)/2.0 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
							outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;

							float up = (pData[m_lHeight*m_lWidth*k+m_lHeight*j+i] - pData[m_lHeight*m_lWidth*(k-1)+m_lHeight*j+i])/2.0*(1.0+b1+b2+b3);;
							float vp = (float)up/(1.0+b1+b2+b3);
							out = (float)M11*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i]-up) + M12*(outData[m_lHeight*m_lWidth*(k-1)+m_lHeight*j+i]-up) + M13*(outData[m_lHeight*m_lWidth*(k-2)+m_lHeight*j+i]-up)+vp;
							outF1 = (float)M21*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i]-up) + M22*(outData[m_lHeight*m_lWidth*(k-1)+m_lHeight*j+i]-up) + M23*(outData[m_lHeight*m_lWidth*(k-2)+m_lHeight*j+i]-up)+vp;
							outF2 = (float)M31*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i]-up) + M32*(outData[m_lHeight*m_lWidth*(k-1)+m_lHeight*j+i]-up) + M33*(outData[m_lHeight*m_lWidth*(k-2)+m_lHeight*j+i]-up)+vp;
							out *= B; outF1 *= B; outF2 *= B;
							outF3 = outF2; outF2 = outF1; outF1 = out;

							outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)out;
						}
						else
						{
							xP1 = (float)pData[m_lHeight*m_lWidth*(k-1)+m_lHeight*j+i];
							xF1 = (float)pData[m_lHeight*m_lWidth*(k+1)+m_lHeight*j+i];
							wC = (float)((xF1- xP1)/2.0 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
							outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;
						}
					}
					for(k=m_lDepth-1-1;k>=0;k--)
					{
						wC = (float)outData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
						out = (float)(B*wC - b1*outF1 - b2*outF2 - b3*outF3)/b0;
						outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)out;
						outF3 = outF2; outF2 = outF1; outF1 = out;
					}
				}
			}
		} break;
		case 2:
		{
			float xC = 0.f;
			float wC = 0.f;
			for(i=0;i<m_lWidth;i++)
			{
				for(j=0;j<m_lHeight;j++)
				{
					for(k=0;k<m_lDepth;k++)
					{
						if(flag)
						{
							xC = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
							xP1 = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
							float xF1 = (float)pData[m_lHeight*m_lWidth*(k+1)+m_lHeight*j+i];
							wP1 = 0;
							wP2 = wP1;
							wP3 = wP1;

							float wC = (float)((xF1 - 2*xC + xP1) - b1*wP1 - b2*wP2 - b3*wP3)/b0;
							outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (short)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;

							flag = false;

							continue;
						}

						if(k == m_lDepth-1)
						{
							flag = true;

							xC = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
							xP1 = (float)pData[m_lHeight*m_lWidth*(k-1)+m_lHeight*j+i];
							xF1 = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
							wC = (float)(xF1 - 2*xC + xP1 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
							outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;

							float up = 0;
							float vp = (float)up/(1.0+b1+b2+b3);
							out = (float)M11*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i]-up) + M12*(outData[m_lHeight*m_lWidth*(k-1)+m_lHeight*j+i]-up) + M13*(outData[m_lHeight*m_lWidth*(k-2)+m_lHeight*j+i]-up)+vp;
							outF1 = (float)M21*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i]-up) + M22*(outData[m_lHeight*m_lWidth*(k-1)+m_lHeight*j+i]-up) + M23*(outData[m_lHeight*m_lWidth*(k-2)+m_lHeight*j+i]-up)+vp;
							outF2 = (float)M31*(outData[m_lHeight*m_lWidth*k+m_lHeight*j+i]-up) + M32*(outData[m_lHeight*m_lWidth*(k-1)+m_lHeight*j+i]-up) + M33*(outData[m_lHeight*m_lWidth*(k-2)+m_lHeight*j+i]-up)+vp;
							out *= B; outF1 *= B; outF2 *= B;
							outF3 = outF2; outF2 = outF1; outF1 = out;

							outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)out;
						}
						else
						{
							xC = (float)pData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
							xP1 = (float)pData[m_lHeight*m_lWidth*(k-1)+m_lHeight*j+i];
							xF1 = (float)pData[m_lHeight*m_lWidth*(k+1)+m_lHeight*j+i];
							wC = (float)(xF1 - 2*xC + xP1 - b1*wP1 - b2*wP2 - b3*wP3)/b0;
							outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)wC;
							wP3 = wP2; wP2 = wP1; wP1 = wC;
						}
					}
					for(k=m_lDepth-1-1;k>=0;k--)
					{
						float wC = (float)outData[m_lHeight*m_lWidth*k+m_lHeight*j+i];
						out = (float)(B*wC - b1*outF1 - b2*outF2 - b3*outF3)/b0;
						outData[m_lHeight*m_lWidth*k+m_lHeight*j+i] = (int)out;
						outF3 = outF2; outF2 = outF1; outF1 = out;
					}
				}
			}
		} break;
	}
	return 0;
}

