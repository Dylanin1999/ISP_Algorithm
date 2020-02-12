#include <iostream>
#include <opencv.hpp>
#include <string>


#define MIN2(a, b) ((a) < (b) ? (a) : (b))
#define MAX2(a, b) ((a) > (b) ? (a) : (b))
#define CLIP3(x, a, b) MIN2(MAX2(a,x), b)

void NL_Means(unsigned char* srcImg, int Height,int Width, int D, int d,int h)
{
	unsigned char* afterProcess = (unsigned char*)malloc(sizeof(unsigned char)* Height*Width);
	int ImgSize = Height*Width;
	int SearchZoneSize = 2*D+1;
	int BlockSize = 2 * d + 1;
	int sumCP = 0;
	int pixelMinus = 0;
	int distances = 0;
	int DisWeight = 0;
	int HH = h*h;
	int pixelSum = 0;
	int Up = 0;
	for (int i{ 0 };i<Height-D;i++)
	{
		if (i<D)
			continue;
		for (int j{ 0 };j<Width-D;j++)
		{
			if (j<D)
				continue;
			for (int block_i = i-D-d;block_i<i+D-d;block_i++)
			{
				for (int block_j = j - D -d; block_j < j + D - d; block_j++)
				{
					for (int zone_i = 0;zone_i<2 * d + 1;zone_i++)
					{
						for (int zone_j = 0; zone_j < 2 * d + 1; zone_j++)
						{
							pixelMinus += pow(srcImg[(block_i + zone_i)*(block_j + zone_j)] - srcImg[i*j],2); //单个zone的所有像素差值相加
						}
					}
					distances = pixelMinus/(BlockSize*BlockSize);
					DisWeight = exp(-(distances / HH));
				}
				sumCP += DisWeight;
				pixelSum += srcImg[i*j] * DisWeight;
			}
			Up = pixelSum / sumCP;
			afterProcess[i*j] = Up;
		}
	}
}

int main()
{
	
	cv::Mat pic = cv::imread("2.jpg");
	cv::Mat gray;
	cv::cvtColor(pic, gray, CV_BGR2GRAY);
	cv::imshow("before", gray);
	NL_Means(gray.data, gray.cols, gray.rows, 10, 5, 2);
	cv::imshow("after", gray);
	cv::waitKey(0);
	return 0;
}