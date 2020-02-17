#include <iostream>
#include <opencv.hpp>
#include <string>
#include "ProgressBar.h"
#include <thread>

#define MIN2(a, b) ((a) < (b) ? (a) : (b))
#define MAX2(a, b) ((a) > (b) ? (a) : (b))
#define Compare(x, a, b) MIN2(MAX2(a,x), b)
void NL_Means(unsigned char* srcImg, int Height,int Width, int D, int d,int h)
{
	if (D<=d)
	{
		std::cout << "Please check your params!"<<std::endl;
	}
	std::cout << "start to process!" << std::endl;
	unsigned char* CopyImg = (unsigned char*)malloc(sizeof(unsigned char)* Height*Width);
	memcpy(CopyImg, srcImg, sizeof(unsigned char)* Height*Width);
	float ImgSize = Height*Width;
	float SearchZoneSize = 2*D+1;
	float BlockSize = 2 * d + 1;
	float BlockArea = BlockSize*BlockSize;
	float sumCP = 0;
	float pixelMinus = 0;
	float distances = 0;
	float DisWeight = 0;
	float pixelSum = 0;
	float Up = 0;
	float HH =h*h;
	for (int i{ D };i<Height-D;i++)
	{
		for (int j{ D };j<Width-D;j++)
		{
			Up = 0;
			sumCP = 0;
			pixelSum = 0;
			//单个像素的search zone求值
			for (int SearchZone_i = -D; SearchZone_i <= D; SearchZone_i++) //blcok的中心纵坐标j的在search zone中的变化范围
			{
				for (int SearchZone_j = -D; SearchZone_j <= D; SearchZone_j++) //blcok的中心横坐标i的在search zone中的变化范围
				{
					//单个block求值
					pixelMinus = 0;
					for (int BlockZone_i = -d; BlockZone_i <=d; BlockZone_i++)
					{
						for (int BlockZone_j = -d; BlockZone_j <=d; BlockZone_j++)
						{
							int block_i = Compare((i + SearchZone_i + BlockZone_i), 0, Height-1);
							int block_j = Compare((j + SearchZone_i + BlockZone_i), 0, Width-1);
							int zone_i = Compare((i - d + BlockZone_i), 0, Height);
							int zone_j = Compare((j - d + BlockZone_j), 0, Width);
							pixelMinus += pow(CopyImg[zone_i*Width + zone_j] - CopyImg[block_i*Width + block_j], 2);
						}
					}
					distances = pixelMinus/(d*d); //将block内的像素差求平均
					DisWeight = exp(-(distances / HH));//与block的欧氏距离
					sumCP += DisWeight; //将search zone求得的所有权重求和
					int certern_x = Compare((i + SearchZone_i), 0, Height - 1);
					int certern_y = Compare((j + SearchZone_j), 0, Width - 1);
					pixelSum += CopyImg[certern_x*Width+ certern_y] * DisWeight; //将像素值和权值相乘   Σui(q)*w(p,q)
				}
			}
			
			if (sumCP!=0)
			{
				Up = pixelSum / sumCP;
				Up = Compare(Up, 0, 255);
				srcImg[i*Width + j] = (int)Up;
			}
		}
	}
}

int main()
{
	
	cv::Mat pic = cv::imread("InputImg.jpg"); //输入图像
	cv::Mat outImg;
	std::vector<cv::Mat> colorSp;
	cv::split(pic, colorSp);
	int Height = pic.rows;
	int Width = pic.cols;
	std::thread t1(NL_Means, colorSp[0].data, Height, Width, 10, 3, 10);//三通道多线程
	std::thread t2(NL_Means, colorSp[1].data, Height, Width, 10, 3, 10);
	std::thread t3(NL_Means, colorSp[2].data, Height, Width, 10, 3, 10);
	
	t1.join();
	t2.join();
	t3.join();
	cv::merge(colorSp, outImg);//通道合并
	cv::imshow("after", outImg);//输出图像
	std::cout << "done" << std::endl;
	cv::imwrite("outImg.jpg", outImg); //保存输出图像
	cv::waitKey(0);
	return 0;
}