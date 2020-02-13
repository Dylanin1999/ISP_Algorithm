#include <iostream>
#include <opencv.hpp>
#include <string>

#define CompZero(dis, sigma) (((dis-(2*sigma*sigma))<0)?0:(dis-(2*sigma*sigma)))

void NL_Means(unsigned char* srcImg, int Height,int Width, int D, int d,int sigma)
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
	float sumCP = 0;
	float pixelMinus = 0;
	float distances = 0;
	float DisWeight = 0;
	float h = 0.4*sigma;
	float HH = h*h;
	float pixelSum = 0;
	float Up = 0;
	for (int i{ 0 };i<Height-D;i++)
	{
		if (i<D)
			continue;
		for (int j{ 0 };j<Width-D;j++)
		{
			if (j<D)
				continue;
			//单个像素的search zone求值
			for (int block_i = i-D+d;block_i<i+D-d;block_i++) //blcok的中心纵坐标j的在search zone中的变化范围
			{
				for (int block_j = j - D + d; block_j < j + D - d; block_j++) //blcok的中心横坐标i的在search zone中的变化范围
				{
					//单个block求值
					for (int zone_i = -d;zone_i<=d;zone_i++) 
					{
						for (int zone_j = -d; zone_j <=d; zone_j++)
						{
							pixelMinus += pow(CopyImg[(block_i + zone_i)*(block_j + zone_j)] - CopyImg[i*j],2); //单个zone的所有像素差值相加
						}
					}
					distances = pixelMinus/(BlockSize*BlockSize); //将block内的像素差求平均
					pixelMinus = 0;
					DisWeight = exp(-(CompZero(distances,sigma) / HH));//与block的欧氏距离
					sumCP += DisWeight; //将search zone求得的所有权重求和
					pixelSum += CopyImg[i*j] * DisWeight; //将像素值和权值相乘   Σui(q)*w(p,q)
				}
				std::cout << "processing : " << double(i*j) / (Height*Width) * 100 << "%" << std::endl;
			}
			Up = pixelSum / sumCP;
			srcImg[i*j] = Up;
			Up = 0;
			sumCP = 0;
		}
	}
}

int main()
{
	
	cv::Mat pic = cv::imread("2.jpg");
	cv::Mat gray;
	cv::cvtColor(pic, gray, CV_BGR2GRAY);
	cv::imshow("before", gray);
	std::cout << "cols: " << gray.cols << "rows: " << gray.rows << std::endl;
	NL_Means(gray.data, gray.cols, gray.rows, 10, 1,10);
	cv::imshow("after", gray);
	cv::waitKey(0);
	return 0;
}