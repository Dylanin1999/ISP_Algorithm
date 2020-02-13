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
			//�������ص�search zone��ֵ
			for (int block_i = i-D+d;block_i<i+D-d;block_i++) //blcok������������j����search zone�еı仯��Χ
			{
				for (int block_j = j - D + d; block_j < j + D - d; block_j++) //blcok�����ĺ�����i����search zone�еı仯��Χ
				{
					//����block��ֵ
					for (int zone_i = -d;zone_i<=d;zone_i++) 
					{
						for (int zone_j = -d; zone_j <=d; zone_j++)
						{
							pixelMinus += pow(CopyImg[(block_i + zone_i)*(block_j + zone_j)] - CopyImg[i*j],2); //����zone���������ز�ֵ���
						}
					}
					distances = pixelMinus/(BlockSize*BlockSize); //��block�ڵ����ز���ƽ��
					pixelMinus = 0;
					DisWeight = exp(-(CompZero(distances,sigma) / HH));//��block��ŷ�Ͼ���
					sumCP += DisWeight; //��search zone��õ�����Ȩ�����
					pixelSum += CopyImg[i*j] * DisWeight; //������ֵ��Ȩֵ���   ��ui(q)*w(p,q)
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