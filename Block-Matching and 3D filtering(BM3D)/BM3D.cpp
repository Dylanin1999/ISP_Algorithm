#include <opencv.hpp>
#include <iostream>
#include <map>
#include <math.h>
#include "IntegralImg.h"
#include "Config.h"
#include "Hadamard.h"
#include "HardThreshold.h"
int GroupNum = 10;

class BlockMessage
{
public:
	int TopLeftX;
	int TopLeftY;
	int width;
	int height;
	float distance;
	cv::Mat Block;
};


std::vector<BlockMessage> Group(cv::Mat src,int SearchWindowSize,int BlockSize,int stride,int PointX,int PointY, cv::Mat pic, std::vector<std::vector<float>> &TD2)
{
	std::vector<BlockMessage> storage;

	//一整个搜索框
	for (int i = (BlockSize - 1) / 2; i < SearchWindowSize - (BlockSize - 1) / 2; i+=stride)
	{
		for (int j = (BlockSize - 1) / 2; j < SearchWindowSize - (BlockSize - 1) / 2; j += stride)
		{
			//block center point (i,j)
			//每一个小的滑动框
			int distance{ 0 };
			for (int K = -(BlockSize - 1) / 2; K <= (BlockSize - 1) / 2; K++)
			{
				for (int L = -(BlockSize - 1) / 2; L <= (BlockSize - 1) / 2; L++)
				{
					auto blockElement = src.at<uchar>(i + K, j + L);
					auto centerElement = src.at<uchar>(SearchWindowSize / 2 + K, SearchWindowSize / 2 + L);

					distance += pow((centerElement-blockElement), 2);
				}
			}
			//将符合条件的block放入到vector中储存起来。
			BlockMessage BMess;
			//if (distance > threshold)
			//{
				BMess.TopLeftX = PointX + i;
				BMess.TopLeftY = PointY + j;
				BMess.height = BlockSize;
				BMess.width = BlockSize;
				BMess.distance = distance;
				cv::Mat dctPic;
				cv::Mat tmp = pic(cv::Rect(BMess.TopLeftX, BMess.TopLeftY, BMess.height + 1, BMess.width + 1));
				tmp.convertTo(tmp, CV_32F);
				cv::dct(tmp, dctPic);
				storage.emplace_back(BMess);
				for (int i = 0; i < dctPic.cols; i++)
				{
					if(TD2.size()<dctPic.cols)
						TD2.emplace_back(std::vector<float>());
					for (int j = 0; j < dctPic.rows; j++)
					{
						TD2[i].emplace_back(dctPic.at<float>(i, j));
					}
				}
			//}
		}
	}
	return storage;
}


bool checkRange(int PointX, int PointY, int SearchWindowSize, int Width, int Height)
{
	if ((PointX + SearchWindowSize) < Width && (PointY + SearchWindowSize) < Height)
		return true;
	else
		return false;
		
}

void Vec2Mat(std::vector<std::vector<float>> Vec, cv::Mat& MatPic)
{
	
}


int main()
{
	cv::Mat pic = cv::imread("lena.jpg",0);
	int Height = pic.rows;
	int Width = pic.cols;

	std::vector<std::vector<BlockMessage>> Block;
	for (int i = 0; i < (Height-SearchWindowSize); i++)
	{
		for (int j = 0; j < (Width - SearchWindowSize); j++)
		{
			//center point is （i,j)
			std::vector<std::vector<float>> TD2;
			std::vector<int> HadamardMatrix;
			if (checkRange)
			{
				cv::Rect Rect(i, j, SearchWindowSize, SearchWindowSize);
			//	std::cout << "right bottom: " << i + SearchWindowSize << ", " << j + SearchWindowSize << std::endl;;
				cv::Mat temp = pic(Rect);
				auto BlockRes = Group(temp, SearchWindowSize, BlockSize, stride, i, j,pic, TD2);
				Block.emplace_back(BlockRes);
			}
			std::vector<std::vector<float>> Res1;
			for (int K = 0; K < TD2.size(); K++)
			{
				Hadmard(HadamardMatrix, TD2[K].size() / 4, 4);
				for (int L = 0; L < HadamardMatrix.size(); L++)
				{
					if (Res1.size() < BlockSize)
					{
						Res1.emplace_back(std::vector<float>());
					}
					TD2[K][L] = TD2[K][L] * HadamardMatrix[L];
					HardThreshold(TD2[K][L], HardThresholdValue);
					//开始将第三维转换为正常的图像
					Res1[(K * L) / 4].emplace_back(TD2[K][L]);
				}
			}
			std::cout << "TD2 size: " << TD2.size() << ", " << TD2[0].size();
		}
	}
	


}

