#include <opencv.hpp>
#include <iostream>
#include <map>
#include <math.h>
#include "IntegralImg.h"
#include "Config.h"
#include "Hadamard.h"
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



void CheckRange(int& x, int& y, int ds, int Height, int Width)
{
	if (x < 0) x = 0;
	if (y < 0) y = 0;
	if ((x + ds) > Width) x = Width - ds;
	if ((y + ds) > Height) y = Height - ds;
}


std::vector<BlockMessage> Group(cv::Mat src,int SearchWindowSize,int BlockSize,int stride,int PointX,int PointY)
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
				cv::Mat tmp = src(cv::Rect(BMess.TopLeftX, BMess.TopLeftY, BMess.height + 1, BMess.width + 1));
				tmp.convertTo(tmp, CV_32F);
				cv::dct(tmp, dctPic);
				storage.emplace_back(BMess);
			//}
		}
	}
	return storage;
}

int main()
{
	//cv::Mat pic = cv::imread("lena.jpg",0);
	//int Height = pic.rows;
	//int Width = pic.cols;

	//std::vector<std::vector<BlockMessage>> Block;

	//for (int i = 0; i < (Height-SearchWindowSize); i++)
	//{
	//	for (int j = 0; j < (Width - SearchWindowSize); j++)
	//	{
	//		//center point is （i,j)
	//		cv::Rect Rect(i, j, SearchWindowSize, SearchWindowSize);
	//		cv::Mat temp = pic(Rect);
	//		auto BlockRes = Group(temp, SearchWindowSize, BlockSize, stride, i, j);
	//		Block.emplace_back(BlockRes);
	//		//int num = BlockRes.size();
	//	}
	//}
	std::vector<std::vector<int>> HadamardMatrix;
	Hadmard(HadamardMatrix);
	for (int i = 0; i < HadamardMatrix.size(); i++)
	{
		for (int j = 0; j < HadamardMatrix[0].size(); j++)
		{
			std::cout << HadamardMatrix[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

