#include <opencv.hpp>
#include <iostream>
#include <map>
#include <math.h>
#include "IntegralImg.h"
#include "Config.h"
#include "Hadamard.h"
#include "HardThreshold.h"
#include "addnoise.h"
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


std::vector<std::vector<::Mat>> GetAllBlock(cv::Mat src, int blocksize, int stride, std::vector<cv::Point2d>& blockLUPoint)
{
	std::vector<std::vector<::Mat>> block;
	for (int i = 0; i < src.rows - blocksize;i+=stride)
	{
		block.push_back(std::vector<cv::Mat>());
		for (int j=0;j<src.cols - blocksize;j+= stride)
		{
			cv::Mat Tmp = src(cv::Rect(i, j, blocksize, blocksize));
			block[i].emplace_back(Tmp);
			blockLUPoint.emplace_back(::Point2d(i, j));
		}
	}
	return block;
}

void BlockDctTrans(std::vector<cv::Mat>& block)
{
	for (int i = 0; i < block.size();i++)
	{
		cv::dct(block[i], block[i]);
	}
}

float calDistance(cv::Mat Center,cv::Mat Block)
{
	if (Center.rows!=Block.rows|| Center.cols != Block.cols)
	{
		return 0;
	}
	float distance = 0;
	for (int i = 0; i < Center.rows; i++)
	{
		for (int j = 0; j < Center.cols; j++)
		{
			distance += pow((Center.at<float>(i, j) - Block.at<float>(i, j)), 2);
		}
	}
	return distance / (Center.rows * Center.cols);
}


std::multimap<float, cv::Mat>  Group(std::vector<std::vector<cv::Mat>> blockVector, int BlockIndexX, int BlockIndexY,int SearchWindowSize,int BlockSize,int stride)
{
	int halfNum = ((SearchWindowSize - BlockSize) / stride) / 2;
	int minIndexX = std::max(0, BlockIndexX - halfNum);
	int maxIndexX = std::min(int(blockVector.size() - 1), BlockIndexX + halfNum);

	int minIndexY = std::max(0, BlockIndexY - halfNum);
	int maxIndexY = std::min(int(blockVector.size() - 1), BlockIndexY + halfNum);

	std::multimap<float, cv::Mat> sortGroup;

	//下面是正式的计算Group方式
	for (int i=minIndexX;i<maxIndexX;i++)
	{
		for (int j=minIndexY;i<maxIndexY;j++)
		{
			float distance = calDistance(blockVector[BlockIndexX][BlockIndexY], blockVector[i][j]);
			sortGroup.insert(std::pair<float, cv::Mat>(distance, blockVector[BlockIndexX][BlockIndexY]));
		}
	}
	return sortGroup;
}


int main()
{
	cv::Mat pic = cv::imread("house.png", 0);
	pic = addGaussianNoise(pic);
	cv::Mat src = pic.clone();
	int Height = pic.rows;
	int Width = pic.cols;

	std::vector<std::vector<BlockMessage>> Block;
	return 0;
}

