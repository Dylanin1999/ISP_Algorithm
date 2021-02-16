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
	for (int i = 0; i < src.rows - blocksize-1;i += stride)
	{
		block.emplace_back(std::vector<cv::Mat>());
		for (int j=0;j<src.cols - blocksize-1;j += stride)
		{
			cv::Mat Tmp = src(cv::Rect(i, j, blocksize, blocksize));
			block[block.size()-1].emplace_back(Tmp);
			blockLUPoint.emplace_back(cv::Point2d(i, j));
		}
	}
	return block;
}

void BlockDctTrans(std::vector<std::vector<::Mat>>& block)
{
	for (int i = 0; i < block.size();i++)
	{
		for (int j=0;j<block.size();j++)
		{
			block[i][j].convertTo(block[i][j], CV_32F);
			cv::dct(block[i][j], block[i][j]);
		}
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

int log2(const int N)
{
	int k = 1;
	int n = 0;
	while (k < N)
	{
		k *= 2;
		n++;
	}
	return n;
}

void wavedec(float* input, int length)
{
	int N = log2(length);
	if (N == 0)return;
	float* tmp = new float[length];
	for (int i = 0; i < length; i++) {
		tmp[i] = input[i];
	}
	for (int k = 0; k < length / 2; k++)
	{
		input[k] = (tmp[2 * k] + tmp[2 * k + 1]) / sqrt(2);
		input[k + length / 2] = (tmp[2 * k] - tmp[2 * k + 1]) / sqrt(2);
	}
	delete[] tmp;
	wavedec(input, length / 2);
	return;
}

void waverec(float* input, int length, int N)
{
	if (log2(length) > N) return;
	float* tmp = new float[length];
	for (int i = 0; i < length; i++) {
		tmp[i] = input[i];
	}
	for (int k = 0; k < length / 2; k++)
	{
		input[2 * k] = (tmp[k] + tmp[k + length / 2]) / sqrt(2);
		input[2 * k + 1] = (tmp[k] - tmp[k + length / 2]) / sqrt(2);
	}
	delete tmp;
	waverec(input, length * 2, N);
}





void tran1d(std::vector<std::vector<cv::Mat>>& block,int BlockSize)
{
	int size = block.size() * block[0].size();
	for (int i=0;i<BlockSize;i++)
	{
		for (int j=0;j< BlockSize;j++)
		{
			float *data = new float[size];
			for (int k = 0;k<block.size();k++)
			{
				for (int L = 0; L < block[k].size(); L++)
				{
					data[k*block.size()+L] = block[k][L].at<float>(i, j);
				}
			}
			wavedec(data, size);
			for (int k = 0; k < block.size(); k++)
			{
				for (int L = 0; L < block[k].size(); L++)
				{
					block[k][L].at<float>(i, j) = data[k * block.size() + L];
				}
			}
		}
	}
}

void shrink(std::vector<std::vector<cv::Mat>>& block, float threshold)
{
	for (int i=0;i<block.size();i++)
	{
		for (int j=0;j<block.size();j++)
		{
			for (int rows = 0;rows<block[i][j].rows;rows++)
			{
				for (int cols = 0;cols<block[i][j].cols;cols++)
				{
					block[i][j].at<float>(rows, cols) = fabs(block[i][j].at<float>(rows, cols)) > threshold ? block[i][j].at<float>(rows, cols) : 0;
				}
			}
		}
	}
}

void aggregation(vector<std::vector<cv::Mat>>& block, std::vector<cv::Point2d> blockLUPoint,Mat &numerator, Mat &denominator,
				 float weight, int patchSize, Mat window)
{
	cv::Rect Rect;
	for (int i = 0;i<block.size();i++)
	{
		for (int j =  0;j<block[0].size();j++)
		{
			Rect.x = blockLUPoint[i * block[0].size() + j].x;
			Rect.y = blockLUPoint[i * block[0].size() + j].y;
			Rect.width = patchSize;
			Rect.height = patchSize;
			numerator(Rect) += weight * (block[i][j].mul(window));
			denominator(Rect) += denominator(Rect) + weight * window;
		}
	}
}

float calculate_weight_hd(vector<std::vector<cv::Mat>>& block, int sigma)
{
	int num{ 0 };
	for (int i = 0; i < block.size(); i++)
	{
		for (int j = 0; j < block[0].size(); j++)
		{
			for (int row = 0;row<block[i][j].rows;row++)
			{
				for (int col = 0;col<block[i][j].cols;col++)
				{
					num++;
				}
			}
		}
	}
	if (num==0)
	{
		return 1;
	}
	else
	{
		return 1.0 / (sigma * sigma * num);
	}
}


Mat gen_kaiser(int beta, int length)//How to do this?
{
	if ((beta == 2) && (length == 8))
	{
		Mat window(length, length, CV_32FC1);
		Mat kai1(length, 1, CV_32FC1);
		Mat kai1_T(1, length, CV_32FC1);
		kai1.at<float>(0, 0) = 0.4387;
		kai1.at<float>(1, 0) = 0.6813;
		kai1.at<float>(2, 0) = 0.8768;
		kai1.at<float>(3, 0) = 0.9858;
		for (int i = 0; i < 4; i++)
		{
			kai1.at<float>(7 - i, 0) = kai1.at<float>(i, 0);
			kai1_T.at<float>(0, i) = kai1.at<float>(i, 0);
			kai1_T.at<float>(0, 7 - i) = kai1.at<float>(i, 0);
		}
		window = kai1 * kai1_T;
		return window;
	}
}


int main()
{
	cv::Mat pic = cv::imread("house.png", 0);
	pic = addGaussianNoise(pic);
	std::vector<cv::Point2d> blockLUPoint;

	std::vector<std::vector<::Mat>> blockGroup = GetAllBlock(pic, BlockSize, stride, blockLUPoint);
	BlockDctTrans(blockGroup);
	tran1d(blockGroup, BlockSize);
	shrink(blockGroup, Threshold);

	Mat denominator_hd(pic.size(), CV_32FC1, Scalar::all(0));
	Mat numerator_hd(pic.size(), CV_32FC1, Scalar::all(0));

	float weight = calculate_weight_hd(blockGroup, sigma);
	int beta = 2;
	cv::Mat kaiser = gen_kaiser(beta, BlockSize);
	aggregation(blockGroup, blockLUPoint, numerator_hd, denominator_hd, weight, BlockSize,kaiser);
	cv::Mat src = pic.clone();
	int Height = pic.rows;
	int Width = pic.cols;
	cv::Mat basic = numerator_hd / denominator_hd;
	cv::imshow("basic", basic);
	std::vector<std::vector<BlockMessage>> Block;
	return 0;
}

