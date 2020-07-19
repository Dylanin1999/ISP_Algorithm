#include <opencv.hpp>
#include <iostream>
#include <map>

#include "IntegralImg.h"

void Grouping(cv::Mat IntegralImg, cv::Mat grayImg, int Ds, int ds, int stride)
{	
	float GroupNum = 10;
	cv::Mat dst = grayImg.clone();
	
	int Height = grayImg.rows;
	int Width = grayImg.cols;
	std::cout << "rows: " << Height << ", cols: " << Width << std::endl;
	float h = 5;

	int blockNum = pow((Ds - ds)*2,2);//水平方向的block个数，一共有blockNum平方个block，选取这个前面的几个
	std::cout << "Block Num:　" << blockNum << std::endl;
	double Hardthreshold = 0;


	double kernelValue = 1.0 / ((2 * ds + 1) * (2 * ds + 1));
	std::cout << "ds: " << ds << std::endl;
	std::cout << "kernelValue: " << kernelValue << std::endl;

	std::vector<std::vector<float> > Weight;
	Weight.resize(2 * (Ds - ds) + 1);
	for (int k = 0; k < 2 * (Ds - ds) + 1; k++)
	{
		Weight[k].resize(2 * (Ds - ds) + 1);
	}
	for (int k = 0; k < 2 * (Ds - ds) + 1; k++)
	{
		for (int j = 0; j < 2 * (Ds - ds) + 1; j++)
		{
			Weight[k][j] = 0;
		}
	}


	for (int row = Ds; row < Height - Ds; row++)
	{
		for (int col = Ds; col < Width - Ds; col++)
		{
			std::multimap<double, int> disMap;
			std::vector<cv::Mat> GroupingRoi;
			int counter = 0;
			//确定block的中心
			float max = 0;
			float sum = 0;
			for (int SearchWindowsX = -Ds + ds; SearchWindowsX < Ds - ds; SearchWindowsX += stride)
			{
				for (int SearchWindowsY = -Ds + ds; SearchWindowsY < Ds - ds; SearchWindowsY += stride)
				{
					float distance = 0;
					//slip的中心
					int blockCenterX = row + SearchWindowsX;
					int blockCenterY = col + SearchWindowsY;

					//固定参考块的左上角和右下角
					int blockLUPouintX = row - ds;
					int blockLUPouintY = col - ds;
					int blockRDPouintX = row + ds;
					int blockRDPouintY = col + ds;

					float blockValue = IntegralImg.at<float>(blockRDPouintX, blockRDPouintY) + IntegralImg.at<float>(blockLUPouintX, blockLUPouintY) - IntegralImg.at<float>(blockLUPouintX, blockRDPouintY) - IntegralImg.at<float>(blockRDPouintX, blockLUPouintY);
					float slipValue = IntegralImg.at<float>(blockCenterX + ds, blockCenterY + ds) - IntegralImg.at<float>(blockCenterX + ds, blockCenterY - ds) - IntegralImg.at<float>(blockCenterX - ds, blockCenterY + ds) + IntegralImg.at<float>(blockCenterX - ds, blockCenterY - ds);


					distance = kernelValue * pow((blockValue - slipValue), 2);
					//计算出的一个block中的近似块的distance
					//根据超参数选择前面的近似块组合成三维数组
					disMap.insert(std::pair<double, int>(distance, counter));
					counter++;
				}
			}
			int t = 0;
			for (auto elem : disMap)
			{
				int Block_num = elem.second;
				int len_y = Block_num / 8;
				int len_x = Block_num % 8;
				int CenterX = row - Ds + ds;
				int CenterY = col - Ds + ds;
				int shiftX = CenterX + len_x * stride;
				int shiftY = CenterY + len_y * stride;
				cv::Mat ROI = grayImg(cv::Rect(shiftX - ds, shiftY - ds, 2 * ds + 1, 2 * ds + 1));
				cv::Mat ROIDct;
				cv::dct(ROI, ROIDct);
				GroupingRoi.push_back(ROIDct);
			}
			//std::cout << "disVector size:" << disVector.size() << std::endl;
			//Weight[Ds - ds + 1][Ds - ds + 1] = max;
			////std::cout << "sum: " << sum << std::endl;
			////std::cout << "Weight: " << Weight[Ds - ds + 1][Ds - ds + 1] << std::endl;

			////归一化权重并求得
			//float PixelSum = 0.0;
			//for (int i = 0; i < 2 * (Ds - ds) + 1; i++)
			//{
			//	for (int j = 0; j < 2 * (Ds - ds) + 1; j++)
			//	{
			//		//Weight[i][j] /= sum;
			//		PixelSum += grayImg.at<uchar>(row - Ds + i, col - Ds + j) * Weight[i][j];

			//	}
			//}
			//PixelSum /= sum;
			////std::cout << "PixelSum: " << PixelSum << std::endl;
			////if (PixelSum > 255) std::cout << "PixelSum: " << PixelSum << std::endl;
			////PixelSum > 255 ? PixelSum = 255 : PixelSum=PixelSum;
			//dst.at<uchar>(row, col) = int(PixelSum);
			////std::cout << "row:" << row << ", col: " << col << std::endl;
		}
	}
	cv::imshow("dst", dst);
	cv::waitKey(0);
}
