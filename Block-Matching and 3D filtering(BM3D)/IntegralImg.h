#pragma once
#include <iostream>
#include <opencv.hpp>

using namespace std;

cv::Mat IntegralImg(cv::Mat src)
{
	int Height = src.rows;
	int Width = src.cols;
	src.convertTo(src, CV_32F);
	//计算出积分图
	//横向添加
	for (int i=0;i<Height;i++)
	{
		for (int j = 1; j < Width; j++)
		{
			src.at<float>(i, j) += src.at<float>(i, j - 1);
		}
	}
	//竖向相加
	for (int i = 1; i < Height; i++)
	{
		for (int j = 0; j < Width; j++)
		{
			src.at<float>(i, j) += src.at<float>(i-1, j);
		}
	}

	//cv::imshow("integralImage", src);
	return src;
}