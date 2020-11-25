#include <iostream>
#include <opencv.hpp>

std::vector<cv::Mat> BuildGaussian(cv::Mat srcImage, int K)
{
	//进行高斯滤波操作
	cv::Mat DstImg;
	cv::Mat tmp;
	//先向上采样
	cv::pyrUp(srcImage, DstImg, cv::Size(srcImage.row * 2, srcImage.col * 2));
	std::vector<cv::Mat> pyr;
	for (int i = 0; i < K; i++)
	{
		cv::GaussianBlur(DstImg, tmp, cv::Size(3,3), 0, 0);
		pyr.push_back(tmp);
		DstImg = tmp.clone();
	}
	return pyr;
}