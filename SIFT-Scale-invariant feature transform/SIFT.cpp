#include <iostream>
#include <opencv.hpp>

std::vector<std::vector<cv::Mat>> BuildGaussian(cv::Mat src, int K)
{
	//进行高斯滤波操作
	float sigma = 1.6;
	std::vector<std::vector<cv::Mat>> pyr;
	std::vector<cv::Mat> temp;
	cv::Mat srcImage;
	cv::Mat tmp;
	//先向上采样
	cv::cvtColor(src, srcImage, cv::COLOR_BGR2GRAY);
	cv::pyrUp(srcImage, srcImage, cv::Size(srcImage.cols * 2, srcImage.rows * 2));
	cv::Mat Mid;
	for (int i = 0; i < K; i++)
	{
		for (int j = 0; j < K; j++)
		{
			cv::GaussianBlur(srcImage, tmp, cv::Size(15, 15), j*(sigma), j*(sigma));
			cv::Mat tt;
			tmp.copyTo(tt);
			temp.push_back(tt);
			if (j == 2)
			{
				Mid = tmp.clone();
			}
			srcImage = tmp.clone();
		}
		pyr.push_back(temp);
		temp.clear();
		cv::pyrDown(Mid, srcImage, cv::Size(int(srcImage.cols / 2), int(srcImage.rows / 2)));
	}

	return pyr;
}


std::vector<std::vector<cv::Mat>> BuildDOG(std::vector<std::vector<cv::Mat>> &Pyr)
{
	std::vector<cv::Mat> Tmp;
	std::vector<std::vector<cv::Mat>> DOG;
	for (int i = 0; i < Pyr.size(); i++)
	{
		for (int j = 1; j < Pyr[0].size(); j++)
		{
			cv::Mat tmp = Pyr[i][j] - Pyr[i][j - 1];
			Tmp.push_back(tmp);
		}
		DOG.push_back(Tmp);
	}
	return DOG;
}

int main()
{
	cv::Mat src = cv::imread("timg.jpg");
	std::vector<std::vector<cv::Mat>> Pyr = BuildGaussian(src, 5);
	std::vector<std::vector<cv::Mat>> Res = BuildDOG(Pyr);
	cv::imshow("-1", src);
	cv::imshow("0", Res[0][1]);
	cv::imshow("1", Res[1][1]);
	cv::imshow("2", Res[2][1]);
	cv::waitKey(0);
	return 0;
}