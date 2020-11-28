#include <iostream>
#include <opencv.hpp>

std::vector<std::vector<cv::Mat>> BuildGaussian(cv::Mat srcImage, int K)
{
	//进行高斯滤波操作
	float sigma = 1.6;
	std::vector<std::vector<cv::Mat>> pyr;
	std::vector<cv::Mat> temp;
	cv::Mat tmp;
	//先向上采样
	cv::pyrUp(srcImage, srcImage, cv::Size(srcImage.cols * 2, srcImage.rows * 2));
	//std::vector<std::vector<cv::Mat>> pyr;

	//pyr.push_back(DstImg);
	for (int i = 0; i < K; i++)
	{
		for (int j = 0; j < K; j++)
		{
			cv::GaussianBlur(srcImage, tmp, cv::Size(15, 15), j*(sigma), j*(sigma));
			temp.push_back(tmp);
			srcImage = tmp.clone();
		}
		pyr.push_back(temp);
		temp.clear();
		cv::pyrDown(srcImage, srcImage, cv::Size(int(srcImage.cols / 2), int(srcImage.rows / 2)));
	}
	return pyr;
}

int main()
{
	cv::Mat src = cv::imread("timg.jpg");
	std::vector<std::vector<cv::Mat>> Res = BuildGaussian(src, 3);
	
	cv::imshow("-1", src);
	cv::imshow("0", Res[0][0]);
	cv::imshow("1", Res[1][0]);
	cv::imshow("2", Res[2][0]);
	cv::waitKey(0);
	return 0;
}