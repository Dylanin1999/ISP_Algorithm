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
			//cv::normalize(tt, tt, 0, 255, cv::NORM_MINMAX);
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


std::vector<cv::Point2i> GetFeaturePoint(std::vector<std::vector<cv::Mat>> Pyr)
{
	int arr[8][2] = { { -1, -1 }, { -1, 0 }, { -1, 1 },
	{ 0, -1 }, { 0, 1 },
	{ 1, -1 }, { 1, 0 }, { 1, 1 } };
	std::vector<cv::Point2i> FeaturePoint;

	int kernelSize = 3;
	for (int i = 0; i < Pyr.size(); i++)
	{
		for (int j = 1; j < Pyr[i].size; j++)
		{
			for (int row = kernelSize / 2; row < Pyr[i][j].rows - kernelSize / 2; row++)
			{
				for (int col = kernelSize / 2; col < Pyr[i][j].cols - kernelSize / 2; col++)
				{
					int maxCounter = 0;
					int minCounter = 0;
					for (int k = 0; k < 8; k++)
					{
						Pyr[i][j].ptr<uchar>(row)[col] > Pyr[i][j-1].ptr<uchar>(row + arr[k][0])[col + arr[k][1]] ? maxCounter++ : minCounter++;
						Pyr[i][j].ptr<uchar>(row)[col] > Pyr[i][j].ptr<uchar>(row + arr[k][0])[col + arr[k][1]] ? maxCounter++ : minCounter++;
						Pyr[i][j].ptr<uchar>(row)[col] > Pyr[i][j+1].ptr<uchar>(row + arr[k][0])[col + arr[k][1]] ? maxCounter++ : minCounter++;
						if (!(maxCounter&&minCounter))
							continue;
					}
					FeaturePoint.push_back(cv::Point2i(i, j));
				}
			}
		}
	}
	return FeaturePoint;
}

int main()
{
	cv::Mat src = cv::imread("timg.jpg");
	std::vector<std::vector<cv::Mat>> Pyr = BuildGaussian(src, 5);
	std::vector<std::vector<cv::Mat>> Res = BuildDOG(Pyr);
	cv::imshow("-1", src);
	cv::imshow("0", Res[0][2]);
	cv::imshow("1", Res[1][2]);
	cv::imshow("2", Res[2][2]);
	cv::waitKey(0);
	return 0;
}