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
		Tmp.clear();
	}
	for (int i = 0; i<DOG.size(); i++)
	{
		std::cout << "DDDD: " << DOG[i].size() << std::endl;
	}
	return DOG;
}


std::vector<std::vector<std::vector<cv::Point2i>>> GetFeaturePoint(std::vector<std::vector<cv::Mat>>& Pyr)
{
	int arr[8][2] = { { -1, -1 }, { -1, 0 }, { -1, 1 },
	{ 0, -1 }, { 0, 1 },
	{ 1, -1 }, { 1, 0 }, { 1, 1 } };
	std::vector<std::vector<std::vector<cv::Point2i>>> FeaturePoint;
	std::vector<cv::Point2i> tmp;
	std::vector<std::vector<cv::Point2i>> crewPic;

	int kernelSize = 3;
	std::cout << "i: " << Pyr.size() << std::endl;
	std::cout << "j: " << Pyr[0].size() << std::endl;
	for (int i = 0; i < Pyr.size(); i++)
	{
		for (int j = 1; j < Pyr[i].size()-1; j++)
		{
			std::cout << "j: " << Pyr[i].size() << std::endl;

			std::cout << "row: " << Pyr[i][j].rows << std::endl;
			std::cout << "cols: " << Pyr[i][j].cols << std::endl;
			for (int row = kernelSize / 2; row < Pyr[i][j].rows - kernelSize / 2-1; row++)
			{
				
				for (int col = kernelSize / 2; col < Pyr[i][j].cols - kernelSize / 2-1; col++)
				{

					int maxCounter = 0;
					int minCounter = 0;
					for (int k = 0; k < 8; k++)
					{
						Pyr[i][j].ptr<uchar>(row)[col] > Pyr[i][j-1].ptr<uchar>(row + arr[k][0])[col + arr[k][1]] ? maxCounter++ : minCounter++;
						Pyr[i][j].ptr<uchar>(row)[col] > Pyr[i][j].ptr<uchar>(row + arr[k][0])[col + arr[k][1]] ? maxCounter++ : minCounter++;
						Pyr[i][j].ptr<uchar>(row)[col] > Pyr[i][j+1].ptr<uchar>(row + arr[k][0])[col + arr[k][1]] ? maxCounter++ : minCounter++;
						if (maxCounter&&minCounter)
							break;
					}
					if (maxCounter&&minCounter)
						continue;
				}
				tmp.push_back(cv::Point2i(i, j));
			}
			crewPic.push_back(tmp);
			tmp.clear();
		}
		FeaturePoint.push_back(crewPic);
		crewPic.clear();
	}
	return FeaturePoint;
}

int main()
{
	cv::Mat src = cv::imread("timg.jpg");
	std::vector<std::vector<cv::Mat>> Pyr = BuildGaussian(src, 5);
	std::vector<std::vector<cv::Mat>> Res = BuildDOG(Pyr);
	std::vector<std::vector<std::vector<cv::Point2i>>> Feature = GetFeaturePoint(Res);
	
	for (int i = 0; i < Feature.size(); i++)
	{
		for (int j = 0; j < Feature[i].size(); j++)
		{
			std::cout << "picRow: " << Res[i][j].rows << ",picCol: " << Res[i][j].cols << std::endl;
			std::cout << "sizeA: " << Feature.size() << ", sizeB: " << Feature[i].size() << ", sizeC: " << Feature[i][j].size() << std::endl;
		}
	}
	system("pause");
	cv::waitKey(0);
	return 0;
}