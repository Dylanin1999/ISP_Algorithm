#include <iostream>
#include <opencv.hpp>

struct KeyPoint
{
	int DOGLayer;
	int OctaveLayer;
	int row;
	int col;
};


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
	return DOG;
}


std::vector<KeyPoint> GetFeaturePoint(std::vector<std::vector<cv::Mat>>& Pyr)
{
	int arr[8][2] = { { -1, -1 }, { -1, 0 }, { -1, 1 },
	{ 0, -1 }, { 0, 1 },
	{ 1, -1 }, { 1, 0 }, { 1, 1 } };
	std::vector<KeyPoint> FeaturePointVector;
	std::vector<cv::Point2i> tmp;
	std::vector<std::vector<cv::Point2i>> crewPic;

	int kernelSize = 3;

	for (int i = 0; i < Pyr.size(); i++)
	{
		for (int j = 1; j < Pyr[i].size()-1; j++)
		{
			for (int row = kernelSize / 2; row < Pyr[i][j].rows - kernelSize / 2-1; row++)
			{
				
				for (int col = kernelSize / 2; col < Pyr[i][j].cols - kernelSize / 2-1; col++)
				{
					KeyPoint FeaturePoint;
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

					FeaturePoint.col = col;
					FeaturePoint.row = row;
					FeaturePoint.DOGLayer = i;
					FeaturePoint.OctaveLayer = j;
					FeaturePointVector.push_back(FeaturePoint);
				}
				
			}
		}
	}
	return FeaturePointVector;
}

void FirstDifference(std::vector<KeyPoint> FeatureKeyPoint, std::vector<std::vector<cv::Mat>> Pyr)
{
	for (int i = 0; i < FeatureKeyPoint.size(); i++)
	{
		int col = FeatureKeyPoint[i].col;
		int row = FeatureKeyPoint[i].row;
		int DOGLayer = FeatureKeyPoint[i].DOGLayer;
		int OctaveLayer = FeatureKeyPoint[i].OctaveLayer;
		//计算一阶差分
		float dx = Pyr[DOGLayer][OctaveLayer].ptr<uchar>[row](col + 1) - Pyr[DOGLayer][OctaveLayer].ptr<uchar>[row](col - 1)*0.5f;
		float dy = Pyr[DOGLayer][OctaveLayer].ptr<uchar>[row + 1](col) - Pyr[DOGLayer][OctaveLayer].ptr<uchar>[row - 1](col)*0.5f;
		float ds = Pyr[DOGLayer][OctaveLayer-1].ptr<uchar>[row](col) - Pyr[DOGLayer][OctaveLayer+1].ptr<uchar>[row](col - 1);

		//计算二阶差分
		float PixValue2 = 2.0f*Pyr[DOGLayer][OctaveLayer].ptr<uchar>[row](col);
		float dxx = Pyr[DOGLayer][OctaveLayer].ptr<uchar>[row](col + 1) + Pyr[DOGLayer][OctaveLayer].ptr<uchar>[row](col - 1) - PixValue2;
		float dyy = Pyr[DOGLayer][OctaveLayer].ptr<uchar>[row + 1](col) + Pyr[DOGLayer][OctaveLayer].ptr<uchar>[row - 1](col) - PixValue2;
		float dss = Pyr[DOGLayer][OctaveLayer + 1].ptr<uchar>[row](col) + Pyr[DOGLayer][OctaveLayer - 1].ptr<uchar>[row](col) - PixValue2;

		float dxy = (Pyr[DOGLayer][OctaveLayer].ptr<uchar>[row + 1](col + 1) - Pyr[DOGLayer][OctaveLayer].ptr<uchar>[row + 1](col - 1) -
			Pyr[DOGLayer][OctaveLayer].ptr<uchar>[row - 1](col + 1) + Pyr[DOGLayer][OctaveLayer].ptr<uchar>[row - 1](col - 1))*0.25f;

		float dxs = (Pyr[DOGLayer][OctaveLayer+1].ptr<uchar>[row ](col + 1) - Pyr[DOGLayer][OctaveLayer+1].ptr<uchar>[row](col - 1) -
			Pyr[DOGLayer][OctaveLayer-1].ptr<uchar>[row](col + 1) + Pyr[DOGLayer][OctaveLayer-1].ptr<uchar>[row](col - 1))*0.25f;

		float dys = (Pyr[DOGLayer][OctaveLayer + 1].ptr<uchar>[row+1](col) - Pyr[DOGLayer][OctaveLayer + 1].ptr<uchar>[row-1](col) -
			Pyr[DOGLayer][OctaveLayer - 1].ptr<uchar>[row+1](col) + Pyr[DOGLayer][OctaveLayer - 1].ptr<uchar>[row-1](col))*0.25f;


		float H[3][3] = { { dxx, dxy, dxs }, { dxy, dyy, dys }, { dxs, dys, dss } };
	}
}

void HessianMatrix()
{
	
}

int main()
{
	cv::Mat src = cv::imread("timg.jpg");
	std::vector<std::vector<cv::Mat>> Pyr = BuildGaussian(src, 5);
	std::vector<std::vector<cv::Mat>> Res = BuildDOG(Pyr);
	std::vector<KeyPoint> Feature = GetFeaturePoint(Res);
	
	/*for (int i = 0; i < Feature.size(); i++)
	{
		for (int j = 0; j < Feature[i].size(); j++)
		{
			std::cout << "picRow: " << Res[i][j].rows << ",picCol: " << Res[i][j].cols << std::endl;
			std::cout << "sizeA: " << Feature.size() << ", sizeB: " << Feature[i].size() << ", sizeC: " << Feature[i][j].size() << std::endl;
		}
	}*/
	system("pause");
	cv::waitKey(0);
	return 0;
}