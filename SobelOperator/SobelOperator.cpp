#include <opencv2/opencv.hpp>
#include <iostream>


enum SobelModel
{
	Horizon = 1,
	Vertical,
	Sobel
};


int HorizonFiter[3][3] = { {-1,0,1},
					  {-2,0,2},
					  {-1,0,1}};

int VerticalFilter[3][3] = { {1,2,1},
					  {0,0,0},
					  {-1,-2,-1} };



int Caculate(const cv::Mat ROI,int mode)
{
	int GXvalue{ 0 };
	int GYvalue{ 0 };
	for (int i{ 0 }; i < 3; i++)
	{
		for (int j{ 0 }; j < 3; j++)
		{
			int pixelValue= ROI.at<float>(i,j);
			int HorizonFiterValue = HorizonFiter[i][j];
			int VerticalFilterValue = VerticalFilter[i][j];
			GXvalue += HorizonFiterValue * pixelValue;
			GYvalue += VerticalFilterValue * pixelValue;
		}
	}
	if (Horizon == mode)
	{
		return int(std::sqrt(std::pow(GXvalue,2)));
	}
	else if (Vertical == mode)
	{
		return int(std::sqrt(std::pow(GYvalue, 2)));
	}
	else
	{
		return int(std::sqrt(std::pow(GYvalue, 2))+ std::sqrt(std::pow(GXvalue, 2)));
	}
}



void SobelFilter(const cv::Mat src, cv::Mat& dst, int mode)
{
	cv::Mat temp;
	if (src.channels()==3)
	{
		cv::cvtColor(src, temp, cv::COLOR_RGB2GRAY);
	}
	else if (src.channels()==1)
	{
		temp = src.clone();
	}

	int picHeight = temp.rows;
	int picWidth = temp.cols;
	cv::Mat dstImg = cv::Mat::zeros(cv::Size(picHeight, picWidth), CV_8UC1);
	temp.convertTo(temp, CV_32F);
	for (int i{ 1 };i<picWidth-1;i+=3)
	{
		for (int j{ 1 };j<picHeight-1;j+=3)
		{
			cv::Mat ROI = temp(cv::Rect(i - 1, j - 1, 3, 3)).clone();
			dstImg.at<uchar>(j, i) = Caculate(ROI, 1);
		}
	}
	dst = dstImg.clone();
}


int main()
{
	std::cout << "HorizonFiter: " << HorizonFiter[0][0] << std::endl;
	cv::Mat SrcPic = cv::imread("lena.jpg");
	
	cv::Mat dstPic;
	SobelFilter(SrcPic, dstPic, 3);
	cv::imshow("result", dstPic);
	cv::waitKey(0);
	return 0;
}
