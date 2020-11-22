#include <iostream>
#include <opencv.hpp>
#include <vector>
#include <string>

using namespace cv;

int arr[16][2] = { { -3, -1 }, { -3, 0 }, { -3, 1 },
{ -2, 2 },
{ -1, 3 }, { 0, 3 }, { 1, 3 },
{ 2, 2 },
{ 3, 1 }, { 3, 0 }, { 3, -1 },
{ 2, -2 },
{ 1, -3 }, { 0, -3 }, { -1, -3 },
{ -2, -2 } };
std::string ImgPath = "123.png";

int FASTCornerDetector()
{
	cv::Mat src = cv::imread(ImgPath);
	cv::Mat gray;
	cv::cvtColor(src, gray, COLOR_BGR2GRAY);
	std::vector<cv::Point2d> Corner;
	int height = gray.rows;
	int width = gray.cols;

	int radius = 3;
	int threshold = 32;
	for (int i = radius; i < height - radius; i++)
	{
		uchar* topRow = gray.ptr<uchar>(i - radius);
		uchar* buttomRow = gray.ptr<uchar>(i - radius);
		uchar* midCol = gray.ptr<uchar>(i);
		int pixCounter = 0;
		for (int j = radius; j < width - radius; j++)
		{
			pixCounter = 0;
			uchar center = midCol[j];
			int counter = 0;
			counter += topRow[j] - center > threshold ? 0 : 1;
			counter += buttomRow[j] - center > threshold ? 0 : 1;
			counter += midCol[j - radius] - center > threshold ? 0 : 1;
			counter += midCol[j + radius] - center> threshold ? 0 : 1;

			if (counter < 3)
				continue;

			for (int k = 0; k < 27; k++)
			{
				//std::cout << "K: " << k << std::endl;
				uchar pixelValue = gray.ptr<uchar>(i + arr[k%16][0])[j + arr[k%16][1]];
				bool flag = abs(pixelValue - center) > threshold ? true : false;
				if (flag)
					pixCounter += 1;
				else
				{
					if (pixCounter < 9)
					{
						pixCounter = 0;
						continue;
					}
				}
				if (pixCounter >= 9)
				{
					Corner.push_back(cv::Point2d(i, j));
					break;
				}
			}
		}
	}

	for (int i = 0; i < Corner.size(); i++)
	{
		cv::circle(src, Corner[i], 3, cv::Scalar(0, 255, 255));
	}
	std::cout << "size: " << Corner.size() << std::endl;
	cv::imshow("pic", src);

	cv::waitKey(0);
	return 0;

}


int main()
{
	FASTCornerDetector();
	return 0;
}