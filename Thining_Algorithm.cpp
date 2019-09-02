#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>

using namespace std;
using namespace cv;
int main()
{
	//读入图像
	cv::Mat img = cv::imread("in.png");
	cv::Mat gray;
	cv::Mat Binary;
	//进行二值化
	cv::cvtColor(img, gray, CV_BGR2GRAY);
	cv::inRange(gray, 200, 255, Binary);
	cv::Mat CopyImg;
	int rows = Binary.rows;
	int cols = Binary.cols;


	std::vector<cv::Point2l> PointSaver1;
	//开始第一轮判断
	for (int i{ 1 };i<rows-1;i++)
	{
		uchar* high = Binary.ptr<uchar>(i-1);
		uchar* mid = Binary.ptr<uchar>(i);
		uchar* low = Binary.ptr<uchar>(i+1);

		for (int j{1};j<cols-1;j++)
		{
			int a1 = mid[j];
			int a2 = high[j];
			int a3 = high[j+1];
			int a4 = mid[j+1];
			int a5 = low[j+1];
			int a6 = low[j];
			int a7 = low[j-1];
			int a8 = mid[j-1];
			int a9 = high[j-1];
			int a[9] = { a1,a2,a3,a4,a5,a6,a7,a8,a9 };

			bool req1 = true;
			bool req2 = true;
			bool req3 = true;
			bool req4 = true;

			//条件1 八领域的和
			if (a[0 == 255])
			{
				int req1_sum{ 0 };
				req1_sum = a[1] + a[2] + a[3] + a[4] + a[5] + a[6] + a[7] + a[8];

				//std::cout << "req1_sum: " << req1_sum << std::endl;
				if (req1_sum >= 510 && req1_sum <= 1530)
				{
					//std::cout << "66666" << std::endl;
					req1 = true;
				}
				else req1 = false;

				//条件2  01的模式
				int req2_sum{ 0 };
				//	std::cout << "req2_sum: " << req2_sum << std::endl;
				for (int k = 2; k < 9; k++)
				{
					if (a[k] == 255 && a[k - 1] == 0)
					{
						req2_sum += 1;
					}
				}
				if (req2_sum == 1) req2 = true;
				else req2 = false;
				//条件三
				int req3_sum = a[1] * a[3] * a[5];
				//std::cout << "req3_sum: " << req3_sum << std::endl;

				if (req3_sum == 0) req3 = true;
				else req3 = false;
				//条件四
				int req4_sum = a[3] * a[5] * a[7];
				//	std::cout << "req3_sum: " << req3_sum << std::endl;

				if (req4_sum == 0) req4 = true;
				else req4 = false;
				if (req1 && req2 && req3 && req4)
				{
					PointSaver1.push_back(Point2l(i, j));
				}
			}
		}
	}
	for (int l=0;l<PointSaver1.size();l++)
	{
		uchar* ptr = Binary.ptr<uchar>(PointSaver1[l].x);
		std::cout << "PointSaver1[l].x: " << PointSaver1[l].x << "PointSaver1[l].y: " << PointSaver1[l].y << std::endl;
		ptr[PointSaver1[l].y] = 0;
	}

	//第二轮判断
	std::vector<cv::Point2l> PointSaver2;
	for (int i{ 1 }; i < rows - 1; i++)
	{
		uchar* high = Binary.ptr<uchar>(i - 1);
		uchar* mid = Binary.ptr<uchar>(i);
		uchar* low = Binary.ptr<uchar>(i + 1);

		for (int j{ 1 }; j < cols - 1; j++)
		{
			int a1 = mid[j];
			int a2 = high[j];
			int a3 = high[j + 1];
			int a4 = mid[j + 1];
			int a5 = low[j + 1];
			int a6 = low[j];
			int a7 = low[j - 1];
			int a8 = mid[j - 1];
			int a9 = high[j - 1];
			int a[9] = { a1,a2,a3,a4,a5,a6,a7,a8,a9 };

			bool req1 = true;
			bool req2 = true;
			bool req3 = true;
			bool req4 = true;

			//条件1 八领域的和
			if (a[0 == 255])
			{
				int req1_sum{ 0 };
				req1_sum = a[1] + a[2] + a[3] + a[4] + a[5] + a[6] + a[7] + a[8];

				//std::cout << "req1_sum: " << req1_sum << std::endl;
				if (req1_sum >= 510 && req1_sum <= 1530)
				{
					//std::cout << "66666" << std::endl;
					req1 = true;
				}
				else req1 = false;

				//条件2  01的模式
				int req2_sum{ 0 };
				//	std::cout << "req2_sum: " << req2_sum << std::endl;
				for (int k = 2; k < 9; k++)
				{
					if (a[k] == 255 && a[k - 1] == 0)
					{
						req2_sum += 1;
					}
				}
				if (req2_sum == 1) req2 = true;
				else req2 = false;
				//条件三
				int req3_sum = a[1] * a[3] * a[5];
				//std::cout << "req3_sum: " << req3_sum << std::endl;

				if (req3_sum == 0) req3 = true;
				else req3 = false;
				//条件四
				int req4_sum = a[1] * a[3] * a[7];
				//	std::cout << "req3_sum: " << req3_sum << std::endl;

				if (req4_sum == 0) req4 = true;
				else req4 = false;
				if (req1 && req2 && req3 && req4)
				{
					PointSaver2.push_back(Point2l(i, j));
				}
			}
		}
	}
	for (int l = 0; l < PointSaver2.size(); l++)
	{
		uchar* ptr = Binary.ptr<uchar>(PointSaver2[l].x);
		std::cout << "PointSaver[l].x: " << PointSaver2[l].x << "PointSaver[l].y: " << PointSaver2[l].y << std::endl;
		ptr[PointSaver2[l].y] = 0;
	}

	cv::imshow("SrcImg", Binary);
	std::cout << "saver_size: " << PointSaver2.size() << std::endl;
	cv::waitKey(0);
}
