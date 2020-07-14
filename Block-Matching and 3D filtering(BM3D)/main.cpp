#include "BM3D.h"
#include "IntegralImg.h"

int main()
{
	//NLM nonmeans;
	//nonmeans.NonLocalMeans();
	cv::Mat pic = cv::imread("lena.jpg");
	cv::Mat gray;
	cv::cvtColor(pic, gray, cv::COLOR_BGR2GRAY);
	cv::Mat res = IntegralImg(gray);
	Grouping(res, gray, 9, 5, 1);
	return 0;
}
