#include "Denoise.h"
#include "IntegralImg.h"

int main()
{
	//NLM nonmeans;
	//nonmeans.NonLocalMeans();
	cv::Mat pic = cv::imread("lena.jpg");
	cv::Mat gray;
	cv::cvtColor(pic, gray, cv::COLOR_BGR2GRAY);
	cv::Mat res = IntegralImg(gray);
	
	 NLM nonmeans;
	//nonmeans.NonLocalMeans()
	nonmeans.NLMWithIntegralImg(res, gray,5,2);
	return 0;
}
