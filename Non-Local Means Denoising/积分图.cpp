#include <opencv.hpp>
#include <iostream>

void PrintResult(unsigned char* out, int Height, int Width)
{
	std::cout << "size of input: " << sizeof(out) << std::endl;
	for (int i = 0; i < Height; i++)
	{
		for (int j = 0; j < Width; j++)
		{
			std::cout << (int)out[i*Width + j];
			std::cout << " ";
		}
		std::cout << std::endl;
	}
}


void Intergal(unsigned char* Inputdata, int Height, int Width)
{
	unsigned char* Outputdata = (unsigned char*)malloc((Height + 1)*(Width + 1) * sizeof(unsigned char));
	int outWidth = Width + 1;
	int outHeight = Height + 1;
	for (int i{ 0 }; i < outWidth; i++)
	{
		//第一行和第一列清零
		Outputdata[i] = 0;
		Outputdata[i*outWidth] = 0;
	}
	for (int i = 1;i<outHeight;i++)
	{
		for (int j = 1;j<outWidth;j++)
		{
			Outputdata[i*outWidth + j] = Outputdata[(i - 1)*(outWidth)+j] + Outputdata[i*outWidth + (j - 1)] - Outputdata[(i - 1)*outWidth + (j - 1)] + Inputdata[(i - 1)*Width + (j - 1)];
		}
	}
	PrintResult(Outputdata, outHeight, outWidth);
	 return ;
}


int main()
{
	uchar mat[3][3] = { {1,2,3},{4,5,6},{7,8,9} };
	cv::Mat pic(cv::Size(3, 3), CV_8UC1, mat);
	Intergal(pic.data, 3, 3);
	return 0;
}
