#include <opencv2\highgui\highgui.hpp>
#include <opencv2\opencv.hpp>
#include <opencv2\core\core.hpp>
#include <vector>
#include <math.h>
#include "bm3d.h"

using namespace cv;
using namespace std;

void addGuassianNoise(const int sigma, const Mat origin, Mat& noisy)
{
	Mat noise(origin.size(), CV_32FC1);
	randn(noise, Scalar::all(0), Scalar::all(sigma));
	for (int i = 0; i < noise.rows; i++)
	{
		const float* Mx = origin.ptr<float>(i);
		float* Mn = noise.ptr<float>(i);
		float* My = noisy.ptr<float>(i);
		for (int j = 0; j < origin.cols; j++)
		{
			My[j] = Mx[j] + Mn[j];
		}
	}
}

float cal_psnr(const Mat x, const Mat y)
{
	float RMS = 0;
	for (int i = 0; i < x.rows; i++)
	{
		const float* Mx = x.ptr<float>(i);
		const float* My = y.ptr<float>(i);
		for (int j = 0; j < x.cols; j++)
		{
			RMS += (My[j] - Mx[j]) * (My[j] - Mx[j]);
		}
	}
	RMS = sqrtf(RMS / (x.rows * x.cols));
	return 20 * log10f(255.0 / RMS);
}

int runBm3d(
	const Mat image_noisy,
	Mat& image_basic,
	Mat& image_denoised
)
{
	

	int Height = image_noisy.rows;
	int Width = image_noisy.cols;
	int Channels = image_noisy.channels();

	
	vector<Mat> block_noisy;//store the patch
	vector<int>row_idx;//patch idx along the row direction
	vector<int>col_idx;
	GetAllBlock(image_noisy, Width, Height, Channels, kHard, pHard, block_noisy, row_idx, col_idx);
	int bn_r = row_idx.size();
	int bn_c = col_idx.size();
	tran2d(block_noisy, kHard);

	vector<int> sim_num;//index number for the selected similar patch in the block vector
	vector<int> sim_idx_row;//index number for the selected similar patch in the original Mat
	vector<int> sim_idx_col;

	vector<Mat>data;//store the data during transforming and shrinking

	Mat kaiser = gen_kaiser(beta, kHard);//2-D kaiser window 
	float weight_hd = 1.0;//weights used for current relevent patch
	Mat denominator_hd(image_noisy.size(), CV_32FC1, Scalar::all(0));
	Mat numerator_hd(image_noisy.size(), CV_32FC1, Scalar::all(0));
	for (int i = 0; i < bn_r; i++)
	{
		for (int j = 0; j < bn_c; j++)
		{
			//for each pack in the block
			sim_num.clear();
			sim_idx_row.clear();
			sim_idx_col.clear();
			data.clear();

			getSimilarPatch(block_noisy, data, sim_num,
				i, j, bn_r, bn_c, int((nHard - kHard) / pHard) + 1, NHard, tao_hard);//block matching

			for (int k = 0; k < sim_num.size(); k++)//calculate idx in the left-top corner
			{
				sim_idx_row.push_back(row_idx[sim_num[k] / bn_c]);
				sim_idx_col.push_back(col_idx[sim_num[k] % bn_c]);
			}

			tran1d(data, kHard);//3-D transforming

			DetectZero(data, lambda3d * sigma);//shrink the cofficient

			weight_hd = calculate_weight_hd(data, sigma);

			Inver3Dtrans(data,kHard);//3-D inverse transforming

			aggregation(numerator_hd, denominator_hd, sim_idx_row, sim_idx_col, data, weight_hd, kHard, kaiser);//aggregation using weigths
		}
	}
	image_basic = numerator_hd / denominator_hd;

	//step 2 wiena filtering
	vector<Mat> block_basic;
	row_idx.clear();
	col_idx.clear();

	GetAllBlock(image_basic, Width, Height, Channels, kHard, pHard, block_basic, row_idx, col_idx);
	bn_r = row_idx.size();
	bn_c = col_idx.size();

	vector<Mat> data_noisy;
	float weight_wien = 1.0;//weights used for current relevent patch
	Mat denominator_wien(image_noisy.size(), CV_32FC1, Scalar::all(0));
	Mat numerator_wien(image_noisy.size(), CV_32FC1, Scalar::all(0));
	for (int i = 0; i < bn_r; i++)
	{
		for (int j = 0; j < bn_c; j++)
		{
			//for each pack in the basic estimate
			sim_num.clear();
			sim_idx_row.clear();
			sim_idx_col.clear();
			data.clear();
			data_noisy.clear();

			getSimilarPatch(block_basic, data, sim_num, i, j, bn_r, bn_c,
				int((nWien - kWien) / pWien) + 1, NWien, tao_wien);//block matching

			for (int k = 0; k < sim_num.size(); k++)//calculate idx in the left-top corner
			{
				sim_idx_row.push_back(row_idx[sim_num[k] / bn_c]);
				sim_idx_col.push_back(col_idx[sim_num[k] % bn_c]);
				data_noisy.push_back(image_noisy(Rect(sim_idx_col[k], sim_idx_row[k], kWien, kWien)));
			}

			tran2d(data,kWien);
			tran2d(data_noisy, kWien);
			tran1d(data, kWien);
			tran1d(data_noisy,kWien);

			gen_wienFilter(data, sigma);
			weight_wien = calculate_weight_wien(data, sigma);
			wienFiltering(data_noisy, data, kWien);

			Inver3Dtrans(data_noisy, kWien);

			aggregation(numerator_wien, denominator_wien,
				sim_idx_row, sim_idx_col, data_noisy, weight_wien, kWien, kaiser);
		}
	}
	image_denoised = numerator_wien / denominator_wien;

	return EXIT_SUCCESS;
}

void GetAllBlock(const Mat img, const int width, const int height, const int channels,
	const int patchSize, const int step, vector<Mat>& block, vector<int>& row_idx, vector<int>& col_idx)
{
	Mat tmp(patchSize, patchSize, CV_32FC1);
	for (int i = 0; i <= height - patchSize; i += step)
	{
		row_idx.push_back(i);
	}
	if ((height - patchSize) % step != 0)
	{
		row_idx.push_back(height - patchSize);
	}
	for (int j = 0; j <= width - patchSize; j += step)
	{
		col_idx.push_back(j);
	}
	if ((width - patchSize) % step != 0)
	{
		col_idx.push_back(width - patchSize);
	}
	for (int i = 0; i < row_idx.size(); i++)
	{
		for (int j = 0; j < col_idx.size(); j++)
		{
			tmp = img(Rect(col_idx[j], row_idx[i], patchSize, patchSize));
			block.push_back(tmp);
		}
	}
}

void tran2d(vector<Mat>& input, int patchsize)
{

	int length = input.size();
	for (int i = 0; i < length; i++)
	{
		dct(input[i], input[i]);
	}
}

void getSimilarPatch(const vector<Mat> block, vector<Mat>& sim_patch, vector<int>& sim_num,
	int i, int j, int bn_r, int bn_c, int area, int maxNum, int tao)
{
	int row_min = max(0, i - (area - 1) / 2);
	int row_max = min(bn_r - 1, i + (area - 1) / 2);
	int row_length = row_max - row_min + 1;

	int col_min = max(0, j - (area - 1) / 2);
	int col_max = min(bn_c - 1, j + (area - 1) / 2);
	int col_length = col_max - col_min + 1;

	const Mat relevence = block[i * bn_c + j];
	Mat tmp;

	float* distance = new float[row_length * col_length];//计算距离
	int* idx = new int[row_length * col_length];//保存下标便于后续聚类计算
	if (!distance) {
		cout << "allocation failure\n";
		system("pause");
	}
	for (int p = 0; p < row_length; p++)
	{
		for (int q = 0; q < col_length; q++)
		{
			tmp = block[(p + row_min) * bn_c + (q + col_min)];
			distance[p * col_length + q] = cal_distance(relevence, tmp);
			//cout << distance[p*col_length + q] << endl;
			idx[p * col_length + q] = p * col_length + q;
		}
	}
	float value; int l;
	//直接排序算法，有待改进！
	for (int k = 1; k < row_length * col_length; k++)
	{
		value = distance[k];
		for (l = k - 1; value < distance[l] && l >= 0; --l)
		{
			distance[l + 1] = distance[l];
			idx[l + 1] = idx[l];
		}
		distance[l + 1] = value;
		idx[l + 1] = k;
	}
	int selectedNum = maxNum;
	while (row_length * col_length < selectedNum)
	{
		selectedNum /= 2;//确保相似块的个数为2的幂
	}
	while (distance[selectedNum - 1] > tao)
	{
		selectedNum /= 2;
	}
	int Row, Col;
	for (int k = 0; k < selectedNum; k++)
	{
		Row = row_min + idx[k] / col_length;
		Col = col_min + idx[k] % col_length;
		tmp = block[Row * bn_c + Col].clone();
		sim_patch.push_back(tmp);
		sim_num.push_back(Row * bn_c + Col);
	}
}

float cal_distance(Mat a, Mat b)
{
	int sy = a.rows;
	int sx = a.cols;
	float sum = 0;
	for (int i = 0; i < sy; i++)
	{
		const float* M1 = a.ptr<float>(i);
		const float* M2 = b.ptr<float>(i);
		for (int j = 0; j < sx; j++)
		{
			sum += (M1[j] - M2[j]) * (M1[j] - M2[j]);
		}
	}
	return sum / (sy * sx);
}

void tran1d(vector<Mat>& input, int patchSize)
{
	int size = input.size();
	int layer = log2(size);
	float* data = new float[size];
	for (int i = 0; i < patchSize; i++)
	{
		for (int j = 0; j < patchSize; j++)
		{
			for (int k = 0; k < size; k++)
			{
				data[k] = input[k].at<float>(i, j);
				//cout << data[k] << endl;
			}
			wavedec(data, size);
			for (int k = 0; k < size; k++)
			{
				input[k].at<float>(i, j) = data[k];
			}
		}
	}
	delete[] data;
}

int log2(const int N)
{
	int k = 1;
	int n = 0;
	while (k < N)
	{
		k *= 2;
		n++;
	}
	return n;
}

void wavedec(float* input, int length)
{
	int N = log2(length);
	if (N == 0)return;
	float* tmp = new float[length];
	for (int i = 0; i < length; i++) {
		tmp[i] = input[i];
	}
	for (int k = 0; k < length / 2; k++)
	{
		input[k] = (tmp[2 * k] + tmp[2 * k + 1]) / sqrt(2);
		input[k + length / 2] = (tmp[2 * k] - tmp[2 * k + 1]) / sqrt(2);
	}
	delete[] tmp;
	wavedec(input, length / 2);
	return;
}

void waverec(float* input, int length, int N)
{
	if (log2(length) > N) return;
	float* tmp = new float[length];
	for (int i = 0; i < length; i++) {
		tmp[i] = input[i];
	}
	for (int k = 0; k < length / 2; k++)
	{
		input[2 * k] = (tmp[k] + tmp[k + length / 2]) / sqrt(2);
		input[2 * k + 1] = (tmp[k] - tmp[k + length / 2]) / sqrt(2);
	}
	delete tmp;
	waverec(input, length * 2, N);
}

void DetectZero(vector<Mat>& input, float threshold)
{
	for (int k = 0; k < input.size(); k++)
	{
		for (int i = 0; i < input[k].rows; i++)
			for (int j = 0; j < input[k].cols; j++)
			{
				if (fabs(input[k].at<float>(i, j)) < threshold)
				{
					input[k].at<float>(i, j) = 0;
				}
			}
	}
}

float calculate_weight_hd(const vector<Mat>input, int sigma)
{
	int num = 0;
	for (int k = 0; k < input.size(); k++)
		for (int i = 0; i < input[k].rows; i++)
			for (int j = 0; j < input[k].cols; j++)
			{
				if (input[k].at<float>(i, j) != 0)
				{
					num++;
				}
			}
	if (num == 0)
	{
		return 1;
	}
	else
		return 1.0 / (sigma * sigma * num);
}

float calculate_weight_wien(const vector<Mat>input, int sigma)
{
	float sum = 0;
	for (int k = 0; k < input.size(); k++)
		for (int i = 0; i < input[k].rows; i++)
			for (int j = 0; j < input[k].cols; j++)
			{
				sum += (input[k].at<float>(i, j)) * (input[k].at<float>(i, j));
			}
	return 1.0 / (sigma * sigma * sum);
}

void Inver3Dtrans(vector<Mat>& input, int patchSize)
{
	int size = input.size();
	int layer = log2(size);
	float* data = new float[size];
	for (int i = 0; i < patchSize; i++)
		for (int j = 0; j < patchSize; j++)
		{
			for (int k = 0; k < size; k++)
			{
				data[k] = input[k].at<float>(i, j);
			}
			waverec(data, 2, layer);
			for (int k = 0; k < size; k++)
			{
				input[k].at<float>(i, j) = data[k];
			}
		}
	for (int k = 0; k < size; k++)
	{
		input[k] = input[k].clone();
		idct(input[k], input[k]);
	}
}

void aggregation(Mat& numerator, Mat& denominator, vector<int>idx_r, vector<int>idx_c,
	const vector<Mat> input, float weight, int patchSize, Mat window)
{
	Rect rect;
	for (int k = 0; k < input.size(); k++)
	{
		rect.x = idx_c[k];
		rect.y = idx_r[k];
		rect.height = patchSize;
		rect.width = patchSize;
		numerator(rect) = numerator(rect) + weight * (input[k].mul(window));
		denominator(rect) = denominator(rect) + weight * window;
	}
}

Mat gen_kaiser(int beta, int length)//How to do this?
{
	if ((beta == 2) && (length == 8))
	{
		Mat window(length, length, CV_32FC1);
		Mat kai1(length, 1, CV_32FC1);
		Mat kai1_T(1, length, CV_32FC1);
		kai1.at<float>(0, 0) = 0.4387;
		kai1.at<float>(1, 0) = 0.6813;
		kai1.at<float>(2, 0) = 0.8768;
		kai1.at<float>(3, 0) = 0.9858;
		for (int i = 0; i < 4; i++)
		{
			kai1.at<float>(7 - i, 0) = kai1.at<float>(i, 0);
			kai1_T.at<float>(0, i) = kai1.at<float>(i, 0);
			kai1_T.at<float>(0, 7 - i) = kai1.at<float>(i, 0);
		}
		window = kai1 * kai1_T;
		return window;
	}
}

void gen_wienFilter(vector<Mat>& input, int sigma)
{
	Mat tmp;
	Mat Sigma(input[0].size(), CV_32FC1, Scalar::all(sigma * sigma));
	for (int k = 0; k < input.size(); k++)
	{
		tmp = input[k].mul(input[k]) + Sigma;
		input[k] = input[k].mul(input[k]) / (tmp.clone());
	}
}

void wienFiltering(vector<Mat>& input, const vector<Mat>wien, int patchSize)
{
	for (int k = 0; k < input.size(); k++)
	{
		input[k] = input[k].mul(wien[k]);
	}
}
