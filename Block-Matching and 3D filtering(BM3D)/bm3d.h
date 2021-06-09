#include <opencv2\highgui\highgui.hpp>
#include <opencv2\opencv.hpp>
#include <opencv2\core\core.hpp>
#include <vector>

using namespace cv;
using namespace std;

void addGuassianNoise(const int sigma,const Mat origin,Mat &noisy);
float cal_psnr(const Mat x, const Mat y);
int log2(const int N);

int runBm3d(const Mat image_noisy,
	Mat &image_basic,Mat &image_denoised);

void GetAllBlock(const Mat img,const int width,const int height,const int channels,
	const int patchSize,const int step,vector<Mat>&block,vector<int>&row_idx,vector<int>&col_idx);

void tran2d( vector<Mat> &input,int patchsize);

void getSimilarPatch(const vector<Mat> block, vector<Mat>& sim_patch, vector<int>& sim_num,
	int i, int j, int bn_r, int bn_c, int area, int maxNum, int tao);

float cal_distance(const Mat a, const Mat b);

void tran1d(vector<Mat>&input,int patchSize);

void DetectZero(vector<Mat>&input, float threshhold);

float calculate_weight_hd(const vector<Mat>input,int sigma);
float calculate_weight_wien(const vector<Mat>input, int sigma);

void Inver3Dtrans(vector<Mat>&input,int patchSize);

void aggregation(Mat &numerator, Mat &denominator, vector<int>idx_r, vector<int>idx_c, const vector<Mat> input,
	float weight,int patchSize,Mat window);

void gen_wienFilter(vector<Mat>&input, int sigma);

void wienFiltering(vector<Mat>&input, const vector<Mat>wien,int patchSize);

Mat gen_kaiser(int beta,int length);
void wavedec(float *input,int length);
void waverec(float* input, int length,int N);

const unsigned nHard = 39;//search area
const unsigned nWien = 39;
const unsigned kHard = 8;//patch size
const unsigned kWien = 8;
const unsigned NHard = 16;//max number
const unsigned NWien = 32;
const unsigned pHard = 3;//step
const unsigned pWien = 3;
int sigma = 25;

const int tao_hard = 2500;
const int tao_wien = 400;

int beta = 2;
float lambda3d = 2.7;
float lambda2d = 0;