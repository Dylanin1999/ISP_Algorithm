/* 
SIFT 库文件
*/

#ifndef SIFT_H
#define SIFT_H

#include "image.h"
#include <list>
#include <vector>

namespace sift {

/****************************************
 *******        SIFT 参数          *******
 ***************************************/

/**************** 金字塔构建参数 *********************/

// 每个Octave中包含的尺度s: 实际计算使用3-5中的值
static int SIFT_INTVLS = 3;

// 高斯模糊Sigma
static float SIFT_SIGMA = 1.6f;
// assumed gaussian blur for input image
static float SIFT_INIT_SIGMA = 0.5f;

// 高斯核半径参数 (2*radius+1)x(2*radius+1), 通常使用2 或 3
static float SIFT_GAUSSIAN_FILTER_RADIUS = 3.0f;

/**************** 梯度方向参数 ******************/

// 关键点方向直方图的Bins数
static int SIFT_ORI_HIST_BINS = 36;

// 梯度方向直方图参数
static float SIFT_ORI_SIG_FCTR = 1.5f;  // sigma
static float SIFT_ORI_RADIUS = 3 * SIFT_ORI_SIG_FCTR;   //邻域窗口半径

// 辅方向阈值: 为了增强匹配的鲁棒性，只保留峰值大于主方向峰值80％的方向作为该关键点的辅方向
static float SIFT_ORI_PEAK_RATIO = 0.8f;

/**************** SIFT关键点筛选参数 *****************/

// 对比度阈值 |D(x)|
static float SIFT_CONTR_THR = 8.0f;

// 主曲率阈值 r
static float SIFT_CURV_THR = 10.0f;

// 子像素下阈值
static float SIFT_KEYPOINT_SUBPiXEL_THR = 0.6f;

// 边界阈值
static int SIFT_IMG_BORDER = 5;

// 最大插值步长
static int SIFT_MAX_INTERP_STEPS = 5;

/**************** SIFT描述子参数 *********************/

// 关键点邻域子区域划分数
static int SIFT_DESCR_WIDTH = 4;

// 每个子区域的梯度划分方向
static int SIFT_DESCR_HIST_BINS = 8;

// 描述子子区域边长参数
static float SIFT_DESCR_SCL_FCTR = 3.f;

// 大梯度值截断阈值
static float SIFT_DESCR_MAG_THR = 0.2f;

// 再归一化处理参数: convert floating-point to unsigned char
static float SIFT_INT_DESCR_FCTR = 512.f;

/****************** SIFT匹配参数 ********************/

// 用于匹配的近邻距离值阈值, 阈值越小匹配越精确, 0.1无匹配, |DR_nearest|/|DR_2nd_nearest|<SIFT_MATCH_NNDR_THR
static float SIFT_MATCH_NNDR_THR = 0.2f; //0.65f


/****************************************
 *******           定义            *******
 ***************************************/

// SIFT关键点: 128维
#define DEGREE_OF_DESCRIPTORS (128)
struct SiftKeypoint {
    int octave;   // octave数量
    int layer;    // layer数量
    float rlayer; // layer实际数量

    float r;     // 归一化的row坐标
    float c;     // 归一化的col坐标
    float scale; // 归一化的scale

    float ri;          // row坐标(layer)
    float ci;          // column坐标(layer)
    float layer_scale; // scale(layer)

    float ori; // 方向(degrees).
    float mag; // 模值

    float descriptors[DEGREE_OF_DESCRIPTORS];
};

// 关键点匹配对
struct MatchPair {
    int r1;
    int c1;
    int r2;
    int c2;
};

/****************************************
 *              基础函数
 ***************************************/

// 高斯模糊
int gaussian_blur(const Image<float> &in, Image<float> &out, std::vector<float> coef1d);

// 高斯一维模糊(水平)&转置图像矩阵
int row_filter_transpose(float *src, float *dst, int w, int h, float *coef1d, int gR);

// 构建Octaves: 下采样
int build_octaves(const Image<unsigned char> &image, std::vector<Image<unsigned char>> &octaves, int nOctaves);

// 生成高斯模糊核
std::vector<std::vector<float>> compute_gaussian_coefs(int nOctaves, int nGpyrLayers);

// 构建高斯金字塔
int build_gaussian_pyramid(std::vector<Image<unsigned char>> &octaves, std::vector<Image<float>> &gpyr, int nOctaves, int nGpyrLayers);

// 构建DoG金字塔
int build_dog_pyr(std::vector<Image<float>> &gpyr, std::vector<Image<float>> &dogPyr, int nOctaves, int nDogLayers);

// 构建梯度金字塔: 用于关键点方向分配
int build_grd_rot_pyr(std::vector<Image<float>> &gpyr, std::vector<Image<float>> &grdPyr, std::vector<Image<float>> &rotPyr, int nOctaves,
                      int nLayers);

// 优化局部极值点
bool refine_local_extrema(std::vector<Image<float>> &dogPyr, int nOctaves, int nDogLayers, SiftKeypoint &kpt);

// Export keypoint list to a file.
int export_kpt_list_to_file(const char *filename, std::list<SiftKeypoint> &kpt_list, bool bIncludeDescpritor);

/****************************************
 *            SIFT 核心函数
 ***************************************/
// SIFT关键点检测
int detect_keypoints(std::vector<Image<float>> &dogPyr, std::vector<Image<float>> &grdPyr, std::vector<Image<float>> &rotPyr, int nOctaves, int nDogLayers, std::list<SiftKeypoint> &kpt_list);

// 生成SIFT描述子
int extract_descriptor(std::vector<Image<float>> &grdPyr, std::vector<Image<float>> &rotPyr, int nOctaves, int nGpyrLayers, std::list<SiftKeypoint> &kpt_list);

/****************************************
 *          SIFT主功能函数
 ***************************************/
// 检测SIFT关键点 & 生成描述子
int sift_cpu(const Image<unsigned char> &image, std::list<SiftKeypoint> &kpt_list, bool bExtractDescriptors);

// 关键点匹配
int match_keypoints(std::list<SiftKeypoint> &kpt_list1, std::list<SiftKeypoint> &kpt_list2, std::list<MatchPair> &match_list);

} // end namespace sift

#endif
