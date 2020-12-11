/*  
    SIFT算法实现
*/

#include "sift.h"
#include "common.h"
#include "image.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <list>

namespace sift {

// 构建Octaves: 下采样
int build_octaves(const Image<unsigned char> &image, std::vector<Image<unsigned char>> &octaves, int nOctaves)
/*
    Parameters:
        image: 输入图片
        octaves: Octaves vector (len: nOctaves)
        nOctaves: Octaves 数量
*/
{
    for (int i = 0; i < nOctaves; i++) {
        if (i == 0) {octaves[i] = image;}
        else {octaves[i] = octaves[(i - 1)].downsample_2x();}
    }
    return 0;
}

// 高斯模糊
int gaussian_blur(const Image<float> &in, Image<float> &out, std::vector<float> coef1d)
/*
    Parameters:
        in: 输入图片
        out: 输出图片
        coef1d: 一维高斯模糊参数
*/
{
    int w = in.w;
    int h = in.h;
    int gR = static_cast<int>(coef1d.size()) / 2;

    Image<float> img_t(h, w);
    // 图像二维高斯模糊: 连续两次高斯一维模糊(水平)&转置图像矩阵实现
    row_filter_transpose(in.data, img_t.data, w, h, &coef1d[0], gR);
    row_filter_transpose(img_t.data, out.data, h, w, &coef1d[0], gR);

    return 0;
}

// 高斯一维模糊(水平)&转置图像矩阵
int row_filter_transpose(float *src, float *dst, int w, int h, float *coef1d, int gR)
{
    float *row_buf = new float[w + gR * 2];
    float *row_start;
    int elemSize = sizeof(float);

    float *srcData = src;
    float *dstData = dst + w * h - 1;
    float partialSum = 0.0f;
    float *coef = coef1d;
    float *prow;

    float firstData, lastData;
    for (int r = 0; r < h; r++) {
        row_start = srcData + r * w;
        memcpy(row_buf + gR, row_start, elemSize * w);
        firstData = *(row_start);
        lastData = *(row_start + w - 1);
        for (int i = 0; i < gR; i++) {
            row_buf[i] = firstData;
            row_buf[i + w + gR] = lastData;
        }

        prow = row_buf;
        dstData = dstData - w * h + 1;
        for (int c = 0; c < w; c++) {
            partialSum = 0.0f;
            coef = coef1d;

            for (int i = -gR; i <= gR; i++) {
                partialSum += (*coef++) * (*prow++);
            }

            prow -= 2 * gR;
            *dstData = partialSum;
            dstData += h;
        }
    }
    delete[] row_buf;
    row_buf = nullptr;

    return 0;
}

// 构建高斯金字塔
int build_gaussian_pyramid(std::vector<Image<unsigned char>> &octaves, std::vector<Image<float>> &gpyr, int nOctaves, int nGpyrLayers)
/*
    Parameters:
        octaves: 构建好的Octaves
        gpyr: 高斯金字塔Vector (len: nOctaves * nGpyrLayers)
        nOctaves: Octaves数量
        nGpyrLayers: 每个Octave中高斯图像数量
*/
{
    // 每个Octave中包含的尺度s
    int nLayers = nGpyrLayers - 3;
    // 生成高斯核
    std::vector<std::vector<float>> gaussian_coefs = compute_gaussian_coefs(nOctaves, nGpyrLayers);

    int w, h;
    for (int i = 0; i < nOctaves; i++) {
        w = octaves[i].w;
        h = octaves[i].h;
        for (int j = 0; j < nGpyrLayers; j++) {
            if (i == 0 && j == 0) {
                // 首个高斯图像生成
                gpyr[0].init(w, h);
                gaussian_blur(octaves[0].to_float(), gpyr[0], gaussian_coefs[j]);
            }
            else if (i > 0 && j == 0) {
                // 降采样时，高斯金字塔上一组图像的初始图像(底层图像)是由前一组图像的倒数第三张图像隔点采样得到的
                gpyr[i * nGpyrLayers] = gpyr[(i - 1) * nGpyrLayers + nLayers].downsample_2x();
            }
            else {
                // Octave内的高斯图像由上一个高斯图像模糊得到
                gpyr[i * nGpyrLayers + j].init(w, h);
                gaussian_blur(gpyr[i * nGpyrLayers + j - 1], gpyr[i * nGpyrLayers + j], gaussian_coefs[j]);
            }
        }
    }
    // 高斯金字塔构建完毕, Octaves不再使用, 释放内存
    octaves.clear();
    return 0;
}

// 生成高斯模糊核
std::vector<std::vector<float>> compute_gaussian_coefs(int nOctaves, int nGpyrLayers)
/*
    Parameters:
        nOctaves: Octaves数量
        nGpyrLayers: 每个Octave中高斯图像数量
*/
{
    // 计算每个高斯图对应的Sigma
    int nLayers = nGpyrLayers - 3;
    float sigma, sigma_pre;
    float sigma0 = SIFT_SIGMA;
    float k = powf(2.0f, 1.0f / nLayers);

    std::vector<float> sig(nGpyrLayers);
    sigma_pre = SIFT_INIT_SIGMA;
    sig[0] = sqrtf(sigma0 * sigma0 - sigma_pre * sigma_pre);
    // 计算Octave组内每层的尺度坐标
    for (int i = 1; i < nGpyrLayers; i++) {
        sigma_pre = powf(k, (float)(i - 1)) * sigma0;
        sigma = sigma_pre * k;
        sig[i] = sqrtf(sigma * sigma - sigma_pre * sigma_pre);
    }

    std::vector<std::vector<float>> gaussian_coefs(nGpyrLayers);
    for (int i = 0; i < nGpyrLayers; i++) {
        // 计算核半径
        float factor = SIFT_GAUSSIAN_FILTER_RADIUS;
        int gR = (sig[i] * factor > 1.0f) ? (int)ceilf(sig[i] * factor) : 1;
        int gW = gR * 2 + 1;
        // 生成高斯核
        gaussian_coefs[i].resize(gW);
        float sum = 0.0f;
        float tmp;
        for (int j = 0; j < gW; j++) {
            tmp = (float)((j - gR) / sig[i]);
            gaussian_coefs[i][j] = expf(tmp * tmp * -0.5f) * (1 + j / 1000.0f);
            sum += gaussian_coefs[i][j];
        }
        // 高斯核归一化
        for (int j = 0; j < gW; j++) {
            gaussian_coefs[i][j] = gaussian_coefs[i][j] / sum;
        }
    }
    return gaussian_coefs;
}

// 构建DoG金字塔
int build_dog_pyr(std::vector<Image<float>> &gpyr, std::vector<Image<float>> &dogPyr, int nOctaves, int nDogLayers)
/*
    Parameters:
        gpyr: 构建好的高斯金字塔vector
        dogPyr: DoG金字塔vector (len: nOctaves * nDogLayers)
        nOctaves: Octaves数量
        nDogLayers: 每个Octave中的DoG数量
*/
{
    int nGpyrLayers = nDogLayers + 1;

    int w, h;
    float *srcData1;
    float *srcData2;
    float *dstData;
    int index = 0;

    for (int i = 0; i < nOctaves; i++) {
        // 每个Octave的底层高斯图信息
        int row_start = i * nGpyrLayers;
        w = gpyr[row_start].w;
        h = gpyr[row_start].h;

        for (int j = 0; j < nDogLayers; j++) {
            dogPyr[i * nDogLayers + j].init(w, h);
            dstData = dogPyr[i * nDogLayers + j].data;

            srcData1 = gpyr[row_start + j].data;
            srcData2 = gpyr[row_start + j + 1].data;

            // 赋 差值 给对应的DoG位置
            index = 0;
            while (index++ < w * h)
                *(dstData++) = *(srcData2++) - *(srcData1++);
        }
    }

    return 0;
}

// 构建梯度金字塔: 用于关键点方向分配
int build_grd_rot_pyr(std::vector<Image<float>> &gpyr, std::vector<Image<float>> &grdPyr, std::vector<Image<float>> &rotPyr, int nOctaves, int nLayers)
/*
    Parameters:
        gpyr: 构建好的高斯金字塔
        grdPyr: 梯度模值金字塔 (len: nOctaves * nGpyrLayers)
        rotPyr: 梯度方向金字塔 (len: nOctaves * nGpyrLayers)
        nOctaves: Octaves数量
        nLayers: 每个Octave中包含的尺度s
*/
{
    int nGpyrLayers = nLayers + 3;
    int w, h;
    float dr, dc;
    float angle;

    float *srcData;
    float *grdData;
    float *rotData;

    for (int i = 0; i < nOctaves; i++) {

        w = gpyr[i * nGpyrLayers].w;
        h = gpyr[i * nGpyrLayers].h;
        for (int j = 1; j <= nLayers; j++) {
            int layer_index = i * nGpyrLayers + j;
            grdPyr[layer_index].init(w, h);
            rotPyr[layer_index].init(w, h);

            srcData = gpyr[layer_index].data;
            grdData = grdPyr[layer_index].data;
            rotData = rotPyr[layer_index].data;

            for (int r = 0; r < h; r++) {
                for (int c = 0; c < w; c++) {
                    dr = get_pixel_f(srcData, w, h, r + 1, c) - get_pixel_f(srcData, w, h, r - 1, c);
                    dc = get_pixel_f(srcData, w, h, r, c + 1) - get_pixel_f(srcData, w, h, r, c - 1);

                    // 快速计算梯度模值
                    grdData[r * w + c] = fast_sqrt_f(dr * dr + dc * dc);
                    // 快速计算梯度方向
                    angle = fast_atan2_f(dr, dc);

                    rotData[r * w + c] = angle;
                }
            }
        }
    }
    return 0;
}

// 生成方向直方图
float compute_orientation_hist_with_gradient(const Image<float> &grdImage, const Image<float> &rotImage, SiftKeypoint &kpt, float *&hist)
/*
    Parameters:
        grdImage: 梯度图
        rotImage: 方向图
        kpt: 极值点
        hist: 直方图
*/
{
    int nBins = SIFT_ORI_HIST_BINS;

    float kptr = kpt.ri;
    float kptc = kpt.ci;
    float kpt_scale = kpt.layer_scale;

    int kptr_i = (int)(kptr + 0.5f);
    int kptc_i = (int)(kptc + 0.5f);
    float d_kptr = kptr - kptr_i;
    float d_kptc = kptc - kptc_i;

    // 邻域参数
    float sigma = SIFT_ORI_SIG_FCTR * kpt_scale;
    int win_radius = (int)(SIFT_ORI_RADIUS * kpt_scale);
    float exp_factor = -1.0f / (2.0f * sigma * sigma);

    float *grdData = grdImage.data;
    float *rotData = rotImage.data;
    int w = grdImage.w;
    int h = grdImage.h;

    int r, c;
    float magni, angle, weight;
    int bin;
    float fbin;

    float *tmpHist = new float[nBins];
    memset(tmpHist, 0, nBins * sizeof(float));

    for (int i = -win_radius; i <= win_radius; i++) // rows
    {
        r = kptr_i + i;
        if (r <= 0 || r >= h - 1) // Cannot calculate dy
            continue;
        for (int j = -win_radius; j <= win_radius; j++) // columns
        {
            c = kptc_i + j;
            if (c <= 0 || c >= w - 1)
                continue;

            magni = grdData[r * w + c];
            angle = rotData[r * w + c];

            fbin = angle * nBins / _2PI;
            weight = expf(((i - d_kptr) * (i - d_kptr) + (j - d_kptc) * (j - d_kptc)) * exp_factor);

            // 梯度方向插值
            bin = (int)(fbin - 0.5f);
            float d_fbin = fbin - 0.5f - bin;

            float mw = weight * magni;
            float dmw = d_fbin * mw;
            tmpHist[(bin + nBins) % nBins] += mw - dmw;
            tmpHist[(bin + 1) % nBins] += dmw;
        }
    }

#define TMPHIST(idx) (idx < 0 ? tmpHist[0] : (idx >= nBins ? tmpHist[nBins - 1] : tmpHist[idx]))

    // Smooth the histogram. Algorithm comes from OpenCV.
    hist[0] = (tmpHist[0] + tmpHist[2]) * 1.0f / 16.0f +
              (tmpHist[0] + tmpHist[1]) * 4.0f / 16.0f +
              tmpHist[0] * 6.0f / 16.0f;
    hist[1] = (tmpHist[0] + tmpHist[3]) * 1.0f / 16.0f +
              (tmpHist[0] + tmpHist[2]) * 4.0f / 16.0f +
              tmpHist[1] * 6.0f / 16.0f;
    hist[nBins - 2] = (tmpHist[nBins - 4] + tmpHist[nBins - 1]) * 1.0f / 16.0f +
                      (tmpHist[nBins - 3] + tmpHist[nBins - 1]) * 4.0f / 16.0f +
                      tmpHist[nBins - 2] * 6.0f / 16.0f;
    hist[nBins - 1] = (tmpHist[nBins - 3] + tmpHist[nBins - 1]) * 1.0f / 16.0f +
                      (tmpHist[nBins - 2] + tmpHist[nBins - 1]) * 4.0f / 16.0f +
                      tmpHist[nBins - 1] * 6.0f / 16.0f;

    for (int i = 2; i < nBins - 2; i++) {
        hist[i] = (tmpHist[i - 2] + tmpHist[i + 2]) * 1.0f / 16.0f +
                  (tmpHist[i - 1] + tmpHist[i + 1]) * 4.0f / 16.0f +
                  tmpHist[i] * 6.0f / 16.0f;
    }

    // 找直方图峰值
    float maxitem = hist[0];
    int max_i = 0;
    for (int i = 0; i < nBins; i++) {
        if (maxitem < hist[i]) {
            maxitem = hist[i];
            max_i = i;
        }
    }

    kpt.ori = max_i * _2PI / nBins;

    delete[] tmpHist;
    tmpHist = nullptr;
    return maxitem;
}

// SIFT关键点检测
int detect_keypoints(std::vector<Image<float>> &dogPyr, std::vector<Image<float>> &grdPyr, std::vector<Image<float>> &rotPyr, int nOctaves, int nDogLayers, std::list<SiftKeypoint> &kpt_list)
/*
    Parameters:
        dogPyr: 构建好的DoG金字塔
        grdPyr: 构建好的梯度模值金字塔
        rotPyr: 构建好的梯度方向金字塔
        nOctaves: Octave数量
        nDogLayers: 每个Octave中DoG的数量
        kpt_list: SIFT关键点列表
*/
{

    float *currData; // 当前DoG数据
    float *lowData;  // 前一个DoG数据
    float *highData; // 后一个DoG数据

    SiftKeypoint kpt;

    int w, h;
    int layer_index;
    int index;
    float val;

    int nBins = SIFT_ORI_HIST_BINS; // 关键点方向直方图的Bins数 36
    float *hist = new float[nBins];
    int nGpyrLayers = nDogLayers + 1;

    // 对比度阈值
    float contr_thr = 0.8f * SIFT_CONTR_THR; // In Lowe's paper, |D(x)|<0.03 will be rejected.

    for (int i = 0; i < nOctaves; i++) {
        w = dogPyr[i * nDogLayers].w;
        h = dogPyr[i * nDogLayers].h;

        for (int j = 1; j < nDogLayers - 1; j++) {
            layer_index = i * nDogLayers + j;

            highData = dogPyr[layer_index + 1].data;
            currData = dogPyr[layer_index].data;
            lowData = dogPyr[layer_index - 1].data;

            for (int r = SIFT_IMG_BORDER; r < h - SIFT_IMG_BORDER; r++) {
                for (int c = SIFT_IMG_BORDER; c < w - SIFT_IMG_BORDER; c++) {
                    index = r * w + c;
                    val = currData[index];

                    // 局部极值bool: 对比度阈值比较 & 27个邻域比较
                    bool bExtrema =
                        //局部极大值判断
                        (val >= contr_thr && 
                         val > highData[index - w - 1] &&
                         val > highData[index - w] &&
                         val > highData[index - w + 1] &&

                         val > highData[index - 1] && 
                         val > highData[index] &&
                         val > highData[index + 1] &&

                         val > highData[index + w - 1] &&
                         val > highData[index + w] &&
                         val > highData[index + w + 1] &&

                         val > currData[index - w - 1] &&
                         val > currData[index - w] &&
                         val > currData[index - w + 1] &&

                         val > currData[index - 1] &&
                         val > currData[index + 1] &&

                         val > currData[index + w - 1] &&
                         val > currData[index + w] &&
                         val > currData[index + w + 1] &&

                         val > lowData[index - w - 1] &&
                         val > lowData[index - w] &&
                         val > lowData[index - w + 1] &&

                         val > lowData[index - 1] && 
                         val > lowData[index] &&
                         val > lowData[index + 1] &&

                         val > lowData[index + w - 1] &&
                         val > lowData[index + w] &&
                         val > lowData[index + w + 1]) || 
                         // 局部极小值判断
                        (val <= -contr_thr && 
                         val < highData[index - w - 1] &&
                         val < highData[index - w] &&
                         val < highData[index - w + 1] &&

                         val < highData[index - 1] && 
                         val < highData[index] &&
                         val < highData[index + 1] &&

                         val < highData[index + w - 1] &&
                         val < highData[index + w] &&
                         val < highData[index + w + 1] &&

                         val < currData[index - w - 1] &&
                         val < currData[index - w] &&
                         val < currData[index - w + 1] &&

                         val < currData[index - 1] &&
                         val < currData[index + 1] &&

                         val < currData[index + w - 1] &&
                         val < currData[index + w] &&
                         val < currData[index + w + 1] &&

                         val < lowData[index - w - 1] &&
                         val < lowData[index - w] &&
                         val < lowData[index - w + 1] &&

                         val < lowData[index - 1] && 
                         val < lowData[index] &&
                         val < lowData[index + 1] &&

                         val < lowData[index + w - 1] &&
                         val < lowData[index + w] &&
                         val < lowData[index + w + 1]);

                    // 极值点筛选
                    if (bExtrema) {
                        kpt.octave = i;
                        kpt.layer = j;
                        kpt.ri = (float)r;
                        kpt.ci = (float)c;
                        // 优化局部极值点: 若不是则丢弃, 寻找下一个极值点
                        bool bGoodKeypoint = refine_local_extrema(dogPyr, nOctaves, nDogLayers, kpt);
                        if (!bGoodKeypoint)
                            continue;
                        // 计算梯度方向直方图峰值
                        float max_mag = compute_orientation_hist_with_gradient(
                            grdPyr[i * nGpyrLayers + kpt.layer],
                            rotPyr[i * nGpyrLayers + kpt.layer], kpt, hist);
                        
                        // 辅方向
                        float threshold = max_mag * SIFT_ORI_PEAK_RATIO;

                        for (int ii = 0; ii < nBins; ii++) {

                            // Use 3 points to fit a curve and find the accurate
                            // location of a keypoints
                            int left = ii > 0 ? ii - 1 : nBins - 1;
                            int right = ii < (nBins - 1) ? ii + 1 : 0;
                            float currHist = hist[ii];
                            float lhist = hist[left];
                            float rhist = hist[right];
                            if (currHist > lhist && currHist > rhist &&
                                currHist > threshold) {
                                // Refer to here:
                                // http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
                                float accu_ii =
                                    ii + 0.5f * (lhist - rhist) /
                                             (lhist - 2.0f * currHist + rhist);

                                // Since bin index means the starting point of a
                                // bin, so the real orientation should be bin
                                // index plus 0.5. for example, angles in bin 0
                                // should have a mean value of 5 instead of 0;
                                accu_ii += 0.5f;
                                accu_ii = accu_ii < 0 ? (accu_ii + nBins)
                                                      : accu_ii >= nBins
                                                            ? (accu_ii - nBins)
                                                            : accu_ii;
                                // The magnitude should also calculate the max
                                // number based on fitting But since we didn't
                                // actually use it in image matching, we just
                                // lazily use the histogram value.
                                kpt.mag = currHist;
                                kpt.ori = accu_ii * _2PI / nBins;
                                kpt_list.push_back(kpt);
                            }
                        }
                    }
                }
            }
        }
    }

    delete[] hist;
    hist = nullptr;
    return 0;
}

// 优化局部极值点
bool refine_local_extrema(std::vector<Image<float>> &dogPyr, int nOctaves, int nDogLayers, SiftKeypoint &kpt)
/*
    Parameters:
        dogPyr: 构建好的DoG金字塔
        nOctaves: Octaves数量
        nDogLayers: 每个Octave中的DoG图数量
        kpt: 待筛选的极值点
*/
{
    int nGpyrLayers = nDogLayers + 1;

    int w = 0, h = 0;
    int layer_idx = 0;
    int octave = kpt.octave;
    int layer = kpt.layer;
    int r = (int)kpt.ri;
    int c = (int)kpt.ci;

    float *currData = nullptr;
    float *lowData = nullptr;
    float *highData = nullptr;

    int xs_i = 0, xr_i = 0, xc_i = 0;
    float tmp_r = 0.0f, tmp_c = 0.0f, tmp_layer = 0.0f;
    float xr = 0.0f, xc = 0.0f, xs = 0.0f;
    float x_hat[3] = {xc, xr, xs};
    float dx = 0.0f, dy = 0.0f, ds = 0.0f;
    float dxx = 0.0f, dyy = 0.0f, dss = 0.0f, dxs = 0.0f, dys = 0.0f,
          dxy = 0.0f;

    tmp_r = (float)r;
    tmp_c = (float)c;
    tmp_layer = (float)layer;

    // 关键点的精确定位: 利用已知的离散空间点插值得到的连续空间极值点的方法叫做子像素插值
    int i = 0;
    for (; i < SIFT_MAX_INTERP_STEPS; i++) {
        c += xc_i;
        r += xr_i;

        layer_idx = octave * nDogLayers + layer;
        w = dogPyr[layer_idx].w;
        h = dogPyr[layer_idx].h;
        currData = dogPyr[layer_idx].data;
        lowData = dogPyr[layer_idx - 1].data;
        highData = dogPyr[layer_idx + 1].data;
        
        // 一阶差分
        dx = (get_pixel_f(currData, w, h, r, c + 1) - get_pixel_f(currData, w, h, r, c - 1)) * 0.5f;
        dy = (get_pixel_f(currData, w, h, r + 1, c) - get_pixel_f(currData, w, h, r - 1, c)) * 0.5f;
        ds = (get_pixel_f(highData, w, h, r, c) - get_pixel_f(lowData, w, h, r, c)) * 0.5f;
        float dD[3] = {-dx, -dy, -ds};

        // 二阶差分
        float v2 = 2.0f * get_pixel_f(currData, w, h, r, c);
        dxx = (get_pixel_f(currData, w, h, r, c + 1) + get_pixel_f(currData, w, h, r, c - 1) - v2);
        dyy = (get_pixel_f(currData, w, h, r + 1, c) + get_pixel_f(currData, w, h, r - 1, c) - v2);
        dss = (get_pixel_f(highData, w, h, r, c) + get_pixel_f(lowData, w, h, r, c) - v2);
        dxy = (get_pixel_f(currData, w, h, r + 1, c + 1) - get_pixel_f(currData, w, h, r + 1, c - 1) -
               get_pixel_f(currData, w, h, r - 1, c + 1) + get_pixel_f(currData, w, h, r - 1, c - 1)) * 0.25f;
        dxs = (get_pixel_f(highData, w, h, r, c + 1) - get_pixel_f(highData, w, h, r, c - 1) -
               get_pixel_f(lowData, w, h, r, c + 1) + get_pixel_f(lowData, w, h, r, c - 1)) * 0.25f;
        dys = (get_pixel_f(highData, w, h, r + 1, c) - get_pixel_f(highData, w, h, r - 1, c) -
               get_pixel_f(lowData, w, h, r + 1, c) + get_pixel_f(lowData, w, h, r - 1, c)) * 0.25f;

        float H[3][3] = {{dxx, dxy, dxs}, {dxy, dyy, dys}, {dxs, dys, dss}};
        float Hinvert[3][3];
        float det;

        // 三阶矩阵求逆 INVERT_3X3 = DETERMINANT_3X3, then SCALE_ADJOINT_3X3;
        float tmp;
        DETERMINANT_3X3(det, H);
        if (fabsf(det) < (std::numeric_limits<float>::min)())
            break;
        tmp = 1.0f / (det);
        SCALE_ADJOINT_3X3(Hinvert, tmp, H);
        MAT_DOT_VEC_3X3(x_hat, Hinvert, dD);

        xs = x_hat[2];
        xr = x_hat[1];
        xc = x_hat[0];

        // 更新temp
        tmp_r = r + xr;
        tmp_c = c + xc;
        tmp_layer = layer + xs;

        // 循环边界检查
        xc_i = ((xc >= SIFT_KEYPOINT_SUBPiXEL_THR && c < w - 2) ? 1 : 0) + ((xc <= -SIFT_KEYPOINT_SUBPiXEL_THR && c > 1) ? -1 : 0);

        xr_i = ((xr >= SIFT_KEYPOINT_SUBPiXEL_THR && r < h - 2) ? 1 : 0) + ((xr <= -SIFT_KEYPOINT_SUBPiXEL_THR && r > 1) ? -1 : 0);

        if (xc_i == 0 && xr_i == 0 && xs_i == 0)
            break;
    }

    if (i >= SIFT_MAX_INTERP_STEPS)
        return false;
    if (fabsf(xc) >= 1.5 || fabsf(xr) >= 1.5 || fabsf(xs) >= 1.5)
        return false;

    // 如果(r, c, layer)越界, 返回False.
    if (tmp_layer < 0 || tmp_layer > (nGpyrLayers - 1) || tmp_r < 0 || tmp_r > h - 1 || tmp_c < 0 || tmp_c > w - 1)
        return false;

    {
        float value = get_pixel_f(currData, w, h, r, c) + 0.5f * (dx * xc + dy * xr + ds * xs);
        if (fabsf(value) < SIFT_CONTR_THR)
            return false;

        // 消除边缘响应:　一个定义不好的高斯差分算子的极值在横跨边缘的地方有较大的主曲率，而在垂直边缘的方向有较小的主曲率
        float trH = dxx + dyy;
        float detH = dxx * dyy - dxy * dxy;
        float response = (SIFT_CURV_THR + 1) * (SIFT_CURV_THR + 1) / (SIFT_CURV_THR);

        if (detH <= 0 || (trH * trH / detH) >= response)
            return false;
    }

    // 当前层中的坐标
    kpt.ci = tmp_c;
    kpt.ri = tmp_r;
    kpt.layer_scale = SIFT_SIGMA * powf(2.0f, tmp_layer / SIFT_INTVLS);


    float norm = powf(2.0f, (float)(octave));
    // 归一化的坐标.
    kpt.c = tmp_c * norm;
    kpt.r = tmp_r * norm;
    kpt.rlayer = tmp_layer;
    kpt.layer = layer;

    // Scale = sigma0 * 2^octave * 2^(layer/S);
    kpt.scale = kpt.layer_scale * norm;

    return true;
}

// 生成SIFT描述子
int extract_descriptor(std::vector<Image<float>> &grdPyr, std::vector<Image<float>> &rotPyr, int nOctaves, int nGpyrLayers, std::list<SiftKeypoint> &kpt_list)
/*
    Parameters:
        grdPyr: 梯度模值金字塔
        rotPyr: 梯度方向金字塔
        nOctaves: Octave数量
        nGpyrLayers: 每个Octave中高斯图个数
        kpt_list: SIFT特征点列表
*/
{
    // 将关键点邻域划分为4*4个子区域
    int nSubregion = SIFT_DESCR_WIDTH;
    int nHalfSubregion = nSubregion >> 1;

    // 每个子区域的梯度划分为8个方向
    int nBinsPerSubregion = SIFT_DESCR_HIST_BINS;
    float nBinsPerSubregionPerDegree = (float)nBinsPerSubregion / _2PI;

    // (rbin, cbin, obin):(row of hist bin, column of hist bin, orientation bin)

    int nBins = nSubregion * nSubregion * nBinsPerSubregion;    //4*4*8=128维向量
    int nHistBins = (nSubregion + 2) * (nSubregion + 2) * (nBinsPerSubregion + 2);
    int nSliceStep = (nSubregion + 2) * (nBinsPerSubregion + 2);
    int nRowStep = (nBinsPerSubregion + 2);
    float *histBin = new float[nHistBins];

    for (std::list<SiftKeypoint>::iterator kpt = kpt_list.begin();
         kpt != kpt_list.end(); kpt++) {
        // 关键点信息
        int octave = kpt->octave;
        int layer = kpt->layer;

        float kpt_ori = kpt->ori;
        float kptr = kpt->ri;
        float kptc = kpt->ci;
        float kpt_scale = kpt->layer_scale;

        // 关键点的最近邻坐标
        int kptr_i = (int)(kptr + 0.5f);
        int kptc_i = (int)(kptc + 0.5f);
        float d_kptr = kptr_i - kptr;
        float d_kptc = kptc_i - kptc;

        int layer_index = octave * nGpyrLayers + layer;
        int w = grdPyr[layer_index].w;
        int h = grdPyr[layer_index].h;
        float *grdData = grdPyr[layer_index].data;
        float *rotData = rotPyr[layer_index].data;

        // 确定计算描述子所需的图像区域
        // 子区域边长
        float subregion_width = SIFT_DESCR_SCL_FCTR * kpt_scale;
        // 图像区域半径
        int win_size = (int)(SQRT2 * subregion_width * (nSubregion + 1) * 0.5f + 0.5f);

        // 归一化的 cos() 和 sin()
        float sin_t = sinf(kpt_ori) / (float)subregion_width;
        float cos_t = cosf(kpt_ori) / (float)subregion_width;

        // Re-init histBin
        memset(histBin, 0, nHistBins * sizeof(float));

        // 计算区域内直方图
        float rr, cc;
        float mag, angle, gaussian_weight;


        float rrotate, crotate;
        float rbin, cbin, obin;
        float d_rbin, d_cbin, d_obin;

        // 区域边界
        int r, c;
        int left = MAX(-win_size, 1 - kptc_i);
        int right = MIN(win_size, w - 2 - kptc_i);
        int top = MAX(-win_size, 1 - kptr_i);
        int bottom = MIN(win_size, h - 2 - kptr_i);

        for (int i = top; i <= bottom; i++)
        {
            for (int j = left; j <= right; j++)
            {
                // 基于(kptr, kptc)的准确位置
                rr = i + d_kptr;
                cc = j + d_kptc;

                // 将坐标轴旋转为关键点的方向，以确保旋转不变性. (i, j)旋转后邻域内采样点的新坐标
                rrotate = (cos_t * cc + sin_t * rr);
                crotate = (-sin_t * cc + cos_t * rr);

                // 对于4*4的bins, 它的真实中心在(1.5, 1.5)
                rbin = rrotate + nHalfSubregion - 0.5f;
                cbin = crotate + nHalfSubregion - 0.5f;

                // rbin和cbin的区间是(-1, d), 区间外的均已被处理过
                if (rbin <= -1 || rbin >= nSubregion || cbin <= -1 || cbin >= nSubregion)
                    continue;

                // 
                r = kptr_i + i;
                c = kptc_i + j;
                mag = grdData[r * w + c];
                angle = rotData[r * w + c] - kpt_ori;
                float angle1 = (angle < 0) ? (_2PI + angle) : angle; // 将角度调整至区间[0, 2PI)
                obin = angle1 * nBinsPerSubregionPerDegree;

                int x0, y0, z0;
                int x1, y1, z1;
                y0 = (int)floor(rbin);
                x0 = (int)floor(cbin);
                z0 = (int)floor(obin);
                d_rbin = rbin - y0;
                d_cbin = cbin - x0;
                d_obin = obin - z0;
                x1 = x0 + 1;
                y1 = y0 + 1;
                z1 = z0 + 1;

                // 高斯加权计算: Lowe建议子区域的像素的梯度大小高斯加权计算
                float exp_scale = -2.0f / (nSubregion * nSubregion);
                gaussian_weight = expf((rrotate * rrotate + crotate * crotate) * exp_scale);
                float gm = mag * gaussian_weight;

                // 插值计算每个种子点八个方向的梯度
                // Tri-linear interpolation
                float vr1, vr0;
                float vrc11, vrc10, vrc01, vrc00;
                float vrco110, vrco111, vrco100, vrco101, vrco010, vrco011,
                    vrco000, vrco001;

                vr1 = gm * d_rbin;
                vr0 = gm - vr1;
                vrc11 = vr1 * d_cbin;
                vrc10 = vr1 - vrc11;
                vrc01 = vr0 * d_cbin;
                vrc00 = vr0 - vrc01;
                vrco111 = vrc11 * d_obin;
                vrco110 = vrc11 - vrco111;
                vrco101 = vrc10 * d_obin;
                vrco100 = vrc10 - vrco101;
                vrco011 = vrc01 * d_obin;
                vrco010 = vrc01 - vrco011;
                vrco001 = vrc00 * d_obin;
                vrco000 = vrc00 - vrco001;

                // int idx =  y0  * nSliceStep + x0  * nRowStep + z0;
                // All coords are offseted by 1. so x=[1, 4], y=[1, 4];
                // data for -1 coord is stored at position 0;
                // data for 8 coord is stored at position 9.
                // z doesn't need to move.
                int idx = y1 * nSliceStep + x1 * nRowStep + z0;
                histBin[idx] += vrco000;

                idx++;
                histBin[idx] += vrco001;

                idx += nRowStep - 1;
                histBin[idx] += vrco010;

                idx++;
                histBin[idx] += vrco011;

                idx += nSliceStep - nRowStep - 1;
                histBin[idx] += vrco100;

                idx++;
                histBin[idx] += vrco101;

                idx += nRowStep - 1;
                histBin[idx] += vrco110;

                idx++;
                histBin[idx] += vrco111;
            }
        }

        // Discard all the edges for row and column.
        // Only retrive edges for orientation bins.
        float *dstBins = new float[nBins];
        for (int i = 1; i <= nSubregion; i++) // slice
        {
            for (int j = 1; j <= nSubregion; j++) // row
            {
                int idx = i * nSliceStep + j * nRowStep;
                // comments: how this line works.
                // Suppose you want to write w=width, y=1, due to circular
                // buffer, we should write it to w=0, y=1; since we use a
                // circular buffer, it is written into w=width, y=1. Now, we
                // fectch the data back.
                histBin[idx] = histBin[idx + nBinsPerSubregion];

                // comments: how this line works.
                // Suppose you want to write x=-1 y=1, due to circular, it
                // should be at y=1, x=width-1; since we use circular buffer,
                // the value goes to y=0, x=width, now, we need to get it back.
                if (idx != 0)
                    histBin[idx + nBinsPerSubregion + 1] = histBin[idx - 1];

                int idx1 = ((i - 1) * nSubregion + j - 1) * nBinsPerSubregion;
                for (int k = 0; k < nBinsPerSubregion; k++) {
                    dstBins[idx1 + k] = histBin[idx + k];
                }
            }
        }

        // 特征向量归一化处理: 特征向量形成后，为了去除光照变化的影响，需要对它们进行归一化处理
        float sum_square = 0.0f;
        for (int i = 0; i < nBins; i++)
            sum_square += dstBins[i] * dstBins[i];

        float thr = fast_sqrt_f(sum_square) * SIFT_DESCR_MAG_THR;

        float tmp = 0.0;
        sum_square = 0.0;

        // 描述子向量门限: 非线性光照，相机饱和度变化对造成某些方向的梯度值过大，而对方向的影响微弱
        // 大梯度值截断: 因此设置门限值(向量归一化后，一般取0.2)截断较大的梯度值
        // float thr = fast_sqrt_f(sum_square) * SIFT_DESCR_MAG_THR;
        for (int i = 0; i < nBins; i++) {
            tmp = fmin(thr, dstBins[i]);
            dstBins[i] = tmp;
            sum_square += tmp * tmp;
        }

        // 再归一化处理: 再次归一化处理，提高特征的鉴别性。同时由于过小的数值难以存储, 因而采用固定因子放大
        float norm_factor = SIFT_INT_DESCR_FCTR / fast_sqrt_f(sum_square);
        for (int i = 0; i < nBins; i++)
            dstBins[i] = dstBins[i] * norm_factor;

        memcpy(kpt->descriptors, dstBins, nBins * sizeof(float));

        if (dstBins) {
            delete[] dstBins;
            dstBins = nullptr;
        }
    }
    if (histBin) {
        delete[] histBin;
        histBin = nullptr;
    }

    return 0;
}

int sift_cpu(const Image<unsigned char> &image, std::list<SiftKeypoint> &kpt_list, bool bExtractDescriptors)
{
    // 每个Octave中包含的尺度s.
    int nLayers = SIFT_INTVLS;
    // 每个Octave中DoG数量: 为了在每组中检测S个尺度的极值点，DOG金字塔每组需S+2层图像
    int nDogLayers = nLayers + 2;
    // 每个Octave中高斯图像数量: DOG金字塔由高斯金字塔相邻两层相减得到，则高斯金字塔每组需S+3层图像
    int nGpyrLayers = nLayers + 3;
    // 金字塔Octave数量(层数): 层数根据图像的原始大小和塔顶图像的大小共同决定
    int nOctaves = (int)my_log2((float)fmin(image.w, image.h)) - 3;

    // 构建Octaves
    std::vector<Image<unsigned char>> octaves(nOctaves);
    build_octaves(image, octaves, nOctaves);

    // 构建高斯金字塔
    std::vector<Image<float>> gpyr(nOctaves * nGpyrLayers);
    build_gaussian_pyramid(octaves, gpyr, nOctaves, nGpyrLayers);

    // 构建DoG金字塔
    std::vector<Image<float>> dogPyr(nOctaves * nDogLayers);
    build_dog_pyr(gpyr, dogPyr, nOctaves, nDogLayers);

    // 构建梯度金字塔: 用于关键点方向分配
    std::vector<Image<float>> grdPyr(nOctaves * nGpyrLayers);
    std::vector<Image<float>> rotPyr(nOctaves * nGpyrLayers);
    build_grd_rot_pyr(gpyr, grdPyr, rotPyr, nOctaves, nLayers);

    // SIFT关键点检测
    detect_keypoints(dogPyr, grdPyr, rotPyr, nOctaves, nDogLayers, kpt_list);

    // 描述子
    if (bExtractDescriptors)
        extract_descriptor(grdPyr, rotPyr, nOctaves, nGpyrLayers, kpt_list);

    return 0;
}

} // end namespace sift
