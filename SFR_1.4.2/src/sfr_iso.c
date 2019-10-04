/*****************************************************************************
 Sfr_iso.c

  				Notice
                                ------

This software is a modification of the "SFR Measurement Algorithm C-code",
copyright PIMA 1998, appearing in Annex A of ISO standard 12233-2000, 
"Photography - Electronic Still Picture Cameras - Resolution Measurement".
Permission to modify this software, and use and distribute the modified 
software, was granted by I3A (the successor organization of PIMA) to the 
MITRE Corporation in 2006.

This MITRE Corporation-modified SFR software was produced for the U.S. 
Government under Contract numbers J-FBI-12-128 and W15P7T-05-C-F600 and is
subject to the Rights in Data-General clause FAR 52.227-l4, Alt. IV (DEC 2007).
    
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

-Redistribution of source code must retain the above copyright notice,
 this list of conditions, and the following disclaimer.
    
-Redistribution in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.
    
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
POSSIBILITY OF SUCH DAMAGE. 

---------------------
  Code modifications by Margaret A. Lepley (MITRE), mlepley@mitre.org
---------------------

*****************************************************************************/
/*****************************************************************************/
/*                                                                           */
/* Revisions to sfr_iso.c                                                    */
/*                                                                           */
/* V1.1 Added user_angle flag to avoid centroid computation                  */
/*                                                           mal 9/06        */
/*                                                                           */
/* V1.3 Removed low contrast check (being performed in wrapper code)         */
/*      Processing can continue even for low/high angles.  Warning printed   */
/*      if more lines of data needed for accurate processing.                */
/*      If less than one phase rotation included in ROI, continue processing */
/*      with all the lines present.                                          */
/*                                                           mal 11/06       */
/*                                                                           */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern FILE *g_mtfout;
#include "sfr.h"
/*****************************************************************************/

#define ALPHA 4.0
#define ESFFILE 8
#define DER3 4
#define PEAK 2
#define ROUND 1

static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0: sqrarg*sqrarg)

unsigned short ftwos(long, double, double *, long, double, double *);
void apply_hamming_window(  unsigned short, unsigned int, unsigned short,
			    double *, long *);
void locate_max_PSF(unsigned int, double *, long *);
void calculate_derivative(unsigned int, double *, double *, double *, int);
unsigned short bin_to_regular_xgrid(unsigned short, double *, double *, 
				    double *, long *,
				    unsigned short, unsigned short);
unsigned short fit(unsigned long, double *, double *, double *, double *, 
		   double *, double *, double *);
unsigned short locate_centroids(double *, double *, double *,
				unsigned short, unsigned short, double *);
unsigned short check_image_data(double *, unsigned short, unsigned short);
unsigned short check_slope(double, unsigned short *, int *, double, int);

/*****************************************************************************/
/* Data passed to this function is assumed to be radiometrically corrected,  */
/* and oriented vertically, with black on left, white on right. The black to */
/* white orientation doesn't appear to be terribly important in this code.   */
/* Radiometric correct data (or at least 0-1 scaled) is important for the    */
/* contrast check.                                                           */
/*     Parameters:
          Input:  farea  = radiometric image data in ROI
                  size_x = number of columns in ROI
		  nrows  = number of rows in ROI

	  Output: freq   = new array of relative spatial frequencies 
	          sfr    = new array of computed SFR at each freq
		  len    = length of freq & sfr arrays
		  slope  = estimated slope value
		  numcycles = number of full periods included in SFR
		  pcnt2  = location of edge mid-point in oversample space
		  off    = shift to center edge in original rows
		  version = 0 = default ([-1 1] deriv, no rounding, no peak)
		            1 = add rounding
			    2 = add peak
			    4 = [-1 0 1] derivative, rounding, no peak
		  iterate = 0 do just a single run
		            1 means more runs after this, don't change farea
		               and let numcycles go as low as 1.0
		  user_angle = flag to indicate if the line fit has
		               been precomputed
		  
                  farea  = size_x*4 of ESF and len*2 of LSF values 
    		                 (if iterate = 0)                            
                  */
/*****************************************************************************/
short sfrProc (double **freq, double **sfr, 
	       int *len,
	       double *farea,
	       unsigned short size_x, int *nrows,
	       double *slope, int *numcycles, int *pcnt2, double *off, double *R2,
	       int version, int iterate, int user_angle)
{
  unsigned short i, j, col, err = 0;
  long pcnt;
  double dt, dt1, sfrc, tmp, tmp2;
  double *temp=NULL, *shifts=NULL, *edgex=NULL, *Signal=NULL;
  double *AveEdge=NULL, *AveTmp=NULL;
  long *counts=NULL;
  int nzero;
  unsigned short size_y;
  unsigned int bin_len;
  double avar, bvar, offset1, offset2, offset;
  double centroid;
  int start_row, center_row;
  double *farea_old;
  double cyclelimit;
  FILE *fp = NULL;

  size_y = *nrows;

  /* Verify input selection dimensions are EVEN */
  if (size_x%2 != 0) { 
    fprintf(stderr, "ROI width is not even.  Does this really matter???\n");
    return 1;
  }

  /* At least this many cycles required. */
  /* For special iterative versions of the code, it can go lower */
  if (iterate) cyclelimit = 1.0;
  else cyclelimit = 5.0;

  /* Allocate memory */
  shifts = (double *)malloc(size_y*sizeof(double));
  temp = (double *)malloc(size_y*sizeof(double));
  edgex = (double *)malloc(size_y*size_x*sizeof(double));
  Signal = (double *)malloc(size_y*size_x*sizeof(double));

  if( !user_angle ) {
     //定位矩心位置
    err = locate_centroids(farea, temp, shifts, size_x, size_y, &offset1); 
    if (err != 0)     { return 2; }

    /* Calculate the best fit line to the centroids */
    /*计算矩心的拟合曲线*/
    err = fit(size_y, temp, shifts, slope, &offset2, R2, &avar, &bvar);
    if (err != 0)     { return 3; }
    
    if (version)
      MTFPRINT4("\nLinear Fit:  R2 = %.3f SE_intercept = %.2f  SE_angle = %.3f\n",  
	      *R2, avar, atan(bvar)*(double)(180.0/M_PI))
  }

  /* Check slope is OK, and set size_y to be full multiple of cycles */
  //检查刀口斜率，以确保后面超采样的质量
  err = check_slope(*slope, &size_y, numcycles, cyclelimit, 1);

  /* Start image at new location, so that same row is center */
  //调整ROI的行
  center_row = *nrows/2;
  start_row = center_row - size_y/2;
  farea_old = farea;
  farea = farea + start_row*size_x;
  /* On center row how much shift to get edge centered in row. */
  /* offset = 0.;  Original code effectively used this (no centering)*/
  if(user_angle)
    offset = *off - size_x/2;
  else
    offset = offset1 + 0.5 + offset2  - size_x/2; 

  *off = offset;
  if(version & ROUND || version & DER3)
    offset += 0.125;

  if (err != 0)     {   /* Slopes are bad.  But send back enough
			   data, so a diagnostic image has a chance. */
    *pcnt2 = 2*size_x;  /* Ignore derivative peak */
    return 4; 
  }

  /* reference the temp and shifts values to the new y centre */
  /* Instead of using the values in shifts, synthesize new ones based on 
     the best fit line. */
  // 基于拟合结果更新shifts
  col = size_y/2;
  for (i=0; i < size_y; i++) {
    shifts[i] = (*slope) * (double)(i-col) + offset;
  }

  /* Don't normalize the data.  It gets normalized during dft process */
  //不要对数据进行归一化，数据在DFT处理过程中会被归一化
  /* To normalize here, set dt = min and dt1 = max of farea data      */
  //如果要在这里初始化，将dt设置为最小值，dt1设置为最大值
  dt = 0.0;
  dt1 = 1.0;

  if (version & ESFFILE)
    fp = fopen("esf.txt","w");

  /* Calculate a long paired list of x values and signal values */
  pcnt = 0;
  for (j = 0; j < size_y; j++) {
    for (i = 0; i < size_x; i++) {
      edgex[pcnt] = (double)i - shifts[j];//计算每个点到刀口的距离
      Signal[pcnt] = ((farea[((j*(long)size_x)+i)]) - dt)/(dt1-dt); //归一化每个点的亮度
      if ((version & ESFFILE) && edgex[pcnt] < size_x/2 + 3 &&
	  edgex[pcnt] > size_x/2 - 3)
	fprintf(fp, "%f %f\n", edgex[pcnt], Signal[pcnt]);
      pcnt++;
    }
  }

  if (version & ESFFILE)
    fclose(fp);

  /* Allocate more memory */
  bin_len = (unsigned int)(ALPHA*size_x);
  AveEdge = (double *)malloc(bin_len*sizeof(double));
  AveTmp = (double *)malloc(bin_len*sizeof(double));
  counts = (long *)malloc(bin_len*sizeof(long));
 
  /* Project ESF data into supersampled bins */
  //进行超采样，生成长度为size_x*ALPHA(4)的当行图像（ESF），保存在AveEdge中
  nzero = bin_to_regular_xgrid((unsigned short)ALPHA, edgex, Signal, 
			       AveEdge, counts, 
			       size_x, size_y); 
  free(counts);
  free(Signal);
  free(edgex);

  /* Compute LSF from ESF.  Not yet centered or windowed. */
  // 将ESF转换为差分图LSF， 并计算LSF的矩心
  if ( (version&DER3) ) 
    calculate_derivative( bin_len, AveTmp, AveEdge, &centroid, 1);
  else
    calculate_derivative( bin_len, AveTmp, AveEdge, &centroid, 0);

  if (iterate == 0) {
    /* Copy ESF to output area */
    for ( i=0; i<bin_len; i++ )
      farea_old[i] = AveTmp[i];
  }

  /* Find the peak/center of LSF */
  //寻找LSF的最大值
  locate_max_PSF( bin_len, AveEdge, &pcnt);//这里获得了最大值的中心

  free(shifts);
  free(temp);

  if(version)
    MTFPRINT3("Off center distance (1/4 pixel units): Peak %ld  Centroid %.2f\n", 
	      pcnt-bin_len/2, centroid-bin_len/2)

  if((version & PEAK) == 0)//忽略差分后的最大值
    pcnt = bin_len/2;  /* Ignore derivative peak */
  else
    MTFPRINT("Shifting peak to center\n")

  /* Here the array length is shortened to ww_in_pixels*ALPHA,
     and the LSF peak is centered and Hamming windowed. */
  //将LSF的最大值移到中心位置，并加上汉明窗
  apply_hamming_window((unsigned short)ALPHA, bin_len,
		       (unsigned short)size_x, 
		       AveEdge, &pcnt);

  /* From now on this is the length used. */
  *len = bin_len/2;

  if (iterate == 0) {
    /* Copy LSF_w to output area */
    for ( i=0; i<bin_len; i++ )
      farea_old[size_x*(int)ALPHA+i] = AveEdge[i];
  }

  tmp = 1.0;
  tmp2 = 1.0/((double)bin_len) ;

  /* Now perform the DFT on AveEdge */
  /* ftwos ( nx, dx, lsf(x), nf, df, sfr(f) ) */
  //ftwos( number, dx, *lsf, ns, ds, *sfr)
  (void) ftwos(bin_len, tmp, AveEdge, (long)(*len), 
	       tmp2, AveTmp); 

  if(*freq==NULL)
    (*freq) = (double *)malloc((*len)*sizeof(double));
  if(*sfr==NULL)
    (*sfr) = (double *)malloc((*len)*sizeof(double));

  for (i=0; i<(*len); i++) {
    sfrc = AveTmp[i];
    (*freq)[i]= ((double)i/(double)size_x);   //每个点对应的频率
    (*sfr)[i] = (double) (sfrc/AveTmp[0]);    //归一化sfr
  }

  /* Free */
  free(AveEdge);
  free(AveTmp);

  *nrows = size_y;
  *pcnt2 = pcnt;

  return(0);
}

/*****************************************************************************/
unsigned short locate_centroids(double *farea, double *temp, 
				double *shifts,
                                unsigned short size_x, unsigned short size_y,
				double *offset) {
  unsigned long i, j;
  double dt, dt1, dt2;

  /* Compute the first difference on each line. Interpolate to find the 
     centroid of the first derivatives. */
  for (j = 0; j < size_y; j++) {
    dt = 0.0;
    dt1 = 0.0;
    for (i = 0; i < size_x-1; i++) {
      dt2 = farea[(j*(long)size_x)+(i+1)] - farea[(j*(long)size_x)+i]; //进行差分
      dt += dt2 * (double)i;
      dt1 += dt2;
    }
    shifts[j]=dt/dt1;
  }

  /* check again to be sure we aren't too close to an edge on the corners. 
     If the black to white transition is closer than 2 pixels from either 
     side of the data box, return an error of 5; the calling program will 
     display an error message (the same one as if there were not a difference 
     between the left and right sides of the box ) */
  //防止边缘离ROI的角太近起码距离两个像素
  if (shifts[size_y-1] < 2  || size_x - shifts[size_y-1] < 2) {
    fprintf(stderr,"** WARNING: Edge comes too close to the ROI corners.\n");
    return 5;
  }
  if (shifts[0] < 2 || size_x - shifts[0] < 2){
    fprintf(stderr,"** WARNING: Edge comes too close to the ROI corners.\n");
    return 5;
  }

  /* Reference rows to the vertical centre of the data box */
  j = size_y/2;
  dt = shifts[j];
  for (i = 0; i < size_y; i++) {
    temp[i] = (double)i - (double)j;
    shifts[i] -= dt;
  }
  *offset = dt;
  return 0;
}
/***************************************************************************/
//这里的xy是反过来的。。
unsigned short fit(unsigned long ndata, double *x, double *y, double *b, 
		   double *a, double *R2, double *avar, double *bvar)
{
  unsigned long i;
  double t,sxoss,syoss,sx=0.0,sy=0.0,st2=0.0;
  double ss,sst,sigdat,chi2,siga,sigb;

  *b=0.0;
  for ( i=0; i < ndata; i++ ) {
    sx += x[i];
    sy += y[i];
  }
  ss=(double)ndata;
  sxoss=sx/ss;
  syoss=sy/ss;
  for ( i=0; i < ndata; i++ ) {
    t = x[i] - sxoss;
    st2 += t*t;
    *b += t * y[i];  
  }
  *b /= st2;         /* slope  */
  *a =(sy-sx*(*b))/ss; /* intercept */
  siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
  sigb=sqrt(1.0/st2);
  chi2=0.0;
  sst=0.0;
  for (i=0; i < ndata; i++) {
    chi2 += SQR( y[i] - (*a) - (*b) * x[i]); 
    sst += SQR( y[i] - syoss); 
  }
  sigdat=sqrt(chi2/(ndata-2));
  siga *= sigdat;
  sigb *= sigdat;
  *R2 = 1.0 - chi2/sst;
  *avar = siga;
  *bvar = sigb;
  return 0;
}

/****************************************************************************/
unsigned short check_slope( double slope, unsigned short *size_y, int *numcycles, double mincyc, int errflag)
{
  double absslope;

  absslope = fabs(slope);

  if( *numcycles <= 0 )
    (*numcycles) = (int)((*size_y)*absslope);

  /* If the slope is too small not enough rows for mincy (typically 5) 
     full periods, then alert the user. */
  if ( absslope < mincyc/(double)(*size_y) ){
    if(errflag == 1) {
      fprintf(stderr, "WARNING: Edge angle (%f) is so shallow it needs\n",
	    atan(slope)*180/M_PI );
      fprintf(stderr, "  %d lines of data (%.1f cycles) for accurate results\n",
	      (int)ceil(mincyc/absslope), mincyc);
      return 0;
    }
    else
      return 1;
  }

  if( absslope > (double)(1.0/4.0) ) { 
      int rows_per_col;
      double bestcycle, x;

      if( absslope > (double)(5.0/4.0) ) {
	fprintf(stderr, "ERROR: Edge angle (%f) is too large\n",
	      atan(slope)*180/M_PI);
	return 1;
      }

      rows_per_col = (int)floor(1/absslope + 0.5);
      x = fabs(1/(double)rows_per_col - absslope);
      bestcycle = 4*rows_per_col*x*ceil(1.0/x/(double)rows_per_col/(double)rows_per_col - 1.0);
      if((int)ceil(mincyc*bestcycle) > *size_y) {
	if(errflag == 1) {
	  fprintf(stderr, "WARNING: Edge angle (%f) will reduce SFR accuracy\n",
	      atan(slope)*180/M_PI);
	  fprintf(stderr, "   if %g * %f = %d lines of data are not in ROI\n\n",
	      mincyc, bestcycle,
	      (int)ceil(mincyc*bestcycle));
	  return 0;
	}
	else
	  return 1;
      }
  }

  /* Figure out how many lines to use for size_y: new window will start at 
     top and go down that number of lines < size_y such that an integer 
     number of x-transitions are made by the edge; for example, if we had a 
     slope of 10 (the edge goes down 10 lines before jumping over one pixel 
     horizontally), and size_y = 35, the new size_y is going to be 30 (an
     integer multiple of 10, less than 35). */

  if( ((*numcycles)/absslope) <= *size_y)
    *size_y = (unsigned short)((*numcycles)/absslope);
  
  return 0;
}

/*****************************************************************************/
/* Notes: this part gets averages and puts them in a number of bins, equal to 
size_x times alpha.  Next a long check is done in case one bin gets no values 
put into it: if this is the case, it will keep checking previous bins until it
finds one with non-zero counts and will use that value as its current bin 
average. If the first bin has zero counts the program checks bins in the 
forward rather than reverse direction. If, in any case, the end of the array 
of bins is reached before finding a non-zero count, the program starts 
checking in the opposite direction. A bin with zero counts is not allowed, 
since each bin will be divided by counts at the end. */

unsigned short bin_to_regular_xgrid(unsigned short alpha,//alpha->指的是超取样的倍数
				    double *edgex, double *Signal, 
				    double *AveEdge, long *counts,
				    unsigned short size_x,
				    unsigned short size_y)
{
  long i, j, k,bin_number;
  long bin_len;
  int nzeros;

  bin_len = size_x * alpha;  //扩大四倍

  for (i=0; i<bin_len; i++) {
    AveEdge[i] = 0;
    counts[i] = 0;
  }
  for (i=0; i<(size_x*(long)size_y); i++) {
    bin_number = (long)floor((double)alpha*edgex[i]);//向下取整
    if (bin_number >= 0) {//bin_number肯定会有相同的，这是什么操作？
      if (bin_number <= (bin_len - 1) ) {
        AveEdge[bin_number] = AveEdge[bin_number] + Signal[i];//把每一行的距离边缘x轴一样远的信号相加
        counts[bin_number] = (counts[bin_number])+1;//记录下对应位置有多少个信号相加
      }
    }
  }

  nzeros = 0;
  for (i=0; i<bin_len; i++) {
    j = 0;
    k = 1;
    //感觉写的有点复杂
    if (counts[i] == 0) {
      nzeros++; //记录有多少个位置是空的，即没有信号
      //K的作用：因为这里的信号为0，找到后面第一个不为零的值，赋给当前这个零
      if (i == 0) {
        while (!j) { //当j==0时，表示此处的信号为0
          if (counts[i+k] != 0) {//第一行的元素
	    AveEdge[i] = AveEdge[i+k]/((double) counts[i+k]);
            j = 1;//充当flag。。。为啥不用布尔类型
	  }
          else k++;//直到找到第一个不为零的数字才会停止
	}
      } else {
        while (!j && ((i-k) >= 0) ) { //j==0&&i-k>=0  j==0说明 counts[i]是0  i-k>0说明  k在i前面，找前面不为零的数值赋给AveEdge[i]
          if ( counts[i-k] != 0) {
	    AveEdge[i] = AveEdge[i-k];   /* Don't divide by counts since it already happened in previous iteration */
	    j = 1;
	  } 
	  else k++;
	}
        if ( (i-k) < 0 ) {//k>i,其实联合上面那段，就是：此处的counts[i]累加次数为零，所以AveEdge[i]也就是0，所以要找到附近一个不为零的值赋给AveEdge[i]
	  k = 1;
	  while (!j) {
	    if (counts[i+k] != 0) {
	      AveEdge[i] = AveEdge[i+k]/((double) counts[i+k]);
	      j = 1;
	    } 
	    else k++;
	  }
	}
      }
    } 
    else 
      AveEdge[i] = (AveEdge[i])/ ((double) counts[i]);//如果此处不为零，直接就求个平均值
  }

  if (nzeros > 0) {//提示信息
    fprintf(stderr, "\nWARNING: %d Zero counts found during projection binning.\n", nzeros);
    fprintf(stderr, "The edge angle may be large, or you may need more lines of data.\n\n");
  }
  return nzeros;
}
/*****************************************************************************/
/* This has been modified from Annex A, to more closely match Annex D and 
   reduce finite difference errors.  Now allows either [-1 1] derivative
   (when separation = 0) or [-1/2 0 1/2] derivative (when separation=1)

   Inputs:   len          length of ESF array
             AveEdge      array of ESF values
             separation   type of derivative 
	                    0 = [-1 1]
			    1 = [-1/2 0 1/2]

   Outputs:  AveTmp       array of original ESF values
             AveEdge      array of derivative (LSF) values
             centroid     centroid of the derivative 
*/
//这里方法写的是计算导数，但是看起来应该是差分
void calculate_derivative(unsigned int len, double *AveTmp, double *AveEdge,
			  double *centroid,int separation)
{
  unsigned long i;
  double dt, dt1;

  dt = 0.0;
  dt1 = 0.0;

  for (i=0; i< len; i++) 
    AveTmp[i] = AveEdge[i];

  for (i=1; i< len-separation; i++) {
    /* Not wasting time with division by 2 since constant factors don't 
       change SFR computation */
    AveEdge[i] = (AveTmp[i+separation]-AveTmp[i-1]); //这里使用的separation=1，那不就是AveTmp[i+1]-AveTmp[i-1]
    if (separation == 1)
      AveEdge[i] /= 2.0;
    dt += AveEdge[i] * (double)i;
    dt1 += AveEdge[i];
  }

  *centroid = dt/dt1;

  AveEdge[0] = AveEdge[1];
  if (separation == 1) AveEdge[len-1] = AveEdge[len-2];
}
/*****************************************************************************/
//点扩展，找到LSF中的最大值
void locate_max_PSF(unsigned int len, 
                    double *AveEdge, long *pcnt2) //点拓展函数
{
  unsigned long i;
  double dt=0.0,dt_new=0.0;
  long left = -1L,right = -1L;

  /* find maximim value in Point Spread Function array */
  //找到在AveEdge中的最大值，这里应该是LSF的最大值啊。。。为什么变成了PSF array
  //dt->最大值
  for (i=0; i<len; i++) {
    dt_new = fabs(AveEdge[i]);
    if ( dt_new > dt) {
      (*pcnt2) = (long) i;
      dt = dt_new;//dt是最大值
    }
  }
  /* find leftmost and rightmost occurrence of maximum */
  for (i=0; i<len; i++) {
    dt_new = fabs(AveEdge[i]);//AveEdge[i]差分后有正有负
    if ( dt_new == dt ) {                  //   [l l l l l l l max]
      if ( left < 0 ) left = i;            //    ↑     ↑        ↑
      right = i;                           //   left center   right
    }
  }
  /* find centre of maxima */
  (*pcnt2) = (right+left)/2;    //最大值离最左边的中心位置
}

/*****************************************************************************/
void apply_hamming_window(  unsigned short alpha,
                            unsigned int oldlen,  //size_x*4
                          unsigned short newxwidth, //size_x
			  double *AveEdge, long *pcnt2)
{
  long i,j,k, begin, end, edge_offset;
  double sfrc;

  /* Shift the AvgEdge[i] vector to centre the lsf in the transform window */
  // 将Avgedge移到中心位置， 两边由于移动造成的空位由0补齐
  edge_offset = (*pcnt2) - (oldlen/2);//不能理解为什么一定要反着计算，pcnt2肯定是小于oldlen/2吧。。
  if (edge_offset != 0) {                                                   //cer=6
    if (edge_offset < 0 ) {  //这里根据分析的话，应该edge_offset只会小于0啊。。     ↓
      for (i=oldlen-1; i > -edge_offset-1; i--)            //   [l l l l l l l max l l l l l]
                  AveEdge[i] = (AveEdge[i+edge_offset]);   //    ↑     ↑        ↑
      for (i=0; i < -edge_offset; i++)                     //   left center=3 right
                  AveEdge[i] = 0.00; /* last operation */    //offset = center-cer=-3
    } else {                                               //            cer=6
                                                           //              ↓
      for (i=0; i < oldlen-edge_offset; i++)               // [0 0 0 l l l l l l l max l l]
                  AveEdge[i] = (AveEdge[i+edge_offset]);   //  ↑           ↑        ↑
      for (i=oldlen-edge_offset; i < oldlen; i++)          //            center=3
		  AveEdge[i] = 0.00;
    }
  }
  /* Multiply the LSF data by a Hamming window of width NEWXWIDTH*alpha */
  //将begin和end两侧的值用0填充,但是感觉没啥用
  begin = (oldlen/2)-(newxwidth*alpha/2);//上下看来，begin只会等于0，因为oldlen=bin_len, bin_len=size_x*alpha == newxwidth*alpha
  if (begin < 0) begin = 0;
  end = (oldlen/2)+(newxwidth*alpha/2);
  if (end > oldlen )  end = oldlen;
  for (i=0; i< begin; i++) 
    AveEdge[i] = 0.0;
  for (i=end; i< oldlen; i++) 
    AveEdge[i] = 0.0;

  // 给begin和end之间的数据加上汉明窗
  // 汉明窗 W(n,α ) = (1 -α ) - α cos(2*PI*n/(N-1)) ，（0≤n≤N-1）
  // 一般情况下，α取0.46
  // 下面计算方法等于窗发生了平移（故符号发生了变化）， 结果是一样的
  for (i=begin,j = -newxwidth*alpha/2; i < end; i++,j++) {
    sfrc = 0.54 + 0.46*cos( (M_PI*(double)j)/(newxwidth*alpha/2) );
    AveEdge[i] = (AveEdge[i])*sfrc; 
  }
  //将lsfbegin的位置左移到index0的位置
  //但在代码中应该是不会起作用的，
  if (begin != 0) /* Shift LSF to begin at index 0 (rather than begin) */
    for (k=0, i=begin; k<newxwidth*alpha; i++,k++) 
      AveEdge[k] = AveEdge[i];

}

/*****************************************************************************/
/* This is the DFT magnitude code                                            */

unsigned short ftwos(long number, double dx, double *lsf, 
		     long ns, double ds, double *sfr)
{
  double a, b, twopi, g;
  long i,j;

    //                n-1              k
    // DFT ==> X[k] = Σ  x[n]e^(-j2π - n)
    //                n=0              N

  twopi = 2.0 * M_PI;
  for (j = 0; j < ns; j++){//ns=1/2*bin_len  前一半
    g = twopi * dx * ds * (double)j;
    for (i = 0, a = 0, b = 0; i < number; i++) { 
      a += lsf[i] * cos(g * (double)(i));
      b += lsf[i] * sin(g * (double)(i)); 
    }
    sfr[j] = sqrt(a * a + b * b); 
  }
  return 0;
}

