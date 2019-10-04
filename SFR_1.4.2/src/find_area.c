/*****************************************************************************
Find_area.c

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern FILE *g_mtfout;
#include "sfr.h"

/*****************************************************************************/
/* External functions */
unsigned short fit(unsigned long, double *, double *, double *, double *, 
		   double *, double *, double *);
unsigned short locate_centroids(double *, double *, double *,
				unsigned short, unsigned short, double *);
unsigned short check_slope(double, unsigned short *, int *, double, int);
/* Internal function */
void find_10percent(double *, int, double, int *, int *);

/*****************************************************************************/
/* Data passed to this function is assumed to be radiometrically corrected,  */
/* and oriented vertically, with black on left, white on right. The black to */
/* white orientation doesn't appear to be terribly important in this code.   */
/*     Parameters:
          Input:  farea  = radiometric image data in original ROI
                  size_x = number of columns in original ROI
		  size_y = number of rows in original ROI
		  verbose  Whether ROI search messages are desired
		           0 = No printout
			   1 = Print info as better ROIs are found

	  Output: ncols  = number of columns in recommended ROI
	          nrows  = number of rows in recommended ROI
		  center_x = column center of recommended ROI
		  center_y = row center of recommended ROI

*/
/*****************************************************************************/
short find_area (double *farea, unsigned short size_x, unsigned short size_y, 
		 int *ncols, int *nrows, int *center_x, int *center_y, 
		 unsigned char verbose)
{
  int numcycles;
  unsigned short len, i, dummy, err = 0;
  double *temp=NULL, *shifts=NULL;
  double R2, avar, bvar, offset1, offset2;
  double slope, min, expect_edge;
  int minlen, mincenter, center, signflag, cnt;
  int center_row, check_row, half_width;
  int min_x, max_x, first1, first2, last1, last2;
  double *tptr, *sptr;

  *nrows = size_y;
  *ncols = size_x;

  /* Allocate memory */
  shifts = (double *)malloc(size_y*sizeof(double));
  temp = (double *)malloc(size_y*sizeof(double));

  center_row = *nrows/2;

  /* These values should change during correct operation */
  //在正确运行下，这些值会被改变
  mincenter = -1;
  minlen = -1;

  err = locate_centroids(farea, temp, shifts, size_x, size_y, &offset1); 

  min = 0.0;
  slope = -1./10.;
  signflag = 1;
  for (i = 0; i<center_row - 1.0/fabs(slope)/2;) {
    center = center_row + signflag*i;
    len = (center< size_y-center)? center: size_y-center;
    sptr = shifts + center - len;
    tptr = temp + center - len;
    len = 2*len+1;
    numcycles = 5;

    cnt = 0;
    for (; numcycles>4; len--, sptr++, tptr++) {
      /* Calculate the best fit line to the centroids */
      if (len <= size_y) {
	err = fit(len, tptr, sptr, &slope, &offset2, &R2, &avar, &bvar);
    
	/* Check slope is OK, and set size_y to be full multiple of cycles */
	dummy = len;
	numcycles = 0;
	err = check_slope(slope, &dummy, &numcycles,5.0, 0);
	
	if (numcycles > 4 && R2 > min) {
	  if(verbose)
	    MTFPRINT7("Center row: %d Rows: %d  Angle %.3f  R2 = %f #Cycles=%d SE =  %.3f\n",
		      center, len, atan(fabs(slope))*180.0/M_PI, 
		      R2, numcycles, atan(bvar)*180/M_PI)
	  min = R2;
	  minlen = len;
	  mincenter = center;
	}
	if (numcycles > 4) cnt++;
      }
      len--;

      /* Calculate the best fit line to the centroids */
      err = fit(len, tptr, sptr, &slope, &offset2, &R2, &avar, &bvar);
    
      /* Check slope is OK, and set size_y to be full multiple of cycles */
      dummy = len;
      numcycles = 0;
      err = check_slope(slope, &dummy, &numcycles, 5.0, 0);
    

      if (numcycles > 4 && R2 > min) {
	if(verbose)
	  MTFPRINT7("Center row: %d Rows: %d  Angle %.3f  R2 = %f #Cycles=%d SE =  %.3f\n",
		    center, len, atan(fabs(slope))*180.0/M_PI, 
		    R2, numcycles, atan(bvar)*180/M_PI)
	min = R2;
	minlen = len;
	mincenter = center;
      }
	if (numcycles > 4) cnt++;
    }
    signflag *= -1;
    if (signflag == -1) i++;
    if (signflag == -1 && cnt == 0) break;
  }

  sptr = shifts + mincenter - minlen/2;
  tptr = temp + mincenter - minlen/2;
  err = fit(minlen, tptr, sptr, &slope, &offset2, &R2, &avar, &bvar);
                                                   
  /* Check slope is OK, and set size_y to be full multiple of cycles */
  dummy = minlen;
  numcycles = 0;
  err = check_slope(slope, &dummy, &numcycles,5.0, 1);
	
  if(verbose)
   MTFPRINT7("Best Center: %d Hgt: %d  Angle %.3f  R2 = %f #Cycles=%d SE =  %.3f\n",
	      mincenter, minlen, atan(fabs(slope))*180.0/M_PI, 
	     R2, numcycles, atan(bvar)*180/M_PI)

  *center_x = *ncols/2;
  *center_y = mincenter;
  *nrows = minlen;

  /* Figure out where edge is on first/last line, and then find
     where derivative goes consistently (5pixels) below 10% that value */
  check_row = mincenter - dummy/2;
  expect_edge = (double)(check_row-center_row)*slope + offset2 + offset1;
  find_10percent(farea+check_row*size_x,size_x,expect_edge,&first1, &last1);
  check_row = mincenter + dummy/2;
  expect_edge = (double)(check_row-center_row)*slope + offset2 + offset1;
  find_10percent(farea+check_row*size_x,size_x,expect_edge,&first2, &last2);

  max_x = (last1 < last2)? last2: last1;
  min_x = (first1 < first2)? first1: first2;

  if(verbose)
  MTFPRINT4("Needs from: %d to %d = %d cols \n",min_x, max_x, max_x-min_x)

  if (max_x-min_x+31 < *ncols) { /* Try adding 15 pixels beyond ends */
    *ncols = max_x-min_x+30;
    *center_x = (max_x+min_x)/2;
    half_width = (*ncols+1)/2;
    if (*center_x < half_width) *center_x = half_width;
    if (size_x - *center_x < half_width) *center_x = size_x - half_width;
    if ((*ncols)%2 != 0) (*ncols)++;
  }
  else if (max_x-min_x > 32) {
    *ncols = max_x-min_x;
    *center_x = (max_x+min_x)/2;
    if ((*ncols)%2 != 0) (*ncols)++;
  }
  else { /* Can't add 15 at each end, but not as big as 32 yet */
    /* Need to add enough to get up to 32 */
    *center_x = (max_x + min_x)/2;
    if (*center_x < 16) *center_x = 16; 
    else if (*center_x > size_x-16) *center_x = size_x - 16;
    if (size_x>= 32) *ncols = 32;
  }

  free(temp);
  free(shifts);

  return(0);
}


/* Figure out where the ESF becomes essentially flat.找出ESF在哪里变平
   In this case 'flat' means that 5 neighboring pixel derivatives
   have less than 10% of the value at the expected mid-edge location 
*/
void find_10percent(double *farea,int len,double expect_edge,
		    int *first, int *last)
{

  int i,j, edge, cnt;
  double maxdiff, diff;

  edge = (int)floor(expect_edge+0.5);
  if (edge < 2) edge = 2;
  maxdiff = farea[edge-1] - farea[edge-2];
  for (i=edge,j=0; j<2 && i<len; j++, i++)
  {
      diff = farea[i] - farea[i-1];
      if (diff > maxdiff) maxdiff = diff;
  }

  maxdiff *= 0.1;

  cnt = 0;
  for (i=edge; i>0; i--) {
    diff = fabs(farea[i]-farea[i-1]);
    if (diff < maxdiff) cnt++;
    else cnt = 0;
    if (cnt == 5) break;
  }
  *first = i-1;
  if(*first < 0) *first = 0;

  cnt = 0;
  for (i=edge; i<len-1; i++) {
    diff = fabs(farea[i+1]-farea[i]);
    if (diff < maxdiff) cnt++;
    else cnt = 0;
    if (cnt == 5) break;
  }
  *last = i+1;
  if(*last > len) *last = len;
}
