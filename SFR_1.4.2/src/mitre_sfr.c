

/*****************************************************************************
Mitre_sfr.c

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
This file contains wrapper code written and is primarily written by MITRE.
Code predominantly modified from IS 12233:2000 Annex A is in sfr_iso.c.
---------------------
*****************************************************************************/

/*****************************************************************************/
/*                                                                           */
/* Revisions to mitre_sfr.c                                                  */
/*                                                                           */
/* V1.1 Added user defined edge option 'a'                                   */
/*      Fixed bugs in horizontal edge diagnostic image code                  */
/*                                                           mal 9/06        */
/*                                                                           */
/* V1.2 Alter ROI box to force width across edge to be even (1 smaller)      */
/*      Change printout of num pixels on either side of edge used in SFR     */
/*      OECF readin extrapolates to next greylevel at either end             */
/*      PIV spec check added as default, plus option to avoid the check      */
/*                                                           mal 10/06       */
/*                                                                           */
/* V1.3 ROI and edge position checks modified to give better error messages  */
/*      Modulation check only warns does not stop processing                 */
/*      Acute angle slope checks warn, but don't stop processing             */
/*      Wide angle slope checks warn, but don't stop processing              */
/*      Edge type printout uses grey levels rather than energy level         */
/*                                                           mal 11/06       */
/*                                                                           */
/* V1.4 Added error and halt if input tiff is not a single band              */
/*      Added explanation of * marking out of PIV spec frequencies           */
/*      Fixed logic bug on problem report error re < 1 cycle in ROI          */
/*                                                           mal 9/09        */
/*                                                                           */
/* V1.4.1 Recompiled using libtiff.                                          */
/*                                                           mal 11/12       */
/*                                                                           */
/* V1.4.2 Allowed reading of PGM with width/height on separate lines.        */
/*                                                           mal 11/13       */
/*****************************************************************************/


/*****************************************************************************/
/* This is the "MAIN" wrapper for the SFR program                            */
/*****************************************************************************/

#include <string.h>

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <ctype.h>
#include <fcntl.h>
#include <errno.h>
#include <math.h>

#ifdef USE_TIFF
#include <tiffio.h>
#else
typedef char TIFF;
#endif

#if !defined(MSDOS)
#include <unistd.h>    /* chdir() */
#include <libgen.h>	  /* dirname()*/
#endif

#if !defined(O_BINARY)
#define O_BINARY 0
#endif

#include "sfr.h"
#include "xterndef.h"

#define MM_PER_INCH 25.4 /* 25.400 millimeters equal 1.000 inches */
#define IQS 1

/*****************************************************************************/
/*                     prototypes                                            */
/*****************************************************************************/

static void get_switches(int, char **,char*,char*);	
static void get_args(char*,char*,TIFF**);
static void print_data_rights_notice(void);
static void print_help(void);
static void print_header(char*,char*);
static void read_in_image(TIFF*, unsigned char, char*, int *, int *);
static void read_scan_line(TIFF*,unsigned char*,short);
static short raw_readscanline(unsigned char*,short);
static void input_area(TIFF *, unsigned char *, int, int, int, int);
static int clipping (int, int, double, unsigned char *, int );
static double *read_in_ref_lut(char*, int *, int *); 
static int read_pgm_header(char *);
static void write_debug_image(unsigned char *, int, int);
static unsigned char *make_debug_image(TIFF *, unsigned char *, unsigned char);
static void draw_lines(unsigned char *, double, unsigned char, int, int, double, int);
static void put_problem(char *, int);
static void print_problems(void);
void slope_bounds( double slope, int size_y, int numcycles, double mincyc, char *problem_string);
int reverse_lut(double *ref_lut, double val, int lowest_val, int highest_val); 

void wait_to_exit(void);

unsigned char *g_debug_array;
unsigned char *g_start_ptr;
int g_debug_fullwidth;
int g_debug_fullheight;
int g_userangle;
float g_pt1x, g_pt1y, g_pt2x, g_pt2y;
static int too_many_problems;

short find_area (double *, unsigned short, unsigned short, int *, int *, 
		 int *, int *, unsigned char);
int sfrProc(double **,double **, int*, double *, 
	    unsigned short, int *, double *, int *, int *, double*, double*, 
	    int, int, int);

/* Windows compilers without rint (MinGW is OK) */
#if defined (_WIN32) & !defined(__GNUC__)  
double rint(x)
double x;
{
  int y = (int)x; 
  double half = (double)y + 0.5;
  double ans;

  if (x >= half)
    ans = (double)y + 1.0;
  else
    ans = (double)y;
  return (ans);
}
#else
extern double rint(double);
#endif

int main(int argc, char **argv)
{
  char problem_string[82];
  unsigned char rotation;
  int i;
  double *farea;
  double slope, scale, b;
  int size_x, size_y;
  int len, err, bin_len;
  int rgt_side, left_side;
  int center;
  int numcycles=0;
  double *Freq=NULL;
  double *disp=NULL;
  double *ref_lut;
  double off, R2;
  int center_x, center_y, new_x, new_y;
  int lowest_val, highest_val, grey_level, first, last;
  int piv_err;

  TIFF* tif = NULL;

#if !defined(MSDOS)
  char *fnamep;
  fnamep = dirname(argv[0]);
  if (argc==1 && fnamep != NULL && fnamep[0] == '/') {
    printf ("Changed to directory of executable:  %s\n\n", fnamep);
    chdir((const char *)fnamep);
  }
#else
  atexit(wait_to_exit);
#endif

  /* Once-only initializations */
  g_scan_image_file_id = 0;
  g_problem_count = 0;                /* # of problems in our problem report */
  g_IQS_problem_count = 0;


  get_switches(argc, argv, image_filename, data_filename);              

  /* always append output to this file */
  if ((g_mtfout = fopen(MTFOUTNAME,"a")) == NULL) {
    fprintf(stderr,"Can't open %s\n",MTFOUTNAME);
    exit(1);
  }

  /* command line args */
  //命令行参数
  get_args(image_filename,data_filename,&tif);//参数获取源文件并读取文件
 
  /* print the user entered data (corner coords, etc.) in the output */
  /* 时间， 版本号， 以及get_args函数里面获取到的信息*/
  print_header(image_filename,data_filename);

  g_target_res = 500;
  if (abs((int)g_ppi - g_target_res) > 10) {
    sprintf(problem_string,
	    "Resolution outside the %d-%d ppi PIV spec range\n",
	    g_target_res-10,g_target_res+10);
    put_problem(problem_string, IQS);
    put_problem("\n", IQS);
  }

  //从ROI 4个角落各读取4个像素点， 用于判断刀口方向， 并确认ROI区域的合法性
  input_area(tif, &rotation, 
	     g_test_pattern_yul, g_test_pattern_ylr, 
	     g_test_pattern_xul, g_test_pattern_xlr);

  /* read the image in - we store it internally in a standard way */
  // 读取ROI图像， 如果并将水平刀口转换为垂直刀口
  read_in_image(tif,rotation,image_filename,&size_x,&size_y);

  /* get reflectance */
  // 创建反射查找表， 将图像灰阶映射在0~1范围内
  ref_lut = read_in_ref_lut(data_filename, &lowest_val, &highest_val);

  // 将原图ROI映射到反射ROI
  len = size_x*size_y;
  farea = (double *)malloc(len*sizeof(double));
  for(i=0; i<len; i++) {    
    grey_level = g_image_array[i];
    if (grey_level < lowest_val || grey_level > highest_val) {
      MTFPRINT("ERROR: OECF range does not cover the entire image\n")
      MTFPRINT2("       Greylevel %d found in image, but not OECF\n", grey_level)
      exit(-1);
    }
    farea[i] = (double)ref_lut[grey_level];
  }

  if(g_autorefine) {   
    find_area(farea, (unsigned short)size_x, (unsigned short)size_y, &new_x, &new_y, &center_x, &center_y, g_extended);
    if (new_x != size_x || new_y != size_y ||
	center_x != size_x/2  || center_y != size_y/2) {
      center_x += (g_test_pattern_xlr + g_test_pattern_xul)/2 - size_x/2;
      center_y += (g_test_pattern_ylr + g_test_pattern_yul)/2 - size_y/2;
      g_test_pattern_xul = center_x - new_x/2;
      g_test_pattern_yul = center_y - new_y/2;
      g_test_pattern_xlr = g_test_pattern_xul + new_x;
      g_test_pattern_ylr = g_test_pattern_yul + new_y;
      read_in_image(tif,rotation,image_filename,&size_x,&size_y);
      len = size_x*size_y;
      for(i=0; i<len; i++)
	farea[i] = (double)ref_lut[g_image_array[i]];
    }
    MTFPRINT("Refined coordinates of slanted edge area:\n")
      MTFPRINT3("Upper Left Col =\t\t%d\tUpper Left Row =\t%d\n",
            g_test_pattern_xul,g_test_pattern_yul)
      MTFPRINT3("Center Col =\t\t\t%d\tCenter Row =\t\t%d\n",
            (g_test_pattern_xlr+g_test_pattern_xul)/2,
	    (g_test_pattern_ylr+g_test_pattern_yul)/2)
      MTFPRINT3("Width =\t\t\t\t%d\tHeight =\t\t%d\n",
            (g_test_pattern_xlr-g_test_pattern_xul),
	    (g_test_pattern_ylr-g_test_pattern_yul))
    MTFPRINT("\n\n")
  }

  if (g_userangle) {
    double fslope;

    slope = (g_pt1x-g_pt2x)/(double)(g_pt1y-g_pt2y);
    b = g_pt1x - slope*g_pt1y;
    if (rotation == RIGHT)
      off = (g_test_pattern_yul + size_y/2)*slope + b - g_test_pattern_xul;
    else {
      off = (g_test_pattern_xlr - 1 - size_y/2 - b ) / slope - g_test_pattern_yul;
      slope = -1.0/slope;
    }
    fslope = fabs(slope)*size_y/2;
    if(  off-fslope < 2  || off+fslope >= size_x-2 ) {
      MTFPRINT("** ERROR: User specified edge does not work with ROI **\n");
      if (off-fslope < 2)
	MTFPRINT3("Edge Angle = %.3f  degrees    located %d pixels before ROI\n",
              (atan(slope)*(double)(180.0/M_PI)), (int)rint(off))  
      else
	MTFPRINT3("Edge Angle = %.3f  degrees    located %d pixels after ROI\n",
              (atan(slope)*(double)(180.0/M_PI)), (int)rint(off-size_x))  
      if (g_debug) {
	MTFPRINT("See location of edge points in diagnostic image\n")
	write_debug_image(g_debug_array,g_debug_fullwidth, g_debug_fullheight);
      }
      exit(-1);
    }
  }  

  //计算ROI中超出取值范围[0,255]外的像素比例
  clipping (0, 255, 0.02, g_image_array, len);

  //计算sfr结果
  /* calculate the sfr on this area */
  err = sfrProc(&Freq, &disp, &bin_len, farea, (unsigned short)size_x, &size_y, &slope, &numcycles,&center,&off, &R2, g_version, 0, g_userangle);

  /* Add messages to problem report */
  slope_bounds( slope, size_y, numcycles, 5.0, problem_string);

  if(g_debug) {
    draw_lines(g_debug_array, slope, rotation, center, size_y, off, size_x);
    write_debug_image(g_debug_array,g_debug_fullwidth, g_debug_fullheight);
  }

  if (!g_userangle)
    MTFPRINT2("R2 of linear edge fit = %.3f\n", R2)
  else 
    MTFPRINT("User input edge location: No R2 available\n")
  MTFPRINT4("Edge Angle = %.3f  degrees  \nCycle Length = %.3f  \t#Cycles = %d\n",
              (atan(slope)*(double)(180.0/M_PI)), fabs(1.0/slope), numcycles )  

  left_side = size_x/2+off;  rgt_side = size_x/2-off;
  if (left_side > bin_len/4) left_side = bin_len/4;
  if (rgt_side > bin_len/4) rgt_side = bin_len/4;

  MTFPRINT2("SFR computed using ~%d pixels on left or top side of the edge,\n", left_side)
  MTFPRINT2("               and ~%d pixels on right or bottom side of the edge,\n", rgt_side)
  MTFPRINT("              all within the ROI\n")
  if (rotation == RIGHT)
    MTFPRINT3("    over %d rows, centered at row %d\n", size_y, (g_test_pattern_yul + g_test_pattern_ylr)/2 )
  else
    MTFPRINT3("    over %d cols, centered at col %d\n", size_y, (g_test_pattern_xul + g_test_pattern_xlr)/2 )

  if (center/4 < 16) {
    MTFPRINT("\nERROR: Too few pixels across the edge for valid SFR computation\n\n")
    exit(-1);
  }
  if ( left_side < 20 || rgt_side < 20) 
    MTFPRINT("Warning: Low width across the edge. SFR values may be suspect.\n")

  /* log any problems we encountered into output */
  if (g_reversepolarity)
    MTFPRINT("\nNOTE: Original image polarity was reversed by user request.\n")

  if (err) {
    MTFPRINT ("** ERROR in computation. SFR values unknown **\n\n")
      exit(-1);
  }

  scale = g_ppi/MM_PER_INCH;
  MTFPRINT("\n\ncy/mm         \t SFR  ")
  if(rotation == RIGHT)
    MTFPRINT("\t  Vert ")
  else
    MTFPRINT("\t  Horz ")
  if(farea[4*size_x-1] - farea[0] > 0) 
    MTFPRINT("black-to-white** edge\n(obj.plane)\n")
  else
    MTFPRINT("white-to-black** edge\n(obj.plane)\n")

  piv_err=0;
  for( i=0; i<bin_len/2; i++) {
    double f, sfr;
    double freq, fd_scale;
    //对加汉明窗导致的功率损失进行补偿？不确定
    freq = M_PI*Freq[i];
    /* [-1 0 0 0 1]  freq /= 1.0; */
    if (g_version & 4) /* [-1 0 1] */
      freq /= 2.0;
    else freq /= 4.0;  /* [-1 1] */
    if (freq == 0.0) fd_scale = 1.0;
    else             fd_scale = freq / sin(freq); 

    f = Freq[i]*scale;
    sfr = disp[i]*fd_scale;

    MTFPRINT3("%9.6f   \t%f ",f,  sfr) 

    if(Freq[i] > 0.5)
      MTFPRINT("\t#\n")    /* Mark frequencies above Nyquist */
    else {
      double lower = 1.02829 - 1.67473e-1*f + 1.06255e-2*f*f - 2.80874e-4*f*f*f;
      double upper = 1.12;
      if ((f<=10 && f >= 0.998) && (sfr < lower || sfr > upper)) {
	     if (g_nocompare)
	       MTFPRINT("\n")
	     else {
		   piv_err = 1;
	       MTFPRINT("\t*\n")    /* Mark out of PIV spec values */
	    }
        sprintf(problem_string,
                "Computed SFR at %.3f cy/mm is outside PIV spec range of %.3f - %.2f\n",
		  f,lower,upper);
        put_problem(problem_string, IQS);
      }
      else MTFPRINT("\n")
    }
  }

  MTFPRINT("\n#  Frequencies above Nyquist\n") 
  if (piv_err)
    MTFPRINT("\n*  SFR is outside PIV spec.  See problem report for details.\n")
  first = reverse_lut(ref_lut, farea[0], lowest_val, highest_val); 
  last = reverse_lut(ref_lut, farea[4*size_x-1], lowest_val, highest_val); 
  if (rotation == RIGHT) {
    MTFPRINT2("** Edge type based on leftside gray level (%d)\n", first)
    MTFPRINT2("    and rightside gray level (%d)\n", last)
  }
  else {
    MTFPRINT2("** Edge type based on top gray level (%d)\n", first)
    MTFPRINT2("    and bottom gray level (%d)\n", last)
  }
  MTFPRINT("    If type isn't correct, use 'f' option to adjust polarity\n\n")

  print_problems();

  if (g_extended){
    int offset, begin, end;
    int j, cnt;
    double threshold,max_edge,current;

    /* Derivatives below this value is considered insignificant for printing */
    max_edge = fabs(farea[4*size_x+bin_len]);
    threshold = max_edge / 25.0;

    /* Find 5 contiguous insignificant values first occur before center */
    offset = center - bin_len;
    j = bin_len;    cnt = 0;
    while (cnt < 5 && j>=0) {
      current = fabs(farea[4*size_x+j]);
      if ( current < threshold ) cnt ++;
      else cnt = 0;
      if (current > max_edge) max_edge = current;
      j--;
    }
    if (j<0) j=0;
    begin = j + offset;

    /* Find 5 contiguous insignificant values first occur after center */
    j = bin_len;    cnt = 0;
    while (cnt < 5 && j<bin_len*2) {
      current = fabs(farea[4*size_x+j]);
      if ( current < threshold ) cnt ++;
      else cnt = 0;
      if (current > max_edge) max_edge = current;
      j++;
    }
    if (j==2*bin_len) j=2*bin_len-1;
    end = j + offset;

    if ( farea[4*size_x+bin_len] < 0 ) max_edge *= -1;


    MTFPRINT("\n\nPixel+ \t\t LSF_w \t\t AvgESF")
    if (g_version&4)
	MTFPRINT(" \t \n")
    else
	MTFPRINT("(with -0.125 offset)\n")

    for( i=begin, j=begin-offset; i<0; j++, i++)
	MTFPRINT3("% .2f \t\t% f\n", (i-center)/4.0,  farea[4*size_x + j]/max_edge)

    for( ; i<end && i<size_x*4; j++, i++) {
        MTFPRINT4("% .2f \t\t% f \t% f\n", (i-center)/4.0,  farea[4*size_x + j]/max_edge,  farea[i])
    }
    for(; i<end; j++, i++)
        MTFPRINT3("% .2f \t\t% f\n", (i-center)/4.0,  farea[4*size_x + j]/max_edge)

    MTFPRINT("\n")

    MTFPRINT("+ Pixel values are centered on the computed line fit\n")
    MTFPRINT("  LSF_w values are scaled to a maximum of 1\n")
    if (!(g_version&4))
       MTFPRINT("  Actual pixel values for the AvgESF are shifted up +0.125 from values here\n")
    MTFPRINT("\n")
  }

  return(EXIT_SUCCESS);

}                               /* end of main */

/*****************************************************************************/
/*                                                                           */
/*****************************************************************************/

void get_switches(int argc, char **argv, char *image_filename, char *data_filename)
{
  int centerx = 0, centery = 0, width = 0, height = 0;
  char switch_entry = ' ';
  int count = 0;
  int read_value = 0;

  g_streamline_args=g_reversepolarity=g_extended=g_debug=g_autorefine=0;
  g_nocompare=0;
  g_userangle=0;
  g_version=0;
  g_center=0;

  /* This is a COMMAND LINE mode for tif images, i.e., 
     read in all args from command line rather than waiting for questions;
     usage: sfr image.tif image.dat center_left_col center_left_row
            width height
  */
  if (argc > 1) {         
      argv++;
      strcpy(image_filename,argv[0]); argv++;
      strcpy(data_filename,argv[0]);  argv++;
      centerx = atoi(argv[0]); argv++;
      centery = atoi(argv[0]); argv++;
      width = atoi(argv[0]); argv++;
      height = atoi(argv[0]); argv++;
      g_ppi = atof(argv[0]); argv++;
      if(argc > 8)  { 
	g_pt1x = atof(argv[0]); argv++;
	g_pt1y = atof(argv[0]); argv++;
	g_pt2x = atof(argv[0]); argv++;
	g_pt2y = atof(argv[0]); argv++;
		  g_userangle = 1;
      }
      g_streamline_args = 1;
  }

 //输出信息
    fprintf(stderr,"SFR version %s \n\n", VERSION);
  fprintf(stderr,"  a      Compute edge tilt angle from user entered points \n");
  fprintf(stderr,"  b      Auto-refine input region \n");
  fprintf(stderr,"  c      ROI defined by center point instead of UL corner\n");
#ifdef USE_TIFF
  fprintf(stderr,"  d      Create diagnostic image (_box.tif) \n");
#else
  fprintf(stderr,"  d      Create diagnostic image (_box.pgm) \n");
#endif
  fprintf(stderr,"  e      Verbose output \n");
  fprintf(stderr,"  f      Reverse image polarity \n");
  fprintf(stderr,"  h      Help & Program Notice \n");
  fprintf(stderr,"  n      Don't compare output to PIV spec \n\n");
  /* 
  fprintf(stderr,"  3      3pt derivative [-1 0 1] vs default [-1 1] \n");
  fprintf(stderr,"  p      Peak centering \n");
  fprintf(stderr,"  r      Bin rounding \n");
  fprintf(stderr,"  q      Projected ESF before bin averaging (esf.txt)\n"); 
  */
  fprintf(stderr,"Enter run options...or...Press RETURN if none:  ");
  
  do
    {
      count++;
      read_value = scanf("%c", &switch_entry);
      switch (switch_entry) 
        {
          case 'a':
          case 'A':
            g_userangle = 1;
            break;
          case 'b':
          case 'B':
            g_autorefine = 1;
            break;
          case 'c':
          case 'C':
            g_center = 1;
            break;
          case 'd':
          case 'D':
            g_debug = 1;
            break;
          case 'e':
          case 'E':
            /* Force verbose output */
            g_extended = 1;
            break;
          case 'f':
          case 'F':
            /* Force image polarity to be reversed */
            g_reversepolarity = 1;
            break;
          case 'h':
	      case 'H':
	        print_data_rights_notice();
	        print_help();
	        break;
          case 'n':
          case 'N':
            g_nocompare = 1;
            break;
          case 'p':
          case 'P':
            g_version |= 2;
            break;
          case 'r':
          case 'R':
            g_version |= 1;
            break;
          case '3':
            g_version |= 5;
            break;
          case 'q':
          case 'Q':
            g_version |= 8;
            break;
          case ' ':
          case '\n':
          case '\0':
            break;
          default:
            fprintf(stderr," INVALID OPTION INPUT\n");
            exit(1);
            break;
        }
    } while((count < 80) && (switch_entry != '\n') && (switch_entry != '\0')
            && (read_value == 1));      

  if( g_streamline_args ) {
    if (g_center) {
      g_test_pattern_xul = centerx - width/2;
      g_test_pattern_yul = centery - height/2;
    }
    else {
      g_test_pattern_xul = centerx;
      g_test_pattern_yul = centery;
    }
    g_test_pattern_xlr = g_test_pattern_xul + width;
    g_test_pattern_ylr = g_test_pattern_yul + height;
  }

}


/*****************************************************************************/
/*                                                                           */
/*****************************************************************************/

void get_args(char *image_filename,char *data_filename, TIFF** tif)
{
  int centerx, centery, width, height;
  
  while ((*tif == NULL) && (g_scan_image_file_id <= 0)) {
    if (!g_streamline_args) {
      fprintf(stderr,"\n\n Enter image filename:   \t");
      scanf("%s",image_filename);
    }

#ifdef USE_TIFF
        {
            char *tmp;
            if(  ((tmp = strstr (image_filename, ".pgm")) == NULL ||
		 	   tmp != strrchr (image_filename, '.')) &&
			((tmp = strstr (image_filename, ".PGM")) == NULL ||
		 	   tmp != strrchr (image_filename, '.') ) )
                *tif = TIFFOpen((char *)image_filename, "r");
        }
#endif
    if (*tif == NULL ) {

      g_pgm = 0;
      read_pgm_header (image_filename);

      g_scan_image_file_id = open(image_filename, (O_RDONLY | O_BINARY));//用二进制只读形式打开tiff文件
      if (g_scan_image_file_id < 0) {  //如果没打开
	fprintf(stderr, "\nCould not open %s\n", image_filename);
	fprintf(stderr, "%s\n\n", strerror(errno));
      }
    }
  }

  //如果打开的不是tiff文件格式
  if (*tif == NULL ) {
    if ( !g_pgm)
    {
      fprintf(stderr,"Image file is not in TIFF or PGM format\n");
      fprintf(stderr," Enter image header bytes, pixel width, pixel height:\t");
      GET_COORD_PAIR(&g_scan_header,&g_scan_cols,&g_scan_rows);
      g_bps=8;
      g_photometric=1;
    }
  }
#ifdef USE_TIFF
  else {
    unsigned short sampleformat;
	unsigned short samplesperpixel;
    g_bps = 8;
    sampleformat = 1;
	samplesperpixel = 1;
    TIFFGetField( *tif, TIFFTAG_IMAGEWIDTH, &g_scan_cols );
    (void) TIFFGetField( *tif, TIFFTAG_IMAGELENGTH, &g_scan_rows );
    (void) TIFFGetField( *tif, TIFFTAG_BITSPERSAMPLE, &g_bps );

	/* Check to make sure it's a 1 band image */
	(void) TIFFGetField( *tif, TIFFTAG_SAMPLESPERPIXEL, &samplesperpixel);
	if (samplesperpixel != 1) {
		MTFPRINT2("ERROR: Cannot read more than one sample per pixel (%d).\n", samplesperpixel)
		exit(-1);
	}

    /* reading in the tag on photometric interpretation - 7/22/94 - gtk */
    (void) TIFFGetField( *tif, TIFFTAG_PHOTOMETRIC, &g_photometric);
    /*    MTFPRINT2("G_photometric = %d\n", g_photometric); */
    (void) TIFFGetField( *tif, TIFFTAG_SAMPLEFORMAT, &sampleformat);
    if (sampleformat != 1) {
      MTFPRINT2("ERROR: Signed or float imagery. SampleFormat = %d\n", sampleformat);
      exit(-1);
    }
    if( g_bps != 8 ){
      MTFPRINT2("ERROR: Cannot handle %d BitsPerPixel TIFF format. Try 8.\n",g_bps);
      exit(-1);
    }
  }
#endif
  g_bytes_per_pixel = (g_bps+7)/8;
  g_max_pixel_value = (1<<g_bps)-1;


  if (!g_streamline_args) {
      int c;

      while ( getchar() != (int)'\n');
      fprintf(stderr, " OECF: if nonlinear enter filename, if linear press RETURN ");
      while ((c = getchar()) == EOF);
      if((char)c == '\r' || (char)c == '\n')
	strcpy(data_filename, "linear");
      else {
	ungetc((int)c,stdin);
	scanf("%s",data_filename);
      }

      fprintf(stderr, " Enter pixels per inch (PPI):   ");
      scanf("%lf", &g_ppi);

      fprintf(stderr,"\n     (ref: col,row = 0,0 at upper left corner of ENTIRE image,\n");
      fprintf(stderr,"     cols increase left to right, rows increase top to bottom.)\n");
      if (g_center)
	fprintf(stderr,"\n Enter Col, Row for center pixel (~ on edge):\t");
      else
	fprintf(stderr,"\n Enter Col, Row for UL pixel of ROI:\t");
      GET_INT_PAIR(&centerx,&centery)

      fprintf(stderr," Enter Width, Height for region of interest:\t");
      GET_INT_PAIR(&width,&height)

      if( g_center) {
	g_test_pattern_xul = centerx - width/2;
	g_test_pattern_yul = centery - height/2;
      }
      else {
	g_test_pattern_xul = centerx;
	g_test_pattern_yul = centery;
      }
      g_test_pattern_xlr = g_test_pattern_xul + width;
      g_test_pattern_ylr = g_test_pattern_yul + height;

      if( g_userangle ) {
	do {
	  fprintf(stderr,"\nEnter edge endpoints to define edge tilt");
	  fprintf(stderr,"\n(left/right for HorzEdge, top/bot for VertEdge)");
	  fprintf(stderr,"\n Enter Col, Row for 1st Edge Point:\t");
	  GET_FLOAT_PAIR(&g_pt1x,&g_pt1y)
          fprintf(stderr," Enter Col, Row for 2nd Edge Point:\t");
	  GET_FLOAT_PAIR(&g_pt2x,&g_pt2y)
	  if (g_pt2x==g_pt1x || g_pt2y == g_pt1y)
	    fprintf(stderr,"*** Edge points cannot have the same row or col value. Try again. *** \n");
	} while (g_pt2x==g_pt1x || g_pt2y == g_pt1y);
      }
  }

}
void print_help(void)
{
  fprintf(stderr,"*************** HELP ***************\n");
  fprintf(stderr," \n");
  fprintf(stderr,"a \tUser Defined Edge Location\n");
  fprintf(stderr,"    Finds edge slope/offset based upon 2 user defined points\n");
  fprintf(stderr,"    Fractional entries allowed if edge falls between pixels\n");
  fprintf(stderr," \n");
  fprintf(stderr,"b \tAutorefine Input Area\n");
  fprintf(stderr,"    Finds size/location of best line fit (max R2).\n");
  fprintf(stderr," \n");
  fprintf(stderr,"c \tROI defined by center point instead of UL corner\n");
  fprintf(stderr,"    Input ROI location using center point instead of UL corner\n");
  fprintf(stderr," \n");
#ifdef USE_TIFF
  fprintf(stderr,"d \tCreate diagnostic image (imagename_box.tif)\n");
#else
  fprintf(stderr,"d \tCreate diagnostic image (imagename_box.pgm)\n");
#endif
  fprintf(stderr,"    Image shows measurement ROI and best line fit\n");
  fprintf(stderr," \n");
  fprintf(stderr,"e \tVerbose output\n");
  fprintf(stderr,"    Prints ESF and LSF curves.\n");
  fprintf(stderr," \n");
  fprintf(stderr,"f \tReverse image polarity\n");
  fprintf(stderr,"    Reverse image polarity before processing; so program reads pure\n"); 
  fprintf(stderr,"    white as gray=255, not gray=0 (original image file is not changed).\n");
  fprintf(stderr," \n");
  fprintf(stderr,"h \tHelp & Program Notice \n\n");
  fprintf(stderr,"    Prints this notice.\n");
  fprintf(stderr," \n");
  fprintf(stderr,"n \tDon't compare output to PIV spec \n\n");
  fprintf(stderr," \n");
  /*
  fprintf(stderr,"3 \t3-pt derivative to get final LSF\n");
  fprintf(stderr,"    [-1 0 1] computation instead of [-1 1]\n"); 
  fprintf(stderr,"    Finite differences corrects for this as well.\n");
  fprintf(stderr,"r \tRound during binning instead of truncation\n");
  fprintf(stderr,"p \tHamming window centered at LSF peak\n");
  fprintf(stderr,"    Vs default: center Hamming window at line-fit edge postion.\n"); 
  fprintf(stderr,"q \tOutput projected ESF data prior to bin/averaging.\n"); 
  */
  fprintf(stderr," \n");
  fprintf(stderr,"Locating ROI on target image:\n");
  fprintf(stderr,"   Enter UL corner col, row and width/height of region of interest.\n");
  fprintf(stderr,"   This bounds data used in computation.\n");
  fprintf(stderr,"\nOECF: \n");
  fprintf(stderr," If Opto-Electronic Conversion Function is nonlinear, then a data file\n");
  fprintf(stderr," must be constructed containg the input energy - ouput gray level value\n");
  fprintf(stderr," pairs, at least covering the range of the image edge gray levels.\n");
  fprintf(stderr," The data file must have one input energy - output gray level pair per\n");
  fprintf(stderr," line; comment lines starting with the pound sign '#' are allowed at\n");
  fprintf(stderr," the beginning of the file.  The program uses linear interpolation\n");
  fprintf(stderr," between the read-in datafile values.\n");
  fprintf(stderr," [Linear OECF does not require a datafile]\n");
  fprintf(stderr,"\nPPI:\n");
  fprintf(stderr," Pixels per inch scale projected to the object plane, e.g., if there is\n");
  fprintf(stderr," a fiducial mark on each side of the target edge with know distance 'D'\n");
  fprintf(stderr," between, then \n\t\t ppi = pixels between fiducials in image/D \n");
  fprintf(stderr," \n");
}

/*****************************************************************************/
/*                                                                           */
/*****************************************************************************/

void print_header(char *image_filename,char *data_filename)
            
{
  time_t tm_s,*tm_ptr;
  char tod[82];
  double fabs(double);
  
  tm_ptr = &tm_s;

  if (time(tm_ptr) == -1)
    strcpy(tod,"DATE = ????");
  else
    strcpy(tod,ctime(tm_ptr));

  MTFPRINT3("\n\n%.24s\tv%s\n",tod,VERSION)

  if(g_version!=0) {
    MTFPRINT3("   Bin rounding = %d    Peak centering = %d  (0=OFF)\n",
	      g_version&1, g_version&2)
    if (g_version&4)
      MTFPRINT("   Derivative filter = [-1 0 1]\n\n")
    else
      MTFPRINT("   Derivative filter = [-1 1]\n\n")
  }

  MTFPRINT2("Image file name = \t%s\n",image_filename)
  if( strcmp(data_filename,"linear") == 0 )
    MTFPRINT("OECF: linear\n")
  else
    MTFPRINT2("OECF: nonlinear (%s)\n",data_filename)

  if (g_scan_image_file_id == 0)
    MTFPRINT3("\nImage width, height:  %d  %d\n", g_scan_cols, g_scan_rows)
  else
    MTFPRINT4("\nImage hdr, width, height:  %d  \t%d  %d\n",
              g_scan_header, g_scan_cols, g_scan_rows)

  MTFPRINT2("PPI scale at object plane: \t%g ppi\n",g_ppi)

  MTFPRINT("\nSlanted edge ROI:\n")
  MTFPRINT3("Upper Left Col =\t\t%d\tUpper Left Row =\t%d\n",
            g_test_pattern_xul,g_test_pattern_yul)
  MTFPRINT3("Lower Right Col =\t\t%d\tLower Right Row =\t%d\n",
            g_test_pattern_xlr-1,g_test_pattern_ylr-1)
  MTFPRINT3("Center Col =\t\t\t%d\tCenter Row =\t\t%d\n",
            (g_test_pattern_xlr+g_test_pattern_xul)/2,
	    (g_test_pattern_ylr+g_test_pattern_yul)/2)
  MTFPRINT3("Width =\t\t\t\t%d\tHeight =\t\t%d\n",
            (g_test_pattern_xlr-g_test_pattern_xul),
	    (g_test_pattern_ylr-g_test_pattern_yul))

  if (g_userangle) {
    MTFPRINT("\nUser defined Edge points:\n")
    MTFPRINT3("Edge Point1 Col =\t\t%g\t  Edge Pt1 Row =\t%g\n", g_pt1x,g_pt1y)
    MTFPRINT3("Edge Point2 Col =\t\t%g\t  Edge Pt2 Row =\t%g\n", g_pt2x,g_pt2y)
  }

  MTFPRINT("\n\n")
}


/*****************************************************************************/
/*								             */
/*****************************************************************************/

/* Read in the image and store it in a big array.  We also flip the image 
   into a standard, horizontal format with the white patch at the top
   and right. */

void read_in_image(TIFF *tif, unsigned char rotation, char *image_name, 
		   int *size_x, int *size_y)
{
  short input_row;
  short array_column_start=0;
  short array_row_index;
  short array_column_index;
  short array_column_delta=1;
  short array_row_delta=1;
  short array_row_start=0;
  unsigned short actual_rows, actual_cols;

  unsigned char *buf;
  short i;

  /* enough to hold a single scan line */
  buf = NULL;
  if (g_scan_image_file_id) {
    buf = (unsigned char*) calloc((unsigned int)g_scan_cols,sizeof(unsigned char));
  }
#ifdef USE_TIFF
  else {
    buf = (unsigned char*) calloc((unsigned int)TIFFScanlineSize(tif),sizeof(unsigned char));
  }
#endif
  if ( buf == NULL ) {
    fprintf(stderr,
	    "\nCan't allocate memory for scanline buffer\n\n");
    exit(-1);
  }

  if (g_debug) {
    int j, k;

    /* given an input image file name of the form "base.ext", here we create 
       an output file with the form "base8_box", where "base8" was taken from 
       the first 8 characters of the input file name - the first 8 characters 
       of "base".  If "base" has less than 8 characters, then all characters 
       before the "." are used. 9/2/94 - gtk */

    i = strlen(image_name);
    j=i;
    while ((image_name[i] != '/')&&(i>=0)) 
      {
	if (image_name[i] == '.') j=i;
	i--;
      }
    i++;
    j=j-i;
    if (j > 20) j=20;		/* no more than 8 characters */
    for (k=0;k<j;k++) 
      g_debug_file_name[k]=image_name[i+k];
    g_debug_file_name[j]='\0';

#ifdef USE_TIFF
    strncat(g_debug_file_name,"_box.tif",8);
#else
    strncat(g_debug_file_name,"_box.pgm",8);
#endif
    fprintf(stderr,"\n\nDiagnostic output image:  %s\n\n",g_debug_file_name);

    g_debug_array = make_debug_image(tif, buf, rotation);
  }

  actual_cols=g_test_pattern_xlr-g_test_pattern_xul;
  actual_rows=g_test_pattern_ylr-g_test_pattern_yul;


  /* enough to hold the whole image */
  g_image_array = (unsigned char*) calloc((int)(actual_rows*actual_cols)*g_bytes_per_pixel,sizeof(char));
  if ( g_image_array == NULL ) {
    fprintf(stderr, "\nCan't allocate memory for image array\n\n");
    exit(-1);
  }

  /* As we read-in the test pattern area, we rotate it to a standard 
     format, with a vertical edge; 
     The variables below are set to cause the rotation to occur. These 
     variables are used to index into the storage array as we process each 
     scan line of the input image. */
  switch(rotation) {
  case BOTTOM:
    /* Horizontal - make top rows left cols */
    /* The first pixel of the input image becomes the bottom
       pixel of the first column in the output array. */
    /* For each input row we fill in a storage column from
       bottom to top; working our way from left to right. */
    array_column_start = 0;	/* left to right */
    array_column_delta = 1;	/* ditto */
    /* bottom to top */
    array_row_delta = -1;
    /* bottom most row */
    array_row_start=actual_cols-1;
    break;
  case RIGHT:
    /* normal case - white patch at top - read in straight*/
    array_column_start = 0;
    array_column_delta = 1;
    array_row_delta=1;
    array_row_start=0;
    break;
  default:
    fprintf(stderr,"\nUnknown rotation!!\n\n");
    exit(-1);
  }

  switch(rotation) {
    
  case BOTTOM:
    array_column_index = array_column_start;
    switch(g_bytes_per_pixel) {
    case 1:
      for ( input_row = g_test_pattern_yul;
	    input_row < g_test_pattern_ylr; input_row++) {
	read_scan_line(tif,buf,input_row);
	array_row_index=array_row_start;
	for (i=g_test_pattern_xul;
	     i<g_test_pattern_xlr;i++) {
	  g_image_array[array_column_index +
		       array_row_index *
		       actual_rows] = buf[i]; 
	  array_row_index+=array_row_delta;
	}
	array_column_index+= array_column_delta;
      }
      break;
    default:
      fprintf(stderr,"\nCannot handle more than 1 byte per pixel!!\n\n");
      exit(-1);
    }
    break;
  case RIGHT:
    array_row_index = array_row_start;
    switch(g_bytes_per_pixel) {
    case 1:
      for ( input_row = g_test_pattern_yul;
	    input_row < g_test_pattern_ylr; input_row++) {
	read_scan_line(tif,buf,input_row);
	array_column_index=array_column_start;
	for (i=g_test_pattern_xul;
	     i<g_test_pattern_xlr;i++) {
	  g_image_array[array_row_index *
		       actual_cols +
		       array_column_index] =
	    buf[i]; 
	  array_column_index+=array_column_delta;
	}
	array_row_index+=array_row_delta;
      }
      break;
    default:
      fprintf(stderr,"\nCannot handle more than 1 byte per pixel!!\n\n");
      exit(-1);
    }
    break;
  default:
    fprintf(stderr,"\nUnknown rotation!!\n\n");
    exit(-1);
  }

  /* switch the dimensions to match the now rotated image chip */
  switch(rotation) {
  case BOTTOM:
    *size_x = actual_rows;
    *size_y = actual_cols;
    break;
  case RIGHT:
    *size_x = actual_cols;
    *size_y = actual_rows;
    break;
  }

  free(buf);

}


/*****************************************************************************/
/*								      */
/*****************************************************************************/
/* Read a scan line into a buffer */
//把scan line读入到缓存中
void read_scan_line(TIFF *tif,unsigned char *buf,short line_num)
{
  int j;

  if (g_scan_image_file_id) {
    if ( raw_readscanline(buf,line_num) < 0 ) {
      fprintf(stderr,
	      "\nBad data read on imageline %d\n\n",line_num);
      exit(-1);
    }
  }
#ifdef USE_TIFF
  else {
    if ( TIFFReadScanline( tif, buf, line_num, 0L ) < 0 ) {
      fprintf(stderr,
	      "\nBad data read on imageline %d\n\n",line_num);
      exit(-1);
    }
  }
#endif

  /* IF the tif tag: PHOTOMETRICINTERPRETATION was read as 'white is zero'; 
     Then the polarity needs to be reversed.  Also, the user may have 
     requested a forced polarity reversal.  This can be useful when reading
     RAW images on little-endian or big-endian systems, or when a TIFF
     image has an incorrectly marked photometric interpretation. 
     Two polarity reversals (tiff photometric=0 + user request) results in 
     no reversal required. 
  */

  if ( (g_reversepolarity || (g_photometric==0)) &&
      !(g_reversepolarity && (g_photometric==0)) ) {
      
      switch (g_bytes_per_pixel) {
      case 1:
	for (j=0;j<(int)g_scan_cols;j++)
	  buf[j] = g_max_pixel_value-buf[j];
	break;
      default:
	fprintf(stderr,"\nCannot handle more than 1 byte per pixel!!\n\n");
	exit(-1);
      }
  }
}

/*****************************************************************************/
/*								             */
/*****************************************************************************/

/* Read in a bit from the 4 corners of the ROI to determine orientation */
void input_area(TIFF *tif, unsigned char *rotation,
		    int rstart, int rend, int cstart, int cend)
{
  short j;
  unsigned char *buf;
  double ave[4];
  double topdiff, botdiff, lftdiff, rgtdiff;
  int extra_height;
  int flag = 0;

  if (g_test_pattern_xul < 0) { g_test_pattern_xul = 0; flag = 1; }
  if (g_test_pattern_yul < 0) { g_test_pattern_yul = 0; flag = 1; }
  if (g_test_pattern_xlr > g_scan_cols) { g_test_pattern_xlr = g_scan_cols; flag = 1; }
  if (g_test_pattern_ylr > g_scan_rows) { g_test_pattern_ylr = g_scan_rows; flag = 1; }

  extra_height = g_test_pattern_ylr-g_test_pattern_yul - 300;//y轴的大小计算，缩小到300
  if (extra_height > 0) {
    flag = 2;
    //假如大小大于300，就会取大小为300的中间部分
    g_test_pattern_yul += (extra_height+1)/2;
    g_test_pattern_ylr = g_test_pattern_yul + 300;
  }
  extra_height = g_test_pattern_xlr-g_test_pattern_xul - 300;//x轴的大小计算，缩小到300
  if (extra_height > 0) {
    flag = 2;
    g_test_pattern_xul += (extra_height+1)/2;
    g_test_pattern_xlr = g_test_pattern_xul + 300;
  }

  if( flag ) {//返回缩小后的坐标点
    rstart = g_test_pattern_yul;
    rend = g_test_pattern_ylr;
    cstart = g_test_pattern_xul;
    cend = g_test_pattern_xlr;
  }

  if (flag == 1) {
    fprintf(stderr, "ROI center/dims changed to stay within the image.\n");
    fprintf(stderr, "  New center (%d,%d)  dims %d x %d\n",
            (g_test_pattern_xlr+g_test_pattern_xul)/2,   //中心x
            (g_test_pattern_ylr+g_test_pattern_yul)/2,   //中心y
            (g_test_pattern_xlr-g_test_pattern_xul),     //ROI x长度
            (g_test_pattern_ylr-g_test_pattern_yul));    //ROI y长度
  }
  if (flag == 2) {
    fprintf(stderr, "ROI dimension reduced to 300.\n");
  }



  /* enough to hold a single scan line */
  buf = NULL;
  if (g_scan_image_file_id) {//正确读取了文件
    buf = (unsigned char*) calloc((unsigned int)g_scan_cols,sizeof(unsigned char));//分配内存buff
  }
#ifdef USE_TIFF
  else {
    buf = (unsigned char*) calloc((unsigned int)TIFFScanlineSize(tif),sizeof(unsigned char));
  }
#endif
  if ( buf == NULL ) {
    fprintf(stderr,
            "\nCan't allocate memory for scanline buffer\n\n"); //没办法分配内存
    exit(-1);
  }

  for(j=0;j<4;j++) 
    ave[j] = 0.0;

  for(j=0;j<2;j++) {
    read_scan_line(tif, buf, rstart+j);
    switch (g_bytes_per_pixel) {
    case 1:
      ave[0] += (double)(buf[cstart] + buf[cstart+1]);
      ave[1] += (double)(buf[cend-1] + buf[cend-2]);
      break;
    default:
      fprintf(stderr,"\nCannot handle more than 1 byte per pixel!!\n\n");
      exit(-1);
    }
    read_scan_line(tif, buf, rend-j-1);
    switch (g_bytes_per_pixel) {
    case 1:
      ave[2] += (double)(buf[cstart] + buf[cstart+1]);
      ave[3] += (double)(buf[cend-1] + buf[cend-2]);
      break;
    default:
      fprintf(stderr,"\nCannot handle more than 1 byte per pixel!!\n\n");
      exit(-1);
    }
  }
  
	free(buf);
	
  for(j=0;j<4;j++) 
    ave[j] /= 4.0;

  topdiff = (ave[1] - ave[0])/(ave[1]+ave[0]);
  botdiff = (ave[3] - ave[2])/(ave[3]+ave[2]);
  lftdiff = (ave[2] - ave[0])/(ave[2]+ave[0]);
  rgtdiff = (ave[3] - ave[1])/(ave[3]+ave[1]);

  if(cend-cstart <= rend-rstart) {//判断是水平边缘还是垂直边缘
    /* Should be a vertical edge */
    if (fabs(topdiff+botdiff) < fabs(lftdiff+rgtdiff)) {
      MTFPRINT("** ERROR:  Edge does not appear to be 'vertical'.\n")
      MTFPRINT("           Rerun with ROI where width > height\n")
      exit (-1);
    }
    if (fabs(topdiff) < 0.2 || fabs(botdiff) < 0.2) 
      MTFPRINT("Warning: Low edge contrast (< 20%% modulation)\n")
    if (topdiff*botdiff > 0) {
      *rotation = RIGHT; 
      if( (g_test_pattern_xlr-g_test_pattern_xul)%2 != 0) {
	MTFPRINT("** Pixels across edge reduced to an even number **\n")
	g_test_pattern_xlr--;
      }
    }
    else {
      MTFPRINT("** ERROR: Can't determine edge orientation\n")
      exit(-1);
    }
  }
  else {
    /* rotate 90 for horizontal edge */
    if (fabs(topdiff+botdiff) > fabs(lftdiff+rgtdiff)) {
      MTFPRINT("** ERROR:  Edge does not appear to be 'horizontal'.\n")
      MTFPRINT("           Rerun with ROI where width < height\n")
      exit (-1);
    }
    if (fabs(lftdiff) < 0.2 || fabs(rgtdiff) < 0.2) 
      MTFPRINT("Warning: Low edge contrast (< 20%% modulation)\n")
    if (lftdiff*rgtdiff > 0) {
      *rotation = BOTTOM;
      if( extra_height%2 != 0) { 
	MTFPRINT("** Pixels across edge reduced to an even number **\n")
	g_test_pattern_ylr--;
      }
    }
    else {
      MTFPRINT("** ERROR: Can't determine edge orientation\n")
      exit(-1);
    }
  }
}

/*****************************************************************************/
/*									     */
/*****************************************************************************/

/* Read in a scan line from a raw image file */
//从原始图像中读取出一行scan line
short raw_readscanline(unsigned char *buf, short i)
{
  int pos;

  if ((unsigned)i > g_scan_rows-1) return(-1);
  if (g_bytes_per_pixel != 1) return(-1);
  pos = g_scan_header + i*g_scan_cols; //确定读出位置
  if (lseek(g_scan_image_file_id,pos,0) < 0) return(-1);//lseek函数的作用是用来重新定位文件读写的位移
  if (read(g_scan_image_file_id,buf,(int)g_scan_cols) != (int)g_scan_cols)
    return(-1);
  else return(0);
}

/*****************************************************************************/
/* Check image for clipping */
//检查图像的剪辑情况
//计算出图像中像素值大于255，小于0的范围的比例
int clipping (int low, int high, double thresh, unsigned char *img, int len) 
{
  int i;
  unsigned char *buf;
  int nlow=0;
  int nhigh=0;
  int status=0;
 //g_bytes_per_pixel 的作用是干啥？？
  switch(g_bytes_per_pixel) {
    case 1:
      buf = img;
      for ( i=0; i<len; i++, buf++) {
	if((int)(*buf) <= low)
	  nlow++;
	else if ((int)(*buf) >= high)
	  nhigh++;
      }
      break;
    default:
      fprintf(stderr,"\nCannot handle more than 1 byte per pixel!!\n\n");
      exit(-1);
  }

  if ( nlow/(double)len > thresh ) {
    MTFPRINT2("Warning: %.1f%% low clipping\n", nlow/(double)len*100.)
    status = 1;
  }
  if ( nhigh/(double)len > thresh ) {
    MTFPRINT2("Warning: %.1f%% high clipping\n", nhigh/(double)len*100.)
    status = 1;
  }

  return(status);
}

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
typedef struct Grey_Reflectance_Pair//灰度映射对
{ 
  double grey;
  double reflectance;
} grey_reflectance_pair;

//比较灰度大小
static int compare_grey_reflectance_pairs(const void *first, const void *second)
{

  if (((grey_reflectance_pair *)first)->grey < ((grey_reflectance_pair *)second)
      ->grey)
    return -1;

  if (((grey_reflectance_pair *)first)->grey > ((grey_reflectance_pair *)second)
      ->grey)
    return 1;

  return 0;

}


double *read_in_ref_lut(char *data_filename, int *lowest_value, int *highest_value) {
  FILE *fd;
  int i, j, tot;
  double *ref_lut;
  double ref, grey, m;
  char str[502];
  grey_reflectance_pair *grey_ref_pairs;
  int ncflag, comment_lines;

  ref_lut = (double *)malloc((g_max_pixel_value + 1)*sizeof(double));

  if (strcmp(data_filename, "linear") == 0) {
    /* Don't bother trying to read a file, just do the obvious. */

   //直接设置最高值和最低值
    *lowest_value = 0; 
    *highest_value = 255;
    for (i=0; i<=g_max_pixel_value; i++)
      ref_lut[i] = (double)i/(double)g_max_pixel_value;//除以最大的像素值
  }
  else {
    fd = fopen(data_filename, "r"); //打开图片
    if (fd == NULL) {
      fprintf(stderr,"\nCould not open %s\n",data_filename);
      fprintf(stderr,"%s\n\n", strerror(errno));
      return(0);
    }

    //记录注释行数
    i=0; comment_lines=0; ncflag = 0;
    while ( fgets(str, 500, fd) != NULL ) {
      if(ncflag==0 && str[0] == '#') comment_lines++; //当读取的行开头为#时表示这是一个注释行
      else {
	ncflag = 1;
	i++;
      }
    }
    tot = i;
    fseek(fd, 0, SEEK_SET);

    grey_ref_pairs = (grey_reflectance_pair *)calloc(tot,sizeof(grey_reflectance_pair));

    /* Skip past comment lines */
    //跳过注释行
    for(i=0; i<comment_lines; i++) fgets(str, 500, fd);
    i = 0;
    while ( fscanf(fd, "%lf %lf", &ref, &grey) == 2  && i<tot) {
      grey_ref_pairs[i].reflectance = ref;
      grey_ref_pairs[i].grey = grey;
      i++;
    }
    tot = i;
    
    /* Sort grey_reflectance_pairs by grey values, low to high */
    //用灰度值从低到高储存灰度映射对
    qsort(grey_ref_pairs, tot, sizeof(grey_reflectance_pair), 
	  compare_grey_reflectance_pairs);
   
    j=0;
    i = ceil(grey_ref_pairs[0].grey);
    if (i<0) i=0;
    if (i>g_max_pixel_value) i=g_max_pixel_value+1;
    *lowest_value = i;
    /* Extrapolate at lower end if needed */
    if (i>0 && ceil(grey_ref_pairs[0].grey) != floor(grey_ref_pairs[0].grey)) {
      *lowest_value -= 1;
      m = (grey_ref_pairs[1].reflectance-grey_ref_pairs[0].reflectance)
	    /(grey_ref_pairs[1].grey-grey_ref_pairs[0].grey);
      ref_lut[i-1] = ((double)i-1-grey_ref_pairs[0].grey)*m 
	             + grey_ref_pairs[0].reflectance;
    }
    for( ; i<grey_ref_pairs[tot-1].grey && i<=g_max_pixel_value;i++) {
      while((double)i > grey_ref_pairs[j].grey && j<tot) j++;
      if (j>=tot) break;
      if((double)i == grey_ref_pairs[j].grey) 
	ref_lut[i] = grey_ref_pairs[j].reflectance;
      else {
	m = (grey_ref_pairs[j].reflectance-grey_ref_pairs[j-1].reflectance)
	    /(grey_ref_pairs[j].grey-grey_ref_pairs[j-1].grey);
	ref_lut[i] = ((double)i-grey_ref_pairs[j-1].grey)*m 
	             + grey_ref_pairs[j-1].reflectance;
      }
    }
    if ( (i<=g_max_pixel_value)  ) {
      if((double) i == grey_ref_pairs[tot-1].grey) {
	ref_lut[i] = grey_ref_pairs[tot-1].reflectance;
      }
      else {
	m = (grey_ref_pairs[tot-1].reflectance-grey_ref_pairs[tot-2].reflectance)
	    /(grey_ref_pairs[tot-1].grey-grey_ref_pairs[tot-2].grey);
	ref_lut[i] = ((double)i-grey_ref_pairs[tot-1].grey)*m 
	             + grey_ref_pairs[tot-1].reflectance;
      }
      i++;
    }
    *highest_value = i-1;
  }

  return (ref_lut);
}

//翻转LUT得到gray图片
int reverse_lut(double *ref_lut, double val, int lowest_val, int highest_val) {
  /* Assumes monotonic OECF curve */
  int grey;
  int sgn = 1;

  if (ref_lut[lowest_val] > ref_lut[highest_val]) sgn = -1;
  val = sgn*val;
  for (grey=lowest_val; grey<highest_val; grey++) {
    if (val <= sgn*(ref_lut[grey]+ref_lut[grey+1])/2.0) break;
  }

  return grey;
}
 
int read_pgm_header(char *image_filename)
{
  FILE *fp;
  int val, rows, cols;
  char header[256];
  int numread;

  fp = fopen(image_filename, "rb");
  if (fp == NULL) return (-2);
   
  do {
    if (fgets(header,255,fp) == NULL)
      return (-1);
  } while (header[0] == '#');
  
  if ( strncmp(header, "P5", 2) != 0 ) return (-1);
  
  do {
    if (fgets(header,255,fp) == NULL)
      return (-1);
  } while (header[0] == '#');

  numread = sscanf(header,"%d %d", &cols, &rows);
  g_scan_cols = cols;
  if (numread != 2) {  
    /*Only one number on this line, so read the next line to get height. */
    fgets(header,255,fp);
    sscanf(header,"%d", &rows);
  }
  g_scan_rows = rows;

  fgets(header,255,fp);
  sscanf(header,"%d", &val);
  if(val <= 255)  
    g_bps = 8;

  g_scan_header = ftell(fp);

  fclose(fp);
  g_pgm = 1;
  g_photometric = 1;

  return(0);  
}
/*****************************************************************************/
/*                                                                           */
/*****************************************************************************/
static void write_debug_image(unsigned char *array, int width, int height)
{
  int i;

#ifdef USE_TIFF
  TIFF *dtif;
  if((dtif = TIFFOpen(g_debug_file_name, "w")) == NULL) {
    fprintf(stderr, "Could not create %s - try running without 'd' option\n", g_debug_file_name);
    exit(-1);
  }
  fprintf(stderr, "\nWriting out diagnostic image (TIF): %s    ", g_debug_file_name);
#else
  FILE *g_ofp;
  g_ofp = fopen(g_debug_file_name, "wb");
  if (g_ofp == NULL) {
    fprintf(stderr,"\n%s%s%s\n","Could not create ",g_debug_file_name,
	    " - try running without 'd' option \n");
    exit(-1);
  }
  fprintf(stderr, "\nWriting out diagnostic image (PGM): %s    ", g_debug_file_name);
#endif

  fflush(stderr);

#ifdef USE_TIFF
  TIFFSetField(dtif, TIFFTAG_IMAGEWIDTH, width);
  TIFFSetField(dtif, TIFFTAG_IMAGELENGTH, height);
  TIFFSetField(dtif, TIFFTAG_SAMPLESPERPIXEL, 1);
  TIFFSetField(dtif, TIFFTAG_BITSPERSAMPLE, 8);
  TIFFSetField(dtif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
  TIFFSetField(dtif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(dtif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
  TIFFSetField(dtif, TIFFTAG_COMPRESSION, COMPRESSION_NONE); 
  i = 8192 / TIFFScanlineSize(dtif);
  TIFFSetField(dtif, TIFFTAG_ROWSPERSTRIP, i==0?1:i);
#else 
  fprintf(g_ofp,"P5\n%d %d\n255\n",width,height);
#endif

  for(i=0; i<height; i++) {  /* Write each row separately */
    /* Now write the line into the debug image file */
#ifdef USE_TIFF
    if( TIFFWriteScanline( dtif, array+i*width, i, 0L ) < 0 ) 
      fprintf(stderr, "\nBad data write on imageline %d\n\n",i);
#else
    fwrite(array+i*width, 1, width, g_ofp);
#endif
  }


#ifdef USE_TIFF
  TIFFClose(dtif);
#else
  if (fclose(g_ofp) != 0)
    fprintf(stderr, "sfr.c: Could not close g_ofp: %s\n\n", strerror(errno));
#endif

  fprintf(stderr, "Done\n\n");
}

static void draw_lines(unsigned char *array, double slope, 
		      unsigned char rotation, 
		      int center, int size_y, double off, int size_x)
{
  int i, j;
  int full_width;
  int begin, row_inc, col_inc;
  unsigned char *start_ptr;
  int indexa, indexb, indexc, indexd;
  double position;

  full_width = g_debug_fullwidth;

  if (rotation == RIGHT) {
    begin = (g_test_pattern_ylr-g_test_pattern_yul)/2 - size_y/2;
    row_inc = full_width;
    col_inc = 1;
  }
  else {
    begin = (g_test_pattern_xlr-g_test_pattern_xul)/2 - size_y/2;
    row_inc = -1;
    col_inc = full_width;
  }
  
  start_ptr = g_start_ptr;
  start_ptr += begin * row_inc;

  for(i=-size_y/2,j=0; j<size_y; j++,i++) {  
    position = (slope*i + off)*4.0 + center;
    if ((g_version & 4) == 0)
      position -= 0.5;
    /* Add 2 to move to center of 4 duplicate pixels. */
    indexa = (int)( floor(position + 2 + 0.25 + 0.5));
    indexb = (int)( floor(position + 2 - 0.25 + 0.5 ));
    start_ptr[indexa*col_inc] += 128;
    if(indexb!=indexa) start_ptr[indexb*col_inc] += 128;

    indexc = (int)( floor(position + 2 - center) - 1 );
    indexd = (int)( ceil(position + 2 + center) + 1 );
    if(j%4 == 0 || j%4 == 1) { /* Dotted line */
      if(indexc >= 0)
	start_ptr[indexc*col_inc] ^= 255;
      if(indexd < size_x*4)
	start_ptr[indexd*col_inc] ^= 255;
    }
    start_ptr += row_inc;
  }

}

    
static unsigned char *make_debug_image(TIFF *tif, unsigned char *buf, 
				       unsigned char rotation) 
{
  int i,j,k,jj, cnt;
  int height,width, full_height, full_width;
  int offset_x, offset_y;
  unsigned char *debug_array, *ptr;

  width = g_scan_cols;
  height = g_scan_rows;
  if (width > g_ppi*1.6) width = g_ppi*1.6;
  if (height > g_ppi*1.5) height = g_ppi*1.5;

  full_width = width;
  full_height = height;
  if(rotation == RIGHT) 
    full_width += 10+4*(g_test_pattern_xlr-g_test_pattern_xul);
  else
    full_height += 10+4*(g_test_pattern_ylr-g_test_pattern_yul+1);
    
  debug_array = (unsigned char *)calloc(full_height*full_width,
					sizeof(unsigned char));

  offset_x = offset_y = 0;
  if(g_test_pattern_xlr > width) offset_x = g_test_pattern_xlr-width;
  if(g_test_pattern_ylr > height) offset_y = g_test_pattern_ylr-height;

  if(rotation == RIGHT) 
    g_start_ptr = debug_array+(g_test_pattern_yul-offset_y)*full_width
                  + width+10;
  else
    g_start_ptr = debug_array+(height+10)*full_width + (g_test_pattern_xul-offset_x);

  cnt = 0;
  for(i=offset_y,k=0; k<height; i++,k++) {

	  read_scan_line(tif, buf,i);

    if ( i>=g_test_pattern_yul && i<=g_test_pattern_ylr) {
      if (rotation == RIGHT) {
	ptr = g_start_ptr+cnt*full_width;
	for(jj=0, j=g_test_pattern_xul;j<g_test_pattern_xlr;j++,jj+=4) 
	  memset(ptr+jj, (int)(buf[j]), 4);
	cnt++;
      }
      else
		  for(j=0;j<4;j++,cnt++) {
	  memcpy((g_start_ptr + cnt*full_width), 
		  (buf+g_test_pattern_xul), 
		  g_test_pattern_xlr-g_test_pattern_xul);
		}
    }
    /* Modify with box if needed */
    if (i==g_test_pattern_yul || i==(g_test_pattern_ylr-1)) {
      for(j=g_test_pattern_xul; j<g_test_pattern_xlr; j++)
	buf[j] += 128;
    }
    else if ( i>g_test_pattern_yul && i<g_test_pattern_ylr) {
      buf[g_test_pattern_xul] += 128;
      buf[g_test_pattern_xlr-1] += 128;
    }

    if (g_userangle) {
      
      if ( fabs(i-g_pt1y) < .6 ) buf[(int)rint(g_pt1x)] += 128;
      if ( fabs(i-g_pt2y) < .6 ) buf[(int)rint(g_pt2x)] += 128;
    }

    memcpy((debug_array+k*full_width), (buf+offset_x), width);
  }

  g_debug_fullwidth = full_width;
  g_debug_fullheight = full_height;
  
  if(rotation != RIGHT) 
    g_start_ptr += g_test_pattern_xlr - g_test_pattern_xul - 1;


  return debug_array;
}

void print_data_rights_notice(void)
{
  fprintf(stderr, "\n\n");
  fprintf(stderr,"*******************************************************************************\n");
  fprintf(stderr, "                               PROGRAM NOTICE\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "This software is a modification of the \"SFR Measurement Algorithm C-code\",\n");
  fprintf(stderr, "copyright PIMA 1998, appearing in Annex A of ISO standard 12233-2000, \n");
  fprintf(stderr, "\"Photography - Electronic Still Picture Cameras - Resolution Measurement\".\n");
  fprintf(stderr, "Permission to modify this software, and use and distribute the modified \n");
  fprintf(stderr, "software, was granted by I3A (the successor organization of PIMA) to the \n");
  fprintf(stderr, "MITRE Corporation in 2006.\n");
  fprintf(stderr," \n");
  fprintf(stderr, "This MITRE Corporation-modified SFR software was produced for the U.S. \n");
  fprintf(stderr, "Government under Contract numbers J-FBI-12-128 and W15P7T-05-C-F600 and is\n");
  fprintf(stderr, "subject to the Rights in Data-General clause FAR 52.227-l4, Alt. IV (DEC 2007).\n");
  fprintf(stderr," \n");
  fprintf(stderr, "Redistribution and use in source and binary forms, with or without modification,\n");
  fprintf(stderr, "are permitted provided that the following conditions are met:\n");
  fprintf(stderr, "-Redistribution of source code must retain the above copyright notice, this list\n");
  fprintf(stderr, "of conditions, and the following disclaimer.\n");
  fprintf(stderr, "-Redistribution in binary form must reproduce the above copyright notice, this\n");
  fprintf(stderr, "list of conditions and the following disclaimer in the documentation and/or\n");
  fprintf(stderr, "other materials provided with the distribution.\n");
  fprintf(stderr," \n");
  fprintf(stderr, "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\"\n");
  fprintf(stderr, "AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE\n");
  fprintf(stderr, "THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITTNESS FOR A PARTICULAR PURPOSE\n");
  fprintf(stderr, "ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE\n"); 
  fprintf(stderr, "LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR\n");
  fprintf(stderr, "CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF\n");
  fprintf(stderr, "SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS\n");
  fprintf(stderr, "INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN \n");
  fprintf(stderr, "CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) \n");
  fprintf(stderr, "ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE \n");
  fprintf(stderr, "POSSIBILITY OF SUCH DAMAGE. \n");
  fprintf(stderr," \n");
  fprintf(stderr,"*******************************************************************************\n");
  fprintf(stderr, "Code modifications by Margaret A. Lepley (MITRE), mlepley@mitre.org \n");
  fprintf(stderr,"*******************************************************************************\n");
  fprintf(stderr," \n");
}

void wait_to_exit(void) {

  fprintf(stderr, "Output appended to %s\n", MTFOUTNAME );

  while ( getchar() != (int)'\n');
  fprintf(stderr, "Press any key to exit  ...");
  fflush(stderr);
  while (getchar() == EOF);
  fprintf(stderr, "\n");
 
}


/******************************************************************************/
/*                                                                            */
/******************************************************************************/

void put_problem(char *s, int problem_type)
{

  if (g_problem_count < MAX_PROBLEMS) {
    g_problem_type[g_problem_count] = (short)problem_type;
    strcpy(&g_problems[g_problem_count++][0],s);
  }
  else
    {
      fprintf(stderr, "Number of problems exceeds count of %d\n", MAX_PROBLEMS)\
;
      too_many_problems = 1;
    }
  if (problem_type == IQS)
    g_IQS_problem_count++;

}


/****************************************************************************/
/*                                                                          */
/****************************************************************************/
/* Prints problems.  */

void print_problems(void)
{
  short i;

  if (g_problem_count == 0) {
    if(g_nocompare)
      MTFPRINT("\nNo spec check applied.")
    else {
      MTFPRINT("\nPIV spec check applied.")
      MTFPRINT("\nNo problems encountered.\n\n")
    }
  }
  else if ((g_nocompare) && g_problem_count-g_IQS_problem_count == 0)
    {
      /* Do nothing when only IQS problems present, but IQS report is
	 supressed. */
    }
  else {
    i=0;

    /* note, "\f" causes a new page */
    MTFPRINT("\n\t\t\tPROBLEM REPORT\n\n")
 
     if(g_nocompare)
       MTFPRINT("No spec check applied.\n")
     else
       MTFPRINT("PIV spec check applied.\n")

     while (g_problem_count--) {
       /* Print when IQS comparison is requested *or* when non-IQS problem */
       if ((!g_nocompare) || g_problem_type[i] != IQS)
	 MTFPRINT2("%s",g_problems[i])
       i++;
     }

    MTFPRINT("\n")

    if (too_many_problems != 0) {
      MTFPRINT("\nMore problems were found, but number exceeds program limits\
.\n")
      too_many_problems = 0;
    }
  }
}


void slope_bounds( double slope, int size_y, int numcycles, double mincyc, char *problem_string)
{
  double absslope;

  absslope = fabs(slope);

  /* If the slope is too small not enough rows for mincy (typically 5) 
     full periods, then alert the user. */
  if ( absslope < mincyc/(double)(size_y) ){
    sprintf(problem_string, "SFR suspect due to shallow edge angle (%f).  Need\n",
	    atan(slope)*180/M_PI );
    put_problem(problem_string, 0);
    sprintf(problem_string, "  %d lines of data (%.1f cycles) for accurate results\n",
	      (int)ceil(mincyc/absslope), mincyc);
    put_problem(problem_string, 0);
    return;
  }

  if( absslope*(double)(size_y) < 1.0) {
    sprintf(problem_string, "Not even one complete phase cycle included in ROI\n");
    put_problem(problem_string, 0);
    return;
  }

  if( absslope > (double)(1.0/4.0) ) { 
      int rows_per_col;
      double bestcycle, x;

      if( absslope > (double)(5.0/4.0) ) {
	sprintf(problem_string, "ERROR: Edge angle (%f) is too large\n",
	      atan(slope)*180/M_PI);
	put_problem(problem_string, 0);
	return;
      }

      rows_per_col = (int)floor(1/absslope + 0.5);
      x = fabs(1/(double)rows_per_col - absslope);
      bestcycle = 4*rows_per_col*x*ceil(1.0/x/(double)rows_per_col/(double)rows_per_col - 1.0);
      if((int)ceil(mincyc*bestcycle) > size_y) {
	  sprintf(problem_string, "SFR suspect due to large edge angle (%f). Need\n",
	      atan(slope)*180/M_PI);
	  put_problem(problem_string, 0);
	  sprintf(problem_string, "  %g * %f = %d lines of data for accurate results\n\n",
	      mincyc, bestcycle,
	      (int)ceil(mincyc*bestcycle));
	  put_problem(problem_string, 0);
      }
  }
}

