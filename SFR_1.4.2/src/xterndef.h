/*****************************************************************************
xterndef.h

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
This file pertains to the wrapper code and is written by MITRE.
Code predominantly modified from IS 12233:2000 Annex A is in sfr_iso.c.
------------------------------------------
*****************************************************************************/
/******************************************************************************/
/*              arrays and file pointers                                      */
/******************************************************************************/
FILE *g_mtfout;
char g_debug_file_name[80];
char image_filename[80];
char data_filename[80];
int g_scan_image_file_id;
int g_ofd; 

unsigned char *g_image_array;

/******************************************************************************/
/*      image, test pattern and patch locations and dimensions                */
/******************************************************************************/

short g_input_column_start;
short g_input_last_column;
short g_input_row_start;
short g_input_last_row;
short g_input_width;
short g_input_height;

unsigned int g_total_image_width;
unsigned int g_total_image_height;
unsigned int g_bytes_per_pixel;
unsigned int g_max_pixel_value;

int g_test_pattern_xul;
int g_test_pattern_yul;
int g_test_pattern_xlr;
int g_test_pattern_ylr;
int g_target_res;

double g_target_width_in_mm;
double g_target_height_in_mm;

double g_target_width;
double g_target_height;

char g_problems[MAX_PROBLEMS][MAX_LINE_LENGTH];
short g_problem_count;
char g_problem_type[MAX_PROBLEMS];
short g_IQS_problem_count;


/******************************************************************************/
/*              things read in                                                */
/******************************************************************************/

unsigned int g_scan_rows;
unsigned int g_scan_cols;
unsigned int g_scan_header;
unsigned short g_bps;
unsigned short g_photometric;

double g_ppi;

int g_version;
unsigned char g_center;
unsigned char g_debug;
unsigned char g_streamline_args;
unsigned char g_reversepolarity;
unsigned char g_extended;
unsigned char g_pgm;
unsigned char g_autorefine;
unsigned char g_nocompare;

