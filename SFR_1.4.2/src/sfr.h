
/*****************************************************************************
sfr.h

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
---------------------
*****************************************************************************/

/*****************************************************************************/
/*                                                                           */
/* Revisions to sfr.h                                                        */
/*                                                                           */
/* V1.2 PIV spec check added as default, plus option to avoid the check      */
/*                                                           mal 10/06       */
/*                                                                           */
/*****************************************************************************/

#define VERSION "1.4.2"

#define PROG_NAME "MITRE_SFR" /* for DOS */

#ifndef M_PI
#define M_PI  3.14159265358979323846
#endif

extern char *strerror(int);

#define TOP 0
#define BOTTOM 2
#define LEFT 1
#define RIGHT 3
#define UNKNOWN 255

#define GET_COORD_PAIR(a,b,c)                                         \
{                                                                       \
                int i;                                                  \
                scanf("%u",a);                                          \
                while(!isdigit(i=getchar())){}                           \
                ungetc(i,stdin);                                        \
                scanf("%u",b);                                          \
                while(!isdigit(i=getchar())){}                           \
                ungetc(i,stdin);                                        \
                scanf("%u",c);                                          \
}

#define GET_INT_PAIR(a,b)                                             \
{                                                                       \
                int i;                                                  \
                scanf("%d",a);                                          \
                while(!isdigit(i=getchar())){}                           \
                ungetc(i,stdin);                                        \
                scanf("%d",b);                                          \
}

#define GET_FLOAT_PAIR(a,b)                                             \
{                                                                       \
                int i;                                                  \
                scanf("%f",a);                                          \
                while(!isdigit(i=getchar())){}                           \
                ungetc(i,stdin);                                        \
                scanf("%f",b);                                          \
}

#define MTFPRINT(a)                                                     \
{                                                                       \
        fprintf(stdout,a);                                              \
        fprintf(g_mtfout,a);                                            \
}

#define MTFPRINT2(a,b)                                                  \
{                                                                       \
        fprintf(stdout,a,b);                                            \
        fprintf(g_mtfout,a,b);                                          \
}

#define MTFPRINT3(a,b,c)                                                \
{                                                                       \
        fprintf(stdout,a,b,c);                                          \
        fprintf(g_mtfout,a,b,c);                                        \
}

#define MTFPRINT4(a,b,c,d)                                              \
{                                                                       \
        fprintf(stdout,a,b,c,d);                                        \
        fprintf(g_mtfout,a,b,c,d);                                      \
}

#define MTFPRINT5(a,b,c,d,e)                                            \
{                                                                       \
        fprintf(stdout,a,b,c,d,e);                                      \
        fprintf(g_mtfout,a,b,c,d,e);                                    \
}

#define MTFPRINT6(a,b,c,d,e,f)                                          \
{                                                                       \
        fprintf(stdout,a,b,c,d,e,f);                                    \
        fprintf(g_mtfout,a,b,c,d,e,f);                                  \
}

#define MTFPRINT7(a,b,c,d,e,f,g)                                        \
{                                                                       \
        fprintf(stdout,a,b,c,d,e,f,g);                                  \
        fprintf(g_mtfout,a,b,c,d,e,f,g);                                \
}

#define MTFPRINT8(a,b,c,d,e,f,g,h)                                      \
{                                                                       \
        fprintf(stdout,a,b,c,d,e,f,g,h);                                \
        fprintf(g_mtfout,a,b,c,d,e,f,g,h);                              \
}

#define MTFPRINT9(a,b,c,d,e,f,g,h,i)                                    \
{                                                                       \
        fprintf(stdout,a,b,c,d,e,f,g,h,i);                              \
        fprintf(g_mtfout,a,b,c,d,e,f,g,h,i);                            \
}

#define MTFPRINT10(a,b,c,d,e,f,g,h,i,j)                                 \
{                                                                       \
        fprintf(stdout,a,b,c,d,e,f,g,h,i,j);                            \
        fprintf(g_mtfout,a,b,c,d,e,f,g,h,i,j);                          \
}

#define MTFPRINT12(a,b,c,d,e,f,g,h,i,j,k,l)                             \
{                                                                       \
        fprintf(stdout,a,b,c,d,e,f,g,h,i,j,k,l);                        \
        fprintf(g_mtfout,a,b,c,d,e,f,g,h,i,j,k,l);                      \
}


#define MTFOUTNAME "SFROUT.txt"


#define MAX_PROBLEMS 60
#define MAX_LINE_LENGTH 132




