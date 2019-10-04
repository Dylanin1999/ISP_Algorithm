#				 Notice
#                                ------
#
# This software is a modification of the "SFR Measurement Algorithm C-code",
# copyright PIMA 1998, appearing in Annex A of ISO standard 12233-2000, 
# "Photography - Electronic Still Picture Cameras - Resolution Measurement".
# Permission to modify this software, and use and distribute the modified 
# software, was granted by I3A (the successor organization of PIMA) to the 
# MITRE Corporation in 2006.
#
# This MITRE Corporation-modified SFR software was produced for the U.S. 
# Government under Contract numbers J-FBI-12-128 and W15P7T-05-C-F600 and is
# subject to the Rights in Data-General clause FAR 52.227-l4, Alt. IV (DEC 2007).
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# -Redistribution of source code must retain the above copyright notice, this
# list of conditions, and the following disclaimer.
# -Redistribution in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 


# To compile a version of this sfr program with the TIFF software:
# - make sure that -DUSE_TIFF is in the CPPFLAGS
# - type "make"
# - executable file will be named mitre_sfr

# A non-tiff capable version can be compiled using 'make sfrraw'
# but in either case the executable will be names mitre_sfr

# uncomment these lines for gcc
CC = gcc
CFLAGS = -O3 -Wall -pedantic 
#CFLAGS = -g -Wall -pedantic 

# Remove -DUSE_TIFF to not use tiff software
CPPFLAGS = \
	-DUSE_TIFF \
	-DUSE_VARARGS=0 \
	-DUSE_PROTOTYPES=1 \
	-DHAVE_IEEEFP=1 \
	-DCOLORIMETRY_SUPPORT \
	-DYCBCR_SUPPORT 
	
# On Windows, with MinGW
LDFLAGS = -static 
# On Mac OSX 
# LDFLAGS =

default:	mitre_sfr


SFR_OBJS = \
	find_area.o \
	sfr_iso.o \
	mitre_sfr.o


${SFR_OBJS}: sfr.h xterndef.h 

mitre_sfr:	$(SFR_OBJS)
	$(LINK.c) -o mitre_sfr $(SFR_OBJS) -ltiff -ljpeg -lz -lm

sfrraw:	$(SFR_OBJS)
	$(LINK.c) -o mitre_sfr $(SFR_OBJS) -lm

clean:
	rm -f *.o core mitre_sfr mitre_sfr.exe
