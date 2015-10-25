##################################################################################
#                                                                                #
#  Library    : libnrec                                                          #
#                                                                                #
#  Filename   : libnrec.linux.make                                               #
#                                                                                #
#  Created    : September 11th 2001                                              #
#                                                                                #
#  Purpose    : LINUX makefile for libnrec.                                      #
#                                                                                #
#  Author(s)  : Peter Baker                                                      #
#                                                                                #
#  This file is the confidential and proprietary product of The Regents of       #
#  the University of California.  Any unauthorized use, reproduction or          #
#  transfer of this file is strictly prohibited.                                 #
#                                                                                #
#  Copyright (2001-2013) The Regents of the University of California.            #
#                                                                                #
#  All rights reserved.                                                          #
#                                                                                #
##################################################################################

LIBRARY_NAME = ../lib/libnrec.a
COMPILER=g++
ARCHIVER=ar
RANLIB=ranlib
OPTIONS= -O2 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE
ADD_OPTIONS=
INCLUDEDIRS=-I../include

ARCHIVE_COMMAND = $(ARCHIVER) r $(LIBRARY_NAME)

OBJ=bico.o \
	factln.o \
	factrl.o \
	fit.o \
	gammln.o \
	gammq.o \
	gcf.o \
	gser.o \
	ln_covsrt.o \
	ln_dot.o \
	ln_gaussj.o \
	ln_mrqcof.o \
	moment.o \
	nrutil.o \
	spline.o \
	zbrent.o

.SUFFIXES:	.cpp

.cpp.o:	
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) $(INCLUDEDIRS) -c $<

libnrec.a:$(OBJ)
	$(ARCHIVE_COMMAND) $(OBJ)
	$(RANLIB) $(LIBRARY_NAME)
clean:	
	/bin/rm -f *.o 	
