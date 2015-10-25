##################################################################################
#                                                                                #
#  Program    : msform                                                           #
#                                                                                #
#  Filename   : msform.linux.cl.make                                             #
#                                                                                #
#  Created    : February 25th 2005                                               #
#                                                                                #
#  Purpose    : LINUX makefile for msform.                                       #
#                                                                                #
#  Author(s)  : Peter Baker                                                      #
#                                                                                #
#  This file is the confidential and proprietary product of The Regents of       #
#  the University of California.  Any unauthorized use, reproduction or          #
#  transfer of this file is strictly prohibited.                                 #
#                                                                                #
#  Copyright (2005-2013) The Regents of the University of California.            #
#                                                                                #
#  All rights reserved.                                                          #
#                                                                                #
##################################################################################

COMPILER=g++
OPTIONS=-O2 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE -DBATCHTAG
ADD_OPTIONS=
STATIC=
INCLUDEDIRS=-I../include
LIBDIRS=-L../lib
LIBS=-lucsf -lgen -lnrec -lm -lexpat -lz

all:
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -c $(INCLUDEDIRS) msform_main.cpp -o msform_main.o
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -o ../bin/msform.cgi msform_main.o $(LIBDIRS) $(LIBS) $(STATIC)
