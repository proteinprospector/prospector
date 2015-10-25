##################################################################################
#                                                                                #
#  Program    : faindex                                                          #
#                                                                                #
#  Filename   : faindex.linux.cl.make                                            #
#                                                                                #
#  Created    : September 11th 2001                                              #
#                                                                                #
#  Purpose    : LINUX makefile for faindex.                                      #
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

COMPILER=g++
OPTIONS=-O2 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE
ADD_OPTIONS=
STATIC=
INCLUDEDIRS=-I../include
LIBDIRS=-L../lib
LIBS=-lucsf -lsingle -lgen -lnrec -lm -lexpat -lz

all:
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -c $(INCLUDEDIRS) faindex_main.cpp -o faindex_main.o
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -o ../bin/faindex.cgi faindex_main.o $(LIBDIRS) $(LIBS) $(STATIC)
