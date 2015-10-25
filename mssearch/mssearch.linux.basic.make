##################################################################################
#                                                                                #
#  Program    : mssearch                                                         #
#                                                                                #
#  Filename   : mssearch.linux.basic.make                                        #
#                                                                                #
#  Created    : September 11th 2001                                              #
#                                                                                #
#  Purpose    : LINUX makefile for mssearch.                                     #
#                                                                                #
#  Author(s)  : Peter Baker                                                      #
#                                                                                #
#  This file is the confidential and proprietary product of The Regents of       #
#  the University of California.  Any unauthorized use, reproduction or          #
#  transfer of this file is strictly prohibited.                                 #
#                                                                                #
#  Copyright (2001-2015) The Regents of the University of California.            #
#                                                                                #
#  All rights reserved.                                                          #
#                                                                                #
##################################################################################

COMPILER=g++
OPTIONS=-O2 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE
ADD_OPTIONS=
STATIC=
LLIB=
INCLUDEDIRS=-I../include
LIBDIRS=-L../lib
LIBS=-lucsf -lsingle -lgen -lnrec -lzip -lm -lexpat -lz -lpthread -lsqlite $(LLIB)

all:
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -c $(INCLUDEDIRS) mssearch_main.cpp -o mssearch_main.o
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -o ../bin/mssearch.cgi mssearch_main.o $(LIBDIRS) $(LIBS) $(STATIC)
