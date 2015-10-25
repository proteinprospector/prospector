##################################################################################
#                                                                                #
#  Program    : msbatch                                                          #
#                                                                                #
#  Filename   : msbatch.linux.make                                               #
#                                                                                #
#  Created    : March 6th 2007                                                   #
#                                                                                #
#  Purpose    : LINUX makefile for msbatch.                                      #
#                                                                                #
#  Author(s)  : Peter Baker                                                      #
#                                                                                #
#  This file is the confidential and proprietary product of The Regents of       #
#  the University of California.  Any unauthorized use, reproduction or          #
#  transfer of this file is strictly prohibited.                                 #
#                                                                                #
#  Copyright (2007-2013) The Regents of the University of California.            #
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
LIBS=-lucsf -ldbase -lgen -lnrec -lexpat -lmysqlclient -lz

all:
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -c $(INCLUDEDIRS) msbatch_main.cpp -o msbatch_main.o
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -o ../bin/msbatch.cgi msbatch_main.o $(LIBDIRS) $(LIBS) $(STATIC)
