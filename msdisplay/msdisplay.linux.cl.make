##################################################################################
#                                                                                #
#  Program    : msdisplay                                                        #
#                                                                                #
#  Filename   : msdisplay.linux.cl.make                                          #
#                                                                                #
#  Created    : May 16th 2007                                                    #
#                                                                                #
#  Purpose    : LINUX makefile for msdisplay.                                    #
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
INCLUDEDIRS=-I. -I../include
LIBDIRS=-L../lib
LIBS=-lraw -lucsf -lgen -lnrec -lm -lexpat -lz

all:
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -c $(INCLUDEDIRS) msdisplay_main.cpp -o msdisplay_main.o
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -c $(INCLUDEDIRS) msdisplay_params.cpp -o msdisplay_params.o
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -o ../bin/msdisplay.cgi msdisplay_main.o msdisplay_params.o $(LIBDIRS) $(LIBS) $(STATIC)
