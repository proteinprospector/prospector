##################################################################################
#                                                                                #
#  Program    : jobStatus                                                        #
#                                                                                #
#  Filename   : jobStatus.linux.make                                             #
#                                                                                #
#  Created    : May 16th 2007                                                    #
#                                                                                #
#  Purpose    : LINUX makefile for jobStatus.                                    #
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
LIBS=-lucsf -ldbase -lgen -lnrec -lm -lexpat -lmysqlclient

all:
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -c $(INCLUDEDIRS) job_status_main.cpp -o job_status_main.o
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -o ../bin/jobStatus.cgi job_status_main.o $(LIBDIRS) $(LIBS) $(STATIC)
