##################################################################################
#                                                                                #
#  Program    : kill                                                             #
#                                                                                #
#  Filename   : kill.linux.full.make                                             #
#                                                                                #
#  Created    : November 16th 2001                                               #
#                                                                                #
#  Purpose    : LINUX makefile for kill.                                         #
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
PSILIB=
INCLUDEDIRS=-I../include
LIBDIRS=-L../lib
LIBS=-lucsf -ldbase -lgen -lnrec -lm -lexpat -lmysqlclient $(PSILIB)

all:
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -c $(INCLUDEDIRS) kill_main.cpp -o kill_main.o
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -o ../bin/kill.cgi kill_main.o $(LIBDIRS) $(LIBS) $(STATIC)
