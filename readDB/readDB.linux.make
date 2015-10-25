##################################################################################
#                                                                                #
#  Program    : readDB                                                           #
#                                                                                #
#  Filename   : readDB.linux.make                                                #
#                                                                                #
#  Created    : May 16th 2007                                                    #
#                                                                                #
#  Purpose    : LINUX makefile for readDB.                                       #
#                                                                                #
#  Author(s)  : Peter Baker                                                      #
#                                                                                #
#  This file is the confidential and proprietary product of The Regents of       #
#  the University of California.  Any unauthorized use, reproduction or          #
#  transfer of this file is strictly prohibited.                                 #
#                                                                                #
#  Copyright (2008-2013) The Regents of the University of California.            #
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
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -c $(INCLUDEDIRS) readDB_main.cpp -o readDB_main.o
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -o ../bin/readDB.cgi readDB_main.o $(LIBDIRS) $(LIBS) $(STATIC)
