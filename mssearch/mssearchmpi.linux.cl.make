##################################################################################
#                                                                                #
#  Program    : mssearch                                                         #
#                                                                                #
#  Filename   : mssearchmpi.linux.cl.make                                        #
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
OPTIONS=-O2 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE -DPP_MULTI -DBATCHTAG
ADD_OPTIONS=
STATIC=
MPILIB=-lmpich
LLIB=
PSILIB=
INCLUDEDIRS=-I../include
LIBDIRS=-L../lib
LIBS=-lraw -lucsf -lmulti -lgen -lnrec -lzip -lm -lz $(MPILIB) -lexpat -lraw -lpthread -lsqlite $(LLIB) $(PSILIB)

all:
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -c $(INCLUDEDIRS) mssearch_main.cpp -o mssearch_main.o
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -o ../bin/mssearchmpi.cgi mssearch_main.o $(LIBDIRS) $(LIBS) $(STATIC)
