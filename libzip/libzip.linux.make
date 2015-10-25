##################################################################################
#                                                                                #
#  Library    : libzip                                                           #
#                                                                                #
#  Filename   : libzip.linux.make                                                #
#                                                                                #
#  Created    : February 19th 2015                                               #
#                                                                                #
#  Purpose    : LINUX makefile for libsqlite.                                    #
#                                                                                #
#  Author(s)  : Peter Baker                                                      #
#                                                                                #
#  This file is the confidential and proprietary product of The Regents of       #
#  the University of California.  Any unauthorized use, reproduction or          #
#  transfer of this file is strictly prohibited.                                 #
#                                                                                #
#  Copyright (2015-2015) The Regents of the University of California.            #
#                                                                                #
#  All rights reserved.                                                          #
#                                                                                #
##################################################################################
C_COMPILER=gcc
LIBRARY_NAME = ../lib/libzip.a
ARCHIVER=ar
RANLIB=ranlib
OPTIONS= -O2 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE
ADD_OPTIONS=
INCLUDEDIRS=-I../include

ARCHIVE_COMMAND = $(ARCHIVER) r $(LIBRARY_NAME)

OBJ=fseek.o \
	unzip.o \
	unzipx.o \
	lz_app.o

.SUFFIXES:	.c

.c.o:	
	$(C_COMPILER) $(OPTIONS) $(ADD_OPTIONS) $(INCLUDEDIRS) -c $<

.cpp.o:	
	$(C_COMPILER) $(OPTIONS) $(ADD_OPTIONS) $(INCLUDEDIRS) -c $<

libzip.a:$(OBJ)
	$(ARCHIVE_COMMAND) $(OBJ)
	$(RANLIB) $(LIBRARY_NAME)
clean:	
	/bin/rm -f *.o 	
