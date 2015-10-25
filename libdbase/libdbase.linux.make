##################################################################################
#                                                                                #
#  Library    : libdbase                                                         #
#                                                                                #
#  Filename   : libdbase.linux.make                                              #
#                                                                                #
#  Created    : March 6th 2007                                                   #
#                                                                                #
#  Purpose    : LINUX makefile for libdbase.                                     #
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

LIBRARY_NAME = ../lib/libdbase.a
COMPILER=g++
ARCHIVER=ar
RANLIB=ranlib
OPTIONS= -O2 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE
ADD_OPTIONS=
INCLUDEDIRS=-I../include -I/usr/include/mysql

ARCHIVE_COMMAND = $(ARCHIVER) r $(LIBRARY_NAME)

OBJ=ld_init.o
.SUFFIXES:	.cpp

.cpp.o:	

	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) $(INCLUDEDIRS) -c $<

libdbase.a:$(OBJ)
	$(ARCHIVE_COMMAND) $(OBJ)
	$(RANLIB) $(LIBRARY_NAME)
clean:	
	/bin/rm -f *.o 	
