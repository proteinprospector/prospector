##################################################################################
#                                                                                #
#  Library    : libgen                                                           #
#                                                                                #
#  Filename   : libgen.linux.make                                                #
#                                                                                #
#  Created    : September 11th 2001                                              #
#                                                                                #
#  Purpose    : LINUX makefile for libgen.                                       #
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

LIBRARY_NAME = ../lib/libgen.a
COMPILER=g++
ARCHIVER=ar
RANLIB=ranlib
OPTIONS= -O2 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE
ADD_OPTIONS=
INCLUDEDIRS=-I../include

ARCHIVE_COMMAND = $(ARCHIVER) r $(LIBRARY_NAME)

OBJ=lg_base64.o \
	lg_email.o \
	lg_io.o \
	lg_memory.o \
	lg_new.o \
	lg_stdio.o \
	lg_stdlib.o \
	lg_string.o \
	lg_time.o \
	lgen_error.o \
	lgen_file.o \
	lgen_getch.o \
	lgen_math.o \
	lgen_net.o \
	lgen_process.o \
	lgen_reg_exp2.o \
	lgen_service.o \
	lgen_uncompress.o \
	lgen_xml.o

.SUFFIXES:	.cpp

.cpp.o:	
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) $(INCLUDEDIRS) -c $<

libgen.a:$(OBJ)
	$(ARCHIVE_COMMAND) $(OBJ)
	$(RANLIB) $(LIBRARY_NAME)
clean:	
	/bin/rm -f *.o 	
