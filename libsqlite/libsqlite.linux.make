##################################################################################
#                                                                                #
#  Library    : libsqlite                                                        #
#                                                                                #
#  Filename   : libsqlite.linux.make                                             #
#                                                                                #
#  Created    : September 2nd 2013                                               #
#                                                                                #
#  Purpose    : LINUX makefile for libsqlite.                                    #
#                                                                                #
#  Author(s)  : Peter Baker                                                      #
#                                                                                #
#  This file is the confidential and proprietary product of The Regents of       #
#  the University of California.  Any unauthorized use, reproduction or          #
#  transfer of this file is strictly prohibited.                                 #
#                                                                                #
#  Copyright (2014-2014) The Regents of the University of California.            #
#                                                                                #
#  All rights reserved.                                                          #
#                                                                                #
##################################################################################

C_COMPILER=gcc
LIBRARY_NAME = ../lib/libsqlite.a
ARCHIVER=ar
RANLIB=ranlib
OPTIONS= -O2 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE
ADD_OPTIONS=
INCLUDEDIRS=-I../include

ARCHIVE_COMMAND = $(ARCHIVER) r $(LIBRARY_NAME)

OBJ=sqlite3.o

.SUFFIXES:	.c

.c.o:	
	$(C_COMPILER) $(OPTIONS) $(ADD_OPTIONS) $(INCLUDEDIRS) -c $<

libsqlite.a:$(OBJ)
	$(ARCHIVE_COMMAND) $(OBJ)
	$(RANLIB) $(LIBRARY_NAME)
clean:	
	/bin/rm -f *.o 	
