##################################################################################
#                                                                                #
#  Library    : libexpat                                                         #
#                                                                                #
#  Filename   : libexpat.linux.make                                              #
#                                                                                #
#  Created    : September 2nd 2013                                               #
#                                                                                #
#  Purpose    : LINUX makefile for liboracle.                                    #
#                                                                                #
#  Author(s)  : Peter Baker                                                      #
#                                                                                #
#  This file is the confidential and proprietary product of The Regents of       #
#  the University of California.  Any unauthorized use, reproduction or          #
#  transfer of this file is strictly prohibited.                                 #
#                                                                                #
#  Copyright (2013-2013) The Regents of the University of California.            #
#                                                                                #
#  All rights reserved.                                                          #
#                                                                                #
##################################################################################

LIBRARY_NAME = ../lib/libexpat.a
COMPILER=g++
ARCHIVER=ar
RANLIB=ranlib
OPTIONS= -O2 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE -DHAVE_EXPAT_CONFIG_H
ADD_OPTIONS=
INCLUDEDIRS=-I../include

ARCHIVE_COMMAND = $(ARCHIVER) r $(LIBRARY_NAME)

OBJ=xmlparse.o \
	xmlrole.o \
	xmltok.o \
	xmltok_impl.o \
	xmltok_ns.o

.SUFFIXES:	.c

.c.o:	
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) $(INCLUDEDIRS) -c $<

liboracle.a:$(OBJ)
	$(ARCHIVE_COMMAND) $(OBJ)
	$(RANLIB) $(LIBRARY_NAME)
clean:	
	/bin/rm -f *.o 	
