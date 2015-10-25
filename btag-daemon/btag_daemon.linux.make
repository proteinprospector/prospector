##################################################################################
#                                                                                #
#  Program    : btag_daemon                                                      #
#                                                                                #
#  Filename   : btag_daemon.linux.make                                           #
#                                                                                #
#  Created    : October 5th 2007                                                 #
#                                                                                #
#  Purpose    : LINUX makefile for btag_daemon.                                  #
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
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -c $(INCLUDEDIRS) btag_daemon_main.cpp -o btag_daemon_main.o
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -o ../bin/btag-daemon btag_daemon_main.o $(LIBDIRS) $(LIBS) $(STATIC)
