##################################################################################
#                                                                                #
#  Program    : makeParamDB                                                      #
#                                                                                #
#  Filename   : make_param_db.linux.full.make                                    #
#                                                                                #
#  Created    : February 15th 2015                                               #
#                                                                                #
#  Purpose    : LINUX makefile for makeParamDB.                                  #
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

COMPILER=g++
OPTIONS=-O2 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE
ADD_OPTIONS=
STATIC=
INCLUDEDIRS=-I../include
LIBDIRS=-L../lib
LIBS=-lucsf -ldbase -lgen -lnrec -lm -lexpat -lpthread -lsqlite -ldl -lmysqlclient -lz

all:
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -c $(INCLUDEDIRS) make_param_db_main.cpp -o make_param_db_main.o
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -o ../bin/makeParamDB make_param_db_main.o $(LIBDIRS) $(LIBS) $(STATIC)
