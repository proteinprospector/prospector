##################################################################################
#                                                                                #
#  Program    : search_compare                                                   #
#                                                                                #
#  Filename   : search_compare.linux.cl.make                                     #
#                                                                                #
#  Created    : October 22th 2003                                                #
#                                                                                #
#  Purpose    : LINUX makefile for search_compare.                               #
#                                                                                #
#  Author(s)  : Peter Baker                                                      #
#                                                                                #
#  This file is the confidential and proprietary product of The Regents of       #
#  the University of California.  Any unauthorized use, reproduction or          #
#  transfer of this file is strictly prohibited.                                 #
#                                                                                #
#  Copyright (2001-2014) The Regents of the University of California.            #
#                                                                                #
#  All rights reserved.                                                          #
#                                                                                #
##################################################################################

COMPILER=g++
OPTIONS=-O2 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE -DBATCHTAG -DRAW_DATA
ADD_OPTIONS=
STATIC=
LLIB=
INCLUDEDIRS=-I. -I../include
LIBDIRS=-L../lib
LIBS=-lraw -lucsf -lgen -lnrec -lm -lexpat -lz -lraw -lpthread -lsqlite $(LLIB)

all:
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -c $(INCLUDEDIRS) sc_anum_res.cpp -o sc_anum_res.o
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -c $(INCLUDEDIRS) sc_params.cpp -o sc_params.o
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -c $(INCLUDEDIRS) sc_search_res.cpp -o sc_search_res.o
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -c $(INCLUDEDIRS) sc_sres_link.cpp -o sc_sres_link.o
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -c $(INCLUDEDIRS) sc_sres_rep.cpp -o sc_sres_rep.o
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -c $(INCLUDEDIRS) sc_pep_xml.cpp -o sc_pep_xml.o
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -c $(INCLUDEDIRS) sc_mzidentml.cpp -o sc_mzidentml.o
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -c $(INCLUDEDIRS) sc_quan.cpp -o sc_quan.o
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -c $(INCLUDEDIRS) sc_xlink.cpp -o sc_xlink.o
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -c $(INCLUDEDIRS) search_compare_main.cpp -o search_compare_main.o
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) -o ../bin/searchCompare.cgi search_compare_main.o sc_anum_res.o sc_params.o sc_search_res.o sc_sres_link.o sc_sres_rep.o sc_pep_xml.o sc_mzidentml.o sc_quan.o sc_xlink.o $(LIBDIRS) $(LIBS) $(STATIC)
