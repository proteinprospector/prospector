/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lg_new.cpp                                                    *
*                                                                             *
*  Created    : July 10th 2000                                                *
*                                                                             *
*  Purpose    : Machine independent memory allocation functions based on new. *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2000-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <cstring>
#ifdef VIS_C_6
#include <new.h>
#else
#include <cstdio>
#include <new>
#endif
#include <lgen_error.h>
#ifndef VIS_C_6
using std::set_new_handler;
#endif

static void gen_out_of_memory ();
static int gen_out_of_memory ( size_t size );

char* gen_new_string ( const char* string )
{
	if ( string == NULL ) return ( NULL );

	char* new_string = new char [strlen ( string ) + 1];
	strcpy ( new_string, string );

	return new_string;
}
char* genNewString ( const char* start, const char* end )
{
	int len = end - start - 1;
	char* newString = new char [len+1];
	strncpy ( newString, start, len );

	return newString;
}
void gen_set_new_handler ()
{
#ifdef VIS_C_6
	_set_new_handler ( gen_out_of_memory );
#else
	set_new_handler ( gen_out_of_memory );
#endif
}
static void gen_out_of_memory ()
{
	ErrorHandler::genError ()->error ( "Memory allocation failure.\n" );
}
static int gen_out_of_memory ( size_t size )
{
	sprintf ( gen_error_message, "Memory allocation failure (size = %d).\n", size );
	ErrorHandler::genError ()->error ( gen_error_message );
	return 0;
}
