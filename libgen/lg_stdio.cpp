/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lg_stdio.cpp                                                  *
*                                                                             *
*  Created    : July 8th 1996                                                 *
*                                                                             *
*  Purpose    : Machine independent basic I/O functions.                      *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1996-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifndef VIS_C
#include <unistd.h>
#endif
#include <fcntl.h>
#include <lg_stdio.h>
#include <lg_string.h>
#include <lgen_error.h>
using std::string;

int gen_open ( const string& path, int flags, const string& errorLocation )
{
	int fd = open ( path.c_str (), flags );

	if ( fd == -1 ) {
		ErrorHandler::genError ()->error ( "File open (open) failure.\n" + errorLocation + "\nFilename: " + path + ".\n" );
	}
	return fd;
}
int gen_close ( int fd, const string& errorLocation )
{
	int ret = close ( fd );

	if ( ret == -1 ) {
		ErrorHandler::genError ()->error ( "File close (close) failure.\n" + errorLocation + "\n" );
	}
	return ret;
}
FILE* gen_fopen ( const string& filename, const char* type, const string& errorLocation )
{
	FILE* fp = fopen ( filename.c_str (), type );

	if ( fp == NULL ) {
		ErrorHandler::genError ()->error ( "File open (fopen) failure.\n" + errorLocation + "\nFilename: " + filename + ".\n" );
	}
	return fp;
}
FILE* gen_fopen_text ( const string& filename, const char* type, const string& errorLocation )
{
	FILE* fp;
#ifdef VIS_C
	char t [10];

	strcpy ( t, type );
	strcat ( t, "t" );
	fp = gen_fopen ( filename, t, errorLocation );
#else
	fp = gen_fopen ( filename, type, errorLocation );
#endif
	return fp;
}
FILE* gen_fopen_binary ( const string& filename, const char* type, const string& errorLocation )
{
	FILE* fp;
#ifdef VIS_C
	char t [10];

	strcpy ( t, type );
	if ( type == "w" && isPrefix ( filename, "\\\\" ) ) {	// Check if this is a UNC drive
		strcat ( t, "+" );									// This gets round 64 MByte fwrite limit on UNC drive
	}
	strcat ( t, "b" );
	fp = gen_fopen ( filename, t, errorLocation );
#else
	fp = gen_fopen ( filename, type, errorLocation );
#endif
	return fp;
}
int gen_fclose ( FILE* fp, const string& errorLocation )
{
	int ret = fclose ( fp );

	if ( ret == EOF ) {
		ErrorHandler::genError ()->error ( "File close (fclose) failure.\n" + errorLocation + "\n" );
	}
	return ret;
}
int gen_fread ( char* ptr, size_t size, size_t nitems, FILE* fp, const string& errorLocation )
{
	return fread ( ptr, size, nitems, fp );
}
int gen_fwrite ( char* ptr, size_t size, size_t nitems, FILE* fp, const string& errorLocation )
{
	int ret = fwrite ( ptr, size, nitems, fp );

	if ( ret != nitems ) {
		ErrorHandler::genError ()->error ( "File write (fwrite) failure.\n" + errorLocation + "\n" );
	}
	return ret;
}
void gen_stdin_buffering_off ()
{
	setvbuf ( stdin, NULL, _IONBF, 0 );
}
void gen_stdout_buffering_off ()
{
	setvbuf ( stdout, NULL, _IONBF, 0 );
}
void gen_stderr_buffering_off ()
{
	setvbuf ( stderr, NULL, _IONBF, 0 );
}
void gen_stdin_binary_mode ()
{
#ifdef VIS_C
	_setmode ( _fileno( stdin ), _O_BINARY );
#endif
}
