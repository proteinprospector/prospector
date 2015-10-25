/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lg_stdlib.cpp                                                 *
*                                                                             *
*  Created    : July 10th 1996                                                *
*                                                                             *
*  Purpose    : Common assorted library functions.                            *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1996-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lgen_error.h>
#include <lgen_file.h>
using std::string;

/*
Visual C++ 6.0 manual entry for system
--------------------------------------

Return Value

If command is NULL and the command interpreter is found, the function returns a nonzero value.
If the command interpreter is not found, it returns 0 and sets errno to ENOENT.
If command is not NULL, system returns the value that is returned by the command interpreter.
It returns the value 0 only if the command interpreter returns the value 0. A return value of
1 indicates an error, and errno is set to one of the following values:

E2BIG - Argument list (which is system-dependent) is too big.

ENOENT - Command interpreter cannot be found.

ENOEXEC - Command-interpreter file has invalid format and is not executable.

ENOMEM - Not enough memory is available to execute command; or available memory has been corrupted; or invalid block exists, indicating that process making call was not allocated properly.

LINUX manual entry for system
-----------------------------

The value returned is -1 on error (e.g. fork failed), and the return status of the command otherwise.
This latter return status is in the format specified in wait(2).
Thus, the exit code of the command will be WEXITSTATUS(status).
In case /bin/sh could not be executed, the exit status will be that of a command
that does exit(127). 

wait (2) info  - http://linux.about.com/od/commands/l/blcmdl2_wait.htm
*/
int genSystem ( const string& command, const string& workingDirectory, bool cont )
{
	string previousWorkingDirectory;
	if ( !workingDirectory.empty () ) {
		previousWorkingDirectory = genCurrentWorkingDirectory ();
		genChangeWorkingDirectory ( workingDirectory );
	}
// NT and 2000 maximum command length is 2047 characters
// XP maximum command length is 8191 characters
// If the command and argument need spaces then the whole lot needs to be in quotes
// Eg ""some program.exe -f "some argument""
// ""C:\Program Files (x86)\R\R-2.2.1\bin\R" --vanilla --slave --quiet --args "C:\Program Files (x86)\UCSF\Prospector\web\temp\Nov_14_2013\datas4bc.1.txt" "C:\Program Files (x86)\UCSF\Prospector\web\temp\Nov_14_2013\images4bc.11.png" < msmsErrorScatter.R"
	int ret = system ( command.c_str () );
	if ( !workingDirectory.empty () ) {
		genChangeWorkingDirectory ( previousWorkingDirectory );
	}
	if ( ret != 0 && !cont ) {
		ErrorHandler::genError ()->error ( workingDirectory + "\n" + command + "\n" );
	}
	return ret;
}
char* gen_getenv ( char* name, const string& error_location_string )
{
	char* pointer = getenv ( name );

	if ( pointer == NULL ) {
		ErrorHandler::genError ()->error ( "Get environment variable (getenv) failure.\n" + error_location_string + "\nEnvironment variable: " + name + "\n" );
	}
	return ( pointer );
}
void gen_qsort ( char* base, int num, int width, int ( *compare ) ( const void*, const void* ) )
{
#ifdef VIS_C
	if ( num > 0 ) qsort ( (void*)base, (size_t)num, (size_t)width, (int(*)(const void*, const void*))compare );
#else
	if ( num > 0 ) qsort ( base, num, width, compare );
#endif
}
