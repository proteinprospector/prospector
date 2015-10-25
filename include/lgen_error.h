/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lgen_error.h                                                  *
*                                                                             *
*  Created    : July 7th 1996                                                 *
*                                                                             *
*  Purpose    : Functions for printing error messages.                        *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1996-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lgen_error_h
#define __lgen_error_h

#include <stdexcept>
#include <string>

#ifdef LGEN_ERROR_MAIN
#define LGEN_ERROR_EXTERN
#else
#define LGEN_ERROR_EXTERN extern
#endif

#define MAX_ERROR_MESSAGE_LENGTH 2000

LGEN_ERROR_EXTERN char gen_error_message [MAX_ERROR_MESSAGE_LENGTH];
LGEN_ERROR_EXTERN char gen_output_message [MAX_ERROR_MESSAGE_LENGTH];

#ifdef VIS_C
std::string genGetLastWindowsError ();
#endif

class ErrorHandler {
	virtual void messageDisplay ( const std::string& messageString );
	virtual void errorDisplay ( const std::string& errorString, bool endProgram );
	static ErrorHandler* _instance;
protected:
	ErrorHandler ();
	static void registration ( ErrorHandler* eh );
public:
	static ErrorHandler* genError ();
	static int getErrorNumber ();
	static void resetErrorNumber ();
	void message ( const std::string& messageString );
	void message ( const std::runtime_error& e );
	void warning ( const std::string& errorString, bool endProgram = false );
	void error ( const std::string& errorString );
	void error ( const std::runtime_error& e );
};

#endif /* ! __lgen_error_h */
