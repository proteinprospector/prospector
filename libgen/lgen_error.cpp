/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lgen_error.cpp                                                *
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
*  Copyright (1996-2008) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef VIS_C
#include <process.h>
#include <windows.h>
#else
#include <errno.h>
#include <stdexcept>
#endif
#include <iostream>
#define LGEN_ERROR_MAIN
#include <lgen_error.h>
using std::string;
using std::cout;
using std::runtime_error;

#ifdef VIS_C
string genGetLastWindowsError ()
{
	DWORD dw = GetLastError ();
	LPVOID lpMsgBuf;
	FormatMessage ( FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
        NULL, dw, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), (LPTSTR) &lpMsgBuf, 0, NULL );
	string s = (char*)lpMsgBuf;
	LocalFree ( lpMsgBuf );
	return s;
}
#endif
ErrorHandler* ErrorHandler::_instance = 0;
ErrorHandler* ErrorHandler::genError ()
{
	if ( _instance == 0 ) {
		_instance = new ErrorHandler;
	}
	return _instance;
}
int ErrorHandler::getErrorNumber ()
{
	return errno;
}
void ErrorHandler::resetErrorNumber ()
{
	errno = 0;
}
void ErrorHandler::registration ( ErrorHandler* eh )
{
	if ( _instance != 0 ) {
		delete _instance;
	}
	_instance = eh;
}
ErrorHandler::ErrorHandler ()
{
	errno = 0;
}
void ErrorHandler::message ( const string& messageString )
{
	messageDisplay ( messageString );
}
void ErrorHandler::message ( const runtime_error& e )
{
	messageDisplay ( string ( e.what () ) + "\n" );
}
void ErrorHandler::warning ( const string& warningString, bool endProgram )
{
	string fullErrorString ( warningString );
	/*	These errors are frequently confusing
	if ( errno ) {
		fullErrorString += "\n";
		fullErrorString += "System error message: ";
		fullErrorString += strerror ( errno );
		fullErrorString += ".\n";
	}
	*/
	errorDisplay ( fullErrorString, endProgram );
}
void ErrorHandler::error ( const string& errorString )
{
	warning ( errorString, true );
	throw ( runtime_error ( errorString.c_str () ) );
}
void ErrorHandler::error ( const runtime_error& e )
{
	warning ( string ( e.what () ) + "\n", true );
	throw ( runtime_error ( e.what () ) );
}
void ErrorHandler::errorDisplay ( const string& errorString, bool endProgram )
{
	cout << "\n";
	cout << "ERROR MESSAGE\n";
	cout << "-------------\n";
	cout << "\n";
	cout << errorString;
}
void ErrorHandler::messageDisplay ( const string& messageString )
{
	cout << messageString;
	cout.flush ();
}
