/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_cgi_val.h                                                  *
*                                                                             *
*  Created    : July 11th 2007                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2007-2009) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_cgi_val_h
#define __lu_cgi_val_h

#include <ostream>
#include <string>
#include <lgen_define.h>

void printCGI ( std::ostream& os, const std::string& name, double value, int sigFig );

template <class T>
void printCGI ( std::ostream& os, const std::string& name, T value )
{
	os << name << "=" << value << "&";
}

void printCGIString ( std::ostream& os, const std::string& name, const std::string& value );
std::string getCommandLineNVPair ( const std::string& name, int value );
std::string getCommandLineNVPair ( const std::string& name, const std::string& value );

void printCGIContainer ( std::ostream& os, const std::string& name, const IntVector& value );
void printCGIContainer ( std::ostream& os, const std::string& name, const StringVector& value );

std::string escapeURL ( const std::string& url );

#endif /* ! __lu_cgi_val_h */
