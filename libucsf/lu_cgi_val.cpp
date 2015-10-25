/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_cgi_val.cpp                                                *
*                                                                             *
*  Created    : July 11th 2007                                                *
*                                                                             *
*  Purpose    :                                                               *
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
#include <lg_io.h>
#include <lg_string.h>
#include <lu_cgi_val.h>

using std::ostream;
using std::string;

void printCGI ( ostream& os, const string& name, double value, int sigFig )
{
	os << name;
	os << "=";
	genPrintSigFig ( os, value, sigFig );
	os << "&";
}
void printCGIString ( ostream& os, const string& name, const string& value )
{
	os << name << "=" << escapeURL ( value ) << "&";
}
string getCommandLineNVPair ( const string& name, int value )
{
	string s;
	s += name;
	s += "=";
	s += gen_itoa ( value );
	return s;
}
string getCommandLineNVPair ( const string& name, const string& value )
{
	string s;
	s += name;
	s += "=";
	s += escapeURL ( value );
	return s;
}
string escapeURL ( const string& url )
{
	string url2;
	for ( StringSizeType i = 0 ; i < url.length () ; i++ ) {
		char c = url [i];
		switch ( c ) {
			char escapedChar [3];
			case '~':case '`':case '!':case '#':case '$':
			case '%':case '^':case '&':case '(':case ')':
			case '+':case '=':case '{':case '[':case '}':
			case ']':case ':':case ';':case '"':case '\'':
			case '<':case ',':case '>':case '?':case '/':
			case '|':case ' ':case '\\':
				url2 += '%';
				sprintf ( escapedChar, "%2X", c );
				url2 += escapedChar;
				break;
			case '\r':
				break;
			case '\n':
				url2 += "%0D%0A";
				break;
			default:
				url2 += c;
				break;
		}
	}
	return ( url2 );
}
void printCGIContainer ( ostream& os, const string& name, const IntVector& value )
{
	if ( value.size () ) {
		os << name << "=";
		for ( IntVectorSizeType i = 0 ; i < value.size () ; i++ ) {
			os << value [i] << "%0D%0A";
		}
		os << "&";
	}
}
void printCGIContainer ( ostream& os, const string& name, const StringVector& value )
{
	if ( value.size () ) {
		os << name << "=";
		for ( StringVectorSizeType i = 0 ; i < value.size () ; i++ ) {
			os << value [i] << "%0D%0A";
		}
		os << "&";
	}
}
