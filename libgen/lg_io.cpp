/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lg_io.cpp                                                     *
*                                                                             *
*  Created    : July 25rd 2000                                                *
*                                                                             *
*  Purpose    : Error checking interface to io library.                       *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2000-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <cmath>
#include <iomanip>
#include <lg_io.h>
#include <lg_string.h>
#include <lgen_error.h>
using std::istream;
using std::ostream;
using std::setw;
using std::scientific;
using std::fixed;
using std::setprecision;
using std::ifstream;
using std::ofstream;
using std::string;
using std::istringstream;
using std::ostringstream;
using std::streambuf;

#ifdef VIS_C
#pragma warning( disable : 4503 )	// Disable warnings about decorated name length being too long
#endif

GenIFStream::GenIFStream () :
	ifstream ()
{
}
GenIFStream::GenIFStream ( const string& fullPath, ios_base::openmode mode ) :
	ifstream ( fullPath.c_str (), mode )
{
	if ( !(*this) ) {
		ErrorHandler::genError ()->error ( "File open failure.\nFilename: " + genEscapePath ( fullPath ) + ".\n" );
	}
}
void GenIFStream::open ( const string& fullPath, ios_base::openmode mode )
{
	ifstream::open ( fullPath.c_str (), mode );
}

GenCommentedIFStream::GenCommentedIFStream ( const string& fullPath ) :
	GenIFStream ( fullPath )
{
}
bool GenCommentedIFStream::getUncommentedLine ( string& line )
{
	char buffer [256];
	while ( getline ( buffer, 256 ) ) {
		line = buffer;
		if ( line.length () != 0 ) {
			if ( line [0] != '#' ) {
				return true;
			}
		}
	}
	return false;
}

GenNameValueStream::GenNameValueStream ( const string& fullPath ) :
	GenCommentedIFStream ( fullPath )
{
	string line;
	while ( getUncommentedLine ( line ) ) {
		istringstream ist ( line );
		string name;
		ist >> name;
		string value = gen_strtrim ( line.substr ( name.length () ) );
		params [name].push_back ( value );
	}
}
bool GenNameValueStream::getValue ( const string& name, string& value )
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () ) {
		value = (*cur).second [0];
		return true;
	}
	else return false;
}
bool GenNameValueStream::getValue ( const string& name, int& value )
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () ) {
		value = atoi ( (*cur).second [0].c_str () );
		return true;
	}
	else return false;
}
bool GenNameValueStream::getValue ( const string& name, unsigned int& value )
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () ) {
		value = strtoul ( (*cur).second [0].c_str (), 0, 10 );
		return true;
	}
	else return false;
}
bool GenNameValueStream::getValue ( const string& name, double& value )
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () ) {
		value = atof ( (*cur).second [0].c_str () );
		return true;
	}
	else return false;
}
bool GenNameValueStream::getValue ( const string& name, StringVector& value ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () ) {
		value = (*cur).second;
		return true;
	}
	else return false;
}
bool GenNameValueStream::getBoolValue ( const string& name, bool defaultVal ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () ) {
		const char* c = (*cur).second [0].c_str ();
		return !genStrcasecmp ( c, "true" ) || !strcmp ( c, "1" );
	}
	else
		return defaultVal;
}
int GenNameValueStream::getIntValue ( const string& name, int defaultVal ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () )
		return atoi ( (*cur).second [0].c_str () );
	else
		return defaultVal;
}
unsigned int GenNameValueStream::getUIntValue ( const string& name, unsigned int defaultVal ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () )
		return strtoul ( (*cur).second [0].c_str (), 0, 10 );
	else
		return defaultVal;
}
double GenNameValueStream::getDoubleValue ( const string& name, double defaultVal ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () )
		return atof ( (*cur).second [0].c_str () );
	else
		return defaultVal;
}
string GenNameValueStream::getStringValue ( const string& name, const string& defaultVal ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () )
		return (*cur).second [0];
	else
		return defaultVal;
}
StringVector GenNameValueStream::getStringVectorValue ( const string& name, const StringVector& defaultVal ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () )
		return (*cur).second;
	else
		return defaultVal;
}
StringVector GenNameValueStream::getNameList () const
{
	StringVector sv;
	for ( MapStringToStringVectorConstIterator i = params.begin () ; i != params.end () ; i++ ) {
		sv.push_back ( (*i).first );
	}
	return sv;
}
GenOFStream::GenOFStream ( const string& fullPath, ios_base::openmode mode ) :
	ofstream ( fullPath.c_str (), mode )
{
	if ( !(*this) ) {
		ErrorHandler::genError ()->error ( "File open failure.\nFilename: " + fullPath + ".\n" );
	}
}
void GenOFStream::open ( const string& fullPath, ios_base::openmode mode )
{
	ofstream::open ( fullPath.c_str (), mode );
}
void genPrint ( ostream& os, double number, int precision, int width )
{
	if ( width != -1 ) {
		os << setw ( width );
	}
	os << fixed;
	if ( precision != -1 ) {
		os << setprecision ( precision );
	}
	os << number;
}
void genPrintSigFig ( ostream& os, double number, int precision )
{
	if ( number == 0.0 ) {
		os << setprecision ( 0 ) << fixed << 0;
	}
	else {
		double absNumber = fabs ( number );
		if ( absNumber >= 1000000 || absNumber < 0.001 ) {
			ostringstream ostr;
			ostr << scientific << setprecision ( precision - 1 );
			ostr << number;
			string s = ostr.str ();		// Strip the leading zeros out of the exponent
			int ind1 = s.find_first_of ( "Ee" ) + 2;
			int ind2 = s.find_first_not_of ( "0", ind1 );
			os << s.substr ( 0, ind1 ) + s.substr ( ind2, s.length () - ind2 );
		}
		else {
			int numDigits;
			if ( absNumber >= 100000 ) numDigits = 6;
			else if ( absNumber >= 10000 ) numDigits = 5;
			else if ( absNumber >= 1000 ) numDigits = 4;
			else if ( absNumber >= 100 ) numDigits = 3;
			else if ( absNumber >= 10 ) numDigits = 2;
			else if ( absNumber >= 1 ) numDigits = 1;
			else if ( absNumber >= 0.1 ) numDigits = 0;
			else if ( absNumber >= 0.01 ) numDigits = -1;
			else if ( absNumber >= 0.001 ) numDigits = -2;
			else numDigits = -3;
			int fixedPrecision = precision - numDigits;
			if ( fixedPrecision > 0 ) os << setprecision ( fixedPrecision );
			else os << setprecision ( 0 );
			os << fixed;
			os << number;
		}
	}
}
// Universal getline, should handle CR, CRLF or LF

// sbumpc - Returns the character currently pointed by the get pointer and advances it one position.
// sgetc  - Returns the character currently pointed by the get pointer.
istream& genUniversalGetLine ( istream& is, string& s )
{
	s.erase ();
	is.unsetf(std::ios_base::skipws);
	istream::sentry se ( is );
	streambuf* sb = is.rdbuf ();
	for ( ; ; ) {
		int c = sb->sbumpc ();
		switch ( c ) {
			case '\r':
				c = sb->sgetc ();
				if ( c == '\n' ) {
					sb->sbumpc ();
				}
				return is;
			case '\n':
				return is;
			case EOF:
				is.clear ( is.rdstate () | std::ios_base::eofbit );
				{									// This bit seems necessary for VIS_C so fn returns 0 at EOF
					istream::sentry se2 ( is );
					sb = is.rdbuf ();
					sb->sbumpc ();
				}
				return is;
			default:
				s += (char)c;
		}
	}
}
