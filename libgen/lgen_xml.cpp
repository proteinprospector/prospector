/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lgen_xml.cpp                                                  *
*                                                                             *
*  Created    : July 12th 2005                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2005-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <stdexcept>
#include <lg_io.h>
#include <lg_string.h>
#include <lgen_xml.h>
#include <expat.h>
using std::string;
using std::runtime_error;

PPExpat::PPExpat () :
	parser ( XML_ParserCreate ( NULL ) )
{
	XML_SetUserData ( parser, this );
	XML_SetElementHandler ( parser, startElement, endElement );
	XML_SetCharacterDataHandler ( parser, characterDataHandler );
}
PPExpat::~PPExpat ()
{
	if ( parser ) XML_ParserFree ( parser );
}
void PPExpat::startElement ( void* userData, const char* name, const char** attributes )
{
	static_cast <PPExpat*> (userData)->startElement ( name, attributes );
}
void PPExpat::endElement ( void* userData, const char* name )
{
	static_cast <PPExpat*> (userData)->endElement ( name );
}
void PPExpat::characterDataHandler ( void* userData, const char* str, int len )
{
	static_cast <PPExpat*> (userData)->characterDataHandler ( str, len );
}
void PPExpat::parseXMLFromFile ( const string& filename )
{
	static const int bsiz = 2048;
	char buf [bsiz];
	GenIFStream ist ( filename );
	for ( ; ; ) {
		ist.read ( buf, bsiz );
		int len = ist.gcount ();
		bool done = ( len < bsiz );
		if ( len ) parseXMLFromString ( buf, len, done );
		if ( done ) break;
	}
}
void PPExpat::parseXMLFromString ( const char* s, int length, bool last )
{
	enum XML_Status xStat = XML_Parse ( parser, s, length, last );
	if ( xStat != XML_STATUS_OK ) {
		if ( XML_GetErrorCode ( parser ) != XML_ERROR_INVALID_TOKEN ) {
			string message;
			if ( xStat == XML_STATUS_ERROR ) {
				message += "The XML parser returned the following error:\n";
				message += XML_ErrorString ( XML_GetErrorCode ( parser ) );
				message += " at line ";
				message += gen_itoa ( XML_GetCurrentLineNumber ( parser ) );
				message += ".\n";
			}
			else if ( xStat == XML_STATUS_SUSPENDED ) {
				message += "The XML parser has returned a suspended status.\n";
			}
			else {
				message += "The XML parser has returned an unknown status.\n";
			}
			throw runtime_error ( message );
		}
	}
}
void PPExpat::freeParser ()
{
	XML_ParserFree ( parser );
	parser = 0;
}
bool PPExpat::getAttributeValue ( const char** attributes, const char* name, bool& value )
{
	int i = 0;
	while ( attributes [i] ) {
		if ( !strcmp ( attributes [i++], name ) ) {
			string a = attributes [i];
			if ( a == "1" || !genStrcasecmp ( a, "true" ) ) value = true;
			else value = false;
			return true;
		}
		else i++;
	}
	return false;
}
bool PPExpat::getAttributeValue ( const char** attributes, const char* name, int& value )
{
	int i = 0;
	while ( attributes [i] ) {
		if ( !strcmp ( attributes [i++], name ) ) {
			value = atoi ( attributes [i] );
			return true;
		}
		else i++;
	}
	return false;
}
bool PPExpat::getAttributeValue ( const char** attributes, const char* name, double& value )
{
	int i = 0;
	while ( attributes [i] ) {
		if ( !strcmp ( attributes [i++], name ) ) {
			value = atof ( attributes [i] );
			return true;
		}
		else i++;
	}
	return false;
}
bool PPExpat::getAttributeValue ( const char** attributes, const char* name, string& value )
{
	int i = 0;
	while ( attributes [i] ) {
		if ( !strcmp ( attributes [i++], name ) ) {
			value = attributes [i];
			return true;
		}
		else i++;
	}
	return false;
}
bool PPExpat::getCVParamNameAndValue ( const char* elementName, const char** attributes, string& name, string& value )
{
	if ( !strcmp ( elementName, "cvParam" ) ) {
		getAttributeValue ( attributes, "name", name );
		bool f = getAttributeValue ( attributes, "value", value );
		return f;
	}
	return false;
}
bool PPExpat::getUserParamNameAndValue ( const char* elementName, const char** attributes, string& name, string& value )
{
	if ( !strcmp ( elementName, "userParam" ) ) {
		getAttributeValue ( attributes, "name", name );
		bool f = getAttributeValue ( attributes, "value", value );
		return f;
	}
	return false;
}
PPExpatStringList::PPExpatStringList ( const string& str ) :
	str ( str ),
	flag ( false )
{
}
PPExpatStringList::~PPExpatStringList () {}
void PPExpatStringList::startElement ( const char* name, const char** attributes )
{
	if ( !strcmp ( name, str.c_str () ) ) flag = true;
}
void PPExpatStringList::characterDataHandler ( const char* str, int len )
{
	if ( flag ) s.append ( str, len );
}
void PPExpatStringList::endElement ( const char* name )
{
	if ( flag ) {
		names.push_back ( s );
		s = "";
		flag = false;
	}
}
