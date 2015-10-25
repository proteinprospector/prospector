/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lgen_xml.h                                                    *
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
*  Copyright (2005-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lgen_xml_h
#define __lgen_xml_h

#include <string>
#include <lgen_define.h>

struct XML_ParserStruct;						// forward definitions so expat.h not required
typedef struct XML_ParserStruct *XML_Parser;

class PPExpat {
	XML_Parser parser;

	static void startElement ( void* userData, const char* name, const char** atts );
	static void endElement ( void* userData, const char* name );
	static void characterDataHandler ( void* userData, const char* str, int len );

	virtual void startElement ( const char* name, const char** attributes ) {}
	virtual void endElement ( const char* name ) {}
	virtual void characterDataHandler ( const char* str, int len ) {}
protected:
	std::string s;	// string for holding character data
	static bool getAttributeValue ( const char** attributes, const char* name, bool& value );
	static bool getAttributeValue ( const char** attributes, const char* name, int& value );
	static bool getAttributeValue ( const char** attributes, const char* name, double& value );
	static bool getAttributeValue ( const char** attributes, const char* name, std::string& value );
	template <class T>
	static bool getCVAttributeValue ( const char* elementName, const char** attributes, const char* name, T& value )
	{
		if ( !strcmp ( elementName, "cvParam" ) ) {
			std::string n;
			getAttributeValue ( attributes, "name", n );
			if ( name == n ) {
				bool f = getAttributeValue ( attributes, "value", value );
				return f;
			}
		}
		return false;
	}
	static bool getCVAttributeName ( const char* elementName, const char** attributes, const char* name )
	{
		if ( !strcmp ( elementName, "cvParam" ) ) {
			std::string n;
			getAttributeValue ( attributes, "name", n );
			if ( name == n ) {
				return true;
			}
		}
		return false;
	}
	static bool getCVParamNameAndValue ( const char* elementName, const char** attributes, std::string& name, std::string& value );
	static bool getUserParamNameAndValue ( const char* elementName, const char** attributes, std::string& name, std::string& value );
	void freeParser ();
public:
	PPExpat ();
	virtual ~PPExpat ();
	void parseXMLFromFile ( const std::string& filename );
	void parseXMLFromString ( const std::string& s )
	{
		parseXMLFromString ( s.c_str (), s.length () );
	}
	void parseXMLFromString ( const char* s, int length, bool last = true );
};

class PPExpatStringList : public PPExpat {
	std::string str;
	bool flag;
	StringVector names;

	void startElement ( const char* name, const char** attributes );
	void characterDataHandler ( const char* str, int len );
	void endElement ( const char* name );
public:
	PPExpatStringList ( const std::string& str );
	~PPExpatStringList ();
	StringVector getNames () const { return names; }
	IntVector getNamesAsInt () const
	{
		IntVector iNames;
		for ( StringVectorSizeType i = 0 ; i < names.size () ; i++ ) {
			iNames.push_back ( atoi ( names [i].c_str () ) );
		}
		return iNames;
	}
};

#endif /* ! __lgen_xml_h */
