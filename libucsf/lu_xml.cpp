/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_xml.cpp                                                    *
*                                                                             *
*  Created    : February 3rd 2003                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2003-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lg_time.h>
#include <lgen_error.h>
#include <lu_getfil.h>
#include <lu_xml.h>
using std::string;
using std::istream;
using std::ostream;
using std::endl;
using std::getline;

void printXMLHeader ( ostream& os )
{
	os << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
}
void printXMLStylesheet ( ostream& os, const string& scriptType, const string& scriptName )
{
	os << "<?xml-stylesheet type=\"";
	os << scriptType;
	os << "\" href=\"";
	os << scriptName;
	os << "\"?>";
	os << endl;
}
void printXMLVersion ( ostream& os )
{
	os << "<?";
	os << genCurrentTimeAndDateString ();
	os << ", ProteinProspector Version ";
	os << Version::instance ().getVersion ();
	os << "?>";
	os << endl;
}

string getXMLFileDate ( const string& fname )
{
	GenIFStream ist ( fname );
	string d;
	string line;
	if ( getline ( ist, line ) ) {
		if ( line.find ( "<?xml" ) != string::npos ) {	// first line states this is an xml file
			if ( getline ( ist, line ) ) {
				if ( line.find ( "<?" ) != string::npos ) {	// this line has the date
					int end = line.find ( "," );
					d = line.substr ( 2, end - 2 );
				}
			}
		}
	}
	return d;
}
void openXMLTag ( ostream& os, int ntab, const string& name, const VectorPairStringString& attr, bool close )
{
	size_t siz = attr.size ();
	for ( size_t i = 0 ; i < ntab ; i++ ) os << '\t';
	os << "<" << name;
	for ( size_t j = 0 ; j < siz ; j++ ) {
		os << " " << attr [j].first << "=" << "\"" << attr [j].second << "\"";
	}
	if ( close ) os << " /";
	os << ">" << endl;
}
void openXMLValueTag ( ostream& os, int ntab, const string& name, const VectorPairStringString& attr, bool close )
{
	size_t siz = attr.size ();
	for ( size_t i = 0 ; i < ntab ; i++ ) os << '\t';
	os << "<" << name;
	for ( size_t j = 0 ; j < siz ; j++ ) {
		os << " " << attr [j].first << "=" << "\"" << attr [j].second << "\"";
	}
	if ( close ) os << " />" << endl;
	else os << ">";
}
void closeXMLTag ( ostream& os, int ntab, const string& name )
{
	for ( size_t i = 0 ; i < ntab ; i++ ) os << '\t';
	os << "</";
	os << name;
	os << ">";
	os << endl;
}
void closeXMLTag ( ostream& os, const string& name )
{
	os << "</" << name << ">" << endl;
}

namespace XMLParser {
string getStringValue ( const string& str, const string& tag )
{
	int nameStartIndex = str.find ( '<' + tag + '>' );
	if ( nameStartIndex == string::npos ) {
		ErrorHandler::genError ()->error ( "XML tag \"" + tag + "\" not found.\n" );
	}
	int valueStartIndex = nameStartIndex + tag.length () + 2;
	int valueLength = str.find ( "</" + tag + '>', valueStartIndex ) - valueStartIndex;
	return str.substr ( valueStartIndex, valueLength );
}
string getStringValue ( const string& str, const string& tag, const string& defaultValue )
{
	int nameStartIndex = str.find ( '<' + tag + '>' );
	if ( nameStartIndex == string::npos ) return defaultValue;
	int valueStartIndex = nameStartIndex + tag.length () + 2;
	int valueLength = str.find ( "</" + tag + '>', valueStartIndex ) - valueStartIndex;
	return str.substr ( valueStartIndex, valueLength );
}
bool getBoolValue ( const string& str, const string& tag )
{
	return getStringValue ( str, tag, "0" ) == "1";
}
char getCharValue ( const string& str, const string& tag )
{
	return (getStringValue ( str, tag ).c_str ())[0];
}
int getIntValue ( const string& str, const string& tag )
{
	return atoi ( getStringValue ( str, tag ).c_str () );
}
double getDoubleValue ( const string& str, const string& tag )
{
	return atof ( getStringValue ( str, tag ).c_str () );
}
double getDoubleValue ( const string& str, const string& tag, const string& defaultValue )
{
	return atof ( getStringValue ( str, tag, defaultValue ).c_str () );
}
double getDoubleValue ( const string& str, const string& tag, double defaultValue )
{
	if ( getStringValue ( str, tag, "DUMMY" ) == "DUMMY" ) return defaultValue;
	else return atof ( getStringValue ( str, tag ).c_str () );
}
void getMZI ( const string& str, double& mOverZ, int& charge, double& intensity )
{
	string s ( XMLParser::getStringValue ( str, "mzi" ) );
	int end1 = s.find ( ',' );
	int start2 = end1 + 1;
	int end2 = s.find ( ',', start2 );
	int start3 = end2 + 1;
	int len = s.length ();
	mOverZ = atof ( s.substr ( 0, end1 ).c_str () );
	charge = atoi ( s.substr ( start2, end2 - start2 ).c_str () );
	intensity = atof ( s.substr ( start3, len - start3 ).c_str () );
}
StringVector getStringVectorValue ( const string& str, const string& tag )
{
	StringVector ret;
	string nameStart = '<' + tag + '>';
	int len = nameStart.length ();
	string nameEnd = "</" + tag + '>';
	int start = 0;
	for ( ; ; ) {
		int nameStartIndex = str.find ( nameStart, start );
		if ( nameStartIndex == string::npos ) break;
		int valueStartIndex = nameStartIndex + len;
		int valueEnd = str.find ( nameEnd, valueStartIndex );
		ret.push_back ( str.substr ( valueStartIndex, valueEnd - valueStartIndex ) );
		start = valueEnd + len + 1;
	}
	return ret;
}

}
XMLStringList::XMLStringList ( const string& str, const string& tag ) :
	str ( str ),
	nameStart ( '<' + tag + '>' ),
	len ( nameStart.length () ),
	nameEnd ( "</" + tag + '>' ),
	start ( 0 )
{
}
bool XMLStringList::getNext ( string& s )
{
	int nameStartIndex = str.find ( nameStart, start );
	if ( nameStartIndex == string::npos ) return false;
	int valueStartIndex = nameStartIndex + len;
	while ( str [valueStartIndex] < 31 ) valueStartIndex++;	//Skip end of line characters
	int valueEnd = str.find ( nameEnd, valueStartIndex );
	s = str.substr ( valueStartIndex, valueEnd - valueStartIndex );
	start = valueEnd + len + 1;
	return true;
}

XMLIStreamList::XMLIStreamList ( istream& ist, const string& tag ) :
	ist ( ist ),
	nameStart ( '<' + tag + '>' ),
	nameEnd ( "</" + tag + '>' )
{
}
bool XMLIStreamList::getNext ( string& s )
{
	bool flag = false;
	string line;
	while ( getline ( ist, line ) ) {
		if ( line.find ( nameStart ) != string::npos ) {	// start tag found
			while ( getline ( ist, line ) ) {
				if ( line.find ( nameEnd ) != string::npos ) return true;	// end tag found
				s += line;
			}
		}
	}
	return false;
}
bool XMLIStreamList::getNextBlock ( string& s )
{
	bool flag = false;
	string line;
	while ( getline ( ist, line ) ) {
		if ( line.find ( nameStart ) != string::npos ) {	// start tag found
			while ( getline ( ist, line ) ) {
				if ( line.find ( nameEnd ) != string::npos ) return true;	// end tag found
				s += line;
				s += '\n';
			}
		}
	}
	return false;
}
void XMLOutputItem::printOpenTag ( ostream& os, int ntab ) const
{
	openXMLTag ( os, ntab, tagName, attr, false );
}
void XMLOutputItem::printCloseTag ( ostream& os, int ntab ) const
{
	closeXMLTag ( os, ntab, tagName );
}
void XMLOutputAttrItem::print ( ostream& os, int ntab ) const
{
	openXMLTag ( os, ntab, tagName, attr, true );
}
XMLOutputContainerItem::~XMLOutputContainerItem ()
{
//	size_t siz = subItems.size ();
//	for ( size_t i = 0 ; i < siz ; i++ ) {
//		delete subItems [i];
//	}
}
void XMLOutputContainerItem::print ( ostream& os, int ntab ) const
{
	size_t siz = subItems.size ();
	if ( siz == 0 && attr.size () == 0 ) return;	// Print nothing if no information
	openXMLTag ( os, ntab, tagName, attr, !siz );
	if ( siz ) {
		for ( int i = 0 ; i < siz ; i++ ) {
			subItems [i]->print ( os, ntab+1 );
		}
		closeXMLTag ( os, ntab, tagName );
	}
}
void XMLOutputValueItem::print ( ostream& os, int ntab ) const
{
	size_t len = value.length ();
	openXMLValueTag ( os, ntab, tagName, attr, !len );
	if ( len ) {
		os << value;
		closeXMLTag ( os, tagName );
	}
}
