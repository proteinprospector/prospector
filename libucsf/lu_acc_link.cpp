/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_acc_link.cpp                                               *
*                                                                             *
*  Created    : September 5th 2001                                            *
*                                                                             *
*  Purpose    : Functions for dealing with accession number links.            *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lg_string.h>
#include <lu_getfil.h>
#include <lu_acc_link.h>
#include <lu_check_db.h>
#include <lu_html.h>
using std::string;
using std::ostream;
using std::endl;

int AccessionNumberLinkInfo::num = 0;
StringVector AccessionNumberLinkInfo::aLink;
AccessionNumberLinkInfo::AccessionNumberLinkInfo ()
{
	init ( "acclinks.txt" );
}
AccessionNumberLinkInfo::AccessionNumberLinkInfo ( const string& linksName )
{
	init ( linksName );
}
void AccessionNumberLinkInfo::init ( const string& fName )
{
	num++;
	ind = num;
	int numFilenamePrefixes;
	char* fileInfo = getParamsFileInfo ( fName, &numFilenamePrefixes );
 
	for ( int i = 0 ; i < numFilenamePrefixes ; i++ ) {
		filenamePrefix.push_back ( ( i == 0 ) ? strtok ( fileInfo, "\t " ) : strtok ( NULL, "\t " ) );
		accessionNumberLink.push_back ( strtok ( NULL, "\n" ) );
	}
	delete [] fileInfo;
}
string AccessionNumberLinkInfo::getAccessionNumberLink ( const string& database ) const
{
	if ( isFullyDecoyDatabase ( database ) )
		return "decoy";
	else {
		int index = longestMatchingPrefix ( database );
		if ( index == -1 )
			return ( "" );
		else
			return ( accessionNumberLink [index] );
	}
}
int AccessionNumberLinkInfo::longestMatchingPrefix ( const string& str ) const
{
	int index;
	int maxlen = 0;
	for ( StringVectorSizeType i = 0 ; i < filenamePrefix.size () ; i++ ) {
		if ( strstr ( str.c_str (), filenamePrefix [i].c_str () ) == str.c_str () ) {
			int len = filenamePrefix [i].length ();
			if ( len > maxlen ) {
				maxlen = len;
				index = i;
			}
		}
	}
	if ( maxlen ) return ( index );
	else return ( -1 );
}
void AccessionNumberLinkInfo::write ( ostream& os, const string& aNum, int num, bool link )
{
	os << endl;
	if ( num >= aLink.size () || aLink [num].empty () )	os << aNum << endl;
	else if ( aLink [num] == "decoy" || aNum [0] == '-' )	os << "decoy" << endl;
	else {
		if ( link ) {
			startJavascript ( os );
			os << "document.writeln ( \"<a target=\\\"_blank\\\" href=\\\"\" + accessLink" << num+1 << " + \"";
			os << aNum;
			os << "\\\"> ";
			os << aNum;
			os << " </a>";
			os << "\" );";
			os << endl;
			endJavascript ( os );
		}
		else os << aNum << endl;
	}
}
void AccessionNumberLinkInfo::printHTML ( ostream& os, const string& database ) const
{
	string a = getAccessionNumberLink ( database );
	if ( !a.empty () && a != "decoy" ) os << "accessLink" << ind << "=\"" << a << "\";\n";
	aLink.resize ( ind );
	aLink [ind-1] = a;
}
void AccessionNumberLinkInfo::write2 ( ostream& os, const string& aNum, bool link ) const
{
	write ( os, aNum, ind-1, link );
}
