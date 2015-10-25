/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_spec_id.cpp                                                *
*                                                                             *
*  Created    : December 13th 2003                                            *
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
#include <cstdlib>
#include <iomanip>
#include <lu_delim.h>
#include <lu_table.h>
#include <lu_html_form.h>
#include <lu_spec_id.h>
#include <lu_cgi_val.h>
#include <lu_param_list.h>
#include <lu_sctag_link.h>
using std::string;
using std::ostream;
using std::stringstream;

static string getNextSubstring ( const string& str, int& start, const string& separator )
{
	if ( start != string::npos ) {
		int end = str.find ( separator, start );
		string subStr = str.substr ( start, end - start );
		start = ( end == string::npos ) ? string::npos : end + separator.length ();
		return subStr;
	}
	else return "";
}
static int getNextInt ( const string& str, int& start, const string& separator )
{
	return atoi ( (getNextSubstring ( str, start, separator )).c_str () );
}
SpecID::SpecID ( const ParameterList* params ) :
	fraction( params->getIntValue	( "fraction", -1 ) ),
	spot	( params->getStringValue( "spot_number", "1" ) ),
	run		( params->getIntValue	( "run", 1 ) ),
	specNum	( params->getIntValue	( "spectrum_number", 1 ) ),
	msmsInfo( params->getStringValue( "msms_info", "" ) )
{
}
SpecID::SpecID ( const string& specID )
{
	int start = 0;
	fraction = getNextInt ( specID, start, "-" );
	spot = getNextSubstring ( specID, start, "-" );
	run = getNextInt ( specID, start, "-" );
	specNum = getNextInt ( specID, start, "-" );
	msmsInfo = getNextSubstring ( specID, start, "-" );
}
void SpecID::putCGI ( ostream& os ) const
{
	printCGI ( os, "fraction", fraction );
	printCGIString ( os, "spot_number", spot );
	printCGI ( os, "run", run );
	printCGI ( os, "spectrum_number", specNum );
	if ( !msmsInfo.empty () ) printCGIString ( os, "msms_info", msmsInfo );
}
void SpecID::printHTMLHidden ( ostream& os ) const
{
	printHTMLFORMHidden ( os, "fraction", fraction );
	printHTMLFORMHidden ( os, "spot_number", spot );
	printHTMLFORMHidden ( os, "run", run );
	printHTMLFORMHidden ( os, "spectrum_number", specNum );
	if ( !msmsInfo.empty () ) printHTMLFORMHidden ( os, "msms_info", msmsInfo );
}
string SpecID::getCommandLineNVPair () const
{
	string s;
	s += ::getCommandLineNVPair ( "fraction", fraction );
	s += " ";
	s += ::getCommandLineNVPair ( "spot_number", spot );
	s += " ";
	s += ::getCommandLineNVPair ( "run", run );
	s += " ";
	s += ::getCommandLineNVPair ( "spectrum_number", specNum );
	s += " ";
	if ( !msmsInfo.empty () ) {
		s += ::getCommandLineNVPair ( "msms_info", msmsInfo );
		s += " ";
	}
	return s;
}
void SpecID::printDelimitedHeader ( ostream& os, bool r, bool s )
{
	delimitedHeader ( os, "RT" );
	if ( r )	delimitedHeader ( os, "Run" );
	if ( s )	delimitedHeader ( os, "Spectrum" );
}
void SpecID::printDelimitedMSMSInfoHeader ( ostream& os )
{
	delimitedHeader ( os, "MSMS Info" );
}
void SpecID::printDelimitedCell ( ostream& os, bool r, bool s ) const
{
	delimitedCell ( os, convertTime ( spot ) );
	if ( r )	delimitedCell ( os, run );
	if ( s )	delimitedCell ( os, specNum );
}
void SpecID::printDelimitedMSMSInfoCell ( ostream& os ) const
{
	delimitedCell ( os, msmsInfo );
}
void SpecID::printDelimitedEmptyCell ( ostream& os, bool r, bool s )
{
	delimitedEmptyCell ( os );
	if ( r )	delimitedEmptyCell ( os );
	if ( s )	delimitedEmptyCell ( os );
}
void SpecID::printTableHeader ( ostream& os, bool r, bool s, const string& styleID, int rowspan )
{
	tableHeader ( os, "RT", styleID, "", false, 0, rowspan );
	if ( r )	tableHeader ( os, "R", styleID, "", false, 0, rowspan );
	if ( s )	tableHeader ( os, "#", styleID, "", false, 0, rowspan );
}
void SpecID::printTableMSMSInfoHeader ( ostream& os, const string& styleID, int rowspan )
{
	tableHeader ( os, "MSMS Info", styleID, "", false, 0, rowspan );
}
void SpecID::printTableCell ( ostream& os, bool r, bool s, const string& styleID, int colspan, int rowspan ) const
{
	tableCell ( os, convertTime ( spot ), true, false, styleID, colspan, rowspan );
	if ( r )	tableCell ( os, run, false, styleID, colspan, rowspan );
	if ( s )	tableCell ( os, specNum, false, styleID, colspan, rowspan );
}
void SpecID::printTableMSMSInfoCell ( ostream& os, const string& styleID, int colspan, int rowspan ) const
{
	tableCell ( os, msmsInfo, true, false, styleID, colspan, rowspan );
}
#ifdef BATCHTAG
void SpecID::printTableCell2 ( ostream& os, bool r, bool s, const string& styleID, const SCMSTagLink& smtl, const string& searchKey, int colspan, int rowspan ) const
{
	tableCellStart ( os, styleID, "", false, colspan, rowspan );
	smtl.write ( os, searchKey, getSpecID (), convertTime ( spot ) );
	tableCellEnd ( os );
	if ( r )	tableCell ( os, run, false, styleID, colspan, rowspan );
	if ( s )	tableCell ( os, specNum, false, styleID, colspan, rowspan );
}
#endif
void SpecID::printTableEmptyCell ( ostream& os, bool r, bool s, const string& styleID )
{
	tableEmptyCell ( os, styleID );
	if ( r )	tableEmptyCell ( os, styleID );
	if ( s )	tableEmptyCell ( os, styleID );
}
string SpecID::convertTime ( const string& oldTime )
{
	if ( oldTime.length () > 4 ) {
		if ( oldTime.substr ( 0, 2 ) == "PT" && oldTime [oldTime.length () - 1] == 'S' ) {
			string seconds = oldTime.substr ( 2, oldTime.length () - 3 );
			int iseconds = atoi ( seconds.substr ( 0, seconds.find_first_of ( "." ) ).c_str () );
			stringstream sstr;
			sstr << iseconds / 60 << '.';
			sstr << std::setw ( 2 ) << std::setfill ( '0' ) << iseconds % 60;
			return sstr.str ();
		}
	}
	return oldTime;
}
