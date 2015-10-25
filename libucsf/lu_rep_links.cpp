/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_rep_links.cpp                                              *
*                                                                             *
*  Created    : April 29th 2014                                               *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2014-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lg_string.h>
#include <lu_param_list.h>
#include <lu_prog.h>
#include <lu_rep_links.h>

using std::string;
using std::ostream;
using std::endl;


ReportLinkProgram::ReportLinkProgram ( const string& linkProgram, const string& linkName ) :
	linkProgram ( linkProgram ),
	linkName ( linkName )
{
}
void ReportLinkProgram::printHTML ( ostream& os, const ParameterList* pList ) const
{
	os << linkName << "=\"";
	os << ProgramLink::getURLStart ( linkProgram );
	os << "?";
	const_cast <ParameterList*> (pList)->removeName ( "page" );
	pList->copyToCGI ( os );
	os << "\";\n";
}
void ReportLinkProgram::write ( ostream& os, const string& p, const string& l ) const
{
	ProgramLink::openLink ( os, linkName, -1, false );
	os << p;
	os << "\\\">";
	os << l;
	ProgramLink::closeLink ( os );
}

ReportLinks::ReportLinks ( const string& url, int rowsPerPage, int pageNumber, int numDataLines ) :
	url ( url ),
	rowsPerPage ( rowsPerPage ),
	pageNumber ( pageNumber )
{
	int numPages = numDataLines / rowsPerPage;
	if ( numDataLines % rowsPerPage ) numPages++;

	int numDisplayed = 10;
	int range = ( (pageNumber-1) / numDisplayed ) + 1;
	if ( range > 1 ) {
		vpss.push_back ( getPair ( "1", "First" ) );
		string sPreviousPage = gen_itoa ( numDisplayed * (range-1) );
		vpss.push_back ( getPair ( sPreviousPage, "Previous" ) );
	}
	int endRow = range*numDisplayed;
	int startRow = endRow+1-numDisplayed;
	for ( int i = startRow ; i <= endRow ; i++ ) {
		if ( i > numPages ) break;
		string sPageNum = gen_itoa ( i );
		vpss.push_back ( getPair ( sPageNum, sPageNum ) );
	}
	int topRange = ( (numPages-1) / numDisplayed ) + 1;
	if ( range != topRange ) {
		string sNextPage = gen_itoa ( endRow + 1 );
		vpss.push_back ( getPair ( sNextPage, "Next" ) );
		string sLastPage = gen_itoa ( numPages );
		vpss.push_back ( getPair ( sLastPage, "Last" ) );
	}
}
PairStringString ReportLinks::getPair ( const string& page, const string& name ) const
{
	string s = url;
	s += "page=";
	s += page;
	return PairStringString ( s, name );
}
void ReportLinks::printHTML ( ostream& os, const ReportLinkProgram& link ) const
{
	for ( int i = 0 ; i < vpss.size () ; i++ ) {
		if ( i == 0 ) {
			os << "<p>" << endl;
		}
		link.write ( os, vpss [i].first, vpss [i].second );
		os << endl;
		if ( i == vpss.size ()-1 ) {
			os << "</p>" << endl;
		}
	}
}
