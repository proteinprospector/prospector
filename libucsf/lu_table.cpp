/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_table.cpp                                                  *
*                                                                             *
*  Created    : March 25rd 2003                                               *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2003-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lg_io.h>
#include <lu_table.h>
using std::ostream;
using std::string;
using std::endl;

void tableStart ( ostream& os, bool border, const string& align, const string& cellspacing, const string& id )
{
	os << "<table";
	if ( border ) {
		os << " border=\"border\"";
		if ( cellspacing.empty () ) {
			os << " cellspacing=\"3\"";
		}
	}
	if ( !align.empty () ) {
		os << " align=\"" << align << "\"";
	}
	if ( !cellspacing.empty () ) {
		os << " cellspacing=\"" << cellspacing << "\"";
	}
	if ( !id.empty () ) {
		os << " id=\"" << id << "\"";
	}
	os << ">" << endl;
}
void tableEnd ( ostream& os )
{
	os << "</table>" << endl;
}
void tableRowStart ( ostream& os )
{
	os << "<tr>" << endl;
}
void tableRowStartDiv ( ostream& os, const string& id, bool open )
{
	os << "<tr id=\"";
	os << id;
	os << "\" style=\"";
	if ( open )	os << "display:block";
	else		os << "display:none";
	os << "\">" << endl;
}
void tableRowEnd ( ostream& os )
{
	os << "</tr>" << endl;
}
void tableEmptyRow ( ostream& os )
{
	os << "<tr><td>&nbsp;</td></tr>" << endl;
}
void tableEmptyNRows ( ostream& os, int n )
{
	for ( int i = 0 ; i < n ; i++ ) tableEmptyRow ( os );
}
void tableHeader ( ostream& os, const string& str, const string& styleID, const string& align, bool nowrap, int colspan, int rowspan )
{
	tableHeaderStart ( os, styleID, align, nowrap, colspan, rowspan );
	os << str << endl;
	tableHeaderEnd ( os );
}
void tableHeaderStart ( ostream& os, const string& styleID, const string& align, bool nowrap, int colspan, int rowspan )
{
	os << "<th";
	if ( !styleID.empty ()  ) {
		os << " class=\"" << styleID << "\"";
	}
	if ( !align.empty () ) {
		os << " align=\"" << align << "\"";
	}
	if ( nowrap ) {
		os << " nowrap=\"nowrap\"";
	}
	if ( colspan ) {
		os << " colspan=\"" << colspan << "\"";
	}
	if ( rowspan ) {
		os << " rowspan=\"" << rowspan << "\"";
	}
	os << ">" << endl;
}
void tableHeaderEnd ( ostream& os )
{
	os << "</th>" << endl;
}
void tableData ( ostream& os, const string& str, const string& styleID )
{
	os << "<td";
	if ( !styleID.empty ()  ) {
		os << " class=\"" << styleID << "\"";
	}
	os << ">";
	os << str;
	os << "</td>" << endl;
}
void tableDataStart ( ostream& os, const string& styleID, const string& align, bool nowrap, int colspan, int rowspan )
{
	os << "<td";
	if ( !styleID.empty ()  ) {
		os << " class=\"" << styleID << "\"";
	}
	if ( !align.empty () ) {
		os << " align=\"" << align << "\"";
	}
	if ( nowrap ) {
		os << " nowrap=\"nowrap\"";
	}
	if ( colspan ) {
		os << " colspan=\"" << colspan << "\"";
	}
	if ( rowspan ) {
		os << " rowspan=\"" << rowspan << "\"";
	}
	os << ">" << endl;
}
void tableDataEnd ( ostream& os )
{
	os << "</td>" << endl;
}
void tableCellStart ( ostream& os, const string& styleID, const string& align, bool nowrap, int colspan, int rowspan )
{
	os << "<td";
	if ( !styleID.empty ()  ) {
		os << " class=\"" << styleID << "\"";
	}
	if ( !align.empty () ) {
		os << " align=\"" << align << "\"";
	}
	if ( nowrap ) {
		os << " nowrap=\"nowrap\"";
	}
	if ( colspan ) {
		os << " colspan=\"" << colspan << "\"";
	}
	if ( rowspan ) {
		os << " rowspan=\"" << rowspan << "\"";
	}
	os << ">" << endl;
}
void tableCellEnd ( ostream& os )
{
	os << "</td>" << endl;
}
void tableCellSigFig ( ostream& os, double number, int sigFig, bool bold, const string& styleID, int colspan, int rowspan )
{
	if ( bold )	os << "<th";
	else		os << "<td";
	if ( !styleID.empty ()  ) {
		os << " class=\"" << styleID << "\"";
	}
	if ( number < 0.0 ) os << " nowrap=\"nowrap\"";
	if ( colspan ) {
		os << " colspan=\"" << colspan << "\"";
	}
	if ( rowspan ) {
		os << " rowspan=\"" << rowspan << "\"";
	}
	os << ">";

	genPrintSigFig ( os, number, sigFig );

	if ( bold )	os << "</th>";
	else		os << "</td>";
	os << endl;
}
void tableCell ( ostream& os, double number, int precision, bool bold, const string& styleID, int colspan, int rowspan )
{
	if ( bold )	os << "<th";
	else		os << "<td";
	if ( !styleID.empty ()  ) {
		os << " class=\"" << styleID << "\"";
	}
	if ( number < 0.0 ) os << " nowrap=\"nowrap\"";
	if ( colspan ) {
		os << " colspan=\"" << colspan << "\"";
	}
	if ( rowspan ) {
		os << " rowspan=\"" << rowspan << "\"";
	}
	os << ">";

	genPrint ( os, number, precision );

	if ( bold )	os << "</th>";
	else		os << "</td>";
	os << endl;
}
void tableCellRange ( ostream& os, double number1, double number2, int precision, bool bold, const string& styleID )
{
	if ( bold )	os << "<th";
	else		os << "<td";
	if ( !styleID.empty ()  ) {
		os << " class=\"" << styleID << "\"";
	}
	os << " nowrap=\"nowrap\"";
	os << ">";

	genPrint ( os, number1, precision );
	os << " - ";
	genPrint ( os, number2, precision );

	if ( bold )	os << "</th>";
	else		os << "</td>";
	os << endl;
}
void tableCell ( ostream& os, const string& str, bool nobr, bool bold, const string& styleID, int colspan, int rowspan )
{
	if ( bold )	os << "<th";
	else		os << "<td";
	if ( !styleID.empty ()  ) {
		os << " class=\"" << styleID << "\"";
	}
	if ( nobr ) os << " nowrap=\"nowrap\"";
	if ( colspan ) {
		os << " colspan=\"" << colspan << "\"";
	}
	if ( rowspan ) {
		os << " rowspan=\"" << rowspan << "\"";
	}
	os << ">";

	os << str;

	if ( bold )	os << "</th>";
	else		os << "</td>";
	os << endl;
}
void tableCell ( ostream& os, char number, bool bold, const string& styleID, int colspan, int rowspan )
{
	if ( bold )	os << "<th";
	else		os << "<td";
	if ( !styleID.empty ()  ) {
		os << " class=\"" << styleID << "\"";
	}
	if ( colspan ) {
		os << " colspan=\"" << colspan << "\"";
	}
	if ( rowspan ) {
		os << " rowspan=\"" << rowspan << "\"";
	}
	os << ">";

	os << number;

	if ( bold )	os << "</th>";
	else		os << "</td>";
	os << endl;
}
void tableCell ( ostream& os, int number, bool bold, const string& styleID, int colspan, int rowspan )
{
	if ( bold )	os << "<th";
	else		os << "<td";
	if ( !styleID.empty ()  ) {
		os << " class=\"" << styleID << "\"";
	}
	if ( colspan ) {
		os << " colspan=\"" << colspan << "\"";
	}
	if ( rowspan ) {
		os << " rowspan=\"" << rowspan << "\"";
	}
	os << ">";

	os << number;

	if ( bold )	os << "</th>";
	else		os << "</td>";
	os << endl;
}
void tableCell ( ostream& os, unsigned int number, bool bold, const string& styleID, int colspan, int rowspan )
{
	if ( bold )	os << "<th";
	else		os << "<td";
	if ( !styleID.empty ()  ) {
		os << " class=\"" << styleID << "\"";
	}
	if ( colspan ) {
		os << " colspan=\"" << colspan << "\"";
	}
	if ( rowspan ) {
		os << " rowspan=\"" << rowspan << "\"";
	}
	os << ">";

	os << number;

	if ( bold )	os << "</th>";
	else		os << "</td>";
	os << endl;
}
void tableEmptyCell ( ostream& os, const string& styleID, int colspan, int rowspan )
{
	os << "<td";
	if ( !styleID.empty ()  ) {
		os << " class=\"" << styleID << "\"";
	}
	if ( colspan ) {
		os << " colspan=\"" << colspan << "\"";
	}
	if ( rowspan ) {
		os << " rowspan=\"" << rowspan << "\"";
	}
	os << " />";
	os << endl;
}
void tableEmptyNCells ( ostream& os, int n, const string& styleID )
{
	for ( int i = 0 ; i < n ; i++ ) tableEmptyCell ( os, styleID );
}
void divStart ( ostream& os, const string& name, bool open )
{
	os << "<div id=\"";
	os << name;
	os << "\" style=\"";
	if ( open )	os << "display:block";
	else		os << "display:none";
	os << "\">" << endl;
}
void divEnd ( ostream& os )
{
	os << "</div>" << endl;
}
void spanStart ( ostream& os, const string& name, bool open )
{
	os << "<span id=\"";
	os << name;
	os << "\" style=\"";
	if ( open )	os << "display:inline";
	else		os << "display:none";
	os << "\">" << endl;
}
void spanEnd ( ostream& os )
{
	os << "</span>" << endl;
}
