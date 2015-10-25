/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_prog.cpp                                                   *
*                                                                             *
*  Created    : September 23rd 1997                                           *
*                                                                             *
*  Purpose    : Functions to define argruments for msdigest, msproduct and    *
*               msisotope when called from a results page.                    *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1997-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lu_getfil.h>
#include <lu_html.h>
#include <lu_prog.h>
using std::string;
using std::ostream;
using std::endl;

const ParameterList* ProgramLink::params = 0;

void ProgramLink::openLink ( ostream& os, const string& linkName, int num, bool newWindow )
{
	startJavascript ( os );
	if ( newWindow )
		os << "document.write ( \"<a target=\\\"_blank\\\" href= \\\"\" + ";
	else
		os << "document.write ( \"<a href= \\\"\" + ";
	os << linkName;
	if ( num != -1 ) os << num;
	os << " + \"";
}
void ProgramLink::closeLink ( ostream& os )
{
	os << "</a>";
	os << "\" );";
	os << endl;
	os << "//-->" << endl;
	os << "</script>";
}
string ProgramLink::getURLStart ( const string& program, const string& type )
{
	return BinDir::instance ().getVirtualBinDir () + program + type;
}
void ProgramLink::putForm ( ostream& os, const string& variableName, const string& script )
{
	os << "<form";
	os << " ";
	os << "target=\"_blank\"";
	os << " ";
	os << "id=\"" << variableName << "\"";
	os << " ";
	os << "method=\"post\"";
	os << " ";
	os << "action=\"";
	os << BinDir::instance ().getVirtualBinDir ();
	os << script;
	os << "\"";
	os << ">" << endl;
}
void ProgramLink::putFormLink ( ostream& os, const string& variableName, const string& text )
{
	os << "<a href=\"javascript:void\"";
	os << " ";
	os << "onclick=\"document.getElementById('" << variableName << "').submit();\"";
	os << ">";
	os << text;
	os << "</a>" << endl;
}
