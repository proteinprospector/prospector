/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_version.cpp                                                *
*                                                                             *
*  Created    : February 18th 1997                                            *
*                                                                             *
*  Purpose    : Function to define and print version numbers on program       *
*               outputs.                                                      *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1997-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lg_time.h>
#include <lgen_error.h>
#include <lu_getfil.h>
using std::string;
using std::ostream;
using std::endl;
using std::cout;
using std::getline;
/* 
	The system and program version number is made up of three parts:
	1). Major revisions.
	2). Minor revisions.
	3). Any change including recompilation of the underlying libraries.
	The numbers only need to be changed if a new release is made. We might also want to signify
	that we aren't working with a release version using something like "1.0.0 debug". 
*/
static string packageName = "ProteinProspector";
static string packageVersionNumber = Version::instance ().getVersion ();
static string packageVersionModifier = "";
static string institutionName = "The Regents of the University of California.";
static string currentYear = genCurrentYear ();

static void printStartYear ( ostream& os, const string& programName );

string getVersionFromPPXMLFile ( const string& filename )
{
	GenIFStream ist ( filename );
	string line;
	getline ( ist, line );
	getline ( ist, line );	// Version is stored on the second line
	int start = line.find ( "Version" ) + 8;
	if ( start != string::npos ) {
		int end = line.find ( "?", start );
		if ( end != string::npos ) {
			return line.substr ( start, end - start );
		}
	}
	return "";
}
void printProgramInformationHTML ( ostream& os, const string& programName )
{
	os << endl;
	os << endl;

	os << "</div>" << endl;		// centerbody
	os << "<div id=\"footer\">" << endl;
		os << programName;
		os << " in " << packageName << " " << packageVersionNumber << packageVersionModifier;
		os << endl;
		os << "<br />" << endl;

		os << "&#169; Copyright (";
		printStartYear ( os, programName );
		os << "-";
		os << currentYear;
		os << ") ";
		os << institutionName;
		os << endl;
	os << "</div>" << endl;

	os << "</div>" << endl;		// content
	os << "</div>" << endl;		// prop
	os << "</body>" << endl;
	os << "</html>" << endl;
}
void printProgramInformation ( const string& programName )
{
	cout << programName;
	cout << " " << packageName << " " << packageVersionNumber << packageVersionModifier;
	cout << "\n";

	cout << "Copyright (";
	printStartYear ( cout, programName );
	cout << "-";
	cout << currentYear;
	cout << ") ";
	cout << institutionName;
	cout << "\n";
}
static void printStartYear ( ostream& os, const string& programName )
{
	if ( programName == "FA-Index" )			{os << "1995";}
	else if ( programName == "MS-Digest" ) 		{os << "1995";}
	else if ( programName == "MS-Fit" )			{os << "1995";}
	else if ( programName == "MS-Product" )		{os << "1995";}
	else if ( programName == "MS-Seq" )			{os << "1995";}
	else if ( programName == "MS-Tag" )			{os << "1995";}
	else if ( programName == "DB-Stat" )		{os << "1996";}
	else if ( programName == "MS-Pattern" )		{os << "1996";}
	else if ( programName == "MS-Comp" )		{os << "1996";}
	else if ( programName == "MS-Isotope" )		{os << "1998";}
	else if ( programName == "MS-Homology" )	{os << "1999";}
	else if ( programName == "MS-Bridge" )		{os << "1999";}
	else if ( programName == "MS-NonSpecific" )	{os << "2002";}
	else if ( programName == "MS-MultiSample" )	{os << "2002";}
	else if ( programName == "BatchTag" )		{os << "2002";}
	else if ( programName == "MS-Display" )		{os << "2002";}
	else if ( programName == "Search Compare" )	{os << "2003";}
	else if ( programName == "MS-Viewer" )		{os << "2011";}
	else if ( programName == "MS-Filter" )		{os << "2012";}
	else ErrorHandler::genError ()->error ( "Function version_number: invalid program name." );
}
