/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_tofil_par.cpp                                              *
*                                                                             *
*  Created    : October 18th 2001                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lgen_file.h>
#include <lu_getfil.h>
#include <lu_tofil_par.h>
#include <lu_param_list.h>
using std::string;
using std::ostream;
using std::endl;

ToFileParameters::ToFileParameters ( const ParameterList* params, const string& programName ) :
	toFile			( params->getBoolValue		( "results_to_file", false ) ),
	outputFilename	( params->getStringValue	( "output_filename", "" ) ),
	scriptName		( params->getStringValue	( "script", "" ) ),
	scriptType		( params->getStringValue	( "script_type", "" ) ),
	script			( scriptName != "" ),
	outputPath		( string ( "results" ) + SLASH + string ( programName ) )
{
}
void ToFileParameters::printOutputFileLink ( ostream& os, const string& linkText, const string& fileSuffix ) const
{
	string filename = outputFilename + fileSuffix;
	os << "<a href=\"";
	os << ParameterList::getVirtualDir ();
	os << genTranslateSlash ( outputPath );
	os << "/";
	os << filename;
	os << "\">";
	os << linkText;
	os << " (" << filename << ")";
	os << "</a>";
	os << "<br />";
	os << endl;
}
bool ToFileParameters::isRandomSearch () const
{
	string s = getOutputFilename ();
	return s.find ( ".exp." ) != string::npos;
}
string ToFileParameters::getOutputPath ( const string& suffix ) const
{
	return adjustPPOutputPath ( outputPath ) + outputFilename + suffix;
}
