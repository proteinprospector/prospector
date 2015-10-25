/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_script.cpp                                                 *
*                                                                             *
*  Created    : January 28th 2003                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2003-2008) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lgen_file.h>
#include <lu_getfil.h>
#include <lu_html.h>
#include <lu_script.h>
#include <lu_param_list.h>
using std::ostream;
using std::string;
using std::cout;
using std::endl;

static void perlChdir ( ostream& os, const string& dir );
static string perlEscapePath ( const string& path );
static void dosCD ( ostream& os, const string& dir );

void writeScript ( ParameterList* paramList )
{
	string filename = paramList->getStringValue ( "script_filename", "" );
	paramList->setValue ( "create_script", "0" );
	paramList->setValue ( "script_filename", "" );
	writePerlScript ( paramList, filename );
	init_html ( cout, "Created Script" );
}
void writePerlScript ( const ParameterList* paramList, const string& filename )
{
	GenOFStream ost ( adjustPPOutputPath ( "scripts" ) + filename + ".pl" );
	ost << "$start = time;" << endl;
	perlChdir ( ost, adjustedPPCurrentWorkingDirectory ( genCurrentWorkingDirectory () ) );
	string path = BinDir::instance ().getBinDir () + SLASH + paramList->getProgramBinaryName ();
	ost << "$c = \"" << perlEscapePath ( path ) << "\";" << endl;
	ost << "$c .= \" - \";" << endl;
	paramList->perlFileParameters ( ost );
	ost << "system \"$c\";" << endl;
	ost << "use integer;" << endl;
	ost << "$seconds = time - $start;" << endl;
	ost << "$days = $seconds / 86400;" << endl;
	ost << "$rem = $seconds % 86400;" << endl;
	ost << "$hours = $rem / 3600;" << endl;
	ost << "$rem = $rem % 3600;" << endl;
	ost << "$minutes = $rem / 60;" << endl;
	ost << "$sec = $rem % 60;" << endl;
	ost << "print \"Process Time = $seconds secs\";" << endl;
	ost << "if ( $days || $hours || $minutes ) {" << endl;
	ost << "print \" (\";" << endl;
	ost << "if ( $days ) { print \"$days days \"; }" << endl;
	ost << "if ( $hours ) { print \"$hours hours \"; }" << endl;
	ost << "if ( $minutes ) { print \"$minutes min \"; }" << endl;
	ost << "if ( $sec ) { print \"$sec sec\"; }" << endl;
	ost << "print \")\";" << endl;
	ost << "}" << endl;
	ost << "print \"\\n\";" << endl;
}
static void perlChdir ( ostream& os, const string& dir )
{
	os << "chdir";
	os << " ";
	os << "\"";
	os << perlEscapePath ( dir );
	os << "\"";
	os << ";";
	os << endl;
}
static string perlEscapePath ( const string& path )
{
#ifdef VIS_C
	string ret;

	for ( StringSizeType i = 0 ; i < path.length () ; i++ ) {
		ret += path [i];
		if ( path [i] == '\\' ) ret += '\\';
	}
	return ret;
#else
	return path;
#endif
}
#ifdef VIS_C
void writeBatchFile ( const ParameterList* paramList, const string& filename )
{
	GenOFStream ost ( filename );
	ost << genCurrentDOSDrive () << endl;
	dosCD ( ost, genCurrentWorkingDirectory () );
	ost << BinDir::instance ().getBinDir () + SLASH + paramList->getProgramBinaryName ();
	ost << " - ";
	paramList->batchFileParameters ( ost );
	ost << endl;
	//ost << "pause" << endl;  // Useful if using at /interactive mode
}
static void dosCD ( ostream& os, const string& dir )
{
	os << "cd";
	os << " ";
	os << dir;
	os << endl;
}
#endif
void writeParamsXML ( ParameterList* paramList, const string& programName )
{
	string filename = paramList->getStringValue ( "script_filename", "" );
	paramList->setValue ( "create_params", "0" );
	paramList->setValue ( "script_filename", "" );
	paramList->XMLParameterFile ( MsparamsDir::instance ().getParamPath ( programName + SLASH + filename + ".xml" ) );
	init_html ( cout, "Created Script" );
}
