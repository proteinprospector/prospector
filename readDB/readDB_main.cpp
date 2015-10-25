/******************************************************************************
*                                                                             *
*  Program    : readDB                                                        *
*                                                                             *
*  Filename   : readDB_main.cpp                                               *
*                                                                             *
*  Created    : March 12th 2007                                               *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2008-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifndef VIS_C
#include <stdexcept>
#endif
#include <ld_init.h>
#include <lu_getfil.h>
#include <lu_param_list.h>
#include <lu_prog_par.h>
#include <lu_prog.h>
#include <lu_html.h>
#include <lu_html_form.h>
using std::ostream;
using std::cout;
using std::string;
using std::runtime_error;

int main ( int argc, char** argv )
{
	initialiseProspector ();

	ParameterList paramList ( argc, argv );
	try {
		ostream& os = cout;
		ProgramLink::setParams ( &paramList );
		init_html ( os, "Database Table" );
		string type = paramList.getStringValue ( "output_type", "HTML" );
		string s = paramList.getStringValue ( "sql_string" );
		if ( s.empty () ) s = paramList.getStringValue ( "sql" );
		try {
			if ( type == "HTML" ) {
				if ( s == "User Info" )
					MySQLPPSDDBase::instance ().printUserInfoHTML ( os );
				else
					MySQLPPSDDBase::instance ().printSQLResultsHTML ( os, s );
			}
			else {
				if ( s == "User Info" )
					MySQLPPSDDBase::instance ().printUserInfoTabDelimitedText ( os );
				else
					MySQLPPSDDBase::instance ().printSQLResultsTabDelimitedText ( os, s );
			}
		}
		catch ( runtime_error e ) {		// Catch database login problems
			ErrorHandler::genError ()->error ( e );
		}
		bool refresh = paramList.getBoolValue ( "refresh" );
		if ( refresh ) {
			int jobStatusRefresh = paramList.getIntValue ( "refresh_time", 5 );
			ParameterList pList ( "readDB", false, false, false, false, false );
			pList.addName ( "refresh", paramList.getStringValue ( "refresh" ) );
			pList.addName ( "refresh_time", paramList.getStringValue ( "refresh_time" ) );
			pList.addName ( "sql", paramList.getStringValue ( "sql" ) );
			pList.addName ( "sql_string", paramList.getStringValue ( "sql_string" ) );
			pList.addName ( "output_type", paramList.getStringValue ( "output_type" ) );
			refreshJavascript ( os, jobStatusRefresh*1000, pList.getURL (), false, true );
		}
		printHTMLFooter ( os, "2008" );
	}
	catch ( runtime_error e ) {
		paramList.writeLogError ( e.what () );
	}
	return 0;
}
