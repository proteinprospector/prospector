/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_btag_submit.cpp                                            *
*                                                                             *
*  Created    : September 10th 2007                                           *
*                                                                             *
*  Purpose    : Functions for dealing with multiply charged data.             *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2007-2011) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef MYSQL_DATABASE
#include <lgen_process.h>
#include <ld_init.h>
#include <lu_getfil.h>

using std::string;

bool btagSubmit ( const string& searchKey, const string& host )
{
	string executable;
	string params;
	string msSearchParams = "-k " + searchKey;
	if ( InfoParams::instance ().getBoolValue ( "multi_process" ) ) {
#ifdef VIS_C
		executable = InfoParams::instance ().getStringValue ( "mpi_run" );
		params = InfoParams::instance ().getStringValue ( "mpi_args" );
		params += " ";
		params += getSystemCall ( "mssearchmpi.cgi" );
		params += " ";
		params += msSearchParams;
#else
		executable = getSystemCall ( "mssearchmpi.pl" );
		params = searchKey;
#endif
	}
	else {
#ifdef VIS_C
		executable = getSystemCall ( "mssearch.cgi" );
		params = msSearchParams;
#else
		executable = getSystemCall ( "mssearch.pl" );
		params = searchKey;
#endif
	}
	int pid = genCreateProcess ( executable, params );
	if ( pid == 0 ) {
		return false;
	}
	else {
		MySQLPPSDDBase::instance ().setJobSearching ( searchKey, host, pid );
		return true;
	}
}
#endif
