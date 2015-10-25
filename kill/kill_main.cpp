/******************************************************************************
*                                                                             *
*  Program    : kill                                                          *
*                                                                             *
*  Filename   : kill_main.cpp                                                 *
*                                                                             *
*  Created    : December 9th 1999                                             *
*                                                                             *
*  Purpose    : Form entry program for killing processes.                     *
*                                                                             *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1999-2010) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifndef VIS_C
#include <stdexcept>
#endif
#include <algorithm>
#include <lg_stdlib.h>
#include <lg_string.h>
#include <lg_time.h>
#include <lgen_file.h>
#include <lgen_process.h>
#include <lu_getfil.h>
#include <lu_html.h>
#include <lu_prog_par.h>
#include <lu_param_list.h>
using std::sort;
using std::string;
using std::cout;
using std::endl;
using std::runtime_error;

namespace {
struct SortByReverseTime
{
	bool operator () ( const string& strA, const string& strB )
	{
		int a = atoi ( strA.substr ( strA.find_first_of ( "_" ) + 1, 6 ).c_str () );
		int b = atoi ( strB.substr ( strB.find_first_of ( "_" ) + 1, 6 ).c_str () );
		return b < a;
	}
};
void cleanUpLogFile ( int pid )
{
	string baseDir = ppBaseDir () + SLASH + "logs";
	StringVector fList = FileList ( baseDir ).getNameList ();
	StringVector dirList;
	for ( StringVectorSizeType i = 0 ; i < fList.size () ; i++ ) {
		if ( genIsDirectory ( baseDir + SLASH + fList [i] ) ) {
			dirList.push_back ( fList [i] );
		}
	}
	sort ( dirList.rbegin (), dirList.rend () ); // Search newest first
	for ( StringVectorConstIterator j = dirList.begin () ; j != dirList.end () ; j++ ) {
		string curDir = baseDir + SLASH + *j;
		string suffix = "_" + gen_itoa ( pid ) + ".txt";
		StringVector fList2 = FileList ( curDir, "", suffix, false ).getNameList ();
		if ( !fList2.empty () ) {
			sort ( fList2.begin (), fList2.end (), SortByReverseTime () );
			string fName = curDir + SLASH + fList2.front ();
			GenOFStream ost ( fName, std::ios_base::out | std::ios_base::app );
			ParameterList::printXML ( ost, "error_message", "Program terminated by abort search button." );
			ParameterList::printXML ( ost, "end_time", genCurrentTimeAndDateString () );
			ost << "</program_log>" << endl;
			ost.close ();
			genRename ( fName, fName.substr ( 0, fName.length () - 4 ) + ".xml" );
			return;
		}
	}
}
}
int main ( int argc, char *argv[] )
{
	initialiseProspector ();
	ParameterList pList ( argc, argv );
	try {
		MSProgramParameters params ( &pList );
		genUnlink ( pList.getStringVectorValue ( "delete" ) );
		int pid = pList.getIntValue ( "pid" );
		string terminationInfo = "Search Terminated";
#ifdef VIS_C
		bool running = isProcessRunning ( pid );
		if ( running ) {
			string command ( getSystemCall ( "kill.exe" ) );
			command += " ";
			command += gen_itoa ( pid );
			genSystem ( command, "", true );
			genSleep ( 500 );	// Wait 1/2 sec
			if ( isProcessRunning ( pid ) ) {
				terminationInfo = "Search Still Running";
			}
		}
#else
		genSystem ( string ( "kill -9 " ) + gen_itoa ( pid ), "", true );
#endif
		cleanUpLogFile ( pid );
		init_html ( cout, terminationInfo );
		cout << "<input type=\"button\" value=\"Search Form\" onclick=\"history.go(-2)\">" << endl;
	}
	catch ( runtime_error e ) {
		pList.writeLogError ( e.what () );
	}
	return 0;
}
