/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lgen_process.cpp                                              *
*                                                                             *
*  Created    : June 23rd 2006                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2006-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <cstdio>
#include <lgen_error.h>
#include <lgen_file.h>
#include <lgen_process.h>
#ifdef VIS_C
#include <windows.h>
#include <Psapi.h>
#else
#include <signal.h>
#include <unistd.h>
#include <algorithm>
#include <lg_string.h>
#include <lg_stdlib.h>
#include <lu_getfil.h>
#endif
using std::string;
#ifndef VIS_C
using std::istringstream;
using std::find;
#endif

#ifdef VIS_C
namespace {
bool SetPrivilege ( HANDLE hToken, LPCTSTR Privilege, BOOL bEnablePrivilege ) 
{ 
	TOKEN_PRIVILEGES tp = { 0 };	// Initialize everything to zero
	 
	LUID luid; 
	DWORD cb = sizeof(TOKEN_PRIVILEGES); 
	if ( !LookupPrivilegeValue( NULL, Privilege, &luid ) ) {
		ErrorHandler::genError ()->error ( "LookupPrivilegeValue failed.\n" );
		return false;
	}
	tp.PrivilegeCount = 1; 
	tp.Privileges[0].Luid = luid; 
	if ( bEnablePrivilege ) { 
		tp.Privileges[0].Attributes = SE_PRIVILEGE_ENABLED; 
	}
	else { 
		tp.Privileges[0].Attributes = 0; 
	} 
	AdjustTokenPrivileges ( hToken, FALSE, &tp, cb, NULL, NULL ); 
	DWORD lastErr = GetLastError ();
    if ( lastErr != ERROR_SUCCESS ) {
		sprintf ( gen_error_message, "AdjustTokenPrivileges 1st pass failed (error = %d).\n", lastErr ); 
		ErrorHandler::genError ()->error ( gen_error_message );
		return false;
	}
	return true;
}
/*		2nd function that apparently does the same thing
bool SetPrivilege (
    HANDLE hToken,          // token handle
    LPCTSTR Privilege,      // Privilege to enable/disable
    BOOL bEnablePrivilege   // TRUE to enable.  FALSE to disable
    )
{
    TOKEN_PRIVILEGES tp;
    LUID luid;
    TOKEN_PRIVILEGES tpPrevious;
    DWORD cbPrevious=sizeof(TOKEN_PRIVILEGES);

    if ( !LookupPrivilegeValue( NULL, Privilege, &luid ) ) {
		ErrorHandler::genError ()->error ( "LookupPrivilegeValue failed.\n" );
		return false;
	}
    
    tp.PrivilegeCount           = 1;	// first pass.  get current privilege setting
    tp.Privileges[0].Luid       = luid;
    tp.Privileges[0].Attributes = 0;
    AdjustTokenPrivileges ( hToken, FALSE, &tp, sizeof(TOKEN_PRIVILEGES), &tpPrevious, &cbPrevious );
	DWORD lastErr = GetLastError ();
    if ( lastErr != ERROR_SUCCESS ) {
		sprintf ( gen_error_message, "AdjustTokenPrivileges 1st pass failed (error = %d).\n", lastErr ); 
		ErrorHandler::genError ()->error ( gen_error_message );
		return false;
	}

    tpPrevious.PrivilegeCount       = 1;	// second pass.  set privilege based on previous setting
    tpPrevious.Privileges[0].Luid   = luid;

    if ( bEnablePrivilege ) {
        tpPrevious.Privileges[0].Attributes |= (SE_PRIVILEGE_ENABLED);
    }
    else {
        tpPrevious.Privileges[0].Attributes ^= (SE_PRIVILEGE_ENABLED & tpPrevious.Privileges[0].Attributes);
    }

    AdjustTokenPrivileges ( hToken, FALSE, &tpPrevious, cbPrevious, NULL, NULL );

	lastErr = GetLastError ();
    if ( lastErr != ERROR_SUCCESS ) {
		sprintf ( gen_error_message, "AdjustTokenPrivileges 2nd pass failed (error = %d).\n", lastErr ); 
		ErrorHandler::genError ()->error ( gen_error_message );
		return false;
	}
    return true;
}
*/
PairStringInt getProcessNameAndID ( DWORD processID )
{
    TCHAR szProcessName [MAX_PATH] = TEXT("");

    // Get a handle to the process.

    HANDLE hProcess = OpenProcess ( PROCESS_QUERY_INFORMATION | PROCESS_VM_READ, FALSE, processID );

    if ( hProcess != NULL ) {		// Get the process name.
        HMODULE hMod;
        DWORD cbNeeded;

        if ( EnumProcessModules ( hProcess, &hMod, sizeof(hMod), &cbNeeded) ) {
            GetModuleBaseName ( hProcess, hMod, szProcessName, sizeof(szProcessName)/sizeof(TCHAR) );
        }
    }
    // Print the process name and identifier.

	PairStringInt p;
	p.first = szProcessName;
	p.second = processID;

    CloseHandle( hProcess );
	return p;
}
}
#endif

// External functions start here

bool killProcess ( int pid )
{
	bool flag = false;
#ifdef VIS_C
	HANDLE hToken;
    if ( !OpenThreadToken ( GetCurrentThread(), TOKEN_ADJUST_PRIVILEGES | TOKEN_QUERY, FALSE, &hToken ) ) {
        if ( GetLastError () == ERROR_NO_TOKEN ) {
            if ( !ImpersonateSelf ( SecurityImpersonation ) ) {
				ErrorHandler::genError ()->error ( "ImpersonateSelf failed.\n" );
			}
            if(!OpenThreadToken(GetCurrentThread(), TOKEN_ADJUST_PRIVILEGES | TOKEN_QUERY, FALSE, &hToken)){
				ErrorHandler::genError ()->error ( "OpenThreadToken failed.\n" );
            }
		}
        else {
			ErrorHandler::genError ()->error ( "OpenThreadToken failed.\n" );
		}
	}
    if ( !SetPrivilege(hToken, SE_DEBUG_NAME, TRUE ) ) {
		ErrorHandler::genError ()->error ( "SetPrivilege failed.\n" );
    }
	SetLastError ( 0 );
	HANDLE h = OpenProcess ( PROCESS_TERMINATE, FALSE, pid );
	if ( h == NULL ) {
		//DWORD errorCode = GetLastError ();
		//sprintf ( gen_error_message, "OpenProcess failed. System error code = %d id = %d.\n", errorCode, pid );
		//ErrorHandler::genError ()->error ( gen_error_message );
		return false;	// Probably means the process isn't running
	}
	else {
		flag = true;
	}
	if ( TerminateProcess ( h, -1 ) == 0 ) {
		CloseHandle ( h );
		ErrorHandler::genError ()->error ( "TerminateProcess failed.\n" );
	}
	CloseHandle ( h );
#else
	string command = "kill -TERM ";
	command += gen_itoa ( pid );
	genSystem ( command, "", true );
#endif
	return flag;
}
void genReduceProcessPriority ()
{
#ifdef VIS_C
	SetPriorityClass ( GetCurrentProcess (), 0x00004000 );
#else
	nice ( 19 );
#endif
}

IntVector getProcessNumberList ()		// Get the list of running pids.
{
	IntVector iv;
#ifdef VIS_C
	DWORD aProcesses[1024];
	DWORD cbNeeded;

    if ( !EnumProcesses( aProcesses, sizeof(aProcesses), &cbNeeded ) )
        return iv;

	DWORD cProcesses = cbNeeded / sizeof(DWORD);	// Calculate how many process identifiers were returned.

    for ( unsigned int i = 0; i < cProcesses ; i++ ) {
		iv.push_back ( aProcesses[i] );
	}
#else
	PPTempFile pidTempFile ( "pid", ".txt" );
	string filename = pidTempFile.getAdjustedPath ();
	string command;
	command += "ps -e > ";
	command += filename;
	genSystem ( command, "", false );
	GenIFStream fromFile ( filename );
	string line;
	getline ( fromFile, line );			// Read and discard the header
	while ( getline ( fromFile, line ) ) {
		istringstream istr ( line );
		int pid;
		string tty;
		string time;
		string process;
		istr >> pid;
		istr >> tty;
		istr >> time;
		istr >> process;
		iv.push_back ( pid );
	}
	fromFile.close ();
	genUnlink ( filename );
#endif
	return iv;
}
bool isProcessRunning ( const string& processName )
{
#ifdef VIS_C
	DWORD aProcesses[1024];
	DWORD cbNeeded;

    if ( !EnumProcesses( aProcesses, sizeof(aProcesses), &cbNeeded ) )
        return false;

	DWORD cProcesses = cbNeeded / sizeof(DWORD);	// Calculate how many process identifiers were returned.

    for ( unsigned int i = 0; i < cProcesses ; i++ ) {
		PairStringInt p = getProcessNameAndID ( aProcesses[i] );
		if ( p.first == processName ) return true;
	}
#else
	PPTempFile pidTempFile ( "pid", ".txt" );
	string filename = pidTempFile.getAdjustedPath ();
	string command;
	command += "ps -e > ";
	command += filename;
	genSystem ( command, "", false );
	GenIFStream fromFile ( filename );
	string line;
	getline ( fromFile, line );			// Read and discard the header
	while ( getline ( fromFile, line ) ) {
		istringstream istr ( line );
		int pid;
		string tty;
		string time;
		string process;
		istr >> pid;
		istr >> tty;
		istr >> time;
		istr >> process;
		if ( process == processName ) return true;
	}
	fromFile.close ();
	genUnlink ( filename );
#endif
	return false;
}
bool isProcessRunning ( int pid )
{
#ifdef VIS_C
	DWORD aProcesses [1024];
	DWORD cbNeeded;

    if ( !EnumProcesses ( aProcesses, sizeof(aProcesses), &cbNeeded ) ) {
  		ErrorHandler::genError ()->error ( "EnumProcesses failed.\n" );
		return false;
	}
	DWORD cProcesses = cbNeeded / sizeof(DWORD);	// Calculate how many process identifiers were returned.

    for ( unsigned int i = 0 ; i < cProcesses ; i++ ) {
		if ( aProcesses[i] == pid ) return true;
	}
	return false;
#else
	IntVector iv = getProcessNumberList ();
	IntVectorIterator iter = find ( iv.begin (), iv.end (), pid );
	return iter != iv.end ();
#endif
}
BoolDeque areProcessesRunning ( const IntVector& pid )
{
	BoolDeque bd;
#ifdef VIS_C
	DWORD aProcesses [1024];
	DWORD cbNeeded;

    if ( !EnumProcesses ( aProcesses, sizeof(aProcesses), &cbNeeded ) ) {
  		ErrorHandler::genError ()->error ( "EnumProcesses failed.\n" );
		return bd;
	}
	DWORD cProcesses = cbNeeded / sizeof(DWORD);	// Calculate how many process identifiers were returned.

	for ( IntVectorSizeType i = 0 ; i < pid.size () ; i++ ) {
		bool flag = false;
		for ( unsigned int j = 0 ; j < cProcesses ; j++ ) {
			if ( aProcesses [j] == pid [i] ) flag = true;
		}
		bd.push_back ( flag );
	}
#else
	IntVector iv = getProcessNumberList ();
	for ( IntVectorSizeType i = 0 ; i < pid.size () ; i++ ) {
		IntVectorIterator iter = find ( iv.begin (), iv.end (), pid [i] );
		bd.push_back ( iter != iv.end () );
	}
#endif
	return bd;
}
void genSleep ( int millisecs )
{
#ifdef VIS_C
	Sleep ( millisecs );
#else
	struct timespec req = { 0 };
	time_t sec = (int)(millisecs / 1000);
	millisecs = millisecs - (sec * 1000);
	req.tv_sec = sec;
	req.tv_nsec = millisecs * 1000000L;
	nanosleep ( &req, &req );
#endif
}
#ifndef VIS_C
static StringVector parseArgs ( const string& str )
{
	string::size_type startInd = 0;
	StringVector sv;
	for ( ; ; ) {
		string::size_type ind = str.find_first_of ( "\" ", startInd );
		if ( ind != string::npos ) {
			if ( str [ind] == '"' ) {
				startInd = ind + 1;
				string::size_type ind = str.find ( "\"", startInd );
				string s = str.substr ( startInd, ind-startInd );
				sv.push_back ( s );								// Empty argument is possible
				if ( ind == str.length () - 1 ) break;
				startInd = ind + 2;
			}
			else {
				string s = str.substr ( startInd, ind-startInd );
				if ( !s.empty () ) sv.push_back ( s );
				startInd = ind + 1;
			}
		}
		else {
			string s = str.substr ( startInd );
			if ( !s.empty () ) sv.push_back ( s );
			break;
		}
	}
	return sv;
}
#endif
int genCreateProcess ( const string& fullExePath, const string& params, const string& dir )
{
	int pid = 0;
	string runDir = dir;
	if ( runDir.empty () ) runDir = genCurrentWorkingDirectory ();
#ifdef VIS_C
	string fullParams = genFilenameFromPath ( fullExePath );
	if ( !params.empty () ) fullParams += " " + params;
	int waitSec = 0;
	DWORD dwExitCode;
	int elapsedMS = 0;
	STARTUPINFO siStartupInfo;
	PROCESS_INFORMATION piProcessInfo;
	memset ( &siStartupInfo, 0, sizeof(siStartupInfo) );
	memset ( &piProcessInfo, 0, sizeof(piProcessInfo) );
	siStartupInfo.cb = sizeof(siStartupInfo);
	if ( CreateProcess ( fullExePath.c_str (), const_cast<char*>(fullParams.c_str()), 0, 0,
//		false, CREATE_DEFAULT_ERROR_MODE, 0, runDir.c_str (), &siStartupInfo, &piProcessInfo ) != 0 ) {
		false, CREATE_NEW_CONSOLE, 0, runDir.c_str (), &siStartupInfo, &piProcessInfo ) != 0 ) {

		GetExitCodeProcess ( piProcessInfo.hProcess, &dwExitCode );

		while ( dwExitCode == STILL_ACTIVE && waitSec ) {
			GetExitCodeProcess ( piProcessInfo.hProcess, &dwExitCode );
			Sleep ( 500 );
			elapsedMS += 500;

			if ( elapsedMS > (waitSec * 1000) ) {
				dwExitCode = 0;
			}
		}
		pid = piProcessInfo.dwProcessId;
	}
	CloseHandle ( piProcessInfo.hProcess );
	CloseHandle ( piProcessInfo.hThread );
#else
	switch ( pid = vfork () ) {
		case -1:	// Fork fails
			return 0;			// Error to be caught by the application
		case 0:
		{
			genChangeWorkingDirectory ( runDir.c_str () );
			StringVector sv = parseArgs ( params );
			string::size_type siz = sv.size ();
			char** c = new char* [siz+2];
			string bin = genFilenameFromPath ( fullExePath );
			c [0] = const_cast<char*>(bin.c_str ());
			for ( int i = 1 ; i <= siz ; i++ ) c [i] = const_cast<char*>(sv [i-1].c_str ());
			c[siz+1] = 0;
			execv ( fullExePath.c_str (), c );
			_exit ( 0 );
		}
			break;
		default:
			break;
	}
#endif
	return pid;
}
#ifndef VIS_C
void genInitSigchld ( void ( *handler ) ( int sigNum ) )
{
	struct sigaction act;
	act.sa_handler = handler;
	act.sa_flags = SA_NOCLDSTOP | SA_RESTART;
	sigfillset ( &act.sa_mask );	// Blocks all signal whilst handler is called
	int ret = sigaction ( SIGCHLD, &act, NULL );
}
void genInitSigterm ( void ( *handler ) ( int sigNum ) )
{
	struct sigaction act;
	act.sa_handler = handler;
	sigemptyset (&act.sa_mask);
	act.sa_flags = 0;
	int ret = sigaction ( SIGTERM, &act, NULL );
}
void genInitSighup ( void ( *handler ) ( int sigNum ) )
{
	struct sigaction act;
	act.sa_handler = handler;
	sigemptyset (&act.sa_mask);
	act.sa_flags = 0;
	int ret = sigaction ( SIGHUP, &act, NULL );
}
#endif
