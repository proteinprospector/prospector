/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lgen_service.cpp                                              *
*                                                                             *
*  Created    : October 22nd 2007                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2007-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef VIS_C
#include <windows.h>
#include <string>
#include <stdexcept>

using std::string;
using std::runtime_error;

/*	This is how certain events are caught in Windows:

Declare handler:

SetConsoleCtrlHandler ( (PHANDLER_ROUTINE) CtrlHandler, true );

Handler function:

bool CtrlHandler( DWORD fdwCtrlType )	// If true returned then interrupt is handled.
{										// If false is returned the next interrupt handler (generally ExitProcess) is called.
	switch ( fdwCtrlType ) {
		case CTRL_C_EVENT:				// Someone presses Ctrl-C.
			return false;
		case CTRL_BREAK_EVENT:			// Someone presses Ctrl-Break.
			return false;
		case CTRL_CLOSE_EVENT:			// Console closed or End Task in Task Manager. This does NOT catch End Process in Task Manager.
			return false;
		case CTRL_SHUTDOWN_EVENT:		// System is shutting down.
			return false;
		case CTRL_LOGOFF_EVENT:			// This means someone (anyone) has logged off. Only services get this signal.
										// We just ignore this.
			return false;
		default:
			return false;
	}
}
*/
class ManualResetEvent {
	HANDLE event;
public:
	enum EventState {
		NonsignaledState,
		SignaledState
	};
	explicit ManualResetEvent ( EventState initialState = NonsignaledState );

	bool signaled () const;
	void signal ();
	void reset ();
	bool wait ( DWORD timeout = INFINITE );
};
ManualResetEvent::ManualResetEvent ( EventState initialState ) :
    event ( CreateEvent (	0,								// Pointer to security attributes structure, 0 for default.				
							true,							// This is true for a manual reset event.
							initialState == SignaledState,	// Initial state, signaled or nonsignaled
							0 ) )							// Name of event object. 0 if not named.
{
	if ( !event ) throw ( runtime_error ( "Can't create an event using CreateEvent." ) );
}
bool ManualResetEvent::signaled () const		// Determines whether the event is currently signaled.
{
    return WaitForSingleObject ( event, 0 ) == WAIT_OBJECT_0;
}
void ManualResetEvent::signal ()				// Sets the state of the event to signaled.
{
    if ( !SetEvent ( event ) ) throw ( runtime_error ( "Can't set an event using SetEvent." ) );
}
void ManualResetEvent::reset ()					// Sets the state of the event to nonsignaled.
{
    if ( !ResetEvent ( event ) ) throw ( runtime_error ( "Can't reset an event using ResetEvent." ) );
}
bool ManualResetEvent::wait ( DWORD timeout )	// Waits for the event to become signaled.
{
    return WaitForSingleObject ( event, timeout ) == WAIT_OBJECT_0;
}

namespace {

SERVICE_STATUS_HANDLE statusHandle;		//The service status handle does not have to be closed.
ManualResetEvent* killServiceEvent = 0;
bool nServiceRunning;
HANDLE hServiceThread;
DWORD nServiceCurrentStatus;

void KillService ();
void ( *daemon ) ();

bool UpdateServiceStatus ( DWORD currentState, DWORD win32ExitCode, DWORD serviceSpecificExitCode, DWORD checkPoint, DWORD waitHint )
{
	SERVICE_STATUS status;
	status.dwServiceType = SERVICE_WIN32_OWN_PROCESS;
	status.dwCurrentState = currentState;
	status.dwControlsAccepted = ( currentState == SERVICE_START_PENDING ) ? 0 : ( SERVICE_ACCEPT_STOP | SERVICE_ACCEPT_SHUTDOWN );
	status.dwWin32ExitCode = ( serviceSpecificExitCode == 0 ) ? win32ExitCode : ERROR_SERVICE_SPECIFIC_ERROR;
	status.dwServiceSpecificExitCode = serviceSpecificExitCode;
	status.dwCheckPoint = checkPoint;
	status.dwWaitHint = waitHint;

	if ( SetServiceStatus ( statusHandle, &status ) == 0 ) {
		KillService ();
		return false;
	}
	else
		return true;
}
void KillService ()
{
	nServiceRunning = false;
	killServiceEvent->signal ();
	UpdateServiceStatus ( SERVICE_STOPPED, NO_ERROR, 0, 0, 0 );
}
void ServiceCtrlHandler ( DWORD nControlCode )
{
	switch ( nControlCode ) {	
		case SERVICE_CONTROL_SHUTDOWN:
		case SERVICE_CONTROL_STOP:
			nServiceCurrentStatus = SERVICE_STOP_PENDING;
			UpdateServiceStatus ( SERVICE_STOP_PENDING, NO_ERROR, 0, 1, 3000 );
			KillService ();		
			return;
		default:
			break;
	}
	UpdateServiceStatus ( nServiceCurrentStatus, NO_ERROR, 0, 0, 0 );
}
DWORD ServiceExecutionThread ( LPDWORD param )
{
	while ( nServiceRunning ) {
		daemon ();
	}
	return 0;
}
bool StartServiceThread ()
{
	DWORD id;
	hServiceThread = CreateThread ( 0, 0, (LPTHREAD_START_ROUTINE)ServiceExecutionThread, 0, 0, &id );
	if ( hServiceThread == 0 ) {
		return false;
	}
	else {
		nServiceRunning = true;
		return true;
	}
}
void WINAPI ServiceMain ( DWORD argc, LPTSTR *argv )
{
	statusHandle = RegisterServiceCtrlHandler ( "", (LPHANDLER_FUNCTION)ServiceCtrlHandler );
	if ( statusHandle == 0 ) {
		throw runtime_error ( "Can't start the service: error in RegisterServiceCtrlHandler" );
	}
	if ( UpdateServiceStatus ( SERVICE_START_PENDING, NO_ERROR, 0, 1, 3000 ) == 0 ) {
		throw runtime_error ( "Can't start the service: error in UpdateServiceStatus - SERVICE_START_PENDING" );
	}
	killServiceEvent = new ManualResetEvent ();
	if ( UpdateServiceStatus ( SERVICE_START_PENDING, NO_ERROR, 0, 2, 1000 ) == 0 ) {
		throw runtime_error ( "Can't start the service: error in UpdateServiceStatus - SERVICE_START_PENDING" );
	}
	if ( StartServiceThread () == false ) {
		throw runtime_error ( "Can't start the service thread." );
	}
	nServiceCurrentStatus = SERVICE_RUNNING;
	if ( UpdateServiceStatus ( SERVICE_RUNNING, NO_ERROR, 0, 0, 0 ) == 0 ) {
		throw runtime_error ( "Can't start the service: error in UpdateServiceStatus - SERVICE_RUNNING" );
	}
	killServiceEvent->wait ();
	delete killServiceEvent;
}
void startService ()
{
	SERVICE_TABLE_ENTRY serviceTable [] = {	{ "", ServiceMain }, { 0, 0 } };
	if ( StartServiceCtrlDispatcher ( serviceTable ) == 0 ) {
		throw runtime_error ( "Can't start the service: error in StartServiceCtrlDispatcher" );
	}
}
}
void genStartDaemon ( void ( *passedDaemon ) () )
{
	daemon = passedDaemon;
	startService ();
}
bool genIsServiceRunning ( const string& serviceName )
{
	SC_HANDLE scm = OpenSCManager ( 0, 0, SC_MANAGER_CONNECT );
	if ( scm == 0 ) {
		throw runtime_error ( "Can't open the Service Control Manager" );
	}
	SC_HANDLE sch = OpenService ( scm, serviceName.c_str (), SERVICE_QUERY_STATUS );
	if ( sch == 0 ) { 
		CloseServiceHandle ( scm );
		throw runtime_error ( "Can't open the service: " + serviceName );
	}
	SERVICE_STATUS serviceStatus;
	if ( !QueryServiceStatus ( sch, &serviceStatus ) ) {	// Error
		CloseServiceHandle ( sch );
		CloseServiceHandle ( scm );
		throw runtime_error ( "Can't query the status for service: " + serviceName );
	}
	return serviceStatus.dwCurrentState == SERVICE_RUNNING ? true : false;
}
bool genStopService ( const string& serviceName )
{
	SC_HANDLE scm = OpenSCManager ( 0, 0, SC_MANAGER_ALL_ACCESS );
	if ( scm == 0 ) {
		throw runtime_error ( "Can't open the Service Control Manager" );
	}
	SC_HANDLE sch = OpenService ( scm, serviceName.c_str (), SERVICE_STOP | SERVICE_QUERY_STATUS );
	if ( sch == 0 ) { 
		CloseServiceHandle ( scm );
		throw runtime_error ( "Can't open the service: " + serviceName );
    }
	SERVICE_STATUS serviceStatus;
	if ( !QueryServiceStatus ( sch, &serviceStatus ) ) {
		CloseServiceHandle ( sch );
		CloseServiceHandle ( scm );
		throw runtime_error ( "Can't query the status for service: " + serviceName );
	}
	if ( serviceStatus.dwCurrentState == SERVICE_STOPPED ) {		// Service has already stopped
		CloseServiceHandle ( sch );
		CloseServiceHandle ( scm );
		return true;
	}
	DWORD dwWaitTime = serviceStatus.dwWaitHint / 10;	// Do not wait longer than the wait hint. A good interval
	if( dwWaitTime < 1000 )			dwWaitTime = 1000;	// is one-tenth the wait hint, but no less than 1 second
	else if ( dwWaitTime > 10000 )  dwWaitTime = 10000;	// and no more than 10 seconds.
	while ( serviceStatus.dwCurrentState == SERVICE_STOP_PENDING ) {
		Sleep ( dwWaitTime );
		if ( !QueryServiceStatus ( sch, &serviceStatus ) ) {			// Check the status again. 
			CloseServiceHandle ( sch ); 
			CloseServiceHandle ( scm );
			throw runtime_error ( "Can't query the status for service: " + serviceName );
		}
		if ( serviceStatus.dwCurrentState == SERVICE_STOPPED ) {		// Service has now stopped
			CloseServiceHandle ( sch );
			CloseServiceHandle ( scm );
			return true;
		}
	}
	if ( !ControlService ( sch, SERVICE_CONTROL_STOP, &serviceStatus ) ) {
		CloseServiceHandle ( sch );
		CloseServiceHandle ( scm );
		throw runtime_error ( "Can't stop the service: " + serviceName );
    }
	while ( serviceStatus.dwCurrentState != SERVICE_STOPPED ) {		// Do not wait longer than the wait hint. A good interval
																	// is one-tenth the wait hint, but no less than 1 second
																	// and no more than 10 seconds. 

		Sleep ( dwWaitTime );
		if ( !QueryServiceStatus ( sch, &serviceStatus ) ) {	// Check the status again. 
			CloseServiceHandle ( sch ); 
			CloseServiceHandle ( scm );
			throw runtime_error ( "Can't query the status for service: " + serviceName );
		}
	}
	CloseServiceHandle ( sch ); 
	CloseServiceHandle ( scm );
	return true;
}
bool genStartService ( const string& serviceName, int argc, const char** argv )
{
	SC_HANDLE scm = OpenSCManager ( 0, 0, SC_MANAGER_CONNECT );
	if ( scm == 0 ) {
		throw runtime_error ( "Can't open the Service Control Manager" );
	}
	SC_HANDLE sch = OpenService ( scm, serviceName.c_str (), SERVICE_START | SERVICE_QUERY_STATUS );
	if ( sch == 0 ) { 
		CloseServiceHandle ( scm );
		throw runtime_error ( "Can't open the service: " + serviceName );
    }
	if ( !StartService ( sch, argc, argv ) ) {
		CloseServiceHandle ( sch );
		CloseServiceHandle ( scm );
		throw runtime_error ( "Can't start the service: " + serviceName );
	}
	SERVICE_STATUS serviceStatus;
	if ( !QueryServiceStatus ( sch, &serviceStatus ) ) {
		CloseServiceHandle ( sch );
		CloseServiceHandle ( scm );
		throw runtime_error ( "Can't query the status for service: " + serviceName );
	}
	DWORD dwWaitTime = serviceStatus.dwWaitHint / 10;	// Do not wait longer than the wait hint. A good interval
	if( dwWaitTime < 1000 )			dwWaitTime = 1000;	// is one-tenth the wait hint, but no less than 1 second
	else if ( dwWaitTime > 10000 )  dwWaitTime = 10000;	// and no more than 10 seconds.
	while ( serviceStatus.dwCurrentState == SERVICE_START_PENDING ) {
		Sleep ( dwWaitTime );

		if ( !QueryServiceStatus ( sch, &serviceStatus ) ) {	// Check the status again. 
			CloseServiceHandle ( sch ); 
			CloseServiceHandle ( scm );
			throw runtime_error ( "Can't query the status for service: " + serviceName );
		}
    } 
    CloseServiceHandle ( sch ); 
    CloseServiceHandle ( scm );

    return serviceStatus.dwCurrentState == SERVICE_RUNNING ? true : false;
}
void genInstallService ( const string& serviceName, const string& serviceDisplayName, const string& serviceBinaryPath, const string& user, const string& password )
{
	SC_HANDLE scm = OpenSCManager (0, 0, SC_MANAGER_CREATE_SERVICE );
	if ( scm == 0 ) {
		throw runtime_error ( "Can't open the Service Control Manager" );
	}
	SC_HANDLE sch = CreateService ( scm,	// SCM handle
		serviceName.c_str (),				// Name of service 
		serviceDisplayName.c_str (),		// Service name to display 
		SERVICE_ALL_ACCESS,					// Desired access 
		SERVICE_WIN32_OWN_PROCESS,			// Service type
		SERVICE_DEMAND_START,				// Start type 
		SERVICE_ERROR_NORMAL,				// Error control type
		serviceBinaryPath.c_str (),			// Path to service's binary 
		0,									// No load ordering group
		0,									// No tag identifier 
		0,									// No dependencies 
		user.c_str (),						// LocalSystem account
		password.c_str () );				// Password
	if ( sch == 0 ) {
		CloseServiceHandle ( scm );
		throw runtime_error ( "Can't create the service: " + serviceName );
    }
	CloseServiceHandle ( sch ); 
	CloseServiceHandle ( scm );
}
void genUninstallService ( const string& serviceName )
{
	if ( genIsServiceRunning ( serviceName ) ) {
		throw runtime_error ( "Can't delete the service as it is running. Please stop it first." );
	}
	SC_HANDLE scm = OpenSCManager ( 0, 0, SC_MANAGER_ALL_ACCESS );
	if ( scm == 0 ) {
		throw runtime_error ( "Can't open the Service Control Manager" );
	}
	SC_HANDLE sch = OpenService (	scm, serviceName.c_str (), DELETE );
	if ( sch == 0 ) { 
		CloseServiceHandle ( scm );
		throw runtime_error ( "Can't open the service: " + serviceName );
    }
	int flag = DeleteService ( sch );
	CloseServiceHandle ( sch ); 
	CloseServiceHandle ( scm );
	if ( !flag ) {
		throw runtime_error ( "Can't delete the service: " + serviceName );
    }
}
#endif
