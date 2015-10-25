/******************************************************************************
*                                                                             *
*  Program    : btag_daemon                                                   *
*                                                                             *
*  Filename   : btag_daemon_main.cpp                                          *
*                                                                             *
*  Created    : September 10th 2007                                           *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2007-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef VIS_C
#include <windows.h>
#include <lgen_service.h>
#else
#include <stdexcept>
#include <sys/types.h>
#include <sys/wait.h>
#endif
#include <algorithm>
#include <lg_io.h>
#include <lg_time.h>
#include <lgen_process.h>
#include <lgen_file.h>
#include <lgen_net.h>
#include <lg_string.h>
#include <ld_init.h>
#include <lu_btag_submit.h>
#include <lu_param_list.h>
#include <lu_getfil.h>
#include <lu_fasta.h>
using std::ostream;
using std::string;
using std::endl;
using std::cout;
using std::set_difference;
using std::inserter;
using std::runtime_error;

namespace {			// Local functions
void stopDaemon ();
void readParameters ();
void logOutput ( const string& message )
{
	static string logFileName = ppBaseDir () + "logs" + SLASH + "btagd.log";
	static GenOFStream ofs ( logFileName, std::ios_base::out | std::ios_base::app );
	ofs << genCurrentTimeAndDateString () << ": " << message << endl;
}
#ifdef VIS_C
void serviceMaintenance ( const string& binaryName, const string& command, const string& user = "", const string& password = "" );
void serviceMaintenance ( const string& binaryName, const string& command, const string& user, const string& password )
{
	try {
		string daemonName = binaryName;
		if ( isSuffix ( daemonName, ".exe" ) ) daemonName = daemonName.substr ( 0, daemonName.length () - 4 );
		if ( command == "install" ) {
			string daemonDisplayName = "Batch-Tag Daemon";
			if ( daemonName != "btag_daemon" ) {
				daemonDisplayName += " " + daemonName;
			}
			string cwd = genCurrentWorkingDirectory ();
			string serviceCommand;
			serviceCommand += "\"" + cwd + SLASH + binaryName + "\"";
			serviceCommand += " ";
			serviceCommand += "run";
			serviceCommand += " ";
			serviceCommand += "\"" + cwd + "\"";
			genInstallService ( daemonName, daemonDisplayName, serviceCommand, user, password );
		}
		else if ( command == "stop" )
			genStopService ( daemonName );
		else
			genUninstallService ( daemonName );
	}
	catch ( runtime_error e ) {
		cout << e.what () << endl;
		exit ( 0 );
	}
}
#else
bool sigtermReceived = false;
void sigchldHandler ( int sigNum )	// This is a signal handler to prevent the child processes becoming zombies
{
	pid_t pid;
	int status;
	while ( ( pid = waitpid ( -1, &status, WNOHANG ) ) > 0 );
}
void sigtermHandler ( int sigNum )
{
	sigtermReceived = true;
}
bool hupSignal = false;
void sighupHandler ( int sigNum )
{
	hupSignal = true;
}
void initSignalHandlers ()
{
	genInitSigchld ( sigchldHandler );
	genInitSigterm ( sigtermHandler );
	genInitSighup ( sighupHandler );
}
#endif
bool sendEmailFlag = false;
void sendEmail ( const DaemonJobItem* dji, const string& url )
{
	if ( sendEmailFlag ) {
		string params;
		params += dji->getEmail ();
		params += " ";
		params += dji->getSearchJobID ();
		params += " ";
		params += "\"" + dji->getProjectName () + "\"";
		params += " ";
		params += "\"" + dji->getResultsName () + "\"";
		params += " ";
		params += "\"" + dji->getJobStatus () + "\"";
		params += " ";
		params += "\"" + url + "\"";
		int pid = genCreateProcess ( getSystemCall ( "mail.pl" ), params );
		if ( pid == 0 ) {
			logOutput ( "Problem running the script mail.pl." );
		}
	}
}
int abortedJobsDeleteDays = 0;
int sessionDeleteDays = 0;
void dailyJobs ()
{
	if ( abortedJobsDeleteDays )	MySQLPPSDDBase::instance ().deleteOldAbortedJobs ( abortedJobsDeleteDays );
	if ( sessionDeleteDays )		MySQLPPSDDBase::instance ().deleteOldSessions ( sessionDeleteDays );
}
void cleanup ( const DaemonJobItem* dji )
{
	static string repHome = InfoParams::instance ().getUserRepository ();
	string resPath = dji->getResultsPath ();
	int end = resPath.rfind ( "results" );
	string dateString = resPath.substr ( end+8 );
	string userDir = resPath.substr ( 0, end );
	string resFile = dji->getResultsFile ();
	string searchKey = dji->getSearchJobKey ();
	string projectPath = dji->getProjectPath ();
	string projectFile = dji->getProjectFile ();
	string projectName = dji->getProjectName ();
	string projectID = dji->getProjectID ();
	logOutput ( "Search " + searchKey + " aborted, status: " + dji->getJobStatus () + ", signal: " + dji->getJobSignal () );
	FileList f1 ( repHome + SLASH + resPath, resFile, "", false );	// Delete the results file(s).
	f1.unlink ();
	genRmdir ( repHome + SLASH + userDir + "results" + SLASH + dateString );	// Try to delete the results date directory. This will fail if the directory isn't empty.

	string dataPath = repHome + SLASH + userDir + "data" + SLASH + dateString + SLASH + searchKey;

	if ( genFileExists ( dataPath ) ) {	// This is a Batch-Tag Web search as there is a data directory with the
											// same name as the search key.
		StringVector res = MySQLPPSDDBase::instance ().getResultsList ( projectID );
		if ( res.size () != 1 ) {	// This is to catch an error that occurs with the LINUX mounter
			throw runtime_error ( "Daemon requested deletion of Batch-Tag Web project with multiple results." );
		}
		genUnlinkDirectory ( dataPath );	// Delete the data directory.

		genUnlink ( repHome + SLASH + projectPath + SLASH + projectFile );// Delete the project file.
		FileList f2 ( repHome + SLASH + projectPath, projectName + ".exp.", "", false );// Delete the .exp files.
		f2.unlink ();
		genRmdir ( repHome + SLASH + userDir + "project" + SLASH + dateString );// Try to delete the project date directory. This will fail if the directory isn't empty.
		MySQLPPSDDBase::instance ().deleteProjectID ( projectID );// Delete the project database entry. This will also delete the search entry.
		genRmdir ( repHome + SLASH + userDir + "data" + SLASH + dateString );	// Try to delete the data date directory. This will fail if the directory isn't empty.
	}
	else {	// This is a Batch-Tag search
		string fullProjectPath = repHome + SLASH + projectPath;
		FileList f3 ( fullProjectPath, projectName + ".exp.", "", false );// Delete the .exp files if appropriate.
		// Check if it is the correct file and is unfinished
		StringVector fList = f3.getNameList ();
		string prefix;
		for ( StringVectorSizeType i = 0 ; i < fList.size () ; i++ ) {
			string f = fList [i];
			if ( isSuffix ( f, "_0" ) ) {	// Only interested in unfinished exp searches
				ParameterList fParams ( fullProjectPath + SLASH + f, false, false, false, false );
				if ( fParams.getStringValue ( "search_key" ) == searchKey ) {
					prefix = f.substr ( 0, f.length () - 2 );
					for ( StringVectorSizeType j = 0 ; j < fList.size () ; j++ ) {
						string f2 = fList [j];
						if ( isPrefix ( f2, prefix ) ) genUnlink ( fullProjectPath + SLASH + f2 );
					}
					break;
				}
			}
		}
		MySQLPPSDDBase::instance ().deleteSearchID ( dji->getSearchJobID () );// Delete the search database entry.
	}
}
int loopTime = 5000;
const int dayInSec = 86400;
int maxSearches;
int maxJobsPerUser;
bool singleServer = false;
GenElapsedTime et;
SetString startedJobs;
SetString activeJobs;
string serverStr;
bool firstLoop = true;

void daemonLoop ()
{
	try {
		static string host = Hostname::instance ().getHostname ();
		readParameters ();
		DaemonJobQueue djq = MySQLPPSDDBase::instance ().getDaemonJobQueue ();
		djq.setActions ( host, maxJobsPerUser );
		VectorDaemonJobItem cleanUpItems = djq.getCleanUpJobItems ();
		for ( VectorDaemonJobItemSizeType i = 0 ; i < cleanUpItems.size () ; i++ ) {	// These searches have been aborted for some reason
			DaemonJobItem* dji = &cleanUpItems [i];
			cleanup ( dji );
			sendEmail ( dji, serverStr + "/jobStatus.cgi?search_key=" + dji->getSearchJobKey () );
			startedJobs.erase ( dji->getSearchJobKey () );
		}
		SetString doneJobs;
		activeJobs = djq.getActiveJobKeys ();
		set_difference ( startedJobs.begin (), startedJobs.end (), activeJobs.begin (), activeJobs.end (), inserter ( doneJobs, doneJobs.begin () ) );
		for ( SetStringConstIterator sKey = doneJobs.begin () ; sKey != doneJobs.end () ; sKey++ ) {
			DaemonJobItem* dji = MySQLPPSDDBase::instance ().getDaemonJobItemByKey ( *sKey );
			sendEmail ( dji, serverStr + "/msform.cgi?form=search_compare&search_key=" + *sKey );
			delete dji;
			startedJobs.erase ( *sKey );
		}
		if ( djq.getNumRunning () < maxSearches ) {		// Only a single search can be started each time round the loop
			string searchKey = djq.getNextSearchJobKey ();
			if ( !searchKey.empty () ) {
				if ( !singleServer ) genSleep ( djq.getNumRunning () * loopTime );
				if ( MySQLPPSDDBase::instance ().setJobStart ( searchKey ) ) {
					if ( !btagSubmit ( searchKey, Hostname::instance ().getHostname () ) ) {
						logOutput ( "Unable to start search." );
						MySQLPPSDDBase::instance ().setJobSubmitted ( searchKey );	// Resubmit the job if it wouldn't start
					}
					else {
						logOutput ( "Started search " + searchKey );
						startedJobs.insert ( searchKey );
					}
				}
			}
		}
		if ( firstLoop || et.getElapsedTime () > dayInSec ) {
			GenElapsedTime newTime;
			et = newTime;						// Reset the time counter
			dailyJobs ();
			firstLoop = false;
		}
#ifdef VIS_C
		genSleep ( loopTime );
#else
		int millisecs = loopTime;
		struct timespec req = { 0 };
		time_t sec = (int)(millisecs / 1000);
		millisecs = millisecs - (sec * 1000);
		req.tv_sec = sec;
		req.tv_nsec = millisecs * 1000000L;
		while ( nanosleep ( &req, &req ) == -1 ) {
			if ( sigtermReceived ) break;
		}
#endif
	}
	catch ( runtime_error e ) {
		logOutput ( e.what () );
		logOutput ( "Daemon finished unexpectedly." );
		exit ( 0 );
	}
}
void stopDaemon ()
{
	for ( SetStringConstIterator sKey = activeJobs.begin () ; sKey != activeJobs.end () ; sKey++ ) {	// Daemon stopped kill active jobs
		DaemonJobItem* dji = MySQLPPSDDBase::instance ().getDaemonJobItemByKey ( *sKey );
		if ( dji ) {
			killProcess ( dji->getPID () );	// Kill the process. Doesn't matter if it isn't running
			MySQLPPSDDBase::instance ().processAbortSignal ( dji->getSearchJobID () );
			MySQLPPSDDBase::instance ().submitAbortedJob ( *sKey, "This search was aborted because the daemon was stopped." );
			cleanup ( dji );				// Clean up the files
			sendEmail ( dji, serverStr + "/jobStatus.cgi?search_key=" + *sKey );
		}
		delete dji;
	}
	logOutput ( "Daemon finished." );
}
void readParameters ()
{
	static bool first = true;
#ifdef VIS_C
	static time_t lastModifiedTime = InfoParams::instance ().getLastModifiedTime ();
	time_t currentModifiedTime = InfoParams::instance ().getLastModifiedTime ();
	if ( first || currentModifiedTime > lastModifiedTime ) {
		lastModifiedTime = currentModifiedTime;
#else
	if ( first || hupSignal ) {
#endif
		logOutput ( "Daemon reading parameter file." );
		sendEmailFlag = InfoParams::instance ( true ).getBoolValue ( "email" );
		string serverName = InfoParams::instance ().getStringValue ( "server_name", "localhost" );
		string serverPort = InfoParams::instance ().getStringValue ( "server_port", "80" );
		serverStr = "http://" + serverName;
		if ( serverPort != "80" ) serverStr += string ( ":" ) + serverPort;
		string virtualDir = InfoParams::instance ().getStringValue ( "virtual_dir", "" );
		if ( !virtualDir.empty () ) serverStr += "/" + virtualDir;
		serverStr += "/" + BinDir::instance ().getBinDir ();

		maxSearches = InfoParams::instance ().getIntValue ( "max_btag_searches", 1 );
		maxJobsPerUser = InfoParams::instance ().getIntValue ( "max_jobs_per_user", 0 );
		singleServer = InfoParams::instance ().getBoolValue ( "single_server", false );
		loopTime = InfoParams::instance ().getIntValue ( "daemon_loop_time", 5 ) * 1000;
		sessionDeleteDays = InfoParams::instance ().getIntValue ( "session_delete_days", 0 );
		abortedJobsDeleteDays = InfoParams::instance ().getIntValue ( "aborted_jobs_delete_days", 0 );
		first = false;
#ifndef VIS_C
		if ( hupSignal ) genSleep ( 1000 ); // If responding to interrupt sleep 1 sec before restarting
		hupSignal = false;
#endif
	}
}
void runService ( const string& currentWorkingDir )
{
	genChangeWorkingDirectory ( currentWorkingDir );
	readParameters ();
#ifndef VIS_C
	initSignalHandlers ();
#endif
	try {			// Check the database connection before starting the daemon
		MySQLPPSDDBase::instance ();
	}
	catch ( runtime_error ) {
		logOutput ( "Daemon couldn't connect to the database." );
		return;
	}
	logOutput ( "Daemon starting." );
	PreloadedDatabases pd;			// Load the preloaded databases into memory
#ifdef VIS_C
	try {
		genStartDaemon ( daemonLoop );
	}
	catch ( runtime_error e ) {
		logOutput ( e.what () );
		return;
	}
#else
	for ( ; ; ) {
		if ( sigtermReceived ) break;
		daemonLoop ();
	}
#endif
	stopDaemon ();
	MySQLPPSDDBase::instance ( false, true );
}

}	// End of namespace
int main ( int argc, char** argv )
{
	bool validArgs = true;
	if ( argc < 2 ) {
		validArgs = false;
	}
	else {
		string type = argv [1];
		if ( type == "run" ) {
			if ( argc == 3 )	runService ( argv [2] );
			else				validArgs = false;
		}
#ifdef VIS_C
		else if ( type == "install" ) {
			if ( argc == 4 )	serviceMaintenance ( argv [0], argv [1], argv [2], argv [3] );
			else				validArgs = false;
		}
		else if ( type == "uninstall" ) {
			if ( argc == 2 )	serviceMaintenance ( argv [0], argv [1] );
			else				validArgs = false;
		}
		else if ( type == "stop" ) {
			if ( argc == 2 )	serviceMaintenance ( argv [0], argv [1] );
			else				validArgs = false;
		}
#endif
		else
			validArgs = false;
	}
	if ( !validArgs ) cout << "Invalid arguments." << endl;
	return 1;
}
