/******************************************************************************
*                                                                             *
*  Program    : msbatch                                                       *
*                                                                             *
*  Filename   : msbatch_main.cpp                                              *
*                                                                             *
*  Created    : February 28th 2007                                            *
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
#ifndef VIS_C
#include <stdexcept>
#endif
#include <lg_string.h>
#include <lgen_file.h>
#include <lgen_process.h>
#include <lgen_service.h>
#include <lgen_uncompress.h>
#include <ld_init.h>
#include <lu_file_type.h>
#include <lu_getfil.h>
#include <lu_param_list.h>
#include <lu_prog_par.h>
#include <lu_html.h>
#include <lu_repository.h>
#include <lu_btag_submit.h>
#include <lu_df_info.h>
#include <lu_mass.h>
using std::runtime_error;
using std::ostream;
using std::cout;
using std::string;
using std::endl;

void startDaemonIfNecessary ()
{
#ifdef VIS_C
	string serviceName = InfoParams::instance ().getStringValue ( "btag_daemon_name", "btag_daemon" );
	if ( !genIsServiceRunning ( serviceName ) ) {
		genStartService ( serviceName );
	}
#else
	string daemonName = InfoParams::instance ().getStringValue ( "btag_daemon_name", "btag-daemon" );
	if ( !isProcessRunning ( daemonName ) ) {		// Start the daemon if it isn't running
		string params = "run";
		params += " ";
		params += genCurrentWorkingDirectory ();
		int pid = genCreateProcess ( getSystemCall ( daemonName ), params );
		if ( pid == 0 ) {
			throw runtime_error ( "Problem starting the Batch-Tag Daemon." );
		}
	}
#endif
}
string getProjectName ( const string& path )
{
	string s = genFilenameFromPath ( path );
	string::size_type len = s.length ();
	string::size_type index;

	if ( isNoCaseSuffix ( s, ".tar.gz" ) )		index = len - strlen ( ".tar.gz" );
	else if ( isNoCaseSuffix ( s, ".tar.z" ) )	index = len - strlen ( ".tar.z" );
	else if ( isNoCaseSuffix ( s, ".tar.bz2" ) )index = len - strlen ( ".tar.bz2" );
	else {
		if ( isNoCaseSuffix ( s, ".gz" ) )		s = s.substr ( 0, len - strlen ( ".gz" ) );
		else if ( isNoCaseSuffix ( s, ".z" ) )	s = s.substr ( 0, len - strlen ( ".z" ) );
		else if ( isNoCaseSuffix ( s, ".bz2" ) )s = s.substr ( 0, len - strlen ( ".bz2" ) );
		index = s.find_last_of ( "." );
	}
	if ( index == string::npos ) return s;
	else return s.substr ( 0, index );
}
void saveParams ( ostream& os, ParameterList& paramList )
{
	init_html ( os, "Saving Batch-Tag Parameters" );
	paramList.removeName ( "save_params" );
	ParameterList cookieParamList ( "", false, false );		// Create empty parameter set
	cookieParamList.appendParameters ( paramList );
	//cookieParamList.removeName ( "accession_nums" );
	cookieParamList.removeName ( "data_source" );
	cookieParamList.removeName ( "max_hits" );
	cookieParamList.removeName ( "output_filename" );
	cookieParamList.removeName ( "password" );
	cookieParamList.removeName ( "project_name" );
	cookieParamList.removeName ( "report_title" );
	cookieParamList.removeName ( "script_filename" );
	cookieParamList.removeName ( "search_name" );
	cookieParamList.removeName ( "upload_data_filename" );
	cookieParamList.removeName ( "upload_data_filepath" );
	cookieParamList.removeName ( "user" );
	cookieParamList.removeName ( "version" );
	if ( paramList.getStringValue ( "data_source" ) == "List of Files" )	// Batch-Tag rather than Batch-Tag Web
		cookieParamList.removeName ( "instrument_name" );
	else
		genUnlink ( paramList.getStringValue ( "upload_data_filepath" ) );	// Delete any uploaded file 
	bool ret = cookieParamList.copyToCookie ( os, "batchtag_params" );
	if ( !ret ) {
		ErrorHandler::genError ()->error ( "Could not save the parameters as their length exceeds the maximum cookie length.\n" );
	}
	os << "<p>Settings saved</p>" << endl;
	os << "<input type=\"button\" value=\"Search Form\" onclick=\"history.go(-1)\">" << endl;
}
void batchSubmission ( ostream& os, ParameterList& paramList )
{
	init_html ( os, "Batch Submission" );
	if ( paramList.empty () ) {
		ErrorHandler::genError ()->error ( "No parameters passed to Prospector program.\n" );
	}
	string user = paramList.getStringValue ( "user" );
	os << "<p>Creating Prospector files.</p>" << endl;
	string searchKey = genRandomString ( 16 );			// Selects a unique string for the search
	string searchName = paramList.getStringValue ( "search_name" );					// Program (eg batchtag)
	string uploadFpath = paramList.getStringValue ( "upload_data_filepath" );	// This is the full path name of the uploaded file as named by PP
	string resultsName = paramList.getStringValue ( "output_filename" );
	bool searchKeyProblem;
	try {
		checkConstantAndVariableModCompatibility ( &paramList );
	}
	catch ( runtime_error e ) {		// Catch database login problems
		genUnlink ( uploadFpath );
		ErrorHandler::genError ()->error ( e );
	}
	try {
		searchKeyProblem = MySQLPPSDDBase::instance ().checkSearchKey ( searchKey );// Checks it is unique
	}
	catch ( runtime_error e ) {		// Catch database login problems
		genUnlink ( uploadFpath );
		ErrorHandler::genError ()->error ( e );
	}
	if ( searchKeyProblem ) {
		genUnlink ( uploadFpath );
		ErrorHandler::genError ()->error ( "Search key not unique.\n" );// The search key has to be unique to carry on.
	}
	if ( resultsName.empty () ) {
		genUnlink ( uploadFpath );
		ErrorHandler::genError ()->error ( "No results name has been chosen.\n" );
	}
	string password;
	bool passwordFlag = paramList.getValue ( "password", password );
	UserInfo* userInfo = MySQLPPSDDBase::instance ().getUserInfo ( user );
	if ( !userInfo || ( passwordFlag && userInfo->getPassword () != password ) ) {													// Get the user information
		genUnlink ( uploadFpath );
		ErrorHandler::genError ()->error ( "Unknown user and/or password.\n" );
	}
	if ( user == "root" ) {													// Get the user information
		genUnlink ( uploadFpath );
		ErrorHandler::genError ()->error ( "Root user can't do searches.\n" );
	}
	if ( userInfo->getIsGuest () ) {
		genUnlink ( uploadFpath );
		ErrorHandler::genError ()->error ( "Guest users can't do searches.\n" );
	}
	if ( !InfoParams::instance ().getBoolValue ( "btag_daemon_remote" ) ) {
		try {
			startDaemonIfNecessary ();
		}
		catch ( runtime_error e ) {
			genUnlink ( uploadFpath );
			ErrorHandler::genError ()->error ( e );
		}
	}
	string userID = userInfo->getUserID ();
	Repository* reposit = new UserRepository ( searchName, userInfo->getDirectoryName () );	// Where to put files given a user. Also creates the directories.
	string projectName;
	if ( uploadFpath.empty () ) {		// The project already exists
		projectName = paramList.getStringValue ( "project_name" );	// From a form with a project name option.
	}
	else {
		string uploadFname = paramList.getStringValue ( "upload_data_filename" );	// This is the original file name.
		string instrument = paramList.getStringValue ( "instrument_name" );
		uploadFname = genFilenameFromPath ( uploadFname );						// IE gives the full path whereas Mozilla give the filename (what we want)
		projectName = getProjectName ( uploadFname );	// The project is named after the uploaded file. This has to be unique.
		if ( projectName.length () > 58 ) {
			genUnlink ( uploadFpath );
			ErrorHandler::genError ()->error ( "The project name " + projectName + " is too long.\n" );
		}
		if ( MySQLPPSDDBase::instance ().checkProject ( userID, projectName ) ) {
			genUnlink ( uploadFpath );
			ErrorHandler::genError ()->error ( "Project already exists.\n" );
		}
		paramList.setValue ( "data_source", "List of Files" );
		paramList.setValue ( "project_name", projectName );
		paramList.removeName ( "upload_data_filename" );
		paramList.removeName ( "upload_data_filepath" );
		string uploadName = genPreprocessFile ( uploadFpath );
		if ( uploadName == uploadFpath || isCompressedUpload ( uploadName ) ) {		// If the file hasn't been preprocessed or it has just been uncompressed it must be a single file
			string shortName = genShortFilenameFromPath ( uploadName );
			string newDir = genDirectoryFromPath ( uploadName ) + SLASH + shortName;
			if ( newDir == uploadName ) {
				newDir += "_1";
			}
			genCreateDirectory ( newDir );
			string newUploadName;
			if ( isCompressedUpload ( uploadName ) )
				newUploadName = newDir + SLASH + genShortFilenameFromPath ( uploadFname );
			else
				newUploadName = newDir + SLASH + uploadFname;
			genRename ( uploadName, newUploadName );
			uploadName = newDir;
		}
		PPProject ppp ( userInfo->getMaxMSMSSpectra (), userInfo->getAllowRaw (), projectName, uploadName, searchKey, reposit, MSMSDataSetInfo::setChargeRange ( &paramList ) );
		if ( ppp.initialised () ) {
			int err = genRename ( uploadName, reposit->getFullDataPath () + SLASH + searchKey );
			if ( err ) {
				genUnlink ( uploadName );
				string e = gen_itoa ( err );
				ErrorHandler::genError ()->error ( "Failed to move uploaded file to repository.\nError code = " + e + "\n" );
			}
			ppp.createProjectFile ( reposit->getFullProjectPath () + SLASH + projectName + ".xml" );
			string projectPath = reposit->getProjectPath ();
			MySQLPPSDDBase::instance ().submitProject ( userID, projectName, projectName + ".xml", reposit->getProjectPath (), instrument );
		}
		else {
			string err = ppp.getErrMessage ();
			if ( ppp.getDeleteFlag () ) {		// Known error - delete the upload directory
				genUnlinkDirectory ( uploadName );
			}
			if ( err.empty () )
				ErrorHandler::genError ()->error ( "Upload file has an invalid format.\n" );
			else
				ErrorHandler::genError ()->error ( err + "\n" );
		}
	}
	paramList.removeName ( "user" );
	paramList.removeName ( "password" );
	paramList.removeName ( "project_name" );
	paramList.removeName ( "output_filename" );					// Search results file.
	paramList.removeName ( "msms_mod_AA_list" );				// These options are used to load the variable mods
	paramList.removeName ( "msms_mod_AA_types" );				// and shouldn't be saved.
	paramList.removeName ( "mod_AA_limit" );
	paramList.removeName ( "mod_AA_file" );
	paramList.removeName ( "motif_offset" );
	paramList.removeName ( "motif" );
	string projectID = MySQLPPSDDBase::instance ().getProjectID ( userID, projectName );
	if ( projectID.empty () ) {
		genUnlink ( uploadFpath );
		ErrorHandler::genError ()->error ( "Project doesn't exist.\n" );
	}
	if ( MySQLPPSDDBase::instance ().checkResults ( userID, projectID, resultsName ) ) {
		genUnlink ( uploadFpath );
		ErrorHandler::genError ()->error ( "Results file of this name already exists for this project.\n" );
	}
	paramList.XMLParameterFile ( reposit->getFullResultsPath () + SLASH + searchKey + ".xml" );	// This is the search parameter file
	MySQLPPSDDBase::instance ().submitSearch ( searchKey, projectID, resultsName, searchKey + ".xml", reposit->getResultsPath () );
	MySQLPPSDDBase::instance ().updateProjectRecordUpdated ( projectID );

	genSleep ( 5000 );
	ParameterList pList ( "jobStatus", false, false, false, false, false );
	pList.addName ( "search_key", searchKey );
	refreshJavascript ( os, 0, pList.getURL (), true );
}
int main ( int argc, char** argv )
{
	initialiseProspector ();

	ParameterList paramList ( argc, argv );							// Process the cgi input including the uploaded file
	try {
		ostream& os = cout;
		if ( paramList.getBoolValue ( "save_params" ) )
			saveParams ( os, paramList );
		else
			batchSubmission ( os, paramList );
	}
	catch ( runtime_error e ) {
		paramList.writeLogError ( e.what () );
	}
	MySQLPPSDDBase::instance ( false, true );
	return 0;
}
