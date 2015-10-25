/******************************************************************************
*                                                                             *
*  Library    : libdbase                                                      *
*                                                                             *
*  Filename   : ld_init.cpp                                                   *
*                                                                             *
*  Created    : February 27th 2007                                            *
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
#include <algorithm>
#include <lg_string.h>
#include <lg_time.h>
#include <lgen_file.h>
#include <lgen_error.h>
#include <lgen_process.h>
#include <lgen_net.h>
#include <lu_getfil.h>
#include <lu_table.h>
#include <lu_delim.h>
#include <lu_param_list.h>
#include <lu_prog.h>
#include <lu_repository.h>
#include <lu_xml.h>
#include <ld_init.h>
#include <my_global.h>
#include <mysql.h>
#include <errmsg.h>
using std::cout;
using std::endl;
using std::sort;
using std::string;
using std::ostream;
using std::ostringstream;
using std::endl;
using std::pair;
using std::make_pair;
using std::find;
using std::map;
using std::runtime_error;

#ifdef VIS_C
#pragma warning( disable : 4503 )	// Disable warnings about decorated name length being too long
#endif

#ifdef VIS_C		// mysql header file config-win.h incorrectly defines bool
#undef bool
#endif

struct DatabaseImportSearchJobInfo {
	string searchStage;
	string searchNumber;
	string numSerial;
	string resultsName;
	string priority;
	string jobStatus;
	string jobSignal;
	string nodeName;
	string percentComplete;
	string searchSubmitted;
	string searchStarted;
	string searchFinished;
	string recordCreated2;

DatabaseImportSearchJobInfo ( const string& searchStage, const string& searchNumber, const string& numSerial, const string& resultsName, const string& priority, const string& jobStatus, const string& jobSignal, const string& nodeName, const string& percentComplete, const string& searchSubmitted, const string& searchStarted, const string& searchFinished, const string& recordCreated2 ) :
	searchStage ( searchStage ),
	searchNumber ( searchNumber ),
	numSerial ( numSerial ),
	resultsName ( resultsName ),
	priority ( priority ),
	jobStatus ( jobStatus ),
	jobSignal ( jobSignal ),
	nodeName ( nodeName ),
	percentComplete ( percentComplete ),
	searchSubmitted ( searchSubmitted ),
	searchStarted ( searchStarted ),
	searchFinished ( searchFinished ),
	recordCreated2 ( recordCreated2 ) {}
};

enum {
	JOB_STATUS_SUBMITTED = 0,
	JOB_STATUS_START = 1,
	JOB_STATUS_SEARCHING = 2,
	JOB_STATUS_MISSING = 3,
	JOB_STATUS_DONE = 4,
	JOB_STATUS_ABORTED_TIMEOUT = 5,
	JOB_STATUS_ABORTED_UNKNOWN = 6,
	JOB_STATUS_ABORTED_USER = 7
};

enum {
	JOB_SIGNAL_RUN = 0,
	JOB_SIGNAL_SUSPEND = 1,
	JOB_SIGNAL_HALT = 2,
	JOB_SIGNAL_ABORT = 3,
	JOB_SIGNAL_RESUME_SUSPENDED = 4,
	JOB_SIGNAL_RESUME_HALTED = 5
};

const string MySQLPPSDDBase::QUOTE = "'";
const string MySQLPPSDDBase::EQUOTE = "',";
const string MySQLPPSDDBase::searchProgram = "batchtag";

char* MySQLPPSDDBase::sock = NULL;
unsigned int MySQLPPSDDBase::flags = 0;

string MySQLPPSDDBase::getJobStatus ( int jobStatus )
{
	switch ( jobStatus ) {
		case JOB_STATUS_SUBMITTED:
			return "submitted";
		case JOB_STATUS_START:
			return "start";
		case JOB_STATUS_SEARCHING:
			return "searching";
		case JOB_STATUS_DONE:
			return "done";
		case JOB_STATUS_ABORTED_TIMEOUT:
			return "aborted - timeout";
		case JOB_STATUS_ABORTED_UNKNOWN:
			return "aborted - unknown";
		case JOB_STATUS_ABORTED_USER:
			return "aborted - user";
	}
	return "illegal job status";
}
string MySQLPPSDDBase::getJobSignal ( int jobSignal )
{
	switch ( jobSignal ) {
		case JOB_SIGNAL_RUN:
			return "run";
		case JOB_SIGNAL_SUSPEND:
			return "suspend";
		case JOB_SIGNAL_HALT:
			return "halt";
		case JOB_SIGNAL_ABORT:
			return "abort";
		case JOB_SIGNAL_RESUME_SUSPENDED:
			return "resume suspended";
		case JOB_SIGNAL_RESUME_HALTED:
			return "resume halted";
	}
	return "illegal job signal";
}
string MySQLPPSDDBase::getCellValue ( const MYSQL_ROW& r, int i )
{
	return r [i] == 0 ? "NULL" : r [i];
}
unsigned int MySQLPPSDDBase::submitSQL ( const string& sql ) const
{
	//pidLogOutput ( "submitSQL1" ); 
	unsigned int numInserted = 0;
	if ( !mysql_query ( conn, sql.c_str () ) ) {
		numInserted = static_cast<unsigned int> ( mysql_affected_rows ( conn ) );
	}
	//pidLogOutput ( "submitSQL2" ); 
	return numInserted;
}
unsigned int MySQLPPSDDBase::deleteQuery ( const string& table, const string& field, const string& value ) const
{
	//pidLogOutput ( "deleteQuery" );
	string sql = "DELETE FROM " + table + " WHERE " + field + " = '" + value + "'";
	return submitSQL ( sql );
}
unsigned int MySQLPPSDDBase::updateQuery ( const string& table, const PairStringString& setClause, const PairStringString& whereClause ) const
{
	//pidLogOutput ( "updateQuery" );
	string sql = "UPDATE " + table + " ";
	sql += "SET " + setClause.first + " = ";
	const string& sec = setClause.second;
	if ( sec == "NULL" || isSuffix ( sec, "()" ) )
		sql += sec;
	else
		sql += "'" + sec + "'";
	sql += " WHERE " + whereClause.first + " = '" + whereClause.second + "'";
	return submitSQL ( sql );
}
unsigned int MySQLPPSDDBase::insertQuery ( const string& table, const StringVector& names, const StringVector& values ) const
{
	//pidLogOutput ( "insertQuery1" );
	string sql = "INSERT INTO " + table + " (";
	for ( StringVectorSizeType i = 0 ; i < names.size () ; i++ ) {
		sql += names [i];
		if ( i != names.size () - 1 ) sql += ",";
	}
	sql += ") ";
	sql += "VALUES ( ";
	for ( StringVectorSizeType j = 0 ; j < values.size () ; j++ ) {
		sql += values [j];
		if ( j != values.size () - 1 ) sql += ",";
	}
	sql += ")";
	unsigned int numInserted = submitSQL ( sql );
	return numInserted;
}
bool MySQLPPSDDBase::checkKeyPresence ( const string& sql ) const		// Check for presence of some key
{
	//pidLogOutput ( "checkKeyPresence1" );
	int rows = 0;
	if ( !mysql_query ( conn, sql.c_str () ) ) {
		MYSQL_RES* res = mysql_store_result ( conn );
		if ( res == NULL ) {
			return true;		// There is a problem so true is returned
		}
		else {
			MYSQL_ROW row;
			while ( row = mysql_fetch_row ( res ) ) {
				rows = atoi ( row [0] );		// Number of rows containing the key (should be zero)
			}
			mysql_free_result ( res );
		}
	}
	//pidLogOutput ( "checkKeyPresence2" );
	return rows != 0;		// Returns true if key is present or there is an error
}
MySQLPPSDDBase::MySQLPPSDDBase ()
{
	//pidLogOutput ( "MySQLPPSDDBase1" );
	mysql_library_init(0, NULL, NULL);		// XP systems sometimes very slow here
	string host	= InfoParams::instance ().getStringValue ( "db_host", "localhost" );
	string user	= InfoParams::instance ().getStringValue ( "db_user", "prospector" );
	string pword = InfoParams::instance ().getStringValue ( "db_password", "pp" );
	string db = InfoParams::instance ().getStringValue ( "db_name", "ppsd" );
	unsigned int port = InfoParams::instance ().getUIntValue ( "db_port", 0 );
	conn = mysql_init ( NULL );
	if ( conn == NULL ) {
		throw runtime_error ( "Problems initiating MYSQL database. Probably a memory problem\n" );
	}
	else {
		if ( mysql_real_connect ( conn, host.c_str (), user.c_str (), pword.c_str (), db.c_str (), port, sock, flags ) == NULL ) {
			mysql_close ( conn );
			ostringstream ost;
			ost << "Problems connecting to the MYSQL database.\n";
			ost << "Error number: " << mysql_errno ( conn ) << '\n';
			ost << "Error: " << mysql_error ( conn ) << '\n';
			ost << "SQL state: " << mysql_sqlstate ( conn ) << '\n';
			throw runtime_error ( ost.str () );
		}
	}
	//pidLogOutput ( "MySQLPPSDDBase2" );
}
MySQLPPSDDBase::MySQLPPSDDBase ( unsigned int conn ) :	// In this state the database connection is closed
	conn ( NULL )
{
}
MySQLPPSDDBase::~MySQLPPSDDBase ()
{
	//pidLogOutput ( "~MySQLPPSDDBase1" );
	if ( conn != NULL ) {
		mysql_close ( conn );
		mysql_library_end();
	}
	//pidLogOutput ( "~MySQLPPSDDBase2" );
}
MySQLPPSDDBase& MySQLPPSDDBase::instance ( bool reset, bool close )
{
	static MySQLPPSDDBase* m = new MySQLPPSDDBase ();
	if ( reset || m->conn == 0 || mysql_ping ( m->conn ) ) {
		delete m;
		m = new MySQLPPSDDBase ();
	}
	if ( close ) {
		delete m;
		m = new MySQLPPSDDBase ( 0 );	// Open a NULL object
	}
	return *m;
}
// Checks if a Search Job Key is unique (not already in the table)
bool MySQLPPSDDBase::checkSearchKey ( const string& searchJobKey ) const
{
	//pidLogOutput ( "checkSearchKey" );
	string sql = "SELECT\
					COUNT(search_job_key)\
				FROM\
					search_jobs\
				WHERE\
					search_job_key = '";
	sql += searchJobKey;
	sql += QUOTE;
	return checkKeyPresence ( sql );
}
// Checks if a project is unique for a given user
bool MySQLPPSDDBase::checkProject ( const string& userID, const string& projectName ) const
{
	//pidLogOutput ( "checkProject" );
	int rows = 0;
	string sql = "SELECT\
					COUNT(project_id)\
				FROM\
					projects\
				WHERE\
					pp_user_id = '" + userID + QUOTE + "\
				AND project_name = '" + projectName + QUOTE;
	return checkKeyPresence ( sql );
}
// Checks if a results name is unique for a given user/project
bool MySQLPPSDDBase::checkResults ( const string& userID, const string& projectID, const string& resultsName ) const
{
	//pidLogOutput ( "checkResults" );
	int rows = 0;
	string sql = "SELECT\
					COUNT(results_name)\
				FROM\
					search_jobs, projects\
				WHERE\
					search_jobs.project_id = '" + projectID + QUOTE + "\
				AND projects.project_id = '" + projectID + QUOTE + "\
				AND search_jobs.results_name = '" + resultsName + QUOTE + "\
				AND search_jobs.job_status <> '" + gen_itoa ( JOB_STATUS_ABORTED_TIMEOUT ) + "'\
				AND search_jobs.job_status <> '" + gen_itoa ( JOB_STATUS_ABORTED_UNKNOWN ) + "'\
				AND search_jobs.job_status <> '" + gen_itoa ( JOB_STATUS_ABORTED_USER ) + "'\
				AND projects.pp_user_id = '" + userID + QUOTE;
	return checkKeyPresence ( sql );
}
// Gets an id number for a given project
string MySQLPPSDDBase::getProjectID ( const string& userID, const string& projectName ) const
{
	//pidLogOutput ( "getProjectID" );
	string sql = "SELECT\
					project_id\
				FROM\
					projects\
				WHERE\
					pp_user_id = '" + userID + QUOTE + "\
				AND project_name = '" + projectName + QUOTE;
	return MySQLPPSDDBase::getSingleString ( conn, sql );
}
// Checks information on a given user
UserInfo* MySQLPPSDDBase::getUserInfo ( const string& userName ) const
{
	//pidLogOutput ( "getUserInfo" );
	UserInfo* userInfo = 0;
	string sql = "SELECT\
					pp_user_id,\
					user_name,\
					secret_phrase,\
					email,\
					directory_name,\
					is_inactive,\
					max_spectra,\
					max_search_time,\
					allow_raw\
				FROM\
					pp_users\
				WHERE\
					user_name = '" + userName + QUOTE;
	if ( !mysql_query ( conn, sql.c_str () ) ) {
		MYSQL_RES* res = mysql_store_result ( conn );
		if ( res == NULL ) {
			return userInfo;
		}
		else {
			MYSQL_ROW row;
			while ( row = mysql_fetch_row ( res ) ) {
				userInfo = new UserInfo ( row );
			}
			mysql_free_result ( res );
		}
	}
	return userInfo;
}
// Checks information on a given project
ProjectInfo* MySQLPPSDDBase::getProjectInfo ( const string& user, const string& project ) const
{
	//pidLogOutput ( "getProjectInfo" );
	ProjectInfo* projectInfo = 0;
	string sql = "SELECT\
					project_id,\
					project_name,\
					project_file,\
					project_path,\
					calibration_index\
				FROM projects, pp_users\
				WHERE\
					projects.project_name = '" + project + QUOTE + "\
				AND pp_users.user_name = '" + user + QUOTE + "\
				AND projects.pp_user_id = pp_users.pp_user_id";
	if ( !mysql_query ( conn, sql.c_str () ) ) {
		MYSQL_RES* res = mysql_store_result ( conn );
		if ( res == NULL ) {
			return projectInfo;
		}
		else {
			MYSQL_ROW row;
			while ( row = mysql_fetch_row ( res ) ) {
				projectInfo = new ProjectInfo ( row );
			}
			mysql_free_result ( res );
		}
	}
	return projectInfo;
}
// Get user id for a given user
string MySQLPPSDDBase::getUserID ( const string& userName ) const
{
	//pidLogOutput ( "getUserID" );
	string sql = "SELECT\
					pp_user_id\
				FROM\
					pp_users\
				WHERE\
					user_name = '" + userName + QUOTE;
	return MySQLPPSDDBase::getSingleString ( conn, sql );
}
// Get user name for a given session key
string MySQLPPSDDBase::getUserName ( const string& key ) const
{
	//pidLogOutput ( "getUserName" );
	string sql = "SELECT\
					user_name\
				FROM\
					pp_users, sessions\
				WHERE\
					sessions.session_key = '" + key + QUOTE + "\
				AND pp_users.pp_user_id = sessions.pp_user_id";
	return MySQLPPSDDBase::getSingleString ( conn, sql );
}
// Submits a user
bool MySQLPPSDDBase::submitUser ( const string& userName, const string& secretPhrase, const string& email, const string& directoryName, const string& firstName, const string& lastName ) const
{
	//pidLogOutput ( "submitUser" );
	StringVector names;
	names.push_back ( "user_name" );
	names.push_back ( "secret_phrase" );
	names.push_back ( "email" );
	names.push_back ( "directory_name" );
	names.push_back ( "first_name" );
	names.push_back ( "last_name" );
	names.push_back ( "record_created" );
	StringVector values;
	values.push_back ( QUOTE + userName + QUOTE );
	values.push_back ( QUOTE + secretPhrase + QUOTE );
	values.push_back ( QUOTE + email + QUOTE );
	values.push_back ( QUOTE + directoryName + QUOTE );
	values.push_back ( QUOTE + firstName + QUOTE );
	values.push_back ( QUOTE + lastName + QUOTE );
	values.push_back ( "NOW()" );
	return insertQuery ( "pp_users", names, values ) != 0;
}
// Submits a project
bool MySQLPPSDDBase::submitProject ( const string& userID, const string& projectName, const string& projectFile, const string& projectPath, const string& instrument ) const
{
	//pidLogOutput ( "submitProject" );
	StringVector names;
	names.push_back ( "pp_user_id" );
	names.push_back ( "project_name" );
	names.push_back ( "project_file" );
	names.push_back ( "project_path" );
	names.push_back ( "instrument" );
	names.push_back ( "record_created" );
	StringVector values;
	values.push_back ( QUOTE + userID + QUOTE );
	values.push_back ( QUOTE + projectName + QUOTE );
	values.push_back ( QUOTE + projectFile + QUOTE );
	values.push_back ( QUOTE + genTranslateSlash ( projectPath ) + QUOTE );
	values.push_back ( QUOTE + instrument + QUOTE );
	values.push_back ( "NOW()" );
	return insertQuery ( "projects", names, values ) != 0;
}
// Submits a project
bool MySQLPPSDDBase::submitCopiedProject ( const string& userID, const string& projectName, const string& projectFile, const string& projectPath, const string& instrument, const string& calibrationIndex, const string& recordCreated ) const
{
	//pidLogOutput ( "submitProject" );
	StringVector names;
	names.push_back ( "pp_user_id" );
	names.push_back ( "project_name" );
	names.push_back ( "project_file" );
	names.push_back ( "project_path" );
	names.push_back ( "instrument" );
	names.push_back ( "calibration_index" );
	names.push_back ( "record_created" );
	StringVector values;
	values.push_back ( QUOTE + userID + QUOTE );
	values.push_back ( QUOTE + projectName + QUOTE );
	values.push_back ( QUOTE + projectFile + QUOTE );
	values.push_back ( QUOTE + genTranslateSlash ( projectPath ) + QUOTE );
	values.push_back ( QUOTE + instrument + QUOTE );
	values.push_back ( QUOTE + calibrationIndex + QUOTE );
	values.push_back ( QUOTE + recordCreated + QUOTE );
	unsigned int numInserted = insertQuery ( "projects", names, values );
	if ( numInserted ) {
		cout << "Project = <b>" << projectName << "</b> added to database.<br /><br />" << endl;
	}
	else {
		cout << "Failed to add project = <b>" << projectName << "</b> to database.<br /><br />" << endl;
	}
	return numInserted != 0;
}
// Submits a search
void MySQLPPSDDBase::submitSearch ( const string& searchJobKey, const string& projectID, const string& resultsName, const string& resultsFile, const string& resultsPath ) const
{
	//pidLogOutput ( "submitSearch" );
	StringVector names;
	names.push_back ( "search_job_key" );
	names.push_back ( "search_program" );
	names.push_back ( "project_id" );
	names.push_back ( "results_name" );
	names.push_back ( "results_file" );
	names.push_back ( "results_path" );
	names.push_back ( "search_submitted" );
	names.push_back ( "record_created" );
	StringVector values;
	values.push_back ( QUOTE + searchJobKey + QUOTE );
	values.push_back ( QUOTE + searchProgram + QUOTE );
	values.push_back ( QUOTE + projectID + QUOTE );
	values.push_back ( QUOTE + resultsName + QUOTE );
	values.push_back ( QUOTE + resultsFile + QUOTE );
	values.push_back ( QUOTE + genTranslateSlash ( resultsPath ) + QUOTE );
	values.push_back ( "NOW()" );
	values.push_back ( "NOW()" );
	unsigned int numInserted = insertQuery ( "search_jobs", names, values );
}
// Submits a finished search
void MySQLPPSDDBase::submitCopiedSearch ( const string& searchJobKey, const string& projectID, const string& resultsFile, const string& resultsPath, const DatabaseImportSearchJobInfo* disji ) const
{
	//pidLogOutput ( "submitFinishedSearch" );
	StringVector names;
	names.push_back ( "search_job_key" );
	names.push_back ( "search_program" );
	names.push_back ( "project_id" );
	names.push_back ( "results_name" );
	names.push_back ( "results_file" );
	names.push_back ( "results_path" );
	names.push_back ( "search_submitted" );
	names.push_back ( "record_created" );
	names.push_back ( "search_finished" );
	names.push_back ( "job_status" );
	names.push_back ( "search_stage" );
	names.push_back ( "search_number" );
	names.push_back ( "num_serial" );
	names.push_back ( "priority" );
	names.push_back ( "job_signal" );
	names.push_back ( "node_name" );
	names.push_back ( "percent_complete" );
	names.push_back ( "search_started" );

	StringVector values;
	values.push_back ( QUOTE + searchJobKey + QUOTE );
	values.push_back ( QUOTE + searchProgram + QUOTE );
	values.push_back ( QUOTE + projectID + QUOTE );
	values.push_back ( QUOTE + disji->resultsName + QUOTE );
	values.push_back ( QUOTE + resultsFile + QUOTE );
	values.push_back ( QUOTE + genTranslateSlash ( resultsPath ) + QUOTE );
	values.push_back ( QUOTE + disji->searchSubmitted + QUOTE );
//	values.push_back ( QUOTE + disji->recordCreated2 + QUOTE );
values.push_back ( "NOW()" );
	values.push_back ( QUOTE + disji->searchFinished + QUOTE );
	values.push_back ( QUOTE + disji->jobStatus + QUOTE );
	values.push_back ( QUOTE + disji->searchStage + QUOTE );
	values.push_back ( QUOTE + disji->searchNumber + QUOTE );
	values.push_back ( QUOTE + disji->numSerial + QUOTE );
	values.push_back ( QUOTE + disji->priority + QUOTE );
	values.push_back ( QUOTE + disji->jobSignal + QUOTE );
	values.push_back ( QUOTE + disji->nodeName + QUOTE );
	values.push_back ( QUOTE + disji->percentComplete + QUOTE );
	values.push_back ( QUOTE + disji->searchStarted + QUOTE );
//for ( int i = 0 ; i < names.size () ; i++ ) {
//	std::cout << names [i] << "||||" << values [i] << "<br />" << std::endl;
//}
	unsigned int numInserted = insertQuery ( "search_jobs", names, values );
	if ( numInserted ) {
		cout << "Results = <b>"  << disji->resultsName << "</b> added to database.<br />" << endl;
	}
	else {
		cout << "Failed to add results = <b>" << disji->resultsName << "</b> to database.<br />" << endl;
	}
//std::cout << numInserted << "<br />" << std::endl;
}
// Submits a session
bool MySQLPPSDDBase::submitSession ( const string& sessionKey, const string& ppUserID ) const
{
	//pidLogOutput ( "submitSession" );
	StringVector names;
	names.push_back ( "session_key" );
	names.push_back ( "pp_user_id" );
	names.push_back ( "record_created" );
	StringVector values;
	values.push_back ( QUOTE + sessionKey + QUOTE );
	values.push_back ( QUOTE + ppUserID + QUOTE );
	values.push_back ( "NOW()" );
	return insertQuery ( "sessions", names, values ) != 0;
}
// Submits an aborted job
bool MySQLPPSDDBase::submitAbortedJob ( const string& searchJobKey, const string& message ) const
{
	string writtenMessage = message;
	if ( message.length () > 510 ) writtenMessage = message.substr ( 0, 510 );	// truncates messages that are too long
	//pidLogOutput ( "submitAbortedJob" );
	string sql = "INSERT INTO aborted_jobs (\
		error_message,\
		search_job_key,\
		search_program,\
		search_stage,\
		search_number,\
		num_serial,\
		results_name,\
		results_file,\
		results_path,\
		priority,\
		job_status,\
		job_signal,\
		job_segment,\
		node_name,\
		node_pid,\
		percent_complete,\
		search_submitted,\
		search_started,\
		project_id,\
		project_name,\
		project_file,\
		project_path,\
		instrument,\
		calibration_index,\
		pp_user_id,\
		user_name,\
		email,\
		directory_name,\
		max_search_time,\
		allow_raw,\
		record_created\
		)\
		SELECT\
				'" + writtenMessage + "',\
				a.search_job_key,\
				a.search_program,\
				a.search_stage,\
				a.search_number,\
				a.num_serial,\
				a.results_name,\
				a.results_file,\
				a.results_path,\
				a.priority,\
				a.job_status,\
				a.job_signal,\
				a.job_segment,\
				a.node_name,\
				a.node_pid,\
				a.percent_complete,\
				a.search_submitted,\
				a.search_started,\
				b.project_id,\
				b.project_name,\
				b.project_file,\
				b.project_path,\
				b.instrument,\
				b.calibration_index,\
				b.pp_user_id,\
				c.user_name,\
				c.email,\
				c.directory_name,\
				c.max_search_time,\
				c.allow_raw,\
				NOW()\
			FROM search_jobs a, projects b, pp_users c\
			WHERE\
				a.search_job_key = '" + searchJobKey + "'\
				AND a.project_id = b.project_id\
				AND b.pp_user_id = c.pp_user_id";
	unsigned int numInserted = submitSQL ( sql );
	return numInserted != 0;
}
// Sets the number of serial searches for a given search
void MySQLPPSDDBase::setNumSerial ( const string& searchJobID, int numSerial ) const
{
	//pidLogOutput ( "setNumSerial" );
	PairStringString setClause = make_pair ( string("num_serial"), gen_itoa ( numSerial ) );
	PairStringString whereClause = make_pair ( string("search_job_id"), searchJobID );
	updateQuery ( "search_jobs", setClause, whereClause );
}
// Updates the search stage for a given search
void MySQLPPSDDBase::updateSearchStage ( const string& searchJobID, int searchStage ) const
{
	//pidLogOutput ( "updateSearchStage" );
	PairStringString setClause = make_pair ( string("search_stage"), gen_itoa ( searchStage ) );
	PairStringString whereClause = make_pair ( string("search_job_id"), searchJobID );
	updateQuery ( "search_jobs", setClause, whereClause );
}
// Updates the search number for a given search
void MySQLPPSDDBase::updateSearchNumber ( const string& searchJobID, int searchNumber ) const
{
	//pidLogOutput ( "updateSearchNumber" );
	string sql = "UPDATE search_jobs ";
	sql += "SET search_number = '" + gen_itoa ( searchNumber ) + "', ";
	sql += "percent_complete = '0' ";
	sql += "WHERE search_job_id = '" + searchJobID + QUOTE;
	unsigned int numInserted = submitSQL ( sql );
}
// Updates the percent complete for a given search
void MySQLPPSDDBase::updatePercentComplete ( const string& searchJobID, int percentComplete ) const
{
	//pidLogOutput ( "updatePercentComplete" );
	PairStringString setClause = make_pair ( string("percent_complete"), gen_itoa ( percentComplete ) );
	PairStringString whereClause = make_pair ( string("search_job_id"), searchJobID );
	updateQuery ( "search_jobs", setClause, whereClause );
}
// Sets the abort signal
void MySQLPPSDDBase::setAbortSignal ( const string& searchJobID ) const
{
	//pidLogOutput ( "setAbortSignal" );
	PairStringString setClause = make_pair ( string("job_signal"), gen_itoa ( JOB_SIGNAL_ABORT ) );
	PairStringString whereClause = make_pair ( string("search_job_id"), searchJobID );
	updateQuery ( "search_jobs", setClause, whereClause );
}
// Changes the password
void MySQLPPSDDBase::setPassword ( const string& user, const string& password ) const
{
	//pidLogOutput ( "setPassword" );
	PairStringString setClause = make_pair ( string("secret_phrase"), password );
	PairStringString whereClause = make_pair ( string("user_name"), user );
	updateQuery ( "pp_users", setClause, whereClause );
}
void MySQLPPSDDBase::setJobSubmitted ( const string& searchKey ) const
{
	//pidLogOutput ( "setJobSubmitted" );
	string sql = "UPDATE search_jobs ";
	sql += "SET job_status = '" + gen_itoa ( JOB_STATUS_SUBMITTED ) + "' ";
	sql += "WHERE search_job_key = '" + searchKey + QUOTE;
	unsigned int numInserted = submitSQL ( sql );
}
bool MySQLPPSDDBase::setJobStart ( const string& searchKey ) const	// Only set to start if previously submitted
{
	//pidLogOutput ( "setJobStart" );
	string sql = "UPDATE search_jobs\
				SET job_status = '" + gen_itoa ( JOB_STATUS_START ) + "'\
				WHERE\
					search_job_key = '" + searchKey + QUOTE + "\
				AND job_status = '" + gen_itoa ( JOB_STATUS_SUBMITTED ) + QUOTE;	// Only set the flag if status previously submitted
	unsigned int numInserted = submitSQL ( sql );
	return numInserted != 0;
}
void MySQLPPSDDBase::setJobSearching ( const string& searchKey, const string& host, int pid ) const
{
	//pidLogOutput ( "setJobSearching" );
	string sql = "UPDATE search_jobs ";
	sql += "SET job_status = '" + gen_itoa ( JOB_STATUS_SEARCHING ) + "', ";
	sql += "node_name = '" + host + "', ";
	sql += "node_pid = '" + gen_itoa ( pid ) + "', ";
	sql += "search_started = NOW() ";
	sql += "WHERE search_job_key = '" + searchKey + QUOTE;
	unsigned int numInserted = submitSQL ( sql );
}
void MySQLPPSDDBase::setJobDone ( const string& searchJobID ) const
{
	//pidLogOutput ( "setJobDone" );
	string sql = "UPDATE search_jobs ";
	sql += "SET job_status = '" + gen_itoa ( JOB_STATUS_DONE ) + "', ";
	sql += "search_finished = NOW() ";
	sql += "WHERE search_job_id = '" + searchJobID + QUOTE;
	unsigned int numInserted = submitSQL ( sql );
}
string MySQLPPSDDBase::getSingleString ( MYSQL* conn, const string& sql )
{
	//pidLogOutput ( "getSingleString" );
	string s;
	if ( !mysql_query ( conn, sql.c_str () ) ) {
		MYSQL_RES* res = mysql_store_result ( conn );
		if ( res == NULL ) {
			return s;
		}
		else {
			MYSQL_ROW row;
			while ( row = mysql_fetch_row ( res ) ) {
				s = MySQLPPSDDBase::getCellValue ( row, 0 );
			}
			mysql_free_result ( res );
		}
	}
	return s;
}
JobItem* MySQLPPSDDBase::getJobItem ( MYSQL* conn, const string& sql )
{
	//pidLogOutput ( "getJobItem" );
	JobItem* jobItem = 0;
	if ( !mysql_query ( conn, sql.c_str () ) ) {
		MYSQL_RES* res = mysql_store_result ( conn );
		if ( res == NULL ) {
			return jobItem;
		}
		else {
			MYSQL_ROW row;
			while ( row = mysql_fetch_row ( res ) ) {
				jobItem = new JobItem ( row );
			}
			mysql_free_result ( res );
		}
	}
	return jobItem;
}
// Gets the search job info given a key
JobItem* MySQLPPSDDBase::getSearchJobByKey ( const string& searchJobKey ) const
{
	//pidLogOutput ( "getSearchJobByKey" );
	string sql = "SELECT\
					search_job_id,\
					search_job_key,\
					priority,\
					num_serial,\
					search_stage,\
					search_number,\
					percent_complete,\
					results_name,\
					job_status,\
					search_submitted,\
					TIMESTAMPDIFF( SECOND, search_started, NOW() )\
				FROM search_jobs\
				WHERE search_job_key = '";
	sql += searchJobKey;
	sql += QUOTE;
	return getJobItem ( conn, sql );
}
// Gets a search id number given a key
string MySQLPPSDDBase::getSearchIDByKey ( const string& searchJobKey ) const
{
	//pidLogOutput ( "getSearchIDByKey" );
	string sql = "SELECT\
					search_job_id\
				FROM search_jobs\
				WHERE search_job_key = '";
	sql += searchJobKey;
	sql += QUOTE;
	return MySQLPPSDDBase::getSingleString ( conn, sql );
}
string MySQLPPSDDBase::getUserNameByKey ( const string& searchJobKey ) const
{
	//pidLogOutput ( "getUserNameByKey" );
	string sql = "SELECT\
					pp_users.user_name\
				FROM search_jobs, projects, pp_users\
				WHERE\
					search_job_key = '" + searchJobKey + QUOTE + "\
				AND search_jobs.project_id = projects.project_id\
				AND projects.pp_user_id = pp_users.pp_user_id";
	return MySQLPPSDDBase::getSingleString ( conn, sql );
}
string MySQLPPSDDBase::getSearchTimeByKey ( const string& searchJobKey ) const
{
	//pidLogOutput ( "getSearchTimeByKey" );
	string sql = "SELECT\
					TIMESTAMPDIFF( SECOND, search_started, search_finished )\
				FROM search_jobs\
				WHERE search_job_key = '";
	sql += searchJobKey;
	sql += QUOTE;
	return GenElapsedTime::getTimeString ( atoi ( MySQLPPSDDBase::getSingleString ( conn, sql ).c_str () ) );
}
string MySQLPPSDDBase::getSearchEndTimeByKey ( const string& searchJobKey ) const
{
	//pidLogOutput ( "getSearchEndTime" );
	string sql = "SELECT\
					search_finished\
				FROM search_jobs\
				WHERE search_job_key = '";
	sql += searchJobKey;
	sql += QUOTE;
	return MySQLPPSDDBase::getSingleString ( conn, sql );
}
// Gets a search files given a key
BatchJobItem* MySQLPPSDDBase::getBatchJobByKey ( const string& searchJobKey ) const
{
	//pidLogOutput ( "getBatchJobByKey" );
	BatchJobItem* batchJobItem = 0;
	string sql = "SELECT\
					pp_user_id,\
					projects.project_id,\
					project_name,\
					project_file,\
					project_path,\
					results_name,\
					results_file,\
					results_path,\
					calibration_index\
				FROM search_jobs, projects\
				WHERE\
					search_job_key = '" + searchJobKey + QUOTE + "\
				AND search_jobs.project_id = projects.project_id";

	if ( !mysql_query ( conn, sql.c_str () ) ) {
		MYSQL_RES* res = mysql_store_result ( conn );
		if ( res == NULL ) {
			return batchJobItem;
		}
		else {
			MYSQL_ROW row;
			while ( row = mysql_fetch_row ( res ) ) {
				batchJobItem = new BatchJobItem ( row );
			}
		}
	}
	return batchJobItem;
}
// Sets a job as aborted
void MySQLPPSDDBase::setJobAborted ( const string& searchJobID ) const
{
	//pidLogOutput ( "setJobAborted" );
	PairStringString setClause = make_pair ( string("job_status"), gen_itoa ( JOB_STATUS_ABORTED_USER ) );
	PairStringString whereClause = make_pair ( string("search_job_id"), searchJobID );
	updateQuery ( "search_jobs", setClause, whereClause );
}
// Sets a job as aborted unknown
bool MySQLPPSDDBase::setJobAbortedUnknown ( const string& searchJobID ) const
{
	//pidLogOutput ( "setJobAbortedUnknown" );
	string sql = "UPDATE search_jobs\
				SET job_status = '" + gen_itoa ( JOB_STATUS_ABORTED_UNKNOWN ) + "'\
				WHERE\
					search_job_id = '" + searchJobID + QUOTE + "\
				AND job_status = '" + gen_itoa ( JOB_STATUS_SEARCHING ) + QUOTE;	// Only set the flag if process is still listed as running
	unsigned int numInserted = submitSQL ( sql );
	return numInserted != 0;
}
bool MySQLPPSDDBase::processAbortSignal ( const string& searchJobID ) const
{
	//pidLogOutput ( "processAbortSignal" );
	string sql = "UPDATE search_jobs\
				SET job_status = '" + gen_itoa ( JOB_STATUS_ABORTED_USER ) + "', \
					job_signal = '" + gen_itoa ( JOB_SIGNAL_RUN ) + "' \
				WHERE\
					search_job_id = '" + searchJobID + QUOTE + "\
				AND job_signal = '" + gen_itoa ( JOB_SIGNAL_ABORT ) + QUOTE;	// Only set the flag if process is still listed as running
	unsigned int numInserted = submitSQL ( sql );
	return numInserted != 0;
}
// Gets the job queue
JobQueue MySQLPPSDDBase::getJobQueue () const
{
	JobQueue jobQueue;
	string sql = "SELECT\
					a.search_job_id,\
					a.search_job_key,\
					a.priority,\
					a.num_serial,\
					a.search_stage,\
					a.search_number,\
					a.percent_complete,\
					a.results_name,\
					a.job_status,\
					a.search_submitted,\
					TIMESTAMPDIFF( SECOND, a.search_started, NOW() )\
				FROM search_jobs a\
				WHERE\
						a.job_status = '" + gen_itoa ( JOB_STATUS_SEARCHING ) + "'\
					OR	a.job_status = '" + gen_itoa ( JOB_STATUS_START ) + "'\
					OR	a.job_status = '" + gen_itoa ( JOB_STATUS_SUBMITTED ) + "'\
				ORDER BY	a.priority,\
							a.search_job_id";
	if ( !mysql_query ( conn, sql.c_str () ) ) {
		MYSQL_RES* res = mysql_store_result ( conn );
		if ( res == NULL ) {
			return jobQueue;
		}
		else {
			MYSQL_ROW row;
			while ( row = mysql_fetch_row ( res ) ) {
				jobQueue.addItem ( row );
			}
			mysql_free_result ( res );
		}
	}
	return jobQueue;
}
// Gets the job queue for the search table
JobQueueForSearchTable MySQLPPSDDBase::getSearchTableJobQueue () const
{
	JobQueueForSearchTable jobQueue;
	string sql = "SELECT\
					a.search_job_id,\
					a.search_job_key,\
					a.num_serial,\
					a.search_stage,\
					a.search_number,\
					a.percent_complete,\
					a.node_name,\
					a.results_name,\
					a.job_status,\
					TIMESTAMPDIFF( SECOND, a.search_started, NOW() ),\
					b.project_name,\
					c.user_name,\
					c.email,\
					c.directory_name\
				FROM search_jobs a, projects b, pp_users c\
				WHERE\
						( a.job_status = '" + gen_itoa ( JOB_STATUS_SEARCHING ) + "'\
					OR	a.job_status = '" + gen_itoa ( JOB_STATUS_START ) + "'\
					OR	a.job_status = '" + gen_itoa ( JOB_STATUS_SUBMITTED ) + "' )\
					AND a.project_id = b.project_id\
					AND b.pp_user_id = c.pp_user_id\
				ORDER BY	a.priority,\
							a.search_job_id";
	if ( !mysql_query ( conn, sql.c_str () ) ) {
		MYSQL_RES* res = mysql_store_result ( conn );
		if ( res == NULL ) {
			return jobQueue;
		}
		else {
			MYSQL_ROW row;
			while ( row = mysql_fetch_row ( res ) ) {
				jobQueue.addItem ( row );
			}
			mysql_free_result ( res );
		}
	}
	return jobQueue;
}
DaemonJobQueue MySQLPPSDDBase::getDaemonJobQueue () const
{
	DaemonJobQueue jobQueue;
	string sql = "SELECT\
					a.search_job_id,\
					a.search_job_key,\
					a.priority,\
					a.results_name,\
					a.results_path,\
					a.results_file,\
					a.node_name,\
					a.node_pid,\
					a.job_status,\
					a.job_signal,\
					b.project_name,\
					b.project_file,\
					b.project_path,\
					b.project_id,\
					c.email\
				FROM search_jobs a, projects b, pp_users c\
				WHERE\
					( a.job_status = '" + gen_itoa ( JOB_STATUS_SEARCHING ) + "' OR a.job_status = '" + gen_itoa ( JOB_STATUS_SUBMITTED ) + "' )\
					AND a.project_id = b.project_id\
					AND b.pp_user_id = c.pp_user_id\
				ORDER BY	a.priority,\
							a.search_job_id";
	if ( !mysql_query ( conn, sql.c_str () ) ) {
		MYSQL_RES* res = mysql_store_result ( conn );
		if ( res == NULL ) {
			return jobQueue;
		}
		else {
			MYSQL_ROW row;
			while ( row = mysql_fetch_row ( res ) ) {
				jobQueue.addItem ( row );
			}
			mysql_free_result ( res );
		}
	}
	return jobQueue;
}
DaemonJobItem* MySQLPPSDDBase::getDaemonJobItemByKey ( const string& searchJobKey ) const
{
	DaemonJobItem* dji = 0;
	string sql = "SELECT\
					a.search_job_id,\
					a.search_job_key,\
					a.priority,\
					a.results_name,\
					a.results_path,\
					a.results_file,\
					a.node_name,\
					a.node_pid,\
					a.job_status,\
					a.job_signal,\
					b.project_name,\
					b.project_file,\
					b.project_path,\
					b.project_id,\
					c.email\
				FROM search_jobs a, projects b, pp_users c\
				WHERE\
					a.search_job_key = '" + searchJobKey + "'\
					AND a.project_id = b.project_id\
					AND b.pp_user_id = c.pp_user_id";
	if ( !mysql_query ( conn, sql.c_str () ) ) {
		MYSQL_RES* res = mysql_store_result ( conn );
		if ( res == NULL ) {
			return dji;
		}
		else {
			MYSQL_ROW row;
			while ( row = mysql_fetch_row ( res ) ) {
				dji = new DaemonJobItem ( row );
			}
		}
	}
	return dji;
}
void JobQueue::addItem ( const MYSQL_ROW& row )
{
	jobItems.push_back ( JobItem ( row ) );
}
int JobQueue::numRunningJobs () const
{
	int num = 0;
	for ( std::vector <JobItem>::size_type i = 0 ; i < jobItems.size () ; i++ ) {
		if ( jobItems [i].isRunning () ) num++;
	}
	return num;
}
int JobQueue::numJobsBeforeJob ( const string& searchJobKey ) const
{
	int num = 0;
	for ( std::vector <JobItem>::size_type i = 0 ; i < jobItems.size () ; i++ ) {
		if ( jobItems [i].getSearchJobKey () != searchJobKey ) num++;
		else break;
	}
	return num;
}
void JobQueue::printJobsBeforeJob ( ostream& os, const string& searchJobKey ) const
{
	int ind = numJobsBeforeJob ( searchJobKey );
	os << "<p>Search job [" << searchJobKey << "] was submitted at " << jobItems [ind].getSearchSubmitted () << ".</p>" << endl;
	os << "<p>There are:" << endl;
	os << "\t<ul>" << endl;
	os << "\t\t<li> " << numRunningJobs () << " Running searches" << endl;
	os << "\t\t<li> " << ind << " Searches in the queue before yours" << endl;
	os << "\t</ul>" << endl;
	os << "</p>" << endl;
}
void JobQueueForSearchTable::addItem ( const MYSQL_ROW& row )
{
	jobItems.push_back ( JobItemForSearchTable ( row ) );
}
void JobQueueForSearchTable::printHTML ( ostream& os, const string& user ) const
{
	if ( jobItems.empty () ) {
		os << "<p>" << endl;
		os << "<b>No searches are running or submitted.</b>" << endl;
		os << "</p>" << endl;
	}
	else {
		tableStart ( os, true );
		for ( std::vector <JobItem>::size_type i = 0 ; i < jobItems.size () ; i++ ) {
			if ( i == 0 ) JobItemForSearchTable::printHTMLHeader ( os, user );
			jobItems [i].printHTML ( os, user );
		}
		tableEnd ( os );
	}
}
JobItem::JobItem ( const MYSQL_ROW& row )
{
	int r = 0;
	searchJobID		= MySQLPPSDDBase::getCellValue ( row, r++ );
	searchJobKey	= MySQLPPSDDBase::getCellValue ( row, r++ );
	priority		= MySQLPPSDDBase::getCellValue ( row, r++ );
	numSerial		= atoi ( MySQLPPSDDBase::getCellValue ( row, r++ ).c_str () );
	searchStage		= atoi ( MySQLPPSDDBase::getCellValue ( row, r++ ).c_str () );
	searchNumber	= atoi ( MySQLPPSDDBase::getCellValue ( row, r++ ).c_str () );
	percentComplete	= atoi ( MySQLPPSDDBase::getCellValue ( row, r++ ).c_str () );
	resultsName		= MySQLPPSDDBase::getCellValue ( row, r++ );
	jobStatus		= atoi ( MySQLPPSDDBase::getCellValue ( row, r++ ).c_str () );
	searchSubmitted	= MySQLPPSDDBase::getCellValue ( row, r++ );
	searchTime		= GenElapsedTime::getTimeString ( atoi ( MySQLPPSDDBase::getCellValue ( row, r++ ).c_str () ) );
}
bool JobItem::isSubmittedOrStart () const
{
	return jobStatus == JOB_STATUS_SUBMITTED || jobStatus == JOB_STATUS_START;
}
bool JobItem::isRunning () const
{
	return jobStatus == JOB_STATUS_SEARCHING;
}
bool JobItem::jobDone () const
{
	return jobStatus == JOB_STATUS_DONE;
}
string JobItem::getJobStatus () const
{
	return MySQLPPSDDBase::getJobStatus ( jobStatus );
}
void JobItem::printJobProgress ( ostream& os ) const
{
	os << "<p>" << endl;
	os << "Search job " << searchJobID << " has been running for " << searchTime;
	os << "<br />";
	if ( searchStage == 1 ) os << "Doing expectation value calibration";
	else os << "Doing search";
	os << "<br />";
	os << "Search number " << searchNumber << '/' << numSerial;
	os << "<br />";
	os << percentComplete;
	os << "% complete" << endl;
	os << "</p>" << endl;
}
JobItemForSearchTable::JobItemForSearchTable ( const MYSQL_ROW& row )
{
	int r = 0;
	searchJobID		= MySQLPPSDDBase::getCellValue ( row, r++ );
	searchJobKey	= MySQLPPSDDBase::getCellValue ( row, r++ );
	numSerial		= atoi ( MySQLPPSDDBase::getCellValue ( row, r++ ).c_str () );
	searchStage		= atoi ( MySQLPPSDDBase::getCellValue ( row, r++ ).c_str () );
	searchNumber	= atoi ( MySQLPPSDDBase::getCellValue ( row, r++ ).c_str () );
	percentComplete	= atoi ( MySQLPPSDDBase::getCellValue ( row, r++ ).c_str () );
	nodeName		= MySQLPPSDDBase::getCellValue ( row, r++ );
	resultsName		= MySQLPPSDDBase::getCellValue ( row, r++ );
	jobStatus		= atoi ( MySQLPPSDDBase::getCellValue ( row, r++ ).c_str () );
	searchTime		= GenElapsedTime::getTimeString ( atoi ( MySQLPPSDDBase::getCellValue ( row, r++ ).c_str () ) );
	projectName		= MySQLPPSDDBase::getCellValue ( row, r++ );
	userName		= MySQLPPSDDBase::getCellValue ( row, r++ );
	email			= MySQLPPSDDBase::getCellValue ( row, r++ );
	directoryName	= MySQLPPSDDBase::getCellValue ( row, r++ );
}
void JobItemForSearchTable::printHTMLHeader ( ostream& os, const string& user )
{
	tableRowStart ( os );
		tableHeader ( os, "ID" );
		tableHeader ( os, "Project" );
		tableHeader ( os, "Results" );
		tableHeader ( os, "Status" );
		tableHeader ( os, "Stage" );
		tableHeader ( os, "Number" );
		tableHeader ( os, "%" );
		if ( user == "root" ) {
			tableHeader ( os, "Node" );
		}
		tableHeader ( os, "Search Time" );
		if ( user == "root" ) {
			tableHeader ( os, "User" );
			tableHeader ( os, "Email" );
			tableHeader ( os, "Directory" );
		}
	tableRowEnd ( os );
}
void JobItemForSearchTable::printHTML ( ostream& os, const string& user ) const
{
	bool rootUser = ( user == "root" );
	bool currentUser = ( user == userName );
	bool rootOrCurrent = ( rootUser || currentUser );
	tableRowStart ( os );
		string s;
		if ( rootOrCurrent )
			s = "<a href=\"" + ProgramLink::getURLStart ( "jobStatus" ) + "?search_key=" + searchJobKey + "\">" + searchJobID + "</a>";
		else
			s = searchJobID;
		tableCell ( os, s );
		if ( rootOrCurrent ) {
			tableCell ( os, projectName );
			tableCell ( os, resultsName );
		}
		else {
			tableEmptyNCells ( os, 2 );
		}
		tableCell ( os, MySQLPPSDDBase::getJobStatus ( jobStatus ) );
		if ( jobStatus == JOB_STATUS_SEARCHING ) {
			if ( searchStage == 1 ) tableCell ( os, "expectation" );
			else tableCell ( os, "search" );
			tableCell ( os, gen_itoa ( searchNumber ) + '/' + gen_itoa ( numSerial ) );
			tableCell ( os, percentComplete );
			if ( user == "root" ) {
				tableCell ( os, nodeName );
			}
			tableCell ( os, searchTime );
		}
		else {
			if ( user == "root" )	tableEmptyNCells ( os, 5 );
			else					tableEmptyNCells ( os, 4 );
		}
		if ( rootUser ) {
			tableCell ( os, userName );
			tableCell ( os, "<a href=\"mailto:" + email + "\">" + email + "</a>" );
			tableCell ( os, directoryName );
		}
	tableRowEnd ( os );
}
DaemonJobItem::DaemonJobItem ( const MYSQL_ROW& row )
{
	int r = 0;
	searchJobID		= MySQLPPSDDBase::getCellValue ( row, r++ );
	searchJobKey	= MySQLPPSDDBase::getCellValue ( row, r++ );
	priority		= atoi ( MySQLPPSDDBase::getCellValue ( row, r++ ).c_str () );
	resultsName		= MySQLPPSDDBase::getCellValue ( row, r++ );
	resultsPath		= MySQLPPSDDBase::getCellValue ( row, r++ );
	resultsFile		= MySQLPPSDDBase::getCellValue ( row, r++ );
	nodeName		= MySQLPPSDDBase::getCellValue ( row, r++ );
	nodePID			= atoi ( MySQLPPSDDBase::getCellValue ( row, r++ ).c_str () );
	jobStatus		= atoi ( MySQLPPSDDBase::getCellValue ( row, r++ ).c_str () );
	jobSignal		= atoi ( MySQLPPSDDBase::getCellValue ( row, r++ ).c_str () );
	projectName		= MySQLPPSDDBase::getCellValue ( row, r++ );
	projectFile		= MySQLPPSDDBase::getCellValue ( row, r++ );
	projectPath		= MySQLPPSDDBase::getCellValue ( row, r++ );
	projectID		= MySQLPPSDDBase::getCellValue ( row, r++ );
	email			= MySQLPPSDDBase::getCellValue ( row, r++ );
}
bool DaemonJobItem::getJobSignalAbort () const
{
	return jobSignal == JOB_SIGNAL_ABORT;
}
bool DaemonJobItem::isSubmitted () const
{
	return jobStatus == JOB_STATUS_SUBMITTED;
}
bool DaemonJobItem::isRunningOrStarted () const
{
	return jobStatus == JOB_STATUS_SEARCHING || jobStatus == JOB_STATUS_START;
}
bool DaemonJobItem::isRunning () const
{
	return jobStatus == JOB_STATUS_SEARCHING;
}
string DaemonJobItem::getJobStatus () const
{
	return MySQLPPSDDBase::getJobStatus ( jobStatus );
}
string DaemonJobItem::getJobSignal () const
{
	return MySQLPPSDDBase::getJobSignal ( jobSignal );
}
DaemonJobQueue::DaemonJobQueue () :
	numRunning ( 0 )
{
}
void DaemonJobQueue::addItem ( const MYSQL_ROW& row )
{
	jobItems.push_back ( DaemonJobItem ( row ) );
}
void DaemonJobQueue::setActions ( const string& host, int maxJobsPerUser )
{
	numRunning = 0;
	nextSearchJobKey = "";
	SetString projRunning;
	for ( VectorDaemonJobItem::size_type ii = 0 ; ii < jobItems.size () ; ii++ ) {	// Get a list of projects with running jobs. Only 1 job per project
																					// can run at the same time
		const DaemonJobItem& ji = jobItems [ii];
		if ( ji.isRunningOrStarted () ) projRunning.insert ( ji.getProjectID () );
	}
	if ( jobItems.empty () ) return;
	IntVector runningPIDs = getProcessNumberList ();			// Get a list of all running processes
	MapPairStringStringToInt mpssi;
	for ( VectorDaemonJobItem::size_type i = 0 ; i < jobItems.size () ; i++ ) {
		DaemonJobItem& ji = jobItems [i];
		string searchJobID = ji.getSearchJobID ();
		string searchJobKey = ji.getSearchJobKey ();
		string nodeName = ji.getNodeName ();
		bool localJob = ( nodeName == host );
		bool notStarted = ( nodeName == "NULL" );
		int pid = ji.getPID ();
		bool limitPerUser = false;
		string projDir = ji.getProjDir ();
		if ( maxJobsPerUser ) {
			if ( nodeName != "NULL" ) {
				PairStringString pss = make_pair ( nodeName, projDir );
				MapPairStringStringToIntIterator cur = mpssi.find ( pss );
				if ( cur != mpssi.end () )	cur->second++;
				else						mpssi [pss] = 1;
			}
			else {
				PairStringString pss = make_pair ( host, projDir );
				MapPairStringStringToIntIterator cur = mpssi.find ( pss );
				if ( cur != mpssi.end () && cur->second >= maxJobsPerUser ) {
					limitPerUser = true;
				}
			}

		}
		if ( ji.getJobSignalAbort () ) {
			if ( localJob ) {			// Only abort searches on the local daemon
				killProcess ( pid );	// Kill the process. Doesn't matter if it isn't running
				MySQLPPSDDBase::instance ().processAbortSignal ( searchJobID );
				MySQLPPSDDBase::instance ().submitAbortedJob ( searchJobKey, "This search was aborted by the user." );
				cleanUpJobItems.push_back ( ji );
			}
			if ( notStarted ) {			// Abort pressed before search started
				if ( MySQLPPSDDBase::instance ().processAbortSignal ( searchJobID ) ) {	// There may be multiple daemons running. Only one should process the abort.
					MySQLPPSDDBase::instance ().submitAbortedJob ( searchJobKey, "This search was aborted by the user." );
					cleanUpJobItems.push_back ( ji );
				}
			}
		}
		else if ( ji.isRunning () && localJob ) {
			IntVectorIterator res = find ( runningPIDs.begin (), runningPIDs.end (), pid );
			if ( res == runningPIDs.end () ) {
				if ( MySQLPPSDDBase::instance ().setJobAbortedUnknown ( searchJobID ) ) {	// Search job has finished unexpectedly
					string err = "This search was aborted for an unknown reason.";
#ifndef VIS_C
					string tmpErrFile = "/tmp/" + searchJobKey + "_err";
					if ( genFileExists ( tmpErrFile ) ) {
						int numEntries;
						char* ch = getFileInfo ( tmpErrFile, '\n', 1, false, &numEntries );
						if ( numEntries == 1 ) err = "The script mssearchmpi.pl returned the following message: \n" + string ( ch );
						delete [] ch;
						genUnlink ( tmpErrFile );
					}
#endif
					MySQLPPSDDBase::instance ().submitAbortedJob ( searchJobKey, err );
					cleanUpJobItems.push_back ( ji );
				}
				else {		// Search may have very recently finished
					activeJobKeys.insert ( searchJobKey );
				}
			}
			else {			// Search still running
				activeJobKeys.insert ( searchJobKey );
				numRunning++;
			}
		}
		else {																														// Other submitted jobs
			if ( nextSearchJobKey.empty () && ji.isSubmitted () && !projRunning.count ( ji.getProjectID () ) && !limitPerUser ) {	// Next job in queue
				nextSearchJobKey = searchJobKey;
			}
		}
	}
}
BatchJobItem::BatchJobItem ( const MYSQL_ROW& row )
{
	int r = 0;
	ppUserID				= MySQLPPSDDBase::getCellValue ( row, r++ );
	projectID				= MySQLPPSDDBase::getCellValue ( row, r++ );
	projectName				= MySQLPPSDDBase::getCellValue ( row, r++ );
	string projectFile		= MySQLPPSDDBase::getCellValue ( row, r++ );
	projectRelativePath		= MySQLPPSDDBase::getCellValue ( row, r++ );
	resultsName				= MySQLPPSDDBase::getCellValue ( row, r++ );
	string resultsFile		= MySQLPPSDDBase::getCellValue ( row, r++ );
	resultsPath				= MySQLPPSDDBase::getCellValue ( row, r++ );
	calibrationIndex		= atoi ( MySQLPPSDDBase::getCellValue ( row, r++ ).c_str () );

	static string baseDir = InfoParams::instance ().getUserRepository ();
	resultsPath		= baseDir + SLASH + resultsPath;
	resultsFullPath	= resultsPath + SLASH + resultsFile;
	projectPath		= baseDir + SLASH + projectRelativePath;
	projectFullPath	= projectPath + SLASH + projectFile;
}
UserInfo::UserInfo ( const MYSQL_ROW& row )
{
	int r = 0;
	userID			= MySQLPPSDDBase::getCellValue ( row, r++ );
	userName		= MySQLPPSDDBase::getCellValue ( row, r++ );
	password		= MySQLPPSDDBase::getCellValue ( row, r++ );
	email			= MySQLPPSDDBase::getCellValue ( row, r++ );
	directoryName	= MySQLPPSDDBase::getCellValue ( row, r++ );
	isInactive		= MySQLPPSDDBase::getCellValue ( row, r++ );
	maxMSMSSpectra	= strtoul ( MySQLPPSDDBase::getCellValue ( row, r++ ).c_str (), 0, 10 );
	maxSearchTime	= strtoul ( MySQLPPSDDBase::getCellValue ( row, r++ ).c_str (), 0, 10 );
	allowRaw		= MySQLPPSDDBase::getCellValue ( row, r++ );
}
ProjectInfo::ProjectInfo ( const MYSQL_ROW& row )
{
	int r = 0;
	projectID		= MySQLPPSDDBase::getCellValue ( row, r++ );
	projectName		= MySQLPPSDDBase::getCellValue ( row, r++ );
	projectFile		= MySQLPPSDDBase::getCellValue ( row, r++ );
	projectPath		= MySQLPPSDDBase::getCellValue ( row, r++ );
	calibrationIndex= atoi ( MySQLPPSDDBase::getCellValue ( row, r++ ).c_str () );
	projectDirectory= InfoParams::instance ().getUserRepository () + SLASH + projectPath;
	projectFullPath = projectDirectory + SLASH + projectFile;
}
// Gets a list of projects for a given user.
StringVector MySQLPPSDDBase::getProjectList ( const string& user, bool listNoResults, const ProjectDateFilter* pdf, bool compressInfo ) const
{
//pidLogOutput ( "getProjectListA" );
	StringVector projects;
	string sql;
	if ( user == "root" ) {
		sql = "SELECT projects.project_name, pp_users.user_name";
		if ( compressInfo ) sql += ", projects.project_path, projects.project_file";
	}
	else {
		sql = "SELECT project_name";
		if ( compressInfo ) sql += ", project_path, project_file";
	}
	sql += " FROM projects, pp_users WHERE ";
	if ( user == "root" )	sql += "projects.pp_user_id = pp_users.pp_user_id";
	else					sql += "pp_users.user_name = '" + user + QUOTE + "AND projects.pp_user_id = pp_users.pp_user_id";

	if ( pdf ) {
		string show = pdf->getDateFilter ();
		if ( show == "Projects Created Between" )	sql += " AND DATE_FORMAT(projects.record_created,'%Y-%m') BETWEEN '";
		if ( show == "Projects Accessed Between" )	sql += " AND DATE_FORMAT(projects.record_updated,'%Y-%m') BETWEEN '";
		if ( show == "Projects Created Between" || show == "Projects Accessed Between" ) {
			sql += pdf->getStartYear () + "-" + pdf->getStartMonth () + "' AND '" + pdf->getEndYear () + "-" + pdf->getEndMonth () + "'";
		}
	}
	if ( !mysql_query ( conn, sql.c_str () ) ) {
		MYSQL_RES* res = mysql_store_result ( conn );
		if ( res == NULL )	return projects;				// Error condition
		else {
			MYSQL_ROW row;
			while ( row = mysql_fetch_row ( res ) ) {
				int r = 0;
				string p = MySQLPPSDDBase::getCellValue ( row, r++ );
				string u;
				if ( user == "root" )	u = MySQLPPSDDBase::getCellValue ( row, r++ );
				bool compressFlag = false;
				if ( compressInfo ) {
					string path = MySQLPPSDDBase::getCellValue ( row, r++ );
					string file = MySQLPPSDDBase::getCellValue ( row, r++ );
					compressFlag = genFileExists ( InfoParams::instance ().getUserRepository () + SLASH + path + SLASH + file + ".7z" );
				}
				string pName = u.empty () ? p : u + "/" + p;
				if ( compressFlag ) pName += "$";
				if ( listNoResults ) projects.push_back ( pName );
				else {
					StringVector sv = getResultsList ( u.empty () ? user : u, p );
					if ( !sv.empty () ) projects.push_back ( pName );
				}
			}
			mysql_free_result ( res );
		}
	}
	return projects;
}
// Gets a list of completed results for a given user and project.
bool MySQLPPSDDBase::projectHasQueuedOrRunningJobs ( const string& user, const string& project ) const
{
	//pidLogOutput ( "projectHasQueuedOrRunningJobs" );
	StringVector results;
	string sql = "SELECT\
					results_name\
				FROM search_jobs, projects, pp_users\
				WHERE\
					projects.project_name = '" + project + QUOTE + "\
				AND pp_users.user_name = '" + user + QUOTE + "\
				AND search_jobs.project_id = projects.project_id\
				AND ( search_jobs.job_status = '" + gen_itoa ( JOB_STATUS_SEARCHING ) + "' OR search_jobs.job_status = '" + gen_itoa ( JOB_STATUS_START ) + "' OR search_jobs.job_status = '" + gen_itoa ( JOB_STATUS_SUBMITTED ) + "' )\
				AND projects.pp_user_id = pp_users.pp_user_id";

	if ( !mysql_query ( conn, sql.c_str () ) ) {
		MYSQL_RES* res = mysql_store_result ( conn );
		if ( res == NULL ) {
			return false;				// Error condition
		}
		else {
			MYSQL_ROW row;
			while ( row = mysql_fetch_row ( res ) ) {
				results.push_back ( MySQLPPSDDBase::getCellValue ( row, 0 ) );
			}
			mysql_free_result ( res );
		}
	}
	return !results.empty ();
}
// Gets a list of completed results for a given project ID.
StringVector MySQLPPSDDBase::getResultsList ( const string& projectID ) const
{
	//pidLogOutput ( "getResultsListA" );
	StringVector results;
	string sql = "SELECT\
					results_name\
				FROM search_jobs\
				WHERE\
					search_jobs.project_id = '" + projectID + QUOTE;

	if ( !mysql_query ( conn, sql.c_str () ) ) {
		MYSQL_RES* res = mysql_store_result ( conn );
		if ( res == NULL ) {
			return results;				// Error condition
		}
		else {
			MYSQL_ROW row;
			while ( row = mysql_fetch_row ( res ) ) {
				results.push_back ( MySQLPPSDDBase::getCellValue ( row, 0 ) );
			}
			mysql_free_result ( res );
		}
	}
	return results;
}
// Gets a list of completed results for a given user and project.
StringVector MySQLPPSDDBase::getResultsList ( const string& user, const string& project ) const
{
	//pidLogOutput ( "getResultsListB" );
	StringVector results;
	string sql = "SELECT\
					results_name\
				FROM search_jobs, projects, pp_users\
				WHERE\
					projects.project_name = '" + project + QUOTE + "\
				AND pp_users.user_name = '" + user + QUOTE + "\
				AND search_jobs.project_id = projects.project_id\
				AND search_jobs.job_status = '" + gen_itoa ( JOB_STATUS_DONE ) + "'\
				AND projects.pp_user_id = pp_users.pp_user_id";

	if ( !mysql_query ( conn, sql.c_str () ) ) {
		MYSQL_RES* res = mysql_store_result ( conn );
		if ( res == NULL ) {
			return results;				// Error condition
		}
		else {
			MYSQL_ROW row;
			while ( row = mysql_fetch_row ( res ) ) {
				results.push_back ( MySQLPPSDDBase::getCellValue ( row, 0 ) );
			}
			mysql_free_result ( res );
		}
	}
	return results;
}
string MySQLPPSDDBase::getFullPath ( MYSQL* conn, const string& sql )
{
	//pidLogOutput ( "getFullPath" );
	string fpath;
	if ( !mysql_query ( conn, sql.c_str () ) ) {
		MYSQL_RES* res = mysql_store_result ( conn );
		if ( res == NULL ) {
			return fpath;				// Error condition
		}
		else {
			MYSQL_ROW row;
			while ( row = mysql_fetch_row ( res ) ) {
				int r = 0;
				string path = MySQLPPSDDBase::getCellValue ( row, r++ );
				string file = MySQLPPSDDBase::getCellValue ( row, r++ );
				fpath = InfoParams::instance ().getUserRepository () + SLASH + path + SLASH + file;
			}
			mysql_free_result ( res );
		}
	}
	return fpath;
}
string MySQLPPSDDBase::getResultsFullPath ( const string& user, const string& project, const string& results ) const
{
	//pidLogOutput ( "getResultsFullPathA" );
	string sql = "SELECT\
					results_path,\
					results_file\
				FROM search_jobs, projects, pp_users\
				WHERE\
					projects.project_name = '" + project + QUOTE + "\
				AND pp_users.user_name = '" + user + QUOTE + "\
				AND search_jobs.results_name = '" + results + QUOTE + "\
				AND search_jobs.project_id = projects.project_id\
				AND projects.pp_user_id = pp_users.pp_user_id";

	return MySQLPPSDDBase::getFullPath ( conn, sql );
}
string MySQLPPSDDBase::getResultsFullPath ( const string& searchKey ) const
{
	//pidLogOutput ( "getResultsFullPathB" );
	string sql = "SELECT\
					results_path,\
					results_file\
				FROM search_jobs\
				WHERE\
					search_job_key = '" + searchKey + QUOTE;

	return MySQLPPSDDBase::getFullPath ( conn, sql );
}
string MySQLPPSDDBase::getProjectFullPath ( const string& user, const string& project ) const
{
	//pidLogOutput ( "getProjectFullPath" );
	string sql = "SELECT\
					project_path,\
					project_file\
				FROM projects, pp_users\
				WHERE\
					projects.project_name = '" + project + QUOTE + "\
				AND pp_users.user_name = '" + user + QUOTE + "\
				AND projects.pp_user_id = pp_users.pp_user_id";

	return MySQLPPSDDBase::getFullPath ( conn, sql );
}
string MySQLPPSDDBase::getInstrument ( const string& user, const string& project ) const
{
	//pidLogOutput ( "getInstrument" );
	string sql = "SELECT\
					instrument\
				FROM projects, pp_users\
				WHERE\
					projects.project_name = '" + project + QUOTE + "\
				AND pp_users.user_name = '" + user + QUOTE + "\
				AND projects.pp_user_id = pp_users.pp_user_id";
	return MySQLPPSDDBase::getSingleString ( conn, sql );
}
string MySQLPPSDDBase::getErrorMessage ( const string& searchJobKey ) const
{
	//pidLogOutput ( "getErrorMessage" );
	string sql = "SELECT\
					error_message\
				FROM aborted_jobs\
				WHERE\
					search_job_key = '" + searchJobKey + QUOTE;
	return MySQLPPSDDBase::getSingleString ( conn, sql );
}
// Updates the project record updated
void MySQLPPSDDBase::updateProjectRecordUpdated ( const string& projectID )
{
	//pidLogOutput ( "updateProjectRecordUpdated" );
	PairStringString setClause = make_pair ( string("record_updated"), string("NOW()") );
	PairStringString whereClause = make_pair ( string("project_id"), projectID );
	updateQuery ( "projects", setClause, whereClause );
}
// Updates the project calibration index
void MySQLPPSDDBase::updateCalIndex ( const string& projectID, int calIndex )
{
	//pidLogOutput ( "updateCalIndex" );
	PairStringString setClause = make_pair ( string("calibration_index"), gen_itoa ( calIndex ) );
	PairStringString whereClause = make_pair ( string("project_id"), projectID );
	updateQuery ( "projects", setClause, whereClause );
}
unsigned int MySQLPPSDDBase::deleteSearchKey ( const string& searchKey ) const
{
	//pidLogOutput ( "deleteSearchKey" );
	return deleteQuery ( "search_jobs", "search_job_key", searchKey );
}
unsigned int MySQLPPSDDBase::deleteSearchID ( const string& searchID ) const
{
	//pidLogOutput ( "deleteSearchID" );
	return deleteQuery ( "search_jobs", "search_job_id", searchID );
}
unsigned int MySQLPPSDDBase::deleteProjectID ( const string& projectID ) const
{
	//pidLogOutput ( "deleteProjectID" );
	return deleteQuery ( "projects", "project_id", projectID );
}
unsigned int MySQLPPSDDBase::deleteUserID ( const string& userID ) const
{
	//pidLogOutput ( "deleteUserID" );
	return deleteQuery ( "pp_users", "pp_user_id", userID );
}
unsigned int MySQLPPSDDBase::deleteUserName ( const string& userName ) const
{
	//pidLogOutput ( "deleteUserName" );
	return deleteQuery ( "pp_users", "user_name", userName );
}
unsigned int MySQLPPSDDBase::deleteSession ( const string& sessionKey ) const
{
	//pidLogOutput ( "deleteSession" );
	return deleteQuery ( "sessions", "session_key", sessionKey );
}
unsigned int MySQLPPSDDBase::deleteOldSessions ( int maxDays ) const
{
	//pidLogOutput ( "deleteOldSessions" );
	return submitSQL ( "DELETE FROM sessions WHERE DATEDIFF(NOW(), record_created) > " + gen_itoa (  maxDays ) );
}
unsigned int MySQLPPSDDBase::deleteOldAbortedJobs ( int maxDays ) const
{
	//pidLogOutput ( "deleteOldAbortedJobs" );
	return submitSQL ( "DELETE FROM aborted_jobs WHERE DATEDIFF(NOW(), record_created) > " + gen_itoa (  maxDays ) );
}
void MySQLPPSDDBase::printSQLResultsHTML ( ostream& os, const string& sql ) const
{
	//pidLogOutput ( "printSQLResultsHTML" );
	pair <StringVector, StringVectorVector> p = getSQLResults ( sql );
	printSQLResultsHTML ( os, p );
}
class SortUserInfo {
	int col;
public:
	SortUserInfo ( int col ) :
		col ( col )
	{
	}
	bool operator () ( const StringVector& a, const StringVector& b ) const
	{
		return atof ( a [col].c_str () ) > atof ( b [col].c_str () );
	}
};
pair <StringVector, StringVectorVector> MySQLPPSDDBase::getUserInfoTable () const
{
	string sql = "select user_name, email, directory_name from pp_users";
	pair <StringVector, StringVectorVector> p = getSQLResults ( sql );
	p.first.push_back ( "GB" );
	string baseDir = InfoParams::instance ().getUserRepository ();
	for ( int i = 0 ; i < p.second.size () ; i++ ) {
		string directory = p.second [i][2];
		string fullPath = baseDir + SLASH + directory [0] + SLASH + directory [1] + SLASH + directory;
		double siz = (double) genDirectorySize ( fullPath ) / 1073741824;
		string aSiz = gen_ftoa ( siz, "%.3f" );
		p.second [i].push_back ( aSiz );
	}
	sort ( p.second.begin (), p.second.end (), SortUserInfo ( 3 ) );
	return p;
}
void MySQLPPSDDBase::printUserInfoHTML ( ostream& os ) const
{
	os << "<p>" << endl;
	os << "<b>";
	os << "Total disk space: ";
	os << gen_ftoa ( genVolumeSize ( InfoParams::instance ().getUserRepository () ) / 1073741824, "%.3f" );
	os << " GB";
	os << "</b>" << endl;
	os << "<br/>" << endl;
	os << "<b>";
	os << "Remaining disk space: ";
	os << gen_ftoa ( genFreeDiskSpace ( InfoParams::instance ().getUserRepository () ) / 1073741824, "%.3f" );
	os << " GB";
	os << "</b>" << endl;
	os << endl;
	os << "</p>" << endl;
	printSQLResultsHTML ( os, getUserInfoTable () );
}
void MySQLPPSDDBase::printUserInfoTabDelimitedText ( ostream& os ) const
{
	os << "<pre>" << endl;
	os << "Total disk space: ";
	os << gen_ftoa ( genVolumeSize ( InfoParams::instance ().getUserRepository () ) / 1073741824, "%.3f" );
	os << " GB";
	os << endl;
	os << "Remaining disk space: ";
	os << gen_ftoa ( genFreeDiskSpace ( InfoParams::instance ().getUserRepository () ) / 1073741824, "%.3f" );
	os << " GB";
	os << endl;
	os << "</pre>" << endl;
	printSQLResultsTabDelimitedText ( os, getUserInfoTable () );
}
void MySQLPPSDDBase::printSQLResultsHTML ( ostream& os, const pair <StringVector, StringVectorVector>& p ) const
{
	//pidLogOutput ( "printSQLResultsHTML" );
	if ( p.first.empty () || p.second.empty () ) {
		ErrorHandler::genError ()->message ( "No results to report.\n" );
	}
	else {
		tableStart ( os, true );
			tableRowStart ( os );
				for ( StringVectorSizeType i = 0 ; i < p.first.size () ; i++ ) {
					tableCell ( os, p.first [i], true, true );
				}
			tableRowEnd ( os );
			for ( StringVectorVectorSizeType j = 0 ; j < p.second.size () ; j++ ) {
				tableRowStart ( os );
					for ( StringVectorSizeType k = 0 ; k < p.second [j].size () ; k++ ) {
						tableCell ( os, p.second [j][k], true, false );
					}
				tableRowEnd ( os );
			}
		tableEnd ( os );
	}
}
void MySQLPPSDDBase::printSQLResultsTabDelimitedText ( ostream& os, const string& sql ) const
{
	//pidLogOutput ( "printSQLResultsTabDelimitedText" );
	pair <StringVector, StringVectorVector> p = getSQLResults ( sql );
	printSQLResultsTabDelimitedText ( os, p );
}
void MySQLPPSDDBase::printSQLResultsTabDelimitedText ( ostream& os, const pair <StringVector, StringVectorVector>& p ) const
{
	//pidLogOutput ( "printSQLResultsTabDelimitedText" );
	if ( p.first.empty () || p.second.empty () ) {
		ErrorHandler::genError ()->message ( "No results to report.\n" );
	}
	else {
		os << "<pre>" << endl;
		delimitedRowStart ( os );
			for ( StringVectorSizeType i = 0 ; i < p.first.size () ; i++ ) {
				delimitedCell ( os, p.first [i] );
			}
		delimitedRowEnd ( os );
		for ( StringVectorVectorSizeType j = 0 ; j < p.second.size () ; j++ ) {
			delimitedRowStart ( os );
				for ( StringVectorSizeType k = 0 ; k < p.second [j].size () ; k++ ) {
					delimitedCell ( os, p.second [j][k] );
				}
			delimitedRowEnd ( os );
		}
		os << "</pre>" << endl;
	}
}
void MySQLPPSDDBase::printSQLResultsXML ( ostream& os, const string& sql ) const
{
	//pidLogOutput ( "printSQLResultsXML" );
	pair <StringVector, StringVectorVector> p = getSQLResults ( sql );
	printSQLResultsXML ( os, p );
}
void MySQLPPSDDBase::printSQLResultsXML ( ostream& os, const pair <StringVector, StringVectorVector>& p ) const
{
	//pidLogOutput ( "printSQLResultsXML" );
	for ( StringVectorVectorSizeType i = 0 ; i < p.second.size () ; i++ ) {
		os << "\t<record>" << endl;
		for ( StringVectorSizeType j = 0 ; j < p.second [i].size () ; j++ ) {
			os << "\t\t";
			ParameterList::printXML ( os, p.first [j], p.second [i][j] );
		}
		os << "\t</record>" << endl;
	}
}
pair <StringVector, StringVectorVector> MySQLPPSDDBase::getSQLResults ( const string& sql ) const
{
	//pidLogOutput ( "getSQLResults" );
	StringVectorVector svv;
	StringVector svh;
	if ( !mysql_query ( conn, sql.c_str () ) ) {
		MYSQL_RES* res = mysql_store_result ( conn );
		if ( res == NULL ) {
			throw runtime_error ( "Problem with SQL query.\n" );
		}
		else {
			unsigned int numFields = mysql_num_fields ( res );
			MYSQL_FIELD* fields = mysql_fetch_fields ( res );
			for ( unsigned int i = 0 ; i < numFields ; i++ ) {
				 svh.push_back ( fields [i].name );
			}
			MYSQL_ROW row;
			while ( row = mysql_fetch_row ( res ) ) {
				StringVector sv;
				for ( unsigned int i = 0 ; i < numFields ; i++ ) {
					sv.push_back ( MySQLPPSDDBase::getCellValue ( row, i ) );
				}
				svv.push_back ( sv );
			}
			mysql_free_result ( res );
		}
	}
	return make_pair ( svh, svv );
}

DatabaseProjectImport::DatabaseProjectImport ( const string& dir )
{
	readProjectExportFile ( dir + SLASH + "project.xml" );
	readSearchJobsExportFile ( dir + SLASH + "search_jobs.xml" );
}
DatabaseProjectImport::~DatabaseProjectImport ()
{
	for ( MapStringToDatabaseImportSearchJobInfoPtr::iterator i = msdisji.begin () ; i != msdisji.end () ; i++ ) {
		delete i->second;
	}
}
void DatabaseProjectImport::readProjectExportFile ( const string& file )
{
	if ( !genFileExists ( file ) ) throw runtime_error ( "The file project.xml is missing from the uploaded archive.\n" );
	char* info = getFileInfo ( file, '\n', 1, false );
	instrument		= XMLParser::getStringVectorValue ( info, "instrument" );
	calibrationIndex= XMLParser::getStringVectorValue ( info, "calibration_index" );
	recordCreated	= XMLParser::getStringVectorValue ( info, "record_created" );
	delete [] info;
}
void DatabaseProjectImport::readSearchJobsExportFile ( const string& file )
{
	if ( !genFileExists ( file ) ) throw runtime_error ( "The file search_jobs.xml is missing from the uploaded archive.\n" );
	char* info = getFileInfo ( file, '\n', 1, false );
	StringVector searchJobKey	= XMLParser::getStringVectorValue ( info, "search_job_key" );
	StringVector searchStage	= XMLParser::getStringVectorValue ( info, "search_stage" );
	StringVector searchNumber	= XMLParser::getStringVectorValue ( info, "search_number" );
	StringVector numSerial		= XMLParser::getStringVectorValue ( info, "num_serial" );
	StringVector resultsName	= XMLParser::getStringVectorValue ( info, "results_name" );
	StringVector priority		= XMLParser::getStringVectorValue ( info, "priority" );
	StringVector jobStatus		= XMLParser::getStringVectorValue ( info, "job_status" );
	StringVector jobSignal		= XMLParser::getStringVectorValue ( info, "job_signal" );
	StringVector nodeName		= XMLParser::getStringVectorValue ( info, "node_name" );
	StringVector percentComplete= XMLParser::getStringVectorValue ( info, "percent_complete" );
	StringVector searchSubmitted= XMLParser::getStringVectorValue ( info, "search_submitted" );
	StringVector searchStarted	= XMLParser::getStringVectorValue ( info, "search_started" );
	StringVector searchFinished	= XMLParser::getStringVectorValue ( info, "search_finished" );
	StringVector recordCreated2	= XMLParser::getStringVectorValue ( info, "record_created" );

	for ( StringVectorSizeType i = 0 ; i < searchJobKey.size () ; i++ ) {
		msdisji [searchJobKey [i]] = new DatabaseImportSearchJobInfo ( searchStage [i], searchNumber [i], numSerial [i], resultsName [i], priority [i], jobStatus [i], jobSignal [i], nodeName [i], percentComplete [i], searchSubmitted [i], searchStarted [i], searchFinished [i], recordCreated2 [i] );
	}
	delete [] info;
}
void DatabaseProjectImport::updateSearchKey ( const string& oldSearchKey, const string& newSearchKey )
{
	MapStringToDatabaseImportSearchJobInfoPtr::iterator cur = msdisji.find ( oldSearchKey );
	if ( cur != msdisji.end () ) {
		DatabaseImportSearchJobInfo* disji = cur->second;
		msdisji.erase ( cur );
		msdisji[newSearchKey] = disji;
	}
}
void DatabaseProjectImport::submit ( const string& userID, const string& projName, const string& projFile, const Repository* reposit, const StringVector& resFiles ) const
{
	MySQLPPSDDBase::instance ().submitCopiedProject ( userID, projName, projFile, reposit->getProjectPath (), instrument [0], calibrationIndex [0], recordCreated [0] );
	string projectID = MySQLPPSDDBase::instance ().getProjectID ( userID, projName );
	for ( StringVectorSizeType i = 0 ; i < resFiles.size () ; i++ ) {
		string key = genShortFilenameFromPath ( resFiles [i] );
		MapStringToDatabaseImportSearchJobInfoPtr::const_iterator cur = msdisji.find ( key );
		if ( cur != msdisji.end () ) {
			DatabaseImportSearchJobInfo* disji = cur->second;
			MySQLPPSDDBase::instance ().submitCopiedSearch ( key, projectID, resFiles [i], reposit->getResultsPath (), disji );
		}
	}
}
ProjectDateFilter::ProjectDateFilter ( const string& dateFilter, const string& startYear, const string& endYear, const string& startMonth, const string& endMonth ) :
	dateFilter ( dateFilter ),
	startYear ( startYear ),
	endYear ( endYear ),
	startMonth ( startMonth ),
	endMonth ( endMonth )
{
}
