/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_param_list.cpp                                             *
*                                                                             *
*  Created    : January 31st 2001                                             *
*                                                                             *
*  Purpose    : Functions to read in parameters.                              *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef VIS_C
#pragma warning( disable : 4503 )	// Disable warnings about decorated name length being too long
#endif
#ifndef VIS_C
#include <stdexcept>
#include <lgen_process.h>
#endif
#ifdef VIS_C
#include <process.h>
#else
#include <sys/types.h>
#include <unistd.h>
#endif
#include <algorithm>
#include <lg_stdio.h>
#include <lg_string.h>
#include <lg_time.h>
#include <lgen_error.h>
#include <lgen_file.h>
#ifndef VIS_C
#include <lg_new.h>
#endif
#include <lu_html.h>
#include <lu_param_list.h>
#include <lu_getfil.h>
#include <lu_xml.h>
#include <lu_html_form.h>
#include <lu_pq_vector.h>
#include <lu_prog.h>
#include <lu_cgi_val.h>
#include <lu_mzidentml.h>
#ifdef MYSQL_DATABASE
#include <ld_init.h>
#endif
using std::string;
using std::istream;
using std::ostream;
using std::istringstream;
using std::stringstream;
using std::getline;
using std::endl;
using std::find;
using std::runtime_error;
using std::pair;
using std::sort;
using std::unique;

class CgiEnvironmentVariables {
	string gatewayInterface;
	string serverName;
	string serverSoftware;
	string serverProtocol;
	string serverPort;
	string requestMethod;
	string pathInfo;
	string pathTranslated;
	string scriptName;
	string documentRoot;
	string queryString;
	string remoteHost;
	string remoteAddr;
	string authType;
	string remoteUser;
	string remoteIdent;
	string contentType;
	string contentLength;
	string httpFrom;
	string httpAccept;
	string httpUserAgent;
	string httpReferer;
	string httpCookie;
public:
	CgiEnvironmentVariables ();
	string getGatewayInterface ()	const { return gatewayInterface; }
	string getServerName ()			const { return serverName; }
	string getServerSoftware ()		const { return serverSoftware; }
	string getServerProtocol ()		const { return serverProtocol; }
	string getServerPort ()			const { return serverPort; }
	string getRequestMethod ()		const { return requestMethod; }
	string getPathInfo ()			const { return pathInfo; }
	string getPathTranslated ()		const { return pathTranslated; }
	string getScriptName ()			const { return scriptName; }
	string getDocumentRoot ()		const { return documentRoot; }
	string getQueryString ()		const { return queryString; }
	string getRemoteHost ()			const { return remoteHost; }
	string getRemoteAddr ()			const { return remoteAddr; }
	string getAuthType ()			const { return authType; }
	string getRemoteUser ()			const { return remoteUser; }
	string getRemoteIdent ()		const { return remoteIdent; }
	string getContentType ()		const { return contentType; }
	string getContentLength ()		const { return contentLength; }
	string getHttpFrom ()			const { return httpFrom; }
	string getHttpAccept ()			const { return httpAccept; }
	string getHttpUserAgent ()		const { return httpUserAgent; }
	string getHttpReferer ()		const { return httpReferer; }
	string getHttpCookie ()			const { return httpCookie; }
	string getEnvVariable ( const char* variable );
	void printHTML ( ostream& os ) const;
	void printHTMLEnv ( ostream& os ) const;
};
CgiEnvironmentVariables::CgiEnvironmentVariables () :

	gatewayInterface( getEnvVariable ( "GATEWAY_INTERFACE" ) ),
	serverName		( getEnvVariable ( "SERVER_NAME" ) ),
	serverSoftware	( getEnvVariable ( "SERVER_SOFTWARE" ) ),
	serverProtocol	( getEnvVariable ( "SERVER_PROTOCOL" ) ),
	serverPort		( getEnvVariable ( "SERVER_PORT" ) ),
	requestMethod	( getEnvVariable ( "REQUEST_METHOD" ) ),
	pathInfo		( getEnvVariable ( "PATH_INFO" ) ),
	pathTranslated	( getEnvVariable ( "PATH_TRANSLATED" ) ),
	scriptName		( getEnvVariable ( "SCRIPT_NAME" ) ),
	documentRoot	( getEnvVariable ( "DOCUMENT_ROOT" ) ),
	queryString		( getEnvVariable ( "QUERY_STRING" ) ),
	remoteHost		( getEnvVariable ( "REMOTE_HOST" ) ),
	remoteAddr		( getEnvVariable ( "REMOTE_ADDR" ) ),
	authType		( getEnvVariable ( "AUTH_TYPE" ) ),
	remoteUser		( getEnvVariable ( "REMOTE_USER" ) ),
	remoteIdent		( getEnvVariable ( "REMOTE_IDENT" ) ),
	contentType		( getEnvVariable ( "CONTENT_TYPE" ) ),
	contentLength	( getEnvVariable ( "CONTENT_LENGTH" ) ),
	httpFrom		( getEnvVariable ( "HTTP_FROM" ) ),
	httpAccept		( getEnvVariable ( "HTTP_ACCEPT" ) ),
	httpUserAgent	( getEnvVariable ( "HTTP_USER_AGENT" ) ),
	httpReferer		( getEnvVariable ( "HTTP_REFERER" ) ),
	httpCookie		( getEnvVariable ( "HTTP_COOKIE" ) )
{
}
void CgiEnvironmentVariables::printHTML ( ostream& os ) const
{
	ParameterList::printHTML ( os, "Gateway Interface", gatewayInterface );
	ParameterList::printHTML ( os, "Server Name", serverName );
	ParameterList::printHTML ( os, "Server Software", serverSoftware );
	ParameterList::printHTML ( os, "Server Protocol", serverProtocol );
	ParameterList::printHTML ( os, "Server Port", serverPort );
	ParameterList::printHTML ( os, "Request Method", requestMethod );
	ParameterList::printHTML ( os, "Path Info", pathInfo );
	ParameterList::printHTML ( os, "Path Translated", pathTranslated );
	ParameterList::printHTML ( os, "Script Name", scriptName );
	ParameterList::printHTML ( os, "Document Root", documentRoot );
	ParameterList::printHTML ( os, "Query String", queryString );
	ParameterList::printHTML ( os, "Remote Host", remoteHost );
	ParameterList::printHTML ( os, "Remote Addr", remoteAddr );
	ParameterList::printHTML ( os, "Auth Type", authType );
	ParameterList::printHTML ( os, "Remote User", remoteUser );
	ParameterList::printHTML ( os, "Remote Ident", remoteIdent );
	ParameterList::printHTML ( os, "Content Type", contentType );
	ParameterList::printHTML ( os, "Content Length", contentLength );
	ParameterList::printHTML ( os, "HTTP From", httpFrom );
	ParameterList::printHTML ( os, "HTTP Accept", httpAccept );
	ParameterList::printHTML ( os, "HTTP User Agent", httpUserAgent );
	ParameterList::printHTML ( os, "HTTP Referer", httpReferer );
	ParameterList::printHTML ( os, "HTTP Cookie", httpCookie );
}
void CgiEnvironmentVariables::printHTMLEnv ( ostream& os ) const
{
	ParameterList::printHTML ( os, "GATEWAY_INTERFACE", gatewayInterface );
	ParameterList::printHTML ( os, "SERVER_NAME", serverName );
	ParameterList::printHTML ( os, "SERVER_SOFTWARE", serverSoftware );
	ParameterList::printHTML ( os, "SERVER_PROTOCOL", serverProtocol );
	ParameterList::printHTML ( os, "SERVER_PORT", serverPort );
	ParameterList::printHTML ( os, "REQUEST_METHOD", requestMethod );
	ParameterList::printHTML ( os, "PATH_INFO", pathInfo );
	ParameterList::printHTML ( os, "PATH_TRANSLATED", pathTranslated );
	ParameterList::printHTML ( os, "SCRIPT_NAME", scriptName );
	ParameterList::printHTML ( os, "DOCUMENT_ROOT", documentRoot );
	ParameterList::printHTML ( os, "QUERY_STRING", queryString );
	ParameterList::printHTML ( os, "REMOTE_HOST", remoteHost );
	ParameterList::printHTML ( os, "REMOTE_ADDR", remoteAddr );
	ParameterList::printHTML ( os, "AUTH_TYPE", authType );
	ParameterList::printHTML ( os, "REMOTE_USER", remoteUser );
	ParameterList::printHTML ( os, "REMOTE_IDENT", remoteIdent );
	ParameterList::printHTML ( os, "CONTENT_TYPE", contentType );
	ParameterList::printHTML ( os, "CONTENT_LENGTH", contentLength );
	ParameterList::printHTML ( os, "HTTP_FROM", httpFrom );
	ParameterList::printHTML ( os, "HTTP_ACCEPT", httpAccept );
	ParameterList::printHTML ( os, "HTTP_USER_AGENT", httpUserAgent );
	ParameterList::printHTML ( os, "HTTP_REFERER", httpReferer );
	ParameterList::printHTML ( os, "HTTP_COOKIE", httpCookie );
}
string CgiEnvironmentVariables::getEnvVariable ( const char* variable )
{
	char* value = getenv ( variable );
	if ( value == NULL ) return ( "" );
	else return value;
}
string ParameterList::virtualDir = "/";
string ParameterList::serverName = "";
string ParameterList::serverPort = "";
string ParameterList::server = "";
string ParameterList::remoteAddr = "";
string ParameterList::method = "";
MapStringToString ParameterList::cookie;
const size_t ParameterList::BLOCK_SIZE = 1048576;	// 1 MByte output block

ParameterList::ParameterList ( int argc, char** argv ) :
	programBinaryName ( genFilenameFromPath ( argv [0] ) )
{
	CgiEnvironmentVariables cgiEnvVar;
	remoteAddr = cgiEnvVar.getRemoteAddr ();
	string scriptName = cgiEnvVar.getScriptName ();
	virtualDir = scriptName.substr ( 0, scriptName.find ( BinDir::instance ().getBinDir () ) );
	serverName = cgiEnvVar.getServerName ();
	serverPort = cgiEnvVar.getServerPort ();
	if ( cgiEnvVar.getContentType ().find ( "multipart/form-data" ) != string::npos )	method = "multipart";
	else if ( cgiEnvVar.getQueryString ().length () )									method = "get";
	else																				method = "post";
	server = string ( "http://" ) + serverName;
	if ( serverPort != "80" ) server += string ( ":" ) + serverPort;
	if ( argc != 1 && string ( argv [1] ) == "-" ) {
		string queryString = "";
		for ( int i = 2 ; i < argc ; i++ ) {
			queryString += argv [i];
			if ( i != argc - 1 ) queryString += "&";
		}
		parseQueryString ( queryString );
	}
	else if ( argc != 1 && string ( argv [1] ) == "-f" ) {
		if ( argc > 2 ) {
			parseXMLStringFromFile ( argv [2] );
			SetString erasedParams;
			for ( int i = 3 ; i < argc ; i++ ) {	// Get additional command line parameters
				string nameValue = argv [i];
				int len = nameValue.length ();
				for ( int j = 0 ; j < len ; j++ ) {
					if ( nameValue [j] == '=' ) {
						string name = nameValue.substr ( 0, j );
						string value = nameValue.substr ( j+1, len-j );
						pair <SetStringIterator, bool> flag = erasedParams.insert ( name );
						if ( flag.second ) {		// This is a new parameter erase the old one
							params.erase ( name );
						}
						params [name].push_back ( convertValue ( value ) );
						break;
					}
				}
			}
		}
	}
#ifdef BATCHTAG
#ifdef MYSQL_DATABASE
	else if ( argc != 1 && string ( argv [1] ) == "-k" ) {	// Get the parameters from a database
		string searchKey = argv [2];
		BatchJobItem* bji;
		try {
			bji = MySQLPPSDDBase::instance ().getBatchJobByKey ( searchKey );
		}
		catch ( runtime_error e ) {		// Catch database login problems
			exit ( 1 );
		}
		if ( !bji ) {
			ErrorHandler::genError ()->error ( "Unknown search key.\n" );
		}
		parseXMLStringFromFile ( bji->getResultsFullPath () );
		setValue ( "search_key", convertValue ( searchKey ) );
	}
#endif
	else if ( argc != 1 && string ( argv [1] ) == "-c" ) {	// Get the parameters from a command line file
		parseXMLStringFromFile ( argv [2] );
	}
#endif
	else {
		string queryString = cgiEnvVar.getQueryString ();
		string contentType	= cgiEnvVar.getContentType ();
		string contentLength = cgiEnvVar.getContentLength ();
		string serverProtocol = cgiEnvVar.getServerProtocol ();
		parseCookie ( cgiEnvVar.getHttpCookie () );

		if ( contentType.find ( "multipart/form-data" ) != string::npos ) {
			parseMultipartQueryEntry ( contentType, contentLength, serverProtocol );
		}
		else if ( queryString.length () ) {			// Method is get
			parseQueryString ( queryString );
		}
		else if ( contentLength.length () ) {		// Method is post
			int cl = atoi ( contentLength.c_str () );
			gen_stdin_binary_mode ();
			for ( int i = 0 ; i < cl ; i++ ) {
				char c = fgetc ( stdin );
				if ( c != '\r' ) queryString += c;
			}
			parseQueryString ( queryString );
		}
	}
	writeLogStart ( cgiEnvVar );
}
ParameterList::ParameterList ( const string& filename, bool paramFile )
{
	GenIFStream fromFile ( filename );
	string line;
	while ( getline ( fromFile, line ) ) {
		if ( line.length () != 0 ) {
			parseLine ( line );
		}
	}
}
ParameterList::ParameterList ( istream& is )
{
	string line;
	while ( getline ( is, line ) ) {
		if ( line.length () != 0 ) {
			if ( line [0] == '>' ) break;
			parseLine ( line );
		}
	}
}
ParameterList::ParameterList ( const string& filename )
{
	GenIFStream fromFile ( filename );
	string queryString;
	fromFile >> queryString;
	parseQueryString ( queryString );
}
ParameterList::ParameterList ( const string& str, bool flag1, bool flag2 )
{
	parseQueryString ( str );
}
ParameterList::ParameterList ( const string& str, bool flag1, bool flag2, bool flag3 )
{
	parseXMLString ( str );
}
ParameterList::ParameterList ( const string& filename, bool flag1, bool flag2, bool flag3, bool flag4 )
{
	GenIFStream ist ( filename );
	string line;
	while ( getline ( ist, line ) ) {
		if ( line.find ( "<parameters>" ) != string::npos ) break;
	}
	string info;
	while ( getline ( ist, line ) ) {
		if ( line.find ( "</parameters>" ) != string::npos ) break;
		info += line;
	}
	parseXMLString ( info );
}
ParameterList::ParameterList ( const string& programBinaryName, bool flag1, bool flag2, bool flag3, bool flag4, bool flag5 ) :
	programBinaryName ( programBinaryName )
{
}
ParameterList::~ParameterList ()
{
	static bool endWritten = false;
	if ( !logFileName.empty () && !endWritten ) {
		writeLogEnd ();
		endWritten = true;
	}
}
#ifndef VIS_C
namespace {						// This catches a term signal and makes sure the log file is closed.
ParameterList* termPList = 0;
void paramSigtermHandler ( int sigNum )
{
	termPList->writeLogError ( "Terminate signal sent to program." );
	termPList->writeLogEnd ();
}
}
void initParamSigtermHandler ( ParameterList* pList )
{
	termPList = pList;
	genInitSigterm ( paramSigtermHandler );
}
#endif
void ParameterList::writeLogStart ( const CgiEnvironmentVariables& cgiEnvVar )
{
	string sName = cgiEnvVar.getScriptName ();
	int pNameStart = sName.rfind ( '/' ) + 1;
	int pNameEnd = sName.rfind ( '.' );
	string pName = sName.substr ( pNameStart, pNameEnd - pNameStart );
	if ( InfoParams::instance ().getBoolValue ( pName + "_logging" ) ) {
		string baseDir = ppBaseDir () + "logs";
		int deleteLogDays = InfoParams::instance ().getIntValue ( "delete_log_days" );
		if ( deleteLogDays > 0 ) {
			FileList fList ( baseDir );
			StringVector names = fList.getNameList ();		// Get a list of subdirectories
			for ( int i = 0 ; i < names.size () ; i++ ) {
				if ( genMoreThanNDaysOld ( names [i], deleteLogDays ) ) {
					genUnlinkDirectory ( baseDir + SLASH + names [i] );
				}
			}
		}
		string directory = baseDir + SLASH + genCurrentYearIntMonthAndDay ( '_' );
		string pidStr = gen_itoa ( getpid () );
		if ( !genCreateDirectory ( directory ) ) {
			ErrorHandler::genError ()->error ( "Problems creating a directory for logging.\n" );
		}
		logFileName = directory + SLASH + pName + "_" + genCurrentHoursMinSec ( 0 ) + "_" + pidStr;
		elapsedTime = new GenElapsedTime;
		GenOFStream ost ( logFileName + ".txt", std::ios_base::out );
		printXMLHeader ( ost );
		printXMLVersion ( ost );
		ost << "<program_log>" << endl;
		printXML ( ost, "pid", pidStr );
		printXML ( ost, "start_time", genCurrentTimeAndDateString () );
		printXML ( ost, "SCRIPT_NAME", sName );
		printXML ( ost, "REMOTE_HOST", cgiEnvVar.getRemoteHost () );
		printXML ( ost, "REMOTE_ADDR", cgiEnvVar.getRemoteAddr () );
		printXML ( ost, "HTTP_USER_AGENT", cgiEnvVar.getHttpUserAgent () );
		printXML ( ost, "HTTP_REFERER", cgiEnvVar.getHttpReferer () );
		if ( InfoParams::instance ().getBoolValue ( pName + "_parameter_logging" ) ) {
			XMLParameters ( ost );
		}
	}
#ifndef VIS_C
	initParamSigtermHandler ( this );
#endif
}
void ParameterList::writeLogError ( const string& err ) const
{
	if ( !logFileName.empty () ) {
		if ( genFileExists ( logFileName + ".txt" ) ) {
			GenOFStream ost ( logFileName + ".txt", std::ios_base::out | std::ios_base::app );
			printXML ( ost, "error_message", gen_strtrim2 ( err ) );
		}
	}
}
void ParameterList::writeLogEnd () const
{
	GenOFStream ost ( logFileName + ".txt", std::ios_base::out | std::ios_base::app );
	printXML ( ost, "end_time", genCurrentTimeAndDateString () );
	printXML ( ost, "search_time", elapsedTime->getElapsedTimeString () );
	ost << "</program_log>" << endl;
	ost.close ();
	genRename ( logFileName + ".txt", logFileName + ".xml" );
	delete elapsedTime;
}
void ParameterList::parseCookie ( const string& c )
{
	int start = 0;
	for ( ; ; ) {
		int end = c.find ( "=", start );
		if ( end == string::npos ) break;
		string n = c.substr ( start, end - start );
		start = end+1;
		end = c.find ( ";", start );
		if ( end == string::npos ) {
			cookie [n] = c.substr ( start );
			break;
		}
		else cookie [n] = c.substr ( start, end - start );
		start = end + 2;
	}
}
void ParameterList::parseLine ( const string& line )
{
	if ( line [0] != '#' ) {
		istringstream ist ( line );
		string name;
		ist >> name;
		string temp;
		getline ( ist, temp );
		string value = temp.substr ( temp.find_first_not_of ( " \t" ) );	// strip leading space
		params [name].push_back ( value );
	}
}
string ParameterList::getFileContents ( const string& name, const string& defaultVal ) const
{
	MapStringToStringConstIterator cur = fileContents.find ( name );
	if ( cur != fileContents.end () )
		return (*cur).second;
	else
		return defaultVal;
}
StringVector ParameterList::getNameList ()
{
	StringVector sv;
	for ( MapStringToStringVectorConstIterator i = params.begin () ; i != params.end () ; i++ ) {
		sv.push_back ( (*i).first );
	}
	return sv;
}
bool ParameterList::getValue ( const string& name, bool& value ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () ) {
		value = (*cur).second [0] == "1";
		return true;
	}
	else return false;
}
bool ParameterList::getValue ( const string& name1, const string& name2, bool& value ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name1 );
	if ( cur != params.end () ) {
		StringVector sv = (*cur).second;
		StringVectorIterator svi = find ( sv.begin (), sv.end (), name2 );
		if ( svi != sv.end () ) {
			value = true;
			return true;
		}
		return false;
	}
	else return false;
}
bool ParameterList::getValue ( const string& name, char& value ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () ) {
		value = (*cur).second [0][0];
		return true;
	}
	else return false;
}
bool ParameterList::getValue ( const string& name, int& value ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () ) {
		value = atoi ( (*cur).second [0].c_str () );
		return true;
	}
	else return false;
}
bool ParameterList::getValue ( const string& name, double& value ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () ) {
		value = atof ( (*cur).second [0].c_str () );
		return true;
	}
	else return false;
}
bool ParameterList::getValue ( const string& name, string& value ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () ) {
		value = (*cur).second [0];
		return true;
	}
	else return false;
}
bool ParameterList::getValue ( const string& name, const char*& value ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () ) {
#ifndef VIS_C
		value = gen_new_string ( (*cur).second[0].c_str () );
#else
		value = (*cur).second [0].c_str ();
#endif
		return true;
	}
	else return false;
}
bool ParameterList::getValue ( const string& name, CharVector& value ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () ) {
		StringVector sv = (*cur).second;
		value.clear ();
		for ( StringVectorSizeType i = 0 ; i < sv.size () ; i++ )
			value.push_back ( sv [i][0] );
		return true;
	}
	else return false;
}
bool ParameterList::getValue ( const string& name, StringVector& value ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () ) {
		value = (*cur).second;
		return true;
	}
	else return false;
}
bool ParameterList::getBoolValue ( const string& name, bool defaultVal ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () )
		return (*cur).second [0] == "1";
	else
		return defaultVal;
}
int ParameterList::getIntValue ( const string& name, int defaultVal ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () )
		return atoi ( (*cur).second [0].c_str () );
	else
		return defaultVal;
}
unsigned int ParameterList::getUIntValue ( const string& name, unsigned int defaultVal ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () )
		return strtoul ( (*cur).second [0].c_str (), 0, 10 );
	else
		return defaultVal;
}
double ParameterList::getDoubleValue ( const string& name, double defaultVal ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () )
		return atof ( (*cur).second [0].c_str () );
	else
		return defaultVal;
}
char ParameterList::getCharValue ( const string& name, char defaultVal ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () )
		return (*cur).second [0][0];
	else
		return defaultVal;
}
string ParameterList::getFileStringValue ( const string& name, const string& defaultVal ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () ) {
#ifdef VIS_C
		return (*cur).second [0];
#else
		const string& s = (*cur).second [0];
		string t;
		for ( string::size_type i = 0 ; i < s.length () ; i++ ) {
			if ( s [i] != '\r' ) t += s [i];
		}
		return t;
#endif
	}
	else
		return defaultVal;
}
string ParameterList::getStringValue ( const string& name, const string& defaultVal ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () )
		return (*cur).second [0];
	else
		return defaultVal;
}
string ParameterList::getStrippedStringValue ( const string& name, const string& defaultVal ) const
{
	return gen_strtrim ( getStringValue ( name, defaultVal )  );
}
IntVector ParameterList::getIntVectorValue ( const string& name ) const
{
	IntVector iv;
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () ) {
		StringVector sv = (*cur).second;
		for ( StringVectorSizeType i = 0 ; i < sv.size () ; i++ ) {
			iv.push_back ( atoi ( sv [i].c_str () ) );
		}
	}
	return iv;
}
StringVector ParameterList::getUniqueStringVectorValue ( const string& name ) const
{
	StringVector defaultVal;
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () ) {
		StringVector sv = (*cur).second;
		sort ( sv.begin (), sv.end (), genStrcasecmpAscending () );
		sv.erase ( unique ( sv.begin (), sv.end () ), sv.end () );
		return sv;
	}
	else
		return defaultVal;
}
StringVector ParameterList::getStringVectorValue ( const string& name ) const
{
	StringVector defaultVal;
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () )
		return (*cur).second;
	else
		return defaultVal;
}
StringVector ParameterList::getPQStringVectorValue ( const string& name ) const
{
	StringVector value;
	const char* temp;
	if ( getValue ( name, temp ) ) {
		getPostQueryVector ( temp, value, '\n' );
	}
	return value;
}
BoolDeque ParameterList::getBoolDequeValue ( const string& name, int size ) const
{
	BoolDeque val;
	for ( int i = 1 ; i <= size ; i++ ) {
		stringstream sstr;
		sstr << i;
		MapStringToStringVectorConstIterator cur = params.find ( name + sstr.str () );
		if ( cur != params.end () )
			val.push_back ( (*cur).second [0] == "1" );
		else
			val.push_back ( false );
	}
	return val;
}
string ParameterList::getDirectoryValue ( const string& name, const string& defaultVal ) const
{
	string val;
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () )
		val = (*cur).second [0];
	else
		val = defaultVal;
	return adjustPPOutputPath ( val );
}
string ParameterList::getCookieValue ( const string& name, const string& defaultVal )
{
	MapStringToStringConstIterator cur = cookie.find ( name );
	if ( cur != cookie.end () )
		return (*cur).second;
	else
		return defaultVal;
}
void ParameterList::removeName ( const string& name )
{
	params.erase ( name );
}
void ParameterList::addName ( const string& name, const string& value )
{
	params [name].push_back ( value );
}
void ParameterList::addOrReplaceName ( const string& name, int value )
{
	params.erase ( name );
	params [name].push_back ( gen_itoa ( value ) );
}
void ParameterList::addOrReplaceName ( const string& name, const string& value )
{
	params.erase ( name );
	params [name].push_back ( value );
}
void ParameterList::addOrReplaceName ( const string& name, const StringVector& value )
{
	params.erase ( name );
	params [name] = value;
}
void ParameterList::setValue ( const string& name, const string& value )
{
	params [name].clear ();
	params [name].push_back ( value );
}
void ParameterList::copyName ( const ParameterList* plist, const string& name )
{
	StringVector sv = plist->getStringVectorValue ( name );
	if ( !sv.empty () ) params [name] = sv;
}
void ParameterList::copyFileContents ( const ParameterList* plist, const string& name )
{
	fileContents [name] = plist->getFileContents ( name );
}
void ParameterList::appendParameters ( const ParameterList* newParams )
{
	for ( MapStringToStringVectorConstIterator i = newParams->params.begin () ; i != newParams->params.end () ; i++ ) {
		params [(*i).first] = (*i).second;
	}
}
void ParameterList::appendParameters ( const ParameterList& newParams )
{
	for ( MapStringToStringVectorConstIterator i = newParams.params.begin () ; i != newParams.params.end () ; i++ ) {
		params [(*i).first] = (*i).second;
	}
}
bool ParameterList::isEqual ( const ParameterList* plist, const string& name ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () )
		return plist->getStringVectorValue ( name ) == (*cur).second;
	else
		return plist->getStringVectorValue ( name ).empty ();
}
void ParameterList::parseQueryString ( const string& queryString )
{
	int length = queryString.length ();
	for ( int i = 0, start = 0 ; i <= length ; i++ ) {
		if ( i == length || queryString [i] == '&' ) {
			string nameValue = queryString.substr ( start, i - start );
			int len = nameValue.length ();
			if ( len ) {
				bool set = false;
				for ( int j = 0 ; j < len ; j++ ) {
					if ( nameValue [j] == '=' ) {
						string name = nameValue.substr ( 0, j );
						string value = nameValue.substr ( j+1, len-j );
						if ( !value.empty () ) {
							params [name].push_back ( convertValue ( value ) );
						}
						set = true;
						break;
					}
				}
				if ( !set ) {		// No equals sign - incorrect format
					ErrorHandler::genError ()->error ( "Program parameters specified in an incorrect format.\n" );
				}
			}
			start = i+1;
		}
	}
}
void ParameterList::parseXMLStringFromFile ( const string& file )
{
	GenIFStream ist ( file );
	string line;
	string info;
	while ( getline ( ist, line ) ) {
		info += line;
		if ( line.find ( "</parameters>" ) != string::npos ) break;
	}
	parseXMLString ( XMLParser::getStringValue ( info, "parameters" ) );
}
void ParameterList::parseXMLString ( const string& xmlString )
{
	int endIndex = 0;
	for ( ; ; ) {
		int startIndex = xmlString.find ( "<", endIndex );
		if ( startIndex == string::npos ) break;
		startIndex += 1;
		endIndex = xmlString.find ( ">", startIndex );
		string name = xmlString.substr ( startIndex, endIndex - startIndex );
		startIndex = endIndex + 1;
		endIndex = xmlString.find ( "<", startIndex );
		string value = xmlString.substr ( startIndex, endIndex - startIndex );
		params [name].push_back ( convertValue ( value ) );
		endIndex += 1;
	}
}
void ParameterList::parseMultipartQueryEntry ( const string& contentType, const string& contentLength, const string& serverProtocol )
{
	unsigned int maxContentLength = InfoParams::instance ().getUIntValue ( "max_upload_content_length" );
	unsigned int pos = contentType.find ( "=" );
	string boundary = contentType.substr ( pos+1, contentType.length () - pos );

	string dataStart = "--";
	dataStart += boundary;

	string dataEnd = "--";
	dataEnd += boundary;
	dataEnd += "--";

	bool http1 = ( serverProtocol == "HTTP/1.0" );
	unsigned int cl = strtoul ( contentLength.c_str (), 0, 10 );
	gen_stdin_binary_mode ();
	if ( maxContentLength && cl > maxContentLength ) {
		if ( maxContentLength == 1 ) {
			writeStdinToLog ( cl );
			ErrorHandler::genError ()->error ( "CGI content written to log file.\n" );
		}
		else {
			for ( unsigned int i = 0 ; i < cl ; i++ ) fgetc ( stdin );
			ErrorHandler::genError ()->error ( "Content length too long.\n" );
		}
	}
	bool nameLine = false;
	string line;
	string name;
	string value;
	bool fileMode = false;
	bool saveSettings = false;
	string emptyFile;
	string sameFilenames;
	SetString uploadFilenames;
	StringVector uploadFilepaths;
	bool uploadsOptional = false;
	for ( ; ; ) {
		char c = fgetc ( stdin );
		if ( c == '\n' ) {
			if ( line == dataEnd ) {
				params [name].push_back ( value );
				break;
			}
			if ( line == dataStart ) {
				if ( !value.empty () ) {
					params [name].push_back ( value );
				}
				nameLine = true;
			}
			else if ( nameLine ) {
				string::size_type fstart = line.find ( "filename=" );
				unsigned int start = line.find ( "\"" );
				unsigned int len = line.length ();
				if ( fstart != string::npos ) {
					fileMode = true;
					name = line.substr ( start+1, fstart - start - 4 );
					if ( !http1 ) while ( fgetc ( stdin ) != '\n' );	// Next line is the file type
					string filename = line.substr ( fstart+10, len - fstart - 11 );
					while ( fgetc ( stdin ) != '\n' );	// Next line blank
					if ( !filename.empty () ) {
						pair <SetStringIterator, bool> flag = uploadFilenames.insert ( filename );
						if ( !flag.second ) sameFilenames = filename;
						bool tempFlag = name.find ( "temp" ) != string::npos;
						string filepath = getTempFileFullPath ( filename, tempFlag );
						bool flag2 = writeUploadFile ( filepath, dataStart );
						if ( flag2 ) {
							uploadFilepaths.push_back ( filepath );
							params [name + "_filename"].push_back ( filename );
							params [name + "_filepath"].push_back ( filepath );
						}
						else emptyFile = filename;
						bool end = false;
						for ( ; ; ) {		// Skip to end of dataStart or dataEnd
							c = fgetc ( stdin );
							if ( c == '-' ) end = true;
							if ( c == '\n' ) break;
						}
						if ( end ) break;
						else {
							line = "";
							continue;
						}
					}
				}
				else
					name = line.substr ( start+1, len - start - 2 );
				if ( name == "save_params" ) saveSettings = true;
				if ( name == "uploads_optional" ) uploadsOptional = true;
				nameLine = false;
				value = "";
			}
			else {	// value line
				if ( !value.empty () ) {
					value += '\n';
				}
				value += line;
			}
			line = "";
		}
		else if ( c != '\r' ) line += c;
	}
	if ( fileMode && !saveSettings && !uploadsOptional ) {
		if ( uploadFilenames.empty () ) {
			ErrorHandler::genError ()->error ( "Upload file not specified.\n" );
		}
		if ( !emptyFile.empty () || !sameFilenames.empty () ) {
			if ( !uploadFilepaths.empty () ) genUnlink ( uploadFilepaths );
			if ( !emptyFile.empty () ) {
				ErrorHandler::genError ()->error ( "Upload file " + emptyFile + " empty.\n" );
			}
			if ( !sameFilenames.empty () ) {
				ErrorHandler::genError ()->error ( "Upload file " + sameFilenames + " used more than once.\n" );
			}
		}
	}
}
bool ParameterList::writeUploadFile ( const string& filepath, const string& dataStart )
{
	unsigned int dataStartLength = dataStart.length ();
	char* f = new char [BLOCK_SIZE+dataStartLength];
	char* lPtr = f;
	FILE* fp = NULL;
	unsigned int bIdx = 0;
	bool ret = false;
	for ( ; ; ) {
		int c = fgetc ( stdin );
		if ( c == EOF ) {
			delete [] f;
			if ( fp != NULL ) gen_fclose ( fp, "writeUploadFile" );
			genUnlink ( filepath );
			throw runtime_error ( "Unexpected EOF in CGI input." );
		}
		*lPtr++ = c;
		if ( c == dataStart [bIdx] ) {		// Look for the end of the file block
			bIdx++;
			if ( bIdx == dataStartLength ) {
				size_t wLen = lPtr - f - dataStartLength - 2;
				if ( wLen && fp == NULL ) {
					fp = gen_fopen_binary ( filepath, "w", "writeUploadFile" );
				}
				if ( fp != NULL ) {
					if ( wLen ) gen_fwrite ( f, wLen * sizeof (char), 1, fp, "parseMultipartQueryEntry" );
					gen_fclose ( fp, "writeUploadFile" );
					fp = NULL;
					ret = true;					// File written
				}
				else
					ret = false;				// The file is of zero length so not written
				delete [] f;
				break;
			}
		}
		else bIdx = 0;
		if ( lPtr - f >= BLOCK_SIZE && bIdx == 0 ) {
			if ( fp == NULL ) fp = gen_fopen_binary ( filepath, "w", "writeUploadFile" );
			gen_fwrite ( f, ( lPtr - f ) * sizeof (char), 1, fp, "writeUploadFile" );
			lPtr = f;
		}
	}
	return ret;
}
string ParameterList::getTempFileFullPath ( const string& filename, bool tempFlag )
{
	string uploadFileDumpDir = InfoParams::instance ().getUploadTemp ();
	string suffix;
	if ( isFileType ( filename, "tar.gz" ) )
		suffix = "tar.gz";
	else if ( isFileType ( filename, "tar.z" ) )
		suffix = "tar.z";
	else if ( isFileType ( filename, "tar.bz2" ) )
		suffix = "tar.bz2";
	else
		suffix = genSuffixFromPath ( filename );
#ifdef VIS_C
	if ( suffix == "tar" || suffix == "tar.gz" || suffix == "tar.bz2" || suffix == "tar.z" || suffix == "tgz" || suffix == "taz" ) {
		uploadFileDumpDir = "temp";		// Temp directory cannot be on a network drive and work with tar on a Windows platform
	}
#endif
	if ( tempFlag ) uploadFileDumpDir = "temp";
	PPTempFile tempFile ( uploadFileDumpDir, suffix, false );
	return tempFile.getFullPath ();
}
void ParameterList::writeStdinToLog ( unsigned int cl )
{
	string logFileName = ppBaseDir () + "logs" + SLASH + "cgi" + gen_itoa ( getpid () ) + ".txt";
	FILE* fp = gen_fopen_binary ( logFileName, "w", "writeStdinToLog" );
	char* f = new char [BLOCK_SIZE];
	char* lPtr = f;
	int off = 0;
	for ( unsigned int i = 0, j = 0 ; i < cl ; i++ ) {
		*lPtr++ = fgetc ( stdin );
		if ( lPtr - f == BLOCK_SIZE ) {
			gen_fwrite ( f, BLOCK_SIZE * sizeof (char), 1, fp, "writeStdinToLog" );
			lPtr = f;
		}
	}
	size_t rest = lPtr - f;
	if ( rest ) gen_fwrite ( f, rest * sizeof (char), 1, fp, "writeStdinToLog" );
	gen_fclose ( fp, "writeStdinToLog" );
}
string ParameterList::convertValue ( const string& oldValue )
{
	string value;
	for ( StringSizeType i = 0 ; i < oldValue.length () ; i++ ) {
		char first;
		char second;
		char digit;
		switch ( oldValue [i] ) {
			case '%':
				first = oldValue [++i];
				second = oldValue [++i];
				digit = (first >= 'A' ? ((first & 0xdf) - 'A')+10 : (first - '0'));
				digit *= 16;
				digit += (second >= 'A' ? ((second & 0xdf) - 'A')+10 : (second - '0'));
				value += digit;
				break;
			case '+':
				value += ' ';
				break;
			default:
				value += oldValue [i];
				break;
		}
	}
	return ( value );
}
void ParameterList::HTMLParameters ( ostream& os ) const
{
	for ( MapStringToStringVectorConstIterator i = params.begin () ; i != params.end () ; i++ ) {
		StringVector value = (*i).second;
		for ( StringVectorSizeType j = 0 ; j < value.size () ; j++ ) {
			printHTML ( os, (*i).first, value [j] );
		}
	}
}
void ParameterList::XMLParameterFile ( const string& filename ) const
{
	GenOFStream ost ( filename );
	printXMLHeader ( ost );
	printXMLVersion ( ost );
	XMLParameters ( ost );
}
void ParameterList::getParams ( StringVector& sv, StringVectorVector& svv ) const
{
	for ( MapStringToStringVectorConstIterator i = params.begin () ; i != params.end () ; i++ ) {
		sv.push_back ( (*i).first );
		svv.push_back ( (*i).second );
	}
}
void ParameterList::pepXMLParameters ( ostream& os, int ntab ) const
{
	for ( MapStringToStringVectorConstIterator i = params.begin () ; i != params.end () ; i++ ) {
		StringVector value = (*i).second;
		for ( StringVectorSizeType j = 0 ; j < value.size () ; j++ ) {
			for ( size_t k = 0 ; k < ntab ; k++ ) os << '\t';
			os << "<parameter name=\"" << (*i).first << "\" value=\"" << value [j] << "\" />" << endl;
		}
	}
}
VectorXMLOutputItemPtr ParameterList::getMZIdentMLUserParameters () const
{
	VectorXMLOutputItemPtr items;
	for ( MapStringToStringVectorConstIterator i = params.begin () ; i != params.end () ; i++ ) {
		StringVector value = (*i).second;
		for ( StringVectorSizeType j = 0 ; j < value.size () ; j++ ) {
			items.push_back ( new MZIdentML_userParam ( (*i).first, value [j] ) );
		}
	}
	return items;
}
void ParameterList::XMLParameters ( ostream& os ) const
{
	os << "<parameters>" << endl;
	for ( MapStringToStringVectorConstIterator i = params.begin () ; i != params.end () ; i++ ) {
		StringVector value = (*i).second;
		for ( StringVectorSizeType j = 0 ; j < value.size () ; j++ ) {
			printXML ( os, (*i).first, escapeURL ( value [j] ) );
		}
	}
	os << "</parameters>" << endl;
}
void ParameterList::perlFileParameters ( ostream& os ) const
{
	for ( MapStringToStringVectorConstIterator i = params.begin () ; i != params.end () ; i++ ) {
		StringVector value = (*i).second;
		for ( StringVectorSizeType j = 0 ; j < value.size () ; j++ ) {
			os << "$c .= \"";
			os << (*i).first << "=";
			os << escapeURL ( value [j] );
			os << " \";";
			os << endl;
		}
	}
}
void ParameterList::batchFileParameters ( ostream& os ) const
{
	for ( MapStringToStringVectorConstIterator i = params.begin () ; i != params.end () ; i++ ) {
		StringVector value = (*i).second;
		for ( StringVectorSizeType j = 0 ; j < value.size () ; j++ ) {
			os << (*i).first << "=";
			os << batchFileEscapeURL ( escapeURL ( value [j] ) );
			os << " ";
		}
	}
}
string ParameterList::getURL () const
{
	string s = ProgramLink::getURLStart ( programBinaryName );
	s += '?';
	for ( MapStringToStringVectorConstIterator i = params.begin () ; i != params.end () ; i++ ) {
		StringVector value = (*i).second;
		for ( StringVectorSizeType j = 0 ; j < value.size () ; j++ ) {
			s += (*i).first;
			s += "=";
			s += escapeURL ( value [j] );
			s += "&";
		}
	}
	s = s.substr ( 0, s.length () - 1 );	// Delete last &
	return s;
}
bool ParameterList::copyToCookie ( ostream& os, const string& cookieName ) const
{
	string cookieValue;
	for ( MapStringToStringVectorConstIterator i = params.begin () ; i != params.end () ; i++ ) {
		StringVector value = (*i).second;
		for ( StringVectorSizeType j = 0 ; j < value.size () ; j++ ) {
			cookieValue += (*i).first;
			cookieValue += "=";
			cookieValue += escapeURL ( value [j] );
			cookieValue += "&";
		}
	}
	return setCookie ( os, cookieName, cookieValue, true );
}
void ParameterList::copyToCGI ( ostream& os ) const
{
	for ( MapStringToStringVectorConstIterator i = params.begin () ; i != params.end () ; i++ ) {
		StringVector value = (*i).second;
		for ( StringVectorSizeType j = 0 ; j < value.size () ; j++ ) {
			os << (*i).first << "=";
			os << escapeURL ( value [j] );
			os << "&";
		}
	}
}
bool ParameterList::copyToCGI ( ostream& os, const string& name ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () ) {
		StringVector sv = (*cur).second;
		for ( StringVectorSizeType i = 0 ; i < sv.size () ; i++ ) {
			os << name << "=" << escapeURL ( sv [i] ) << "&";
		}
		return true;
	}
	else return false;
}
bool ParameterList::copyToCGI ( ostream& os, const string& name, int num ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () ) {
		StringVector sv = (*cur).second;
		os << name << "=" << escapeURL ( sv [num] ) << "&";
		return true;
	}
	else return false;
}
void ParameterList::copyToHiddenFormEntry ( ostream& os ) const
{
	for ( MapStringToStringVectorConstIterator i = params.begin () ; i != params.end () ; i++ ) {
		StringVector value = (*i).second;
		for ( StringVectorSizeType j = 0 ; j < value.size () ; j++ ) {
			os << "<input type=\"hidden\" name=\"" << (*i).first << "\" value=\"" << value [j] << "\" />" << std::endl;
		}
	}
}
bool ParameterList::copyToHiddenFormEntry ( ostream& os, const string& name ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () ) {
		StringVector sv = (*cur).second;
		for ( StringVectorSizeType i = 0 ; i < sv.size () ; i++ ) {
			os << "<input type=\"hidden\" name=\"" << name << "\" value=\"" << sv [i] << "\" />" << std::endl;
		}
		return true;
	}
	else return false;
}
bool ParameterList::copyToHiddenFormEntry ( ostream& os, const string& name, int num ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () ) {
		StringVector sv = (*cur).second;
		os << "<input type=\"hidden\" name=\"" << name << "\" value=\"" << sv [num] << "\" />" << std::endl;
		return true;
	}
	else return false;
}
string ParameterList::getCommandLineNVPair ( const string& name ) const
{
	string s;
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () ) {
		StringVector sv = (*cur).second;
		for ( StringVectorSizeType i = 0 ; i < sv.size () ; i++ ) {
			s += ::getCommandLineNVPair ( name, sv [i] );
			s += " ";
		}
		return s;
	}
	else return s;
}
bool ParameterList::copyToHiddenFormJavascriptEntry ( ostream& os, const string& name ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () ) {
		StringVector sv = (*cur).second;
		for ( StringVectorSizeType i = 0 ; i < sv.size () ; i++ ) {
			printHTMLFORMJavascriptHidden ( os, name, sv [i] );
		}
		return true;
	}
	else return false;
}
bool ParameterList::copyToHiddenFormJavascriptEntry ( ostream& os, const string& name, int num ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () ) {
		StringVector sv = (*cur).second;
		printHTMLFORMJavascriptHidden ( os, name, sv [num] );
		return true;
	}
	else return false;
}
void ParameterList::printHTMLContainer ( ostream& os, const string& name, const StringVector& value )
{
	os << name << ": ";
	os << "<b>";
	for ( StringVectorSizeType i = 0 ; i < value.size () ; i++ ) {
		os << value [i];
		if ( i < value.size () - 1 ) os << ", ";
	}
	os << "</b><br />" << endl;
}
void ParameterList::printXMLContainer ( ostream& os, const string& name, const IntVector& value )
{
	for ( IntVectorSizeType i = 0 ; i < value.size () ; i++ ) {
		printXML ( os, name, value [i] );
	}
}
void ParameterList::printXMLContainer ( ostream& os, const string& name, const DoubleVector& value )
{
	for ( DoubleVectorSizeType i = 0 ; i < value.size () ; i++ ) {
		printXML ( os, name, value [i] );
	}
}
void ParameterList::printXMLContainer ( ostream& os, const string& name, const StringVector& value )
{
	for ( StringVectorSizeType i = 0 ; i < value.size () ; i++ ) {
		printXML ( os, name, value [i] );
	}
}
string ParameterList::spaceToPlus ( const string& str1 )
{
	string str2 = str1;

	for ( int i = 0 ; str2[i] ; i++ ) if ( str2[i] == ' ' ) str2[i] = '+';

	return ( str2 );
}
string ParameterList::batchFileEscapeURL ( const string& url )
{
	string url2;
	for ( StringSizeType i = 0 ; i < url.length () ; i++ ) {
		char c = url [i];
		switch ( c ) {
			case '%':
				url2 += '%';
				break;
		}
		url2 += c;
	}
	return ( url2 );
}
