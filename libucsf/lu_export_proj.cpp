/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_export_proj.cpp                                            *
*                                                                             *
*  Created    : June 29th 2012                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2012-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef MYSQL_DATABASE
#include <lg_string.h>
#include <lgen_error.h>
#include <lgen_file.h>
#include <lgen_uncompress.h>
#include <ld_init.h>
#include <lu_export_proj.h>
#include <lu_file_type.h>
#include <lu_getfil.h>
#include <lu_html.h>
#include <lu_param_list.h>
#include <lu_proj_file.h>
#include <lu_xml.h>
using std::getline;
using std::string;
using std::cout;
using std::endl;
using namespace FileTypes;

namespace {
void checkCalProject ( const StringVector& projects );
void checkRunningJobs ( const string& user, const StringVector& projects );
StringVector getResultsFiles ( const string& user, const string& project, const string& results );
StringVector getExpFiles ( const string& projectDir, const string& projectName );
StringVector getCentroidAndRawFiles ( const string& dataToExport, const string& projectFullPath, bool uploadOnly = false );
void writeProjectXML ( const string& fullPath, const string& project, const string& user );
void writeSearchJobsXML ( const string& fullPath, const string& projectID );
void copyResultsFiles ( const string& fullPath, const StringVector& resultsFiles );
void copyProjectFiles ( const string& fullPath, const string& projectFile, const StringVector& expFiles );
void copyDataFiles ( const string& fullPath, const StringVector& dataFiles );
void createCompressedFile ( const string& fullPath, const string& project, UpdatingJavascriptMessage* ujm );
void exportResults ( const string& dataToExport, const string& user, const string& project, const string& results = "", UpdatingJavascriptMessage* ujm = 0 );
//void exportResults ( const string& dataToExport, const string& user, const StringVector& results ); // Not currently used.

bool compressZ7 ( const string& f, const UpdatingJavascriptMessage* ujm );
bool compressProjectFiles ( const string& projectFile, const StringVector& expFiles, const UpdatingJavascriptMessage* ujm );
void compressResultsFiles ( const StringVector& resultsFiles, const UpdatingJavascriptMessage* ujm );
void compressDataFiles ( const StringVector& dataFiles, const UpdatingJavascriptMessage* ujm );
void compressResults ( const string& dataToExport, const string& user, const string& project, const string& results, const UpdatingJavascriptMessage* ujm );

void uncompressZ7 ( const string& f, const UpdatingJavascriptMessage* ujm );
void uncompressProjectFiles ( const string& projectFile, const StringVector& expFiles, const UpdatingJavascriptMessage* ujm );
void uncompressResultsFiles ( const StringVector& resultsFiles, const UpdatingJavascriptMessage* ujm );
void uncompressDataFiles ( const StringVector& dataFiles, const UpdatingJavascriptMessage* ujm );
void uncompressResults ( const string& dataToExport, const string& user, const string& project, const string& results, const UpdatingJavascriptMessage* ujm );
}
namespace {
// Export functions
void checkCalProject ( const StringVector& projects )
{
	for ( StringVectorSizeType i = 0 ; i < projects.size () ; i++ ) {
		if ( projects [i].rfind ( ".cal." ) != string::npos ) {
			ErrorHandler::genError ()->error ( "The export/compression of calibrated projects is not currently supported.\n" );
		}
	}
}
void checkRunningJobs ( const string& user, const StringVector& projects )
{
	for ( StringVectorSizeType i = 0 ; i < projects.size () ; i++ ) {
		if ( MySQLPPSDDBase::instance ().projectHasQueuedOrRunningJobs ( user, projects [i] ) ) {
			ErrorHandler::genError ()->error ( "One of the projects selected for export has running jobs.\n" );
		}
	}
}
StringVector getResultsFiles ( const string& user, const string& project, const string& results )
{
//	Get the results files
	StringVector resultsFiles;
	StringVector resultsList;
	if ( !results.empty () )	resultsList.push_back ( results );
	else						resultsList = MySQLPPSDDBase::instance ().getResultsList ( user, project );
	for ( StringVectorSizeType i = 0 ; i < resultsList.size () ; i++ ) {
		string resultsFullPath = MySQLPPSDDBase::instance ().getResultsFullPath ( user, project, resultsList [i] );
		string searchKey = genShortFilenameFromPath ( resultsFullPath );
		string d = genDirectoryFromPath ( resultsFullPath );
		string f = genFilenameFromPath ( resultsFullPath );
		StringVector sv;
		sv.push_back ( searchKey + "*" );
		FileList fileList ( d, sv );
		for ( int j = 0 ; j < fileList.size () ; j++ ) {
			string file = d + SLASH + fileList [j];
			resultsFiles.push_back ( file );
		}
	}
	return resultsFiles;
}
StringVector getExpFiles ( const string& projectDir, const string& projectName )
{
	FileList fileList ( projectDir, projectName + ".exp.", "", false );
	StringVector expFiles;
	for ( int i = 0 ; i < fileList.size () ; i++ ) {
		expFiles.push_back ( projectDir + SLASH + fileList [i] );
	}
	return expFiles;
}
StringVector getCentroidAndRawFiles ( const string& dataToExport, const string& projectFullPath, bool uploadOnly )
{
	ProjectFile pf ( projectFullPath );
	StringVector dataFiles;
	if ( uploadOnly && !pf.isUploadProject () ) return dataFiles;
	if ( dataToExport != "No Data" ) {				// Add centroid files
		for ( StringVectorSizeType i = 0 ; i < pf.getNumFiles () ; i++ ) {
			dataFiles.push_back ( pf.getCentroidPath ( i ) );
		}
	}
	if ( dataToExport == "All Data" ) {				// Add raw files
		for ( StringVectorSizeType j = 0 ; j < pf.getNumFiles () ; j++ ) {
			string path = pf.getAbsoluteRawPath ( j );
			if ( !path.empty () ) {
				dataFiles.push_back ( path );
			}
		}
	}
	return dataFiles;
}
void writeProjectXML ( const string& fullPath, const string& project, const string& user )
{
	string ppUserID = MySQLPPSDDBase::instance ().getUserID ( user );
	GenOFStream ost ( fullPath + SLASH + "project.xml" );
	printXMLHeader ( ost );
	printXMLVersion ( ost );
	ost << "<project>" << endl;
		MySQLPPSDDBase::instance ().printSQLResultsXML ( ost, "select * from projects where project_name = '" + project + "' and pp_user_id = '" + ppUserID + "'" );
	ost << "</project>" << endl;
	ost.close ();
}
void writeSearchJobsXML ( const string& fullPath, const string& projectID )
{
	GenOFStream ost ( fullPath + SLASH + "search_jobs.xml" );
	printXMLHeader ( ost );
	printXMLVersion ( ost );
	ost << "<search_jobs>" << endl;
		MySQLPPSDDBase::instance ().printSQLResultsXML ( ost, "select * from search_jobs where project_id = '" + projectID + "'" );
	ost << "</search_jobs>" << endl;
	ost.close ();
}
void copyResultsFiles ( const string& fullPath, const StringVector& resultsFiles )
{
	string resultsDir = fullPath + SLASH + "results";
	genCreateDirectory ( resultsDir );
	for ( StringVectorSizeType i = 0 ; i < resultsFiles.size () ; i++ ) {
		copyFile ( resultsFiles [i], resultsDir + SLASH + genFilenameFromPath ( resultsFiles [i] ) );
	}
}
void copyProjectFiles ( const string& fullPath, const string& projectFile, const StringVector& expFiles )
{
	string projectDir = fullPath + SLASH + "project";
	genCreateDirectory ( projectDir );
	copyFile ( projectFile, projectDir + SLASH + genFilenameFromPath ( projectFile ) );
	for ( StringVectorSizeType i = 0 ; i < expFiles.size () ; i++ ) {
		copyFile ( expFiles [i], projectDir + SLASH + genFilenameFromPath ( expFiles [i] ) );
	}
}
void copyDataFiles ( const string& fullPath, const StringVector& dataFiles )
{
	if ( !dataFiles.empty () ) {
		string dataDir = fullPath + SLASH + "data";
		genCreateDirectory ( dataDir );
		for ( StringVectorSizeType i = 0 ; i < dataFiles.size () ; i++ ) {
			string path = dataFiles [i];
			string f = genFilenameFromPath ( path );
			if ( genIsDirectory ( path ) ) {
				genCreateDirectory ( dataDir + SLASH + f );
				FileList fList ( path, "", "", false );
				StringVector sv = fList.getNameList ();
				for ( StringVectorSizeType j = 0 ; j < sv.size () ; j++ ) {
					string df2 = sv [j];
					string path2 =  path + SLASH + df2;
					if ( !genIsDirectory ( path2 ) ) {
						copyFile ( path2, dataDir + SLASH + f + SLASH + df2 );
					}
				}
			}
			else
				copyFile ( path, dataDir + SLASH + f );
		}
	}
}
void createCompressedFile ( const string& fullPath, const string& project, UpdatingJavascriptMessage* ujm )
{
	if ( ujm )	ujm->addMessage ( cout, "<b>Starting zip file creation for project " + project + "</b>." );
	PPTempFile pptf2 ( "", "" );
	bool flag = gen7zaCreate ( pptf2.getFullPath () + SLASH + project, fullPath + SLASH + "*", "zip" );
	if ( flag ) {
		cout << "<a href=\"" << pptf2.getURL () + "/" + project + ".zip" << "\">";
		cout << "Archive file (" << project << ".zip)</a><br />" << endl;
	}
	if ( ujm )	ujm->addMessage ( cout, "<b>Finished zip file creation for project " + project + "</b>." );
}
void exportResults ( const string& dataToExport, const string& user, const string& project, const string& results, UpdatingJavascriptMessage* ujm )
{
	ProjectInfo* pi = MySQLPPSDDBase::instance ().getProjectInfo ( user, project );
	string projectPath = pi->getProjectFullPath ();
	bool compressed = false;
	if ( genFileExists ( projectPath + ".7z" ) ) {		// Check that the project isn't compressed
		ujm->startWriteMessage ( cout );
		uncompressResults ( "All Data", user, project, results, ujm );
		compressed = true;
		ujm->endWriteMessage ( cout );
	}
	StringVector resultsFiles = getResultsFiles ( user, project, results );

	StringVector expFiles = getExpFiles ( pi->getProjectDirectory (), pi->getProjectName () );
	StringVector dataFiles = getCentroidAndRawFiles ( dataToExport, projectPath );

	PPTempFile pptf ( "", "" );
	string fullPath = pptf.getFullPath ();
	genCreateDirectory ( fullPath );

	writeProjectXML ( fullPath, project, user );
	writeSearchJobsXML ( fullPath, pi->getProjectID () );
	copyResultsFiles ( fullPath, resultsFiles );
	copyProjectFiles ( fullPath, projectPath, expFiles );
	copyDataFiles ( fullPath, dataFiles );
	if ( compressed ) {
		ujm->startWriteMessage ( cout );
		compressResults ( "All Data", user, project, results, ujm );
		ujm->endWriteMessage ( cout );
	}
	createCompressedFile ( fullPath, project, ujm );
	genUnlinkDirectory ( fullPath );
}
/* Not currently used.
void exportResults ( const string& dataToExport, const string& user, const StringVector& results )
{
	for ( StringVectorSizeType i = 0 ; i < results.size () ; i++ ) {
		StringVector sv = genGetSubstrings ( results [i], '/' );
		if ( sv.size () == 3 )
			exportResults ( dataToExport, sv [0], sv [1], sv [2] );
		else
			exportResults ( dataToExport, user, sv [0], sv [1] );
	}
}
*/

// Compression functions
bool compressZ7 ( const string& f, const UpdatingJavascriptMessage* ujm )
{
	if ( isFileType ( f, Z7 ) ) {
		if ( ujm )	ujm->addMessage ( cout, "<b>" + genShortFilenameFromPath ( f ) + "</b> already compressed.</b>" );
		else		ErrorHandler::genError ()->message ( "<b>" + genShortFilenameFromPath ( f ) + "</b> already compressed.</b>\n" );
		return false;
	}
	else {
		if ( genFileExists ( f + ".7z" ) ) {
			if ( ujm )	ujm->addMessage ( cout, "<b>" + genFilenameFromPath ( f ) + "</b> already compressed." );
			else		ErrorHandler::genError ()->message ( "<b>" + genFilenameFromPath ( f ) + "</b> already compressed.\n" );
			return false;
		}
		else {
			if ( ujm )	ujm->addMessage ( cout, "Compressing <b>" + genFilenameFromPath ( f ) + "</b>" );
			else		ErrorHandler::genError ()->message ( "Compressing <b>" + genFilenameFromPath ( f ) + "</b>\n" );
			gen7zaCreate ( f );
			return true;
		}
	}
}
bool compressProjectFiles ( const string& projectFile, const StringVector& expFiles, const UpdatingJavascriptMessage* ujm )
{
	if ( !compressZ7 ( projectFile, ujm ) ) return false;
	for ( StringVectorSizeType i = 0 ; i < expFiles.size () ; i++ ) {
		if ( !compressZ7 ( expFiles [i], ujm ) ) return false;
	}
	return true;
}
void compressResultsFiles ( const StringVector& resultsFiles, const UpdatingJavascriptMessage* ujm )
{
	for ( StringVectorSizeType i = 0 ; i < resultsFiles.size () ; i++ ) {
		compressZ7 ( resultsFiles [i], ujm );
	}
}
void compressDataFiles ( const StringVector& dataFiles, const UpdatingJavascriptMessage* ujm )
{
	if ( !dataFiles.empty () ) {
		for ( StringVectorSizeType i = 0 ; i < dataFiles.size () ; i++ ) {
			string path = dataFiles [i];
			if ( genIsDirectory ( path ) ) {
				FileList fList ( path, "", "", false );
				StringVector sv = fList.getNameList ();
				for ( StringVectorSizeType j = 0 ; j < sv.size () ; j++ ) {
					string df2 = sv [j];
					string path2 =  path + SLASH + df2;
					if ( !genIsDirectory ( path2 ) ) {
						compressZ7 ( path2, ujm );
					}
				}
			}
			else
				compressZ7 ( path, ujm );
		}
	}
}
void compressResults ( const string& dataToExport, const string& user, const string& project, const string& results, const UpdatingJavascriptMessage* ujm )
{
	StringVector resultsFiles = getResultsFiles ( user, project, results );
	ProjectInfo* pi = MySQLPPSDDBase::instance ().getProjectInfo ( user, project );
	string projectPath = pi->getProjectFullPath ();
	if ( genFileExists ( projectPath + ".7z" ) ) {
		if ( ujm )	ujm->addMessage ( cout, "<b>" + genFilenameFromPath ( projectPath ) + "</b> already compressed." );
		else		ErrorHandler::genError ()->message ( "<b>" + genFilenameFromPath ( projectPath ) + "</b> already compressed.\n" );
		return;
	}
	StringVector expFiles = getExpFiles ( pi->getProjectDirectory (), pi->getProjectName () );
	StringVector dataFiles = getCentroidAndRawFiles ( dataToExport, projectPath, true );
	bool flag = compressProjectFiles ( projectPath, expFiles, ujm );
	if ( flag ) {
		compressResultsFiles ( resultsFiles, ujm );
		compressDataFiles ( dataFiles, ujm );
	}
}
// Uncompression functions
void uncompressZ7 ( const string& f, const UpdatingJavascriptMessage* ujm )
{
	if ( !isFileType ( f, Z7 ) ) {
		if ( ujm )	ujm->addMessage ( cout, "<b>" + genFilenameFromPath ( f ) + "</b> not compressed.</b>" );
		else		ErrorHandler::genError ()->message ( "<b>" + genFilenameFromPath ( f ) + "</b> not compressed.</b>\n" );
	}
	else {
		if ( !genFileExists ( f ) ) {
			if ( ujm )	ujm->addMessage ( cout, "<b>" + genShortFilenameFromPath ( f ) + "</b> not compressed." );
			else		ErrorHandler::genError ()->message ( "<b>" + genShortFilenameFromPath ( f ) + "</b> not compressed.\n" );
		}
		else {
			if ( ujm )	ujm->addMessage ( cout, "Uncompressing <b>" + genFilenameFromPath ( f ) + "</b>" );
			else		ErrorHandler::genError ()->message ( "Uncompressing <b>" + genFilenameFromPath ( f ) + "</b>\n" );
			gen7zaUncompress ( f );
		}
	}
}
void uncompressProjectFiles ( const string& projectFile, const StringVector& expFiles, const UpdatingJavascriptMessage* ujm )
{
	uncompressZ7 ( projectFile + ".7z", ujm );
	for ( StringVectorSizeType i = 0 ; i < expFiles.size () ; i++ ) {
		uncompressZ7 ( expFiles [i], ujm );
	}
}
void uncompressResultsFiles ( const StringVector& resultsFiles, const UpdatingJavascriptMessage* ujm )
{
	for ( StringVectorSizeType i = 0 ; i < resultsFiles.size () ; i++ ) {
		uncompressZ7 ( resultsFiles [i], ujm );
	}
}
void uncompressDataFiles ( const StringVector& dataFiles, const UpdatingJavascriptMessage* ujm )
{
	if ( !dataFiles.empty () ) {
		for ( StringVectorSizeType i = 0 ; i < dataFiles.size () ; i++ ) {
			string path = dataFiles [i];
			if ( genFileExists ( path ) && genIsDirectory ( path ) ) {
				FileList fList ( path, "", "", false );
				StringVector sv = fList.getNameList ();
				for ( StringVectorSizeType j = 0 ; j < sv.size () ; j++ ) {
					string df2 = sv [j];
					string path2 =  path + SLASH + df2;
					if ( !genIsDirectory ( path2 ) ) {
						uncompressZ7 ( path2, ujm );
					}
				}
			}
			else
				uncompressZ7 ( path + ".7z", ujm );
		}
	}
}
void uncompressResults ( const string& dataToExport, const string& user, const string& project, const string& results, const UpdatingJavascriptMessage* ujm )
{
	StringVector resultsFiles = getResultsFiles ( user, project, results );
	ProjectInfo* pi = MySQLPPSDDBase::instance ().getProjectInfo ( user, project );
	StringVector expFiles = getExpFiles ( pi->getProjectDirectory (), pi->getProjectName () );
	uncompressProjectFiles ( pi->getProjectFullPath (), expFiles, ujm );
	StringVector dataFiles = getCentroidAndRawFiles ( dataToExport, pi->getProjectFullPath (), true );
	uncompressResultsFiles ( resultsFiles, ujm );
	uncompressDataFiles ( dataFiles, ujm );
}
string getExpFile ( const string& resultsFullPath )
{
	string eFile;
	GenIFStream ist ( resultsFullPath );
	string line;
	while ( getline ( ist, line ) ) {
		if ( line.find ( "<expect_coeff_file>" ) != string::npos ) {
			eFile = ParameterList::convertValue ( XMLParser::getStringValue ( line, "expect_coeff_file" ) );
		}
		if ( line.find ( "</parameters>" ) != string::npos ) break;
	}
	return eFile;
}
void checkProjects ( const string& user, const string& project, bool verbose )
{
	cout << "<b>User=</b>" << user << " <b>Project=</b>" << project << "<br />" << endl;

	// Check results files

	StringVector resultsList = MySQLPPSDDBase::instance ().getResultsList ( user, project );
	SetString expFiles;
	for ( StringVectorSizeType i = 0 ; i < resultsList.size () ; i++ ) {
		string resultsFullPath = MySQLPPSDDBase::instance ().getResultsFullPath ( user, project, resultsList [i] );
		string actualResultsFullPath = resultsFullPath;
		bool exists = genFileExists ( actualResultsFullPath );
		bool exists7z = false;
		if ( !exists ) {
			if ( genFileExists ( actualResultsFullPath + ".7z" ) ) {
				actualResultsFullPath += ".7z";
				exists7z = true;
			}
		}
		if ( exists7z ) {
			uncompressZ7 ( actualResultsFullPath, 0 );
			exists = true;
		}
		if ( !exists ) {
			cout << actualResultsFullPath << " missing<br />" << endl;
		}
		else {
			if ( verbose ) {
				cout << actualResultsFullPath << "<br />" << endl;
			}
			string eFile = getExpFile ( resultsFullPath );
			if ( !eFile.empty () ) expFiles.insert ( eFile );
			if ( exists7z ) compressZ7 ( resultsFullPath, 0 );
		}
	}
	ProjectInfo* pi = MySQLPPSDDBase::instance ().getProjectInfo ( user, project );

	// Check exp files

	for ( SetStringConstIterator j = expFiles.begin () ; j != expFiles.end () ; j++ ) {
		string fName = pi->getProjectDirectory () + SLASH + *j + ".xml";
		bool exists = genFileExists ( fName );
		if ( !exists ) {
			if ( genFileExists ( fName + ".7z" ) ) {
				fName += ".7z";
				exists = true;
			}
		}
		if ( !exists ) cout << fName << " missing<br />" << endl;
		else {
			if ( verbose ) cout << fName << "<br />" << endl;
		}
	}

	// Check project file

	string projectPath = pi->getProjectFullPath ();
	string actualProjectPath = projectPath;
	bool exists = genFileExists ( actualProjectPath );
	bool exists7z = false;

	if ( !exists ) {
		if ( genFileExists ( actualProjectPath + ".7z" ) ) {
			actualProjectPath += ".7z";
			exists7z = true;
		}
	}
	if ( exists7z ) {
		uncompressZ7 ( actualProjectPath, 0 );
		exists = true;
	}
	if ( !exists ) {
		cout << projectPath << " missing, unable to check data files<br />" << endl;
	}
	else {
		if ( verbose ) cout << actualProjectPath << "<br />" << endl;
		// Check data files

		StringVector dataFiles = getCentroidAndRawFiles ( "All Data", pi->getProjectFullPath (), false );
		for ( StringVectorSizeType k = 0 ; k < dataFiles.size () ; k++ ) {
			string path = dataFiles [k];
			string actualPath = path;
			bool dExists = genFileExists ( actualPath );
			bool dExists7z = false;
			if ( !dExists ) {
				if ( genFileExists ( actualPath + ".7z" ) ) {
					actualPath += ".7z";
					dExists7z = true;
				}
			}
			if ( dExists7z ) {
				dExists = true;
			}
			if ( !dExists ) cout << path << " missing<br />" << endl;
			else {
				if ( verbose ) cout << actualPath << "<br />" << endl;
			}
		}
		if ( exists7z ) compressZ7 ( projectPath, 0 );
	}
}
}

void exportData ( const ParameterList& paramList )
{
	string user = paramList.getStringValue ( "user" );
	StringVector results = paramList.getStringVectorValue ( "results_file" );
	string dataToExport = paramList.getStringValue ( "data_to_export" );
	//exportResults ( dataToExport, user, results );	// Not currently used.
	StringVector projects = paramList.getStringVectorValue ( "project_name" );
	checkCalProject ( projects );
	checkRunningJobs ( user, projects );
	UpdatingJavascriptMessage* ujm = new UpdatingJavascriptMessage;
	for ( StringVectorSizeType i = 0 ; i < projects.size () ; i++ ) {
		StringVector sv = genGetSubstrings ( projects [i], '/' );
		if ( sv.size () == 2 )
			exportResults ( dataToExport, sv [0], sv [1], "", ujm );
		else
			exportResults ( dataToExport, user, projects [i], "", ujm );
	}
	ujm->deletePreviousMessage ( cout );
	delete ujm;
}
void compressProjects ( const ParameterList& paramList )
{
	ErrorHandler::genError ()->message ( "<b>Please let the compression complete. Compression of large files can take a few minutes.</b>\n\n" );
	string user = paramList.getStringValue ( "user" );
	StringVector projects = paramList.getStringVectorValue ( "project_name" );
	checkCalProject ( projects );
	checkRunningJobs ( user, projects );
	UpdatingJavascriptMessage* ujm = new UpdatingJavascriptMessage;
	for ( StringVectorSizeType i = 0 ; i < projects.size () ; i++ ) {
		ujm->startWriteMessage ( cout );
		ujm->addMessage ( cout, "Compressing project <b>" + projects [i] + "</b>" );
		StringVector sv = genGetSubstrings ( projects [i], '/' );
		if ( sv.size () == 2 )
			compressResults ( "All Data", sv [0], sv [1], "", ujm );
		else
			compressResults ( "All Data", user, projects [i], "", ujm );
		ujm->addMessage ( cout, "Project <b>" + projects [i] + " compressed</b>." );
		ujm->endWriteMessage ( cout );
	}
	ujm->deletePreviousMessage ( cout );
	delete ujm;
}
void uncompressProjects ( const ParameterList& paramList )
{
	ErrorHandler::genError ()->message ( "<b>Please let the uncompression complete. Uncompression of large files can take a few minutes.</b>\n\n" );
	string user = paramList.getStringValue ( "user" );
	StringVector projects = paramList.getStringVectorValue ( "project_name" );
	checkCalProject ( projects );
	checkRunningJobs ( user, projects );
	UpdatingJavascriptMessage* ujm = new UpdatingJavascriptMessage;
	for ( StringVectorSizeType i = 0 ; i < projects.size () ; i++ ) {
		ujm->startWriteMessage ( cout );
		ujm->addMessage ( cout, "Uncompressing project <b>" + projects [i] + "</b>" );
		StringVector sv = genGetSubstrings ( projects [i], '/' );
		if ( sv.size () == 2 )
			uncompressResults ( "All Data", sv [0], sv [1], "", ujm );
		else
			uncompressResults ( "All Data", user, projects [i], "", ujm );
		ujm->addMessage ( cout, "Project <b>" + projects [i] + "</b> uncompressed.n" );
		ujm->endWriteMessage ( cout );
	}
	ujm->deletePreviousMessage ( cout );
	delete ujm;
}
bool uncompressProjects2 ( const ParameterList& paramList )
{
	string user = paramList.getStringValue ( "user" );
	StringVector projects = paramList.getStringVectorValue ( "project_name" );
	bool html = false;
	for ( StringVectorSizeType i = 0 ; i < projects.size () ; i++ ) {
		ProjectInfo* pi = MySQLPPSDDBase::instance ().getProjectInfo ( user, projects [i] );
		string projectPath = pi->getProjectFullPath ();
		bool compressed = false;
		if ( genFileExists ( projectPath + ".7z" ) ) {		// Check that the project isn't compressed
			if ( !html ) {
				init_html ( cout, "Uncompress Projects" );
				ErrorHandler::genError ()->message ( "<b>Please let the uncompression complete. Uncompression of large files can take a few minutes.</b>\n\n" );
				html = true;
			}
			uncompressResults ( "All Data", user, projects [i], "", 0 );
			compressed = true;
		}
	}
	return html;
}
void checkProjects ( const ParameterList& paramList )
{
	string user = paramList.getStringValue ( "user" );
	StringVector projects = paramList.getStringVectorValue ( "project_name" );
	checkRunningJobs ( user, projects );
	bool verbose = false;
	for ( StringVectorSizeType i = 0 ; i < projects.size () ; i++ ) {
		StringVector sv = genGetSubstrings ( projects [i], '/' );
		if ( sv.size () == 2 )
			checkProjects ( sv [0], sv [1], verbose );
		else
			checkProjects ( user, projects [i], verbose );
	}
}
#endif
