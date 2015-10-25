/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_import_proj.cpp                                            *
*                                                                             *
*  Created    : June 28th 2012                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2012-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef MYSQL_DATABASE
#ifndef VIS_C
#include <stdexcept>
#endif
#include <algorithm>
#include <lg_string.h>
#include <lgen_file.h>
#include <lgen_uncompress.h>
#include <ld_init.h>
#include <lu_file_type.h>
#include <lu_getfil.h>
#include <lu_import_proj.h>
#include <lu_param_list.h>
#include <lu_repository.h>
#include <lu_version.h>
#include <lu_xml.h>

using std::string;
using std::runtime_error;
using std::make_pair;
using std::endl;
using std::ostream;
using std::find;

using namespace FileTypes;

class ImportProjectFile {
	bool rawDataIncluded;
	char* info;
	int numLines;
	StringVector centroid;
	StringVector centroidFiles;
	StringVector raw;
	StringVector rawFiles;
	void getCentroidFiles ();
	void getRawFiles ();
	static string centroidBaseDir;
	static string rawBaseDir;
	static string userBaseDir;
	bool checkDataFile ( const string& f, bool centroid ) const;
	char getPathType ( const string& f ) const;
public:
	ImportProjectFile ( const string& file, bool rawDataIncluded );
	~ImportProjectFile ();
	void writeNewProjectFile ( const string& newFile, const string& path );
	bool checkDataFiles () const;
};
string ImportProjectFile::centroidBaseDir = adjustPPOutputPath ( InfoParams::instance ().getCentroidDir () );
string ImportProjectFile::rawBaseDir = adjustPPOutputPath ( InfoParams::instance ().getRawDir () );
string ImportProjectFile::userBaseDir = InfoParams::instance ().getUserRepository ();

ImportProjectFile::ImportProjectFile ( const string& file, bool rawDataIncluded ) :
	rawDataIncluded ( rawDataIncluded )
{
	if ( !genFileExists ( file ) ) throw runtime_error ( "The project file is missing from the uploaded archive.\n" );
	info = getFileInfo ( file, '\n', 1, false, &numLines );
	centroid = XMLParser::getStringVectorValue ( info, "centroid" );
	raw = XMLParser::getStringVectorValue ( info, "raw" );
	getCentroidFiles ();
	getRawFiles ();
}
ImportProjectFile::~ImportProjectFile ()
{
	delete [] info;
}
void ImportProjectFile::getCentroidFiles ()
{
	for ( int i = 0 ; i < centroid.size () ; i++ ) {
		centroidFiles.push_back ( genFilenameFromPath ( centroid [i] ) );
	}
}
void ImportProjectFile::getRawFiles ()
{
	for ( int i = 0 ; i < raw.size () ; i++ ) {
		rawFiles.push_back ( genFilenameFromPath ( raw [i] ) );
	}
}
char ImportProjectFile::getPathType ( const string& f ) const
{
	if ( genIsFullPath ( f ) )	return 'F';
	else if ( f[0] == '$' )			return '$';
	else if ( f[0] == '#' )			return '#';
	else							return '0';
}
bool ImportProjectFile::checkDataFile ( const string& f, bool centroid ) const
{
	string filename;
	if ( genIsFullPath ( f ) ) filename = f;
	else if ( f[0] == '$' ) {
		if ( centroid )	filename = centroidBaseDir + f.substr ( 1 );
		else			filename = rawBaseDir + f.substr ( 1 );
	}
	//else if ( f[0] == '#' ) filename = userBaseDir + SLASH + f.substr ( 1 );
	else if ( f[0] == '#' ) return false;	// Imported project with no data can only be from a data repository
	else return false;					// Illegal path
	return genFileExists ( filename );
}
bool ImportProjectFile::checkDataFiles () const
{
	for ( int i = 0 ; i < centroid.size () ; i++ ) {
		if ( !checkDataFile ( centroid [i], true ) ) return false;
	}
	for ( int j = 0 ; j < raw.size () ; j++ ) {
		if ( !checkDataFile ( raw [j], false ) ) return false;
	}
	return true;
}
void ImportProjectFile::writeNewProjectFile ( const string& newFile, const string& path )
{
	GenOFStream ofs ( newFile );
	for ( int i = 0, j = 0, k = 0 ; i < numLines ; i++ ) {
		string line = ( i == 0 ) ? strtok ( info, "\n" ) : strtok ( NULL, "\n" );
		if ( line.find ( "<centroid>" ) != string::npos ) {
			if ( path.empty () )
				ofs << "<centroid>" << centroid [j++] << "</centroid>" << endl;
			else
				ofs << "<centroid>" << path << "/" << centroidFiles [j++] << "</centroid>" << endl;
		}
		else if ( line.find ( "<raw>" ) != string::npos ) {
			if ( rawDataIncluded ) {
				if ( path.empty () )
					ofs << "<raw>" << raw [k++] << "</raw>" << endl;
				else
					ofs << "<raw>" << path << "/" << rawFiles [k++] << "</raw>" << endl;
			}
		}
		else {
			ofs << line << endl;
		}
	}
}
ImportProject::ImportProject ( const ParameterList& paramList ) :
	newProjectName ( paramList.getStringValue ( "new_project_name" ) )
{
	init ( paramList );
	try {												// These are the files with the database entries
		dpi = new DatabaseProjectImport ( uploadName );
	}
	catch ( runtime_error e ) {
		genUnlinkDirectory ( uploadName );
		throw e;
	}
	bool rawDataIncluded = getDataFiles ();
	getProjectFiles ();
	try {
		ipf = new ImportProjectFile ( uploadName + SLASH + "project" + SLASH + projFile, rawDataIncluded );
	}
	catch ( runtime_error e ) {
		genUnlinkDirectory ( uploadName );
		throw e;
	}
	if ( dataFiles.empty () ) {
		if ( !ipf->checkDataFiles () ) {
			genUnlinkDirectory ( uploadName );
			throw runtime_error ( "One or more of the data files in the uploaded project is not present." );
		}
	}
	getResultsFiles ();
}
ImportProject::~ImportProject ()
{
	delete dpi;
	delete ipf;
}
void ImportProject::init ( const ParameterList& paramList )
{
	string fpath = paramList.getStringValue ( "upload_temp_filepath" );
	string user = paramList.getStringValue ( "user" );
	if ( user == "root" ) {
		genUnlink ( fpath );
		throw runtime_error ( "Root user can't create projects." );
	}
	UserInfo* userInfo = MySQLPPSDDBase::instance ().getUserInfo ( user );
	if ( !userInfo ) {
		genUnlink ( fpath );
		throw runtime_error ( "Unknown user." );
	}
	userID = userInfo->getUserID ();
	if ( !newProjectName.empty () && MySQLPPSDDBase::instance ().checkProject ( userID, newProjectName ) ) {
		genUnlink ( fpath );
		throw runtime_error ( "Project already exists." );
	}
	reposit = new UserRepository ( "batchtag", userInfo->getDirectoryName () );	// Where to put files given a user. Also creates the directories.
	uploadName = genPreprocessFile ( fpath );
}
bool ImportProject::getDataFiles ()
{
	bool raw = false;
	FileList dataFileList ( uploadName + SLASH + "data", "", "", false );
	for ( FileList::size_type i = 0 ; i < dataFileList.size () ; i++ ) {
		string f = dataFileList [i];
		if ( f != "." && f != ".." ) {
			dataFiles.push_back ( f );
			if ( isFileType ( f, WIFF ) || isFileType ( f, RAW ) || genIsDirectory ( uploadName + SLASH + "data" + SLASH + f ) ) {
				raw = true;
			}
		}
	}
	if ( dataFiles.empty () )	return true;		// Data in repository
	else						return raw;
	// There don't need to be any data files.
}
void ImportProject::getProjectFiles ()
{
	string pPath = uploadName + SLASH + "project";
	if ( !genFileExists ( pPath ) ) {
		genUnlinkDirectory ( uploadName );
		throw runtime_error ( "The imported archive doesn't contain a project directory." );
	}
	FileList projFileList ( pPath, "", ".xml", false );
	for ( FileList::size_type i = 0 ; i < projFileList.size () ; i++ ) {
		string f = projFileList [i];
		if ( f.find ( ".exp." ) != string::npos ) {
			expFiles.push_back ( f );
			if ( Version::instance ().isNewerVersion ( getVersionFromPPXMLFile ( pPath + SLASH + f ) ) ) {
				genUnlinkDirectory ( uploadName );
				throw runtime_error ( "The imported project contains searches from a newer version than you have installed." );
			}
		}
		else if ( isFileType ( f, XML ) ) {
			if ( isXMLFile ( pPath + SLASH + f, "project" ) ) {
				projFile = f;
				projName = genShortFilenameFromPath ( f );
				if ( Version::instance ().isNewerVersion ( getVersionFromPPXMLFile ( pPath + SLASH + f ) ) ) {
					genUnlinkDirectory ( uploadName );
					throw runtime_error ( "The imported project file is from a newer version than you have installed." );
				}
			}
			else {
				genUnlinkDirectory ( uploadName );
				throw runtime_error ( "The project file in the uploaded archive has an illegal format." );
			}
		}
		else {
			genUnlinkDirectory ( uploadName );
			throw runtime_error ( "The project directory in the uploaded archive contains an illegal file." );
		}
	}
	if ( !newProjectName.empty () ) projName = newProjectName;	// Check if project is to be renamed

	if ( MySQLPPSDDBase::instance ().checkProject ( userID, projName ) ) {
		genUnlinkDirectory ( uploadName );
		throw runtime_error ( "Project already exists." );
	}
	if ( !newProjectName.empty () ) {	// Rename all the files
		PairStringString edit1 = make_pair ( string ("<project_name>"), "<project_name>" + newProjectName + "</project_name>" );
		if ( !importDataStreamEditor ( pPath + SLASH + projFile, pPath + SLASH + projName + ".xml", edit1 ) ) {
			genUnlinkDirectory ( uploadName );
			throw runtime_error ( "Problem encountered when renaming the project." );
		}
		projFile = projName + ".xml";
		for ( StringVectorSizeType i = 0 ; i < expFiles.size () ; i++ ) {
			string eFile = expFiles [i];
			string newExpName = projName + eFile.substr ( eFile.rfind ( ".exp." ) );
			PairStringString edit2;
			edit2 = make_pair ( string ("<output_filename>"), "<output_filename>" + genShortFilenameFromPath ( newExpName ) + "</output_filename>" );
			if ( !importDataStreamEditor ( pPath + SLASH + eFile, pPath + SLASH + newExpName, edit2 ) ) {
				genUnlinkDirectory ( uploadName );
				throw runtime_error ( "Problem encountered when renaming the project." );
			}
			expFiles [i] = newExpName;
		}
	}
}
void ImportProject::getResultsFiles ()
{
	string rPath = uploadName + SLASH + "results";
	FileList resFileList ( rPath, "", "", false );
	for ( FileList::size_type i = 0 ; i < resFileList.size () ; i++ ) {
		string f = resFileList [i];
		if ( isSuffix ( f, ".xml" ) ) {
			resFiles.push_back ( f );
			if ( Version::instance ().isNewerVersion ( getVersionFromPPXMLFile ( rPath + SLASH + f ) ) ) {
				genUnlinkDirectory ( uploadName );
				throw runtime_error ( "The imported project contains searches from a newer version than you have installed." );
			}
		}
		if ( isSuffix ( f, ".disc.txt" ) )	discFiles.push_back ( f );
	}
	parseResultsFiles ();
}
bool ImportProject::importDataStreamEditor ( const string& file, const string& newFile, const PairStringString& edit )
{
	int numLines;
	char* info = getFileInfo ( file, '\n', 1, false, &numLines );
	genUnlink ( file );
	GenOFStream ofs ( newFile );
	bool tagFlag = true;
	for ( int i = 0 ; i < numLines ; i++ ) {
		string line = ( i == 0 ) ? strtok ( info, "\n" ) : strtok ( NULL, "\n" );
		if ( tagFlag && line.find ( edit.first ) != string::npos ) {
			ofs << edit.second << endl;
			tagFlag = false;
		}
		else
			ofs << line << endl;
	}
	if ( tagFlag )	return false;
	else			return true;
}
bool ImportProject::importDataStreamEditor ( const string& file, const string& newFile, const string& searchKey, const string& newProject )
{
	int numLines;
	char* info = getFileInfo ( file, '\n', 1, false, &numLines );
	genUnlink ( file );
	GenOFStream ofs ( newFile );
	bool searchKeyFlag = !searchKey.empty ();
	bool newProjectFlag = !newProject.empty ();
	for ( int i = 0, j = 0, k = 0 ; i < numLines ; i++ ) {
		string line = ( i == 0 ) ? strtok ( info, "\n" ) : strtok ( NULL, "\n" );
		if ( newProjectFlag && line.find ( "<expect_coeff_file>" ) != string::npos ) {
			string eFile = XMLParser::getStringValue ( line, "expect_coeff_file" );
			string newExpectCoeffFile = newProject + eFile.substr ( eFile.rfind ( ".exp." ) );;
			ofs << "<expect_coeff_file>" << newExpectCoeffFile << "</expect_coeff_file>" << endl;
			newProjectFlag = false;
		}
		else if ( searchKeyFlag && line.find ( "<search_key>" ) != string::npos ) {
			ofs << "<search_key>" << searchKey << "</search_key>" << endl;
			searchKeyFlag = false;
		}
		else
			ofs << line << endl;
	}
	if ( !searchKey.empty () && searchKeyFlag ) return false;
	return true;
}
void ImportProject::parseResultsFiles ()
{
	string rPath = uploadName + SLASH + "results";
	for ( StringVectorSizeType i = 0 ; i < resFiles.size () ; i++ ) {
		string origSearchKey = genShortFilenameFromPath ( resFiles [i] );
		string searchKey = origSearchKey;
		int j = 0;
		for ( ; ; j++ ) {
			if ( MySQLPPSDDBase::instance ().checkSearchKey ( searchKey ) ) {
				searchKey = genRandomString ( 16 );
			}
			else break;
			if ( j > 6 ) {										// Some problem getting a unique search key
				genUnlinkDirectory ( uploadName );
				throw runtime_error ( "Search key already used." );
			}
		}
		if ( j > 0 || !newProjectName.empty () ) {	// Search key changed
			resFiles [i] = searchKey + ".xml";
			if ( !importDataStreamEditor ( rPath + SLASH + origSearchKey + ".xml", rPath + SLASH + resFiles [i], searchKey, newProjectName ) ) {
				genUnlinkDirectory ( uploadName );
				throw runtime_error ( "Problems encountered when renaming the project or changing the search key." );
			}
			dpi->updateSearchKey ( origSearchKey, searchKey );

			if ( j > 0 ) {
				string oldDiscFName = origSearchKey + ".disc.txt";
				StringVectorIterator iter = find ( discFiles.begin (), discFiles.end (), oldDiscFName );
				if ( iter != discFiles.end () ) {
					discFiles [iter-discFiles.begin()] = searchKey + ".disc.txt";
					genRename ( rPath + SLASH + oldDiscFName, rPath + SLASH + searchKey + ".disc.txt" );
				}
			}
		}
	}
}
void ImportProject::writeDataFiles ( const string& key ) const
{
	string dataPath = reposit->getFullDataPath () + "/" + key;
	genCreateDirectoryPath ( dataPath );
	for ( int i = 0 ; i < dataFiles.size () ; i++ ) {
		string f = dataFiles [i];
		string path = uploadName + SLASH + "data" + SLASH + f;
		if ( genIsDirectory ( path ) ) {
			genCreateDirectory ( dataPath + SLASH + f );
			FileList fList ( path, "", "", false );
			StringVector sv = fList.getNameList ();
			for ( StringVectorSizeType j = 0 ; j < sv.size () ; j++ ) {
				string df2 = sv [j];
				string path2 =  path + SLASH + df2;
				if ( !genIsDirectory ( path2 ) ) {
					copyFile ( path2, dataPath + SLASH + f + SLASH + df2 );
				}
			}
		}
		else
			copyFile ( path, dataPath + SLASH + f ); 
	}
}
void ImportProject::writeResultsFiles () const
{
	string resPath = reposit->getFullResultsPath ();
	genCreateDirectoryPath ( resPath );
	for ( int i = 0 ; i < resFiles.size () ; i++ ) {
		copyFile ( uploadName + SLASH + "results" + SLASH + resFiles [i], resPath + SLASH + resFiles [i] ); 
	}
	for ( int j = 0 ; j < discFiles.size () ; j++ ) {
		copyFile ( uploadName + SLASH + "results" + SLASH + discFiles [j], resPath + SLASH + discFiles [j] ); 
	}
}
void ImportProject::writeProjectFiles ( const string& key ) const
{
	string projFileDataPath;	// Only set if the data is uploaded
	if ( !key.empty () ) projFileDataPath = reposit->getProjectFileDataPath ( true ) + "/" + key;
	string projPath = reposit->getFullProjectPath ();
	genCreateDirectoryPath ( projPath );
	for ( int i = 0 ; i < expFiles.size () ; i++ ) {
		copyFile ( uploadName + SLASH + "project" + SLASH + expFiles [i], projPath + SLASH + expFiles [i] ); 
	}
	ipf->writeNewProjectFile ( projPath + SLASH + projFile, projFileDataPath );
}
void ImportProject::write ( ostream& os ) const
{
	string key;
	if ( !dataFiles.empty () ) {
		key = genRandomString ( 16 );
		writeDataFiles ( key );
	}
	writeResultsFiles ();
	writeProjectFiles ( key );
	dpi->submit ( userID, projName, projFile, reposit, resFiles );	//	Deal with database
	genUnlinkDirectory ( uploadName );
}
#endif
