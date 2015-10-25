/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_get_file.cpp                                               *
*                                                                             *
*  Created    : July 20th 1996                                                *
*                                                                             *
*  Purpose    : Utility function for reading parameter and database files.    *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1996-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef VIS_C
#include <process.h>
#else
#include <sys/types.h>
#include <unistd.h>
#endif
#include <algorithm>
#include <lg_string.h>
#include <lg_time.h>
#include <lgen_error.h>
#include <lgen_file.h>
#include <lu_getfil.h>
#include <lu_param_list.h>
#include <lu_check_db.h>
using std::string;
using std::copy;
using std::back_inserter;
using std::sort;
using std::endl;

class TempDir {
	string tempDir;
	string fullPathDir;
	TempDir ();
public:
	static TempDir& instance ();
	string getTempDir () const { return tempDir; }
	string getTempFullPathDir () const { return fullPathDir; }
};

class PackageHome {
	string packageHome;
	PackageHome ();
public:
	static PackageHome& instance ();
	string getPackageHome () const { return packageHome; }
};

TempDir::TempDir ()
{
	tempDir = "temp";
	fullPathDir = ppBaseDir () + tempDir;
}
TempDir& TempDir::instance ()
{
	static TempDir d;
	return d;
}
SeqdbDir::SeqdbDir ()
{
	seqdbDir = adjustPPOutputPath ( InfoParams::instance ().getSeqdbDir () );
}
SeqdbDir& SeqdbDir::instance ()
{
	static SeqdbDir d;
	return d;
}
StringVector SeqdbDir::getDatabaseList ( bool dbSearchFlag ) const
{
	StringVector sv = FileList ( seqdbDir, "", ".idc", true ).getNameList ();
	sort ( sv.begin (), sv.end (), genStrcasecmpAscending () );		// Sort the databases into order
	NoCaseSetString ss;
	for ( int i = 0 ; i < sv.size () ; i++ ) {
		string db = sv [i];
		string dbPath = seqdbDir + db;
		if ( !genFileExists ( dbPath + getDatabaseSuffix ( db ) ) ) continue;	// Check for database file
		if ( !genFileExists ( dbPath + ".idc" ) ) continue;
		if ( !genFileExists ( dbPath + ".idi" ) ) continue;
		if ( !genFileExists ( dbPath + ".idp" ) ) continue;
		if ( !genFileExists ( dbPath + ".tax" ) ) continue;
		if ( !genFileExists ( dbPath + ".tl" ) ) continue;
		if ( is_numeric_acc_number_database ( db ) ) {
			if ( !genFileExists ( dbPath + ".acn" ) ) continue;
		}
		else {
			if ( !genFileExists ( dbPath + ".acc" ) ) continue;
		}
		if ( !is_dna_database ( db ) ) {
			if ( !genFileExists ( dbPath + ".mw" ) ) continue;
			if ( !genFileExists ( dbPath + ".pi" ) ) continue;
			if ( genFileSize ( dbPath + ".mw" ) != genFileSize ( dbPath + ".pi" ) ) continue;
		}
		ss.insert ( db );
		if ( dbSearchFlag ) {	// this is a database search program (ie not MS-Bridge/Digest)
			int len = 0;
			if ( isSuffix ( db, ".random" ) ) len = 7;
			if ( isSuffix ( db, ".reverse" ) )len = 8;
			if ( len && ss.find ( db.substr ( 0, db.length () - len ) ) != ss.end () ) ss.insert ( db + ".concat" );
		}
	}
	StringVector sv2;
	copy ( ss.begin (), ss.end (), back_inserter ( sv2 ) );
	return sv2;
}
StringVector SeqdbDir::getUserDatabaseList ( bool dbSearchFlag ) const
{
	StringVector sv;
	sv.push_back ( "User Protein" );
	StringVector nl = getDatabaseList ( dbSearchFlag );
	copy ( nl.begin (), nl.end (), back_inserter ( sv ) );
	return sv;
}
string SeqdbDir::getDatabaseSuffix ( const string& databaseName ) const
{
	static string databaseSuffix = InfoParams::instance ().getStringValue ( "database_suffix" );
	if ( !databaseSuffix.empty () && genFileExists ( seqdbDir + databaseName + databaseSuffix ) ) return databaseSuffix;
	else return "";
}
string SeqdbDir::getDatabasePath ( const string& databaseName ) const
{
	static string databaseSuffix = InfoParams::instance ().getStringValue ( "database_suffix" );
	string path = seqdbDir + databaseName;
	if ( !databaseSuffix.empty () && genFileExists ( path + databaseSuffix ) ) path += databaseSuffix;
	return path;
}
string SeqdbDir::getDatabasePathCreateOrAppend ( const string& databaseName ) const
{
	static string databaseSuffix = InfoParams::instance ().getStringValue ( "database_suffix" );
	string path = seqdbDir + databaseName;
	if ( genFileExists ( path ) ) return path;	// Existing database with no suffix
	path += databaseSuffix;
	return path;									// In other cases the database has a suffix which may be empty
}
Version::Version ()
{
	version = "5.14.3";
#ifdef BATCHTAG
#ifndef MYSQL_DATABASE
	version += " CL";
#endif
#else
	version += " Basic";
#endif
}
Version& Version::instance ()
{
	static Version d;
	return d;
}
void Version::getVersionParts ( const string& v, int& vMain, int& vSub, int& vInc )
{
	int end = v.find ( '.' );
	if ( end == string::npos ) {
		vMain = atoi ( v.c_str () );
		vSub = 0;
		vInc = 0;
	}
	else {
		vMain = atoi ( v.substr ( 0, end ).c_str () );
		int start = end+1;
		end = v.find ( '.', start );
		if ( end == string::npos ) {
			vSub = atoi ( v.substr ( start ).c_str () );
			vInc = 0;
		}
		else {
			vSub = atoi ( v.substr ( start, end-start ).c_str () );
			vInc = atoi ( v.substr ( end+1 ).c_str () );
		}
	}
}
bool Version::isOlderVersion ( const string& v ) const
{
	return isOlderVersion ( v, version );
}
bool Version::isNewerVersion ( const string& v ) const
{
	return isNewerVersion ( v, version );
}
bool Version::isOlderVersion ( const string& v1, const string& v2 )
{
	if ( v1 == v2 ) return false;		// Same versions
	int v1Main, v1Sub, v1Inc;
	int v2Main, v2Sub, v2Inc;
	getVersionParts ( v1, v1Main, v1Sub, v1Inc );
	getVersionParts ( v2, v2Main, v2Sub, v2Inc );
	if ( v1Main < v2Main ) return true;
	if ( v1Main > v2Main ) return false;
	if ( v1Sub < v2Sub ) return true;
	if ( v1Sub > v2Sub ) return false;
	if ( v1Inc < v2Inc ) return true;
	if ( v1Inc > v2Inc ) return false;
	return false;
}
bool Version::isNewerVersion ( const string& v1, const string& v2 )
{
	if ( v1 == v2 ) return false;
	return !isOlderVersion ( v1, v2 );
}
PackageHome::PackageHome ()
{
	packageHome = "";	// Installer version
}
PackageHome& PackageHome::instance ()
{
	static PackageHome d;
	return d;
}
MsparamsDir::MsparamsDir () :
	globalFlag ( false )
{
	msparamsDir = adjustPPOutputPath ( "params" );
	GenNameValueStream nvs ( msparamsDir + "info.txt" );
	if ( nvs.getValue ( "params_dir", globalMSParamsDir ) ) {
		globalMSParamsDir = adjustPPOutputPath ( globalMSParamsDir );
		globalFlag = true;
	}
}
MsparamsDir& MsparamsDir::instance ()
{
	static MsparamsDir d;
	return d;
}
string MsparamsDir::getParamPath ( const string& filename ) const
{
	if ( globalFlag && genFileExists ( globalMSParamsDir + filename ) ) return globalMSParamsDir + filename;
	else return msparamsDir + filename;
}
HTMLDir::HTMLDir ()
{
	htmlDir = "html";
}
HTMLDir& HTMLDir::instance ()
{
	static HTMLDir d;
	return d;
}
string HTMLDir::getVirtualHTMLDir () const
{
	return ParameterList::getVirtualDir () + htmlDir;
}
string HTMLDir::getVirtualHTMLJavaDir () const
{
	return getVirtualHTMLDir () + "/java/";
}
string HTMLDir::getVirtualHTMLImagesDir () const
{
	return getVirtualHTMLDir () + "/images/";
}
string HTMLDir::getHTMLFullJavaDir () const
{
	return ParameterList::getServer () + getVirtualHTMLJavaDir ();
}
BinDir::BinDir ()
{
	binDir = "cgi-bin";
}
BinDir& BinDir::instance ()
{
	static BinDir d;
	return d;
}
string BinDir::getVirtualBinDir () const
{
	return ParameterList::getVirtualDir () + binDir + "/";
}
SystemCallBinDir::SystemCallBinDir ()
{
	string binDir = BinDir::instance ().getBinDir ();
	if ( isSuffix ( genCurrentWorkingDirectory (), binDir ) )			// XP or UNIX
		systemCallBinDir = ".";
	else {
		systemCallBinDir = PackageHome::instance ().getPackageHome ();	// NT OR 2000
		if ( !systemCallBinDir.empty () ) systemCallBinDir += SLASH;
		systemCallBinDir += binDir;
	}
}
SystemCallBinDir& SystemCallBinDir::instance ()
{
	static SystemCallBinDir d;
	return d;
}
string getSystemCall ( const string& executable )
{
	string command ( SystemCallBinDir::instance ().getSystemCallBinDir () );
	command += SLASH;
	command += executable;
	return command;
}
static void removeComments ( char* info, int size, int* newSize )
{
	int i, j;
	bool checkFirstChar = true;
	bool comment = false;
	char character;

	for ( i = 0, j = 0 ; i < size ; i++ ) {
		character = info [i];
		if ( checkFirstChar ) {
			if ( character == '#' ) {
				comment = true;
			}
			checkFirstChar = false;
		}
		if ( !comment ) {
			info [j] = character;
			j++;
		}
		if ( character == '\n' ) {
			checkFirstChar = true;
			comment = false;
		}
	}
	*newSize = j;
}
static int getNumSeparators ( char* pointer, int numBytes, char separator )
{
	int numSeparators = 0;

	for ( int i = 0 ; i < numBytes ; i++ ) {
		if ( *pointer++ == separator ) {
			numSeparators++;
		}
	}
	return ( numSeparators );
}
char* getParamsFileInfo ( const string& filename )
{
	int dummy;
	return getFileInfo ( MsparamsDir::instance ().getParamPath ( filename ), '\n', 1, true, &dummy );
}
char* getParamsFileInfo ( const string& filename, int* numEntries )
{
	return getFileInfo ( MsparamsDir::instance ().getParamPath ( filename ), '\n', 1, true, numEntries );
}
char* getFileInfo ( const string& filename, char separator, int separatorsPerEntry, bool deleteComments )
{
	int dummy;
	return getFileInfo ( filename, separator, separatorsPerEntry, deleteComments, &dummy );
}
char* getFileInfo ( const string& filename, char separator, int separatorsPerEntry, bool deleteComments, int* numEntries )
{
	GenIFStream ist1 ( filename );

	int size = (int) genFileSize ( filename );
	char* info = new char [size+1];
	ist1.read ( (char*) info, size );
	/* size may be different to the physical file size for an microsoft file in
	which the <CR> characters have been stripped out */
	size = ist1.gcount ();

	if ( deleteComments ) {
		int oldSize = size;
		removeComments ( info, oldSize, &size );
	}
	info [size] = 0;
	*numEntries = getNumSeparators ( info, size, separator ) / separatorsPerEntry;

	return ( info );
}
char* getFileAsCharPtr ( const string& filename )
{
	GenIFStream ist1 ( filename );

	int size = (int) genFileSize ( filename );
	char* info = new char [size+1];
	ist1.read ( (char*) info, size );
	/* size may be different to the physical file size for an microsoft file in
	which the <CR> characters have been stripped out */
	size = ist1.gcount ();
	info [size] = 0;
	return info;
}
string adjustPPOutputPath ( const string& path )
{
	string newPath;
	if ( genIsFullPath ( path ) )	// If path is a full path
		newPath = path;
	else {
		string workDir = genCurrentWorkingDirectory ();
		string runDirectory = genFilenameFromPath ( workDir );
		if ( runDirectory == BinDir::instance ().getBinDir () ) {	// If the current working directory is the bin directory
			newPath = "..";
			newPath += SLASH;
			newPath += path;
		}
		else {
			newPath = workDir;
			newPath += SLASH;
			string pHome = PackageHome::instance ().getPackageHome ();
			if ( !pHome.empty () ) {
				newPath += pHome;
				newPath += SLASH;
			}
			newPath += path;
		}
	}
	if ( newPath [newPath.length () - 1] != SLASH ) {
		newPath += SLASH;
	}
	return newPath;
}
string adjustedPPCurrentWorkingDirectory ( const string& path )
{
	string runDirectory = genFilenameFromPath ( genCurrentWorkingDirectory () );
	if ( runDirectory == BinDir::instance ().getBinDir () ) {
		string newPath = path.substr ( 0, path.length () - runDirectory.length () - 1 );
		return newPath;
	}
	else
		return path;
}
string ppBaseDir ()
{
	string baseDir;
	string currentDirectory = genCurrentWorkingDirectory ();
	string runDirectory = genFilenameFromPath ( currentDirectory );
	string binDir = BinDir::instance ().getBinDir (); 
	if ( runDirectory == binDir )
		baseDir = currentDirectory.substr ( 0, currentDirectory.find ( binDir ) );
	else
		baseDir = currentDirectory + SLASH + PackageHome::instance ().getPackageHome ();
	if ( baseDir [baseDir.length () - 1] != SLASH ) {
		baseDir += SLASH;
	}
	return baseDir;
}
PPTempFile::PPTempFile ( const string& directory, const string& suffix, bool dummy ) : // Creates a temp file in a given directory
	directory ( directory ),
	fullPathDir ( genIsFullPath ( directory ) ? directory : ppBaseDir () + directory ),
	adjustedDir ( adjustPPOutputPath ( directory ) )
{
	string s ( tmpnam (NULL) );
#ifdef VIS_C
	if ( s [0] == '\\' ) s = s.substr ( 1 );
#else
	s = genReplaceSubstrings ( s, "/tmp/", "" );
#endif
	if ( !suffix.empty () ) {
		if ( s [s.length ()-1] != '.' ) s += '.';
		s += suffix;
	}
	else {
		if ( s [s.length ()-1] == '.' ) s = s.substr ( 0, s.length ()-1 );	// delete dot if no suffix
	}
	filename = s;
}
PPTempFile::PPTempFile ( const string& prefix, const string& suffix ) :
	directory ( TempDir::instance ().getTempDir () + SLASH + genCurrentDateString ( '_' ) ),
	fullPathDir ( TempDir::instance ().getTempFullPathDir () + SLASH + genCurrentDateString ( '_' ) ),
	adjustedDir ( adjustPPOutputPath ( directory ) )
{
	string actualPrefix = prefix;
	string s ( tmpnam (NULL) );
#ifdef VIS_C
	if ( s [0] == '\\' ) s = s.substr ( 1 );
#else
	s = genReplaceSubstrings ( s, "/tmp/", "" );
#endif
	actualPrefix += s;
	string fullPath = TempDir::instance ().getTempFullPathDir ();
	FileList fList ( fullPath );
	StringVector names = fList.getNameList ();
	for ( int i = 0 ; i < names.size () ; i++ ) {
		if ( genMoreThan2DaysOld ( names [i] ) ) {
			genUnlinkDirectory ( fullPath + SLASH + names [i] );
		}
	}
	if ( !genCreateDirectory ( adjustedDir ) ) {
		ErrorHandler::genError ()->error ( "Unable to create directory to hold data file.\n" );
	}
	filename = genNextFreeFileName ( adjustedDir, actualPrefix, suffix );
}
string PPTempFile::getURL () const
{
	string url;
	string vdp = InfoParams::instance ().getStringValue ( "virtual_dir_proxy" );
	if ( !vdp.empty () )
		url += "/" + vdp + "/";
	else
		url += ParameterList::getVirtualDir ();
	url += genTranslateSlash ( directory + SLASH + filename );
	return url;
}
InfoParams::InfoParams () :
	GenNameValueStream ( MsparamsDir::instance ().getParamPath ( "info.txt" ) )
{
}
InfoParams& InfoParams::instance ( bool reset )
{
	static InfoParams* d = new InfoParams ();
	if ( reset ) {
		delete d;
		d = new InfoParams ();
	}
	return *d;
}
time_t InfoParams::getLastModifiedTime ()
{
	return genLastModifyTime ( MsparamsDir::instance ().getParamPath ( "info.txt" ) );
}
string InfoParams::getPlatformValue ( const string& param ) const
{
	string v;
#ifdef VIS_C
	v = instance ().getStringValue ( param + "_win" );
#else
	v = instance ().getStringValue ( param + "_unix" );
#endif
	if ( v.empty () ) {
		v = instance ().getStringValue ( param );
	}
	return v;
}
void pidLogOutput ( const string& message, bool showTime )
{
	static string logFileName = ppBaseDir () + "logs" + SLASH + gen_itoa ( getpid () ) + ".txt";
	static GenOFStream ofs ( logFileName, std::ios_base::out | std::ios_base::app );
	if ( showTime ) ofs << genCurrentTimeAndDateString () << ": ";
	ofs << message << endl;
}
