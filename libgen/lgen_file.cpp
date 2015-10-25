/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lgen_file.cpp                                                 *
*                                                                             *
*  Created    : June 20th 1996                                                *
*                                                                             *
*  Purpose    : File and directory information functions.                     *
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
#include <sys/stat.h>
#include <sstream>
#include <errno.h>
#include <lg_stdio.h>
#include <lg_string.h>
#include <lgen_error.h>
#include <lgen_file.h>
#ifdef VIS_C
#include <windows.h>
#include <direct.h>
#else
#include <sys/param.h>
#include <unistd.h>
#include <netinet/in.h>
#include <sys/statvfs.h>
#include <dirent.h>
#include <sys/vfs.h>
#include <lgen_reg_exp.h>
#endif
using std::string;
using std::istringstream;

GENINT64 genFileSize ( const string& full_name )
{
#ifdef VIS_C
	struct _stati64 stats;

	if ( _stati64 ( full_name.c_str (), &stats ) == -1 ) {
#else
	struct stat stats;
	if ( stat ( full_name.c_str (), &stats ) == -1 ) {
#endif
		ErrorHandler::genError ()->error ( "File status (stat) failure.\nFunction name: genFileSize.\nFilename: " + full_name + ".\n" );
	}
	return ( stats.st_size );
}
GENINT64 genDirectorySize ( const string& name )	// Calculates the disk space taken up by a directory
{
	FileList fList ( name, "", "", false );
	StringVector sv = fList.getNameList ();
	GENINT64 siz = 0;
	for ( StringVectorSizeType i = 0 ; i < sv.size () ; i++ ) {
		if ( sv [i] != "." && sv [i] != ".." ) {
			string path = name + SLASH + sv [i];
			if ( genIsDirectory ( path ) ) siz += genDirectorySize ( path );	// Called recursively to deal with subdirectories
			else siz += genFileSize ( path );
		}
	}
	return siz;
}
string stripFilenameChars ( const string& filename )
{
	return gen_strstripchars ( filename, "\\/:*?\"<>|" );
}
bool isFileType ( const string& filename, const string& type )
{
	string t = "." + type;
	StringSizeType tLen = t.length ();
	return filename.length () > tLen && !genStrcasecmp ( filename.substr ( filename.length () - tLen ).c_str (), t.c_str () );
}
bool genFileExists ( const string& full_name )
{
	struct stat stats;

	return ( stat ( full_name.c_str (), &stats ) == 0 );
}
time_t genLastModifyTime ( const string& full_name )
{
	struct stat stats;

	if ( stat ( full_name.c_str (), &stats ) == -1 ) {
		ErrorHandler::genError ()->error ( "File status (stat) failure.\nFunction name: genLastModifyTime.\nFilename: " + full_name + ".\n" );
	}
	return ( stats.st_mtime );
}
FileList::FileList ( const string& directory, const string& prefix, const string& extension, bool deleteExtension ) :
	directory ( directory )
{
	string extensionToDelete;
	if ( deleteExtension ) extensionToDelete = extension;
#ifdef VIS_C
	string d = directory;
	if ( !d.empty () ) if ( d [d.length () - 1] != SLASH ) d += SLASH;
	string filespec = d + prefix + "*" + extension;

	addNames ( filespec, extensionToDelete );
#else
	DIR* dirp;
	if ( ( dirp = opendir ( directory.c_str () ) ) == NULL ) {		/*Directory doesn't exist*/
		return;
	}
	struct dirent *dp;
	for ( dp = readdir ( dirp ) ; dp != NULL ; dp = readdir ( dirp ) ) {
		string filename ( dp->d_name );
		if ( filename.length () >= prefix.length () + extension.length () ) {
			if ( !filename.compare ( 0, prefix.length (), prefix ) ) {
				if ( !filename.compare ( filename.length () - extension.length (), extension.length (), extension ) ) {
					nameList.push_back ( filename.substr ( 0, filename.length () - extensionToDelete.length () ) );
				}
			}
		}
	}
	closedir ( dirp );
#endif
}
FileList::FileList ( const string& directory, const string& prefix ) :
	directory ( directory )
{
#ifdef VIS_C
	string d = directory;
	if ( !d.empty () ) if ( d [d.length () - 1] != SLASH ) d += SLASH;
	addSubDirectoryNames ( d + prefix + string ( "*" ) );
#else
	DIR* dirp;
	if ( ( dirp = opendir ( directory.c_str () ) ) == NULL ) {		/*Directory doesn't exist*/
		return;
	}
	struct dirent *dp;
	for ( dp = readdir ( dirp ) ; dp != NULL ; dp = readdir ( dirp ) ) {
		string filename ( dp->d_name );
		if ( genIsDirectory ( directory + string ( "/" ) + filename ) ) {
			if ( filename [0] != '.' ) {
				if ( filename.length () >= prefix.length () ) {
					if ( !filename.compare ( 0, prefix.length (), prefix ) ) {
						nameList.push_back ( filename );
					}
				}
			}
		}
	}
	closedir ( dirp );
#endif
}
FileList::FileList ( const string& directory, const StringVector& regExpList ) :
	directory ( directory )
{
	string dir = directory;
	if ( directory [directory.length () - 1] != SLASH ) {
		dir += SLASH;
	}
#ifdef VIS_C
	for ( StringVectorConstIterator i = regExpList.begin () ; i != regExpList.end () ; i++ ) {
		string filespec = dir + *i;
		addNames ( filespec );
	}
#else
	DIR* dirp;
	if ( ( dirp = opendir ( dir.c_str () ) ) == NULL ) {		/*Directory doesn't exist*/
		return;
	}
	struct dirent *dp;
	for ( dp = readdir ( dirp ) ; dp != NULL ; dp = readdir ( dirp ) ) {
		string filename ( dp->d_name );
		for ( StringVectorConstIterator i = regExpList.begin () ; i != regExpList.end () ; i++ ) {
			RegularExpression re ( (*i).c_str (), false );
			if ( re.isPresent ( filename.c_str () ) ) {
				nameList.push_back ( filename );
			}
		}
	}
	closedir ( dirp );
#endif
}
#ifdef VIS_C
void FileList::addNames ( const string& filespec, const string& extensionToDelete )
{
	struct _finddata_t file_info;
	long hFile;
	if ( ( hFile = _findfirst ( filespec.c_str (), &file_info ) ) != -1L ) {
		do {
			string filename ( file_info.name );
			if ( !extensionToDelete.empty () ) {
				int len = filename.length () - extensionToDelete.length ();
				if ( filename.substr ( len ) == extensionToDelete ) {	// Check it really is the correct extension
					nameList.push_back ( filename.substr ( 0, len ) );
				}
			}
			else nameList.push_back ( filename );
		} while ( _findnext( hFile, &file_info ) == 0 );
		ErrorHandler::genError ()->resetErrorNumber ();	/* Make sure system error message isn't printed */
		_findclose( hFile );
	}
}
void FileList::addSubDirectoryNames ( const string& filespec )
{
	struct _finddata_t file_info;
	long hFile;
	if ( ( hFile = _findfirst ( filespec.c_str (), &file_info ) ) != -1L ) {
		do {
			if ( file_info.attrib & _A_SUBDIR ) {
				if ( file_info.name [0] != '.' ) {
					nameList.push_back ( file_info.name );
				}
			}
		} while ( _findnext( hFile, &file_info ) == 0 );
		ErrorHandler::genError ()->resetErrorNumber ();	/* Make sure system error message isn't printed */
		_findclose( hFile );
	}
}
#endif
void FileList::unlink () const
{
	for ( StringVectorSizeType i = 0 ; i < nameList.size () ; i++ ) {
		genUnlink ( directory + SLASH + nameList [i] );
	}
}
double genFreeDiskSpace ( const string& fileName )
{
#ifdef VIS_C
	DWORD SectorsPerCluster;
	DWORD BytesPerSector;
	DWORD NumberOfFreeClusters;
	DWORD TotalNumberOfClusters;
	char* str = NULL;
	int size;
	char* dummy;

	size = GetFullPathName ( fileName.c_str (), 0, str, &dummy );
	if ( size == 0 ) {
		ErrorHandler::genError ()->error ( "The file " + fileName + " does not exist." );
	}
	str = new char [size];
	GetFullPathName ( fileName.c_str (), size, str, &dummy );
	str [3] = 0;
	GetDiskFreeSpace ( str, &SectorsPerCluster, &BytesPerSector, &NumberOfFreeClusters, &TotalNumberOfClusters );
	delete [] str;

	return ( (double)NumberOfFreeClusters * (double)SectorsPerCluster * (double)BytesPerSector );
#else
	struct statvfs buf;
	if ( statvfs ( fileName.c_str (), &buf ) < 0 ) {
		ErrorHandler::genError ()->error ( "The file " + fileName + " does not exist." );
	}
	return ( buf.f_bavail * (double)buf.f_bsize );
#endif
}
double genVolumeSize ( const string& fileName )
{
#ifdef VIS_C
	DWORD SectorsPerCluster;
	DWORD BytesPerSector;
	DWORD NumberOfFreeClusters;
	DWORD TotalNumberOfClusters;
	char* str = NULL;
	int size;
	char* dummy;

	size = GetFullPathName ( fileName.c_str (), 0, str, &dummy );
	if ( size == 0 ) {
		ErrorHandler::genError ()->error ( "The file " + fileName + " does not exist." );
	}
	str = new char [size];
	GetFullPathName ( fileName.c_str (), size, str, &dummy );
	str [3] = 0;
	GetDiskFreeSpace ( str, &SectorsPerCluster, &BytesPerSector, &NumberOfFreeClusters, &TotalNumberOfClusters );
	delete [] str;

	return ( (double)TotalNumberOfClusters * (double)SectorsPerCluster * (double)BytesPerSector );
#else
	struct statvfs buf;
	if ( statvfs ( fileName.c_str (), &buf ) < 0 ) {
		ErrorHandler::genError ()->error ( "The file " + fileName + " does not exist." );
	}
	return ( buf.f_blocks * (double)buf.f_bsize );
#endif
}
// MAC \r
// DOS \r\n
// UNIX \n
bool isTextFile ( const string& fileName, int maxLen )
{
	FILE* fp1 = gen_fopen_binary ( fileName, "r", "isTextFile" );
	bool ret = false;
	for ( int i = 0 ; i < maxLen ; i++ ) {
		int c = getc ( fp1 );
		if ( c == EOF ) break;
		if ( c == '\r' || c == '\n' ) {
			ret = true;
			break;
		}
	}
	gen_fclose ( fp1, "isTextFile" );
	return ret;
}
bool checkForCR ( const string& file_name )
{
	FILE* fp1 = gen_fopen_binary ( file_name, "r", "checkForCR" );
	bool ret = false;

	for ( ; ; ) {
		int c = getc ( fp1 );
		if ( c == EOF ) break;
		if ( c == '\r' ) {
			ret = true;
			break;
		}
		if ( c == '\n' ) break;
	}
	gen_fclose ( fp1, "checkForCR" );
	return ret;
}
bool checkForCR ( const string& file_name, int maxLen )	// Use this version if you aren't sure it is a text file
{
	FILE* fp1 = gen_fopen_binary ( file_name, "r", "checkForCR" );
	bool ret = false;

	for ( int i = 0 ; i < maxLen ; i++ ) {
		int c = getc ( fp1 );
		if ( c == EOF ) break;
		if ( c == '\r' ) {
			ret = true;
			break;
		}
		if ( c == '\n' ) break;
	}
	gen_fclose ( fp1, "checkForCR" );
	return ret;
}
bool macTextFile ( const string& file_name )
{
	FILE* fp1 = gen_fopen_binary ( file_name, "r", "macTextFile" );
	bool ret = false;

	for ( ; ; ) {
		int c = getc ( fp1 );
		if ( c == EOF ) break;
		if ( c == '\r' ) {
			ret = true;
			continue;
		}
		if ( ret == true ) {
			if ( c == '\n' ) ret = false;
			break;
		}
	}
	gen_fclose ( fp1, "macTextFile" );
	return ret;
}
bool macTextFile ( const string& file_name, int maxLen )	// Use this version if you aren't sure it is a text file
{
	FILE* fp1 = gen_fopen_binary ( file_name, "r", "macTextFile" );
	bool ret = false;

	for ( int i = 0 ; i < maxLen ; i++ ) {
		int c = getc ( fp1 );
		if ( c == EOF ) break;
		if ( c == '\r' ) {
			ret = true;
			continue;
		}
		if ( ret == true ) {
			if ( c == '\n' ) ret = false;
			break;
		}
	}
	gen_fclose ( fp1, "macTextFile" );
	return ret;
}
void dos2Unix ( const string& name1 )
{
	string name2 = name1 + ".new";
	dos2Unix ( name1, name2 );
	genUnlink ( name1 );
	genRename ( name2, name1 );
}
void dos2Unix ( const string& name1, const string& name2 )
{
	FILE* fp1 = gen_fopen_binary ( name1, "r", "dos2Unix" );
	FILE* fp2 = gen_fopen_binary ( name2, "w+", "dos2Unix" );

	for ( ; ; ) {
		int c = getc ( fp1 );
		if ( c == EOF ) break;
		if ( c != '\r' ) putc ( c, fp2 );
	}
	gen_fclose ( fp1, "dos2Unix" );
	gen_fclose ( fp2, "dos2Unix" );
}
void mac2Unix ( const string& name1 )
{
	string name2 = name1 + ".new";
	mac2Unix ( name1, name2 );
	genUnlink ( name1 );
	genRename ( name2, name1 );
}
void mac2Unix ( const string& name1, const string& name2 )
{
	FILE* fp1 = gen_fopen_binary ( name1, "r", "dos2Unix" );
	FILE* fp2 = gen_fopen_binary ( name2, "w+", "dos2Unix" );

	for ( ; ; ) {
		int c = getc ( fp1 );
		if ( c == EOF ) break;
		if ( c != '\r' ) putc ( c, fp2 );
		else putc ( '\n', fp2 );
	}
	gen_fclose ( fp1, "dos2Unix" );
	gen_fclose ( fp2, "dos2Unix" );
}
bool convertMacFiles ( const string& name )
{
	if ( checkForCR ( name, 256 ) && macTextFile ( name ) ) {
		mac2Unix ( name );
		return true;
	}
	return false;
}
void convertToUnixTextWithMessage ( const string& name )
{
	if ( checkForCR ( name, 256 ) ) {
		ErrorHandler::genError ()->message ( "Converting file " + name + " to UNIX text format.\n" );
		if ( macTextFile ( name ) ) mac2Unix ( name );
		else dos2Unix ( name );
	}
}
void convertToUnixText ( const string& name )
{
	if ( checkForCR ( name, 256 ) ) {
		if ( macTextFile ( name ) ) mac2Unix ( name );
		else dos2Unix ( name );
	}
}
void copyFileWithMessage ( const string& name1, const string& name2 )
{
	string fileName = genFilenameFromPath ( name1 );
	ErrorHandler::genError ()->message ( "Copying " + fileName + ".\n" );
	copyFile ( name1, name2 );
}
void copyFile ( const string& name1, const string& name2 )
{
	FILE* fp1 = gen_fopen_binary ( name1, "r", "copyFile" );
	FILE* fp2 = gen_fopen_binary ( name2, "w+", "copyFile" );

	for ( ; ; ) {
		int c = getc ( fp1 );
		if ( c == EOF ) break;
		putc ( c, fp2 );
	}
	gen_fclose ( fp1, "copyFile" );
	gen_fclose ( fp2, "copyFile" );
}
bool genIsFullPath ( const string& path )
{
#ifdef VIS_C
	if ( path.length () < 2 ) return false;
	if ( path [1] == ':' || path [1] == SLASH ) return true;	// SLASH for UNC drive
	else return false;
#else
	if ( path.length () == 0 ) return false;
	if ( path [0] == SLASH ) return true;
	else return false;
#endif
}
string genSuffixFromPath ( const string& path )
{
	string::size_type index = path.find_last_of ( "." );
	if ( index == string::npos ) return "";
	else return path.substr ( index + 1, path.length () - index - 1 );
}
string genFilenameFromPath ( const string& path )
{
	string::size_type index = path.find_last_of ( "\\/" );
	if ( index == string::npos ) return path;
	else return path.substr ( index + 1, path.length () - index - 1 );
}
string genShortFilenameFromPath ( const string& path )
{
	string s = genFilenameFromPath ( path );
	string::size_type len = s.length ();
	string::size_type index;

	if ( isNoCaseSuffix ( s, ".tar.gz" ) )		index = len - strlen ( ".tar.gz" );
	else if ( isNoCaseSuffix ( s, ".tar.z" ) )	index = len - strlen ( ".tar.z" );
	else if ( isNoCaseSuffix ( s, ".tar.bz2" ) )index = len - strlen ( ".tar.bz2" );
	else if ( isNoCaseSuffix ( s, ".7z" ) )		index = len - strlen ( ".7z" );
	else
		index = s.find_last_of ( "." );
	if ( index == string::npos ) return s;
	else return s.substr ( 0, index );
}
string genDirectoryFromPath ( const string& path )
{
	string::size_type index = path.find_last_of ( "\\/" );
	if ( index == string::npos ) return "";
	else return path.substr ( 0, index );
}
string genCurrentWorkingDirectory ()
{
	string ret;
#ifdef VIS_C
	char* buffer = _getcwd ( NULL, 1 );
#else
	char* buffer = getcwd ( NULL, 1024 );
#endif
	if ( buffer == NULL ) {
		ErrorHandler::genError ()->error ( "The current path length is too long.\n" );
	}
	ret = buffer;
	free ( buffer );
	return ret;
}
void genChangeWorkingDirectory ( const string& workingDirectory )
{
#ifdef VIS_C
	int ret = _chdir ( workingDirectory.c_str () );
#else
	int ret = chdir ( workingDirectory.c_str () );
#endif
	if ( ret != 0 ) {
		ErrorHandler::genError ()->error ( "The specified working directory could not be found.\n" );
	}
}
#ifdef VIS_C
string genCurrentDOSDrive ()
{
	string str = genCurrentWorkingDirectory ();
	return str.substr ( 0, 2 );
}
#endif
string genNextFreeFileName ( const string& directory, const string& prefix, const string& suffix )
{
	FileList f ( directory, prefix, suffix, true );
	int nextFile = 0;
	for ( FileList::size_type i = 0 ; i < f.size () ; i++ ) {
		if ( isPrefix ( f [i], prefix ) ) {
			istringstream istr ( f [i].substr ( prefix.length () ) );
			int number;
			istr >> number;
			nextFile = genMax ( nextFile, number );
		}
	}
	nextFile++;
	return prefix + gen_itoa ( nextFile ) + suffix;
}
bool genCreateDirectory ( const string& dir )
{
#ifdef VIS_C
	if( _mkdir( dir.c_str () ) != 0 ) {
#else
	if( mkdir( dir.c_str (), 0xffffffff ) != 0 ) {
#endif
		if ( errno == EEXIST ) return true;
		return false;
	}
	return true;
}
bool genCreateNewDirectory ( const string& dir )
{
#ifdef VIS_C
	if( _mkdir( dir.c_str () ) != 0 ) {
#else
	if( mkdir( dir.c_str (), 0xffffffff ) != 0 ) {
#endif
		return false;
	}
	return true;
}
void genCreateDirectory ( const string& directory, bool errorIfPresent, bool deleteFiles )
{
	if ( !genCreateNewDirectory ( directory ) ) {
		if ( errorIfPresent ) {
			ErrorHandler::genError ()->error ( "Directory exists or problems creating directory.\n" );
		}
		else {
			if ( errno == EEXIST ) {
//				ErrorHandler::genError ()->message ( "Directory " + directory + " exists.\n" );
				if ( deleteFiles ) {
					StringVector regExpList;
					regExpList.push_back ( "*" );
					FileList fList ( directory, regExpList );
					StringVector fListNames = fList.getNameList ();
					for ( StringVectorSizeType i = 0 ; i < fListNames.size () ; i++ ) {
						if ( fListNames [i] != "." && fListNames [i] != ".." ) {
							genUnlink ( directory + SLASH + fListNames [i] );
						}
					}
				}
			}
			else {
				ErrorHandler::genError ()->error ( "Problems creating directory " + directory + ".\n" );
			}
		}
	}
}
string genMakeSubDirectory ( const string& dir, const string& subDir )
{
	string path = dir;
	if ( dir != "" ) {
		if ( dir [dir.length () - 1] != SLASH ) path += SLASH;
	}
	path += subDir;
	return path;
}
namespace {
string getNextString ( const string& s, const string& delim, string::size_type& start, string::size_type& end )
{
	end = s.find_first_of ( delim, start );
	string str = s.substr ( start, end-start );
	start = end + 1;
	return str;
}
}
void genCreateDirectoryPath ( const string& ssName )
{
	string::size_type start = 0;
	string::size_type end;
	string outDir;
	if ( ssName.length () && ssName [0] == '/' ) {
		start = 1;
		outDir = "/";
	}
	else if ( ssName.length () > 2 && ssName [0] == '\\' && ssName [1] == '\\' ) {	// UNC path
		start = 2;
		outDir = "\\\\";
		string s = getNextString ( ssName, "\\/", start, end );	// Get the map name
		outDir += s;
		outDir += SLASH;
	}
	do {
		string s = getNextString ( ssName, "\\/", start, end );
		outDir += s;
		if ( ( outDir.length () != 1 && !(outDir.length () < 4 && outDir [1] == ':' ) ) || outDir [1] == '\\' ) {
			if ( !genFileExists ( outDir ) ) genCreateDirectory ( outDir, false, false );
		}
		outDir += SLASH;
	} while ( end != string::npos );
}
void genRmdir ( const string& dir )
{
	rmdir ( dir.c_str () );
}
/*
unlink
------

VISUAL C++ 6.0
--------------

returns 0 if successful. 

Otherwise, the function returns -1 and sets errno to EACCES (13), which means the path
specifies a read-only file, or to ENOENT (2), which means the file or path is not found
or the path specified a directory.
*/
void genUnlink ( const string& file )
{
	unlink ( file.c_str () );
}
void genUnlink ( const StringVector& files )
{
	for ( StringVectorSizeType i = 0 ; i < files.size () ; i++ ) {
		genUnlink ( files [i].c_str () );
	}
}
void genUnlinkDirectory ( const string& name )	// Deletes a directory and all its contents
{
	FileList fList ( name, "", "", false );
	StringVector sv = fList.getNameList ();
	for ( StringVectorSizeType i = 0 ; i < sv.size () ; i++ ) {
		if ( sv [i] != "." && sv [i] != ".." ) {
			string path = name + SLASH + sv [i];
			if ( genIsDirectory ( path ) ) genUnlinkDirectory ( path );	// Called recursively to get rid of subdirectories
			else genUnlink ( path );
		}
	}
	genRmdir ( name );
}
bool genIsDirectory ( const string& filenameFullPath )
{
	struct stat d_name_stat;
	// Using stat() because dp->d_type returns DT_UNKNOWN or '0' on reiserfs on debian. This appears more universal.
	// if ( dp->d_type & DT_DIR ) { // Previous code
	string path = filenameFullPath;
	if ( path.find_last_of ( "\\/" ) == filenameFullPath.length () - 1 ) { // If the last character is a slash delete it
		path = path.substr ( 0, filenameFullPath.length () - 1 );
	}
	if ( stat ( path.c_str (), &d_name_stat ) != 0 ) {
		sprintf( gen_error_message, "Cannot get stat() for %s.", path.c_str () );
		ErrorHandler::genError ()->error ( gen_error_message );
	}
	if ( d_name_stat.st_mode & S_IFDIR ) return true;
	else return false;
}
string genTranslateSlash ( const string& s )
{
	return genTranslateCharacter ( s, '\\', '/' );
}
void genCopyDirectoryContents ( const string& oldName, const string& newName, bool message )
{
	FileList fList ( oldName, "", "", false );
	StringVector sv = fList.getNameList ();
	for ( StringVectorSizeType i = 0 ; i < sv.size () ; i++ ) {
		string f = sv [i];
		if ( f != "." && f != ".." ) {
			string oldPath = oldName + SLASH + f;
			string newPath = newName + SLASH + f;
			if ( genIsDirectory ( oldPath ) ) {
				genCreateDirectory ( newPath );
				genCopyDirectoryContents ( oldPath, newPath, message );	// Called recursively to deal with subdirectories
			}
			else {
				if ( message )	copyFileWithMessage ( oldPath, newPath );
				else			copyFile ( oldPath, newPath );
			}
		}
	}
}
int genRename ( const string& oldName, const string& newName, bool message )
{
#ifdef VIS_C
	if ( genIsDirectory ( oldName ) ) {		// rename doesn't work on Windows if moving directory to a network drive
		SHFILEOPSTRUCT fileOp;

		fileOp.wFunc = FO_MOVE;
		fileOp.fFlags = FOF_SILENT | FOF_NOCONFIRMATION | FOF_NOERRORUI | FOF_NOCONFIRMMKDIR;

		size_t fromLen = oldName.length ();
		char* from = new char [fromLen + 2];
		from [fromLen+1] = 0;				// must be doubly null terminated
		strcpy ( from, oldName.c_str () );
		fileOp.pFrom = from;

		size_t toLen = newName.length ();
		char* to = new char [toLen + 2];
		to [toLen+1] = 0;					// must be doubly null terminated
		strcpy ( to, newName.c_str () );
		fileOp.pTo = to;

		int err = SHFileOperation ( &fileOp );
		delete [] from;
		delete [] to;
		return err;
	}
	else
		return rename ( oldName.c_str (), newName.c_str () );
#else
	if ( genDifferentPartitions ( oldName, newName ) ) {
		bool isDir = genIsDirectory ( oldName );
		if ( isDir ) {
			genCreateDirectory ( newName );
			genCopyDirectoryContents ( oldName, newName, message );
			genUnlinkDirectory ( oldName );
		}
		else {
			if ( message )	copyFileWithMessage ( oldName, newName );
			else			copyFile ( oldName, newName );
			genUnlink ( oldName );
		}
	}
	else {
		return rename ( oldName.c_str (), newName.c_str () );
	}
#endif
}
bool genDifferentPartitions ( const string& file1, const string& file2 )
{
	struct stat stat1;
	if ( stat ( file1.c_str (), &stat1 ) == -1 ) {
		string parentDir = genDirectoryFromPath ( file1 );		// Check parent directory
		if ( stat ( parentDir.c_str (), &stat1 ) == -1 ) {
			ErrorHandler::genError ()->error ( "File status (stat) failure.\nFunction name: genDifferentPartitions.\nFilename: " + file1 + "\n" );
		}
	}
	struct stat stat2;
	if ( stat ( file2.c_str (), &stat2 ) == -1 ) {
		string parentDir = genDirectoryFromPath ( file2 );		// Check parent directory
		if ( stat ( parentDir.c_str (), &stat2 ) == -1 ) {
			ErrorHandler::genError ()->error ( "File status (stat) failure.\nFunction name: genDifferentPartitions.\nFilename: " + file2 + "\n" );
		}
	}
	return stat1.st_dev != stat2.st_dev;
}
StringVector getRepositoryKeyList ( const string& d, int keyLen )
{
	StringVector sv;
	string dir1 = d;
	FileList fList1 ( dir1 );
	StringVector n1 = fList1.getNameList ();
	for ( StringVectorSizeType i = 0 ; i < n1.size () ; i++ ) {
		if ( n1 [i].length () == 1 ) {
			string dir2 = dir1 + n1 [i] + SLASH;
			FileList fList2 ( dir2 );
			StringVector n2 = fList2.getNameList ();
			for ( StringVectorSizeType j = 0 ; j < n2.size () ; j++ ) {
				if ( n2 [j].length () == 1 ) {
					string dir3 = dir2 + n2 [j] + SLASH;
					FileList fList3 ( dir3 );
					StringVector n3 = fList3.getNameList ();
					for ( StringVectorSizeType k = 0 ; k < n3.size () ; k++ ) {
						if ( n3 [k].length () == keyLen ) {
							sv.push_back ( n3 [k] );
						}
					}
				}
			}
		}
	}
	return sv;
}
