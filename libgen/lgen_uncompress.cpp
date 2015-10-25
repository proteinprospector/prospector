/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lgen_uncompress.cpp                                           *
*                                                                             *
*  Created    : March 2nd 2007                                                *
*                                                                             *
*  Purpose    : Uncompression functions.                                      *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2007-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lgen_error.h>
#include <lg_stdlib.h>
#include <lg_string.h>
#include <lgen_file.h>
#include <lu_getfil.h>
using std::string;

namespace {
string getUnzipError ( int err )
{
	static MapIntToString uzErr;
	if ( uzErr.empty () ) {
		uzErr [0] = "Normal; no errors or warnings detected.";
		uzErr [1] = "One or more warning errors were encountered, but processing completed successfully anyway.\nThis includes zipfiles where one or more files was skipped due to unsupported compression method or encryption with an unknown password.";
		uzErr [2] = "A generic error in the zipfile format was detected.\nProcessing may have completed successfully anyway; some broken zipfiles created by other archivers have simple work-arounds.";
		uzErr [3] = "A severe error in the zipfile format was detected.\nProcessing probably failed immediately.";
		uzErr [4] = "unzip was unable to allocate memory for one or more buffers during program initialization.";
		uzErr [5] = "unzip was unable to allocate memory or unable to obtain a tty to read the decryption password(s).";
		uzErr [6] = "unzip was unable to allocate memory during decompression to disk.";
		uzErr [7] = "unzip was unable to allocate memory during in-memory decompression.";
		uzErr [8] = "[currently not used]";
		uzErr [9] = "The specified zipfiles were not found.";
		uzErr [10] = "Invalid options were specified on the command line.";
		uzErr [11] = "No matching files were found.";
		uzErr [50] = "The disk is (or was) full during extraction.";
		uzErr [51] = "The end of the ZIP archive was encountered prematurely.";
		uzErr [80] = "The user aborted unzip prematurely with control-C (or similar)";
		uzErr [81] = "Testing or extraction of one or more files failed due to unsupported compression methods or unsupported decryption.";
		uzErr [82] = "No files were found due to bad decryption password(s).";	//(If even one file is successfully processed, however, the exit status is 1.)
	}
	MapIntToStringConstIterator cur = uzErr.find ( err );
	if ( cur == uzErr.end () )
		return "Unknown error.";
	else
		return (*cur).second;
}
}

string genUnzip ( const string& filename, bool deleteFile, bool createSubdirectory )
{
	string command;
#ifdef VIS_C
	command += SystemCallBinDir::instance ().getSystemCallBinDir ();
	command += SLASH;
#endif
	command += "unzip";	// Info-ZIP unzip utility 
	command += " -q";		// Quiet mode
	command += " -P pword";	// Dummy password
	string newDir;
	if ( createSubdirectory ) {
		newDir = filename.substr ( 0, filename.length () - 4 );
		command += " -d";		// Subdirectory name
		command += " ";
		command += "\"" + newDir + "\"";		// Remove the .zip suffix
	}
	command += " ";
	command += "\"" + filename + "\"";
	int ret = genSystem ( command, "", true );
#ifndef VIS_C
	if ( WIFEXITED ( ret ) ) { 
		ret = WEXITSTATUS ( ret );
	}
#endif
	if ( deleteFile ) genUnlink ( filename );
	if ( ret != 0 && ret != 1 ) {
		if ( createSubdirectory ) genUnlinkDirectory ( newDir );
		ErrorHandler::genError ()->error ( "\nunzip error message:\n" + getUnzipError ( ret ) + "\n" );
	}
	return newDir;
}
/*

http://www.7-zip.org/download.html

Usage: 7za <command> [<switches>...] <archive_name> [<file_names>...]
       [<@listfiles...>]

<Commands>
  a: Add files to archive
  b: Benchmark
  d: Delete files from archive
  e: Extract files from archive (without using directory names)
  l: List contents of archive
  t: Test integrity of archive
  u: Update files to archive
  x: eXtract files with full paths
<Switches>
  -ai[r[-|0]]{@listfile|!wildcard}: Include archives
  -ax[r[-|0]]{@listfile|!wildcard}: eXclude archives
  -bd: Disable percentage indicator
  -i[r[-|0]]{@listfile|!wildcard}: Include filenames
  -m{Parameters}: set compression Method
  -o{Directory}: set Output directory
  -p{Password}: set Password
  -r[-|0]: Recurse subdirectories
  -scs{UTF-8 | WIN | DOS}: set charset for list files
  -sfx[{name}]: Create SFX archive
  -si[{name}]: read data from stdin
  -slt: show technical information for l (List) command
  -so: write data to stdout
  -ssc[-]: set sensitive case mode
  -ssw: compress shared files
  -t{Type}: Set type of archive
  -v{Size}[b|k|m|g]: Create volumes
  -u[-][p#][q#][r#][x#][y#][z#][!newArchiveName]: Update options
  -w[{path}]: assign Work directory. Empty path means a temporary directory
  -x[r[-|0]]]{@listfile|!wildcard}: eXclude filenames
  -y: assume Yes on all queries
 
0 - Normal (no errors or warnings detected)
1 - Warning (Non fatal error(s)). For example, some files cannot be read during compressing. So they were not compressed
2 - Fatal error
7 - Bad command line parameters
8 - Not enough memory for operation
255 - User stopped the process with control-C (or similar)
*/
namespace {
string get7zaError ( int err )
{
	static MapIntToString uzErr;
	if ( uzErr.empty () ) {
		uzErr [0] = "Normal (no errors or warnings detected).";
		uzErr [1] = "Warning (Non fatal error(s)).\nFor example, some files cannot be read during compressing.\nSo they were not compressed.";
		uzErr [2] = "Fatal error.";
		uzErr [7] = "Bad command line parameters.";
		uzErr [8] = "Not enough memory for operation.";
		uzErr [255] = "User stopped the process with control-C (or similar).";
	}
	MapIntToStringConstIterator cur = uzErr.find ( err );
	if ( cur == uzErr.end () )
		return "Unknown error.";
	else
		return (*cur).second;
}
}
//
// 7za a -t7z files.7z *.txt
//
// a - archive or add
//
// t - type of archive you want to create (7z, gzip, zip, bzip2, tar, iso, udf)
//
// archive name (here files.7z)
//
// list of files (with wildcard)
//
// -o dir (eg -o c:\Doc) destination directory
//
// Instead of *.txt you can specify an entire subdirectory (eg subdir\)
//
// -y - answer yes to everything
//
string gen7zaUncompressStart ()
{
	string command;
#ifdef VIS_C
	command += SystemCallBinDir::instance ().getSystemCallBinDir ();
	command += SLASH;
#endif
	command += "7za";		// 7za utility 
	command += " x";		// Extract
	command += " -y";		// Answer yes to everything
	command += " -ppword";	// Dummy password
	return command;
}
string gen7za ( const string& filename, bool deleteFile, bool createSubdirectory, const string& suffix = ".7z" )
{
	int len = filename.length () - suffix.length ();
	string command = gen7zaUncompressStart ();
	string newDir;
	if ( createSubdirectory ) {
		newDir = filename.substr ( 0, len );		// Remove the .7z suffix
		command += " -o";		// Subdirectory name
		command += "\"" + newDir + "\"";
	}
	command += " ";
	command += "\"" + filename + "\"";
#ifndef VIS_C
	command += " > /dev/null";
#endif
	int ret = genSystem ( command, "", true );
#ifndef VIS_C
	if ( WIFEXITED ( ret ) ) { 
		ret = WEXITSTATUS ( ret );
	}
#endif
	if ( ret != 0 && ret != 1 ) {
		if ( createSubdirectory ) genUnlinkDirectory ( newDir );
		ErrorHandler::genError ()->error ( "\n7za error message:\n" + get7zaError ( ret ) + "\n" );
	}
	if ( deleteFile ) genUnlink ( filename );
	return newDir;
}
void gen7zaUncompress ( const string& filename )
{
	string command = gen7zaUncompressStart ();
	command += " -o";		// Subdirectory name
	command += "\"" + genDirectoryFromPath ( filename ) + "\"";
	command += " ";
	command += "\"" + filename + "\"";
#ifndef VIS_C
	command += " > /dev/null";
#endif
	int ret = genSystem ( command, "", true );
#ifndef VIS_C
	if ( WIFEXITED ( ret ) ) { 
		ret = WEXITSTATUS ( ret );
	}
#endif
	if ( ret != 0 && ret != 1 ) {
		ErrorHandler::genError ()->error ( "\n7za error message:\n" + get7zaError ( ret ) + "\n" );
	}
	genUnlink ( filename );
}
// 7za a -t7z files.7z *.txt
bool gen7zaCreate ( const string& archiveName, const string& filename, const string& type )
{
	string command;
#ifdef VIS_C
	command += SystemCallBinDir::instance ().getSystemCallBinDir ();
	command += SLASH;
#endif
	command += "7za";
	command += " a";
	command += " -y";
	command += " -t";
	command += type;
	command += ' ';
	command += "\"" + archiveName + "." + type + "\"";
	command += ' ';
	command += "\"" + filename + "\"";
#ifndef VIS_C
	command += " > /dev/null";
#endif
	int ret = genSystem ( command, "", true );
#ifndef VIS_C
	if ( WIFEXITED ( ret ) ) { 
		ret = WEXITSTATUS ( ret );
	}
#endif
	if ( ret != 0 && ret != 1 ) {
		ErrorHandler::genError ()->message ( "\n7za error message:\n" + get7zaError ( ret ) + "\n" );
		return false;
	}
	return true;
}
bool gen7zaCreate ( const string& filename )
{
	bool flag = gen7zaCreate ( filename, filename, "7z" );
	if ( flag ) genUnlink ( filename );
	return flag;
}
/*
http://www.rarlab.com/rar_add.htm

UNRAR 3.91 freeware      Copyright (c) 1993-2009 Alexander Roshal

Usage:     unrar <command> -<switch 1> -<switch N> <archive> <files...>
               <@listfiles...> <path_to_extract\>

<Commands>
  e             Extract files to current directory
  l[t,b]        List archive [technical, bare]
  p             Print file to stdout
  t             Test archive files
  v[t,b]        Verbosely list archive [technical,bare]
  x             Extract files with full path

<Switches>
  -             Stop switches scanning
  ac            Clear Archive attribute after compression or extraction
  ad            Append archive name to destination path
  ai            Ignore file attributes
  ap<path>      Set path inside archive
  av-           Disable authenticity verification check
  c-            Disable comments show
  cfg-          Disable read configuration
  cl            Convert names to lower case
  cu            Convert names to upper case
  dh            Open shared files
  ep            Exclude paths from names
  ep3           Expand paths to full including the drive letter
  f             Freshen files
  id[c,d,p,q]   Disable messages
  ierr          Send all messages to stderr
  inul          Disable all messages
  ioff          Turn PC off after completing an operation
  kb            Keep broken extracted files
  n<file>       Include only specified file
  n@            Read file names to include from stdin
  n@<list>      Include files listed in specified list file
  o[+|-]        Set the overwrite mode
  oc            Set NTFS Compressed attribute
  or            Rename files automatically
  ow            Save or restore file owner and group
  p[password]   Set password
  p-            Do not query password
  r             Recurse subdirectories
  ri<P>[:<S>]   Set priority (0-default,1-min..15-max) and sleep time in ms
  sl<size>      Process files with size less than specified
  sm<size>      Process files with size more than specified
  ta<date>      Process files modified after <date> in YYYYMMDDHHMMSS format
  tb<date>      Process files modified before <date> in YYYYMMDDHHMMSS format
  tn<time>      Process files newer than <time>
  to<time>      Process files older than <time>
  ts<m,c,a>[N]  Save or restore file time (modification, creation, access)
  u             Update files
  v             List all volumes
  ver[n]        File version control
  vp            Pause before each volume
  x<file>       Exclude specified file
  x@            Read file names to exclude from stdin
  x@<list>      Exclude files listed in specified list file
  y             Assume Yes on all queries
*/
string genUnrar ( const string& filename, bool deleteFile, bool createSubdirectory )
{
	string newDir;
	if ( createSubdirectory ) {
		newDir = filename.substr ( 0, filename.length () - 4 );		// Remove the .rar suffix
		genCreateNewDirectory ( newDir );							// Directory must exist before program run
	}
	string command;
#ifdef VIS_C
	command += SystemCallBinDir::instance ().getSystemCallBinDir ();
	command += SLASH;
#endif
	command += "unrar";		// unrar utility 
	command += " x";		// Extract
	command += " -inul";	// Disable messages
	command += " -p-";		// Skip password
	command += " -y";		// Always answer yes
	command += " ";
	command += "\"" + filename + "\"";
	if ( createSubdirectory ) {
		command += " ";
		command += "\"" + newDir + "\"";
	}
	int ret = genSystem ( command, "", true );
#ifndef VIS_C
	if ( WIFEXITED ( ret ) ) { 
		ret = WEXITSTATUS ( ret );
	}
#endif
	if ( deleteFile ) genUnlink ( filename );
	if ( ret != 0 ) {
		if ( createSubdirectory ) genUnlinkDirectory ( newDir );
		ErrorHandler::genError ()->error ( "\nunrar error code: " + gen_itoa ( ret ) + "\n" );
	}
	return newDir;
}
string genCommon ( const string& filename, bool deleteFile )
{
	string command;
	string newFilename = filename.substr ( 0, filename.length () - 4 );
#ifdef VIS_C
	command += SystemCallBinDir::instance ().getSystemCallBinDir ();
	command += SLASH;
#endif
	command += "common";		// unrar utility 
	command += " ";
	command += "\"" + filename + "\"";
	command += " -dmgf";		// Convert to mgf
	command += " -s255";		// Maximum spectrum size
	command += " -o";
	command += "\"" + newFilename + "\"";
	int ret = genSystem ( command, "", true );
#ifndef VIS_C
	if ( WIFEXITED ( ret ) ) { 
		ret = WEXITSTATUS ( ret );
	}
#endif
	if ( deleteFile ) genUnlink ( filename );
	if ( ret != 0 ) {
		genUnlink ( newFilename );
		ErrorHandler::genError ()->error ( "\ncommon error code: " + gen_itoa ( ret ) + "\n" );
	}
	return newFilename;
}
string genUntar ( const string& filename, bool deleteFile, bool createSubdirectory )
{
	string newDir;
	if ( createSubdirectory ) {						// To tar into a subdirectory then that has to be the working directory
		newDir = filename.substr ( 0, filename.length () - 4 );		// Remove the .tar suffix
		genCreateDirectory ( newDir, true, false );
	}
	string command;
#ifdef VIS_C
	command += "\"";
	command += "\"";
	command += ppBaseDir ();
	command += "cgi-bin";
	command += SLASH;
#endif
	command += "tar";
#ifdef VIS_C
	command += "\"";
#endif
	command += " -x";
	command += " -f ";
	command += "\"" + filename + "\"";
#ifdef VIS_C
	command += "\"";
#endif
	string previousWorkingDirectory;
	if ( createSubdirectory ) {
		previousWorkingDirectory = genCurrentWorkingDirectory ();
		genChangeWorkingDirectory ( newDir );
	}
	int ret = genSystem ( command, "", true );
#ifndef VIS_C
	if ( WIFEXITED ( ret ) ) { 
		ret = WEXITSTATUS ( ret );
	}
#endif
	if ( createSubdirectory ) {
		genChangeWorkingDirectory ( previousWorkingDirectory );
	}
	if ( deleteFile ) genUnlink ( filename );
	if ( ret != 0 ) {
		if ( createSubdirectory ) genUnlinkDirectory ( newDir );
		ErrorHandler::genError ()->error ( "\ntar error code: " + gen_itoa ( ret ) + "\n" );
	}
	return newDir;
}
string genGunzip ( const string& filename )
{
	string command;
#ifdef VIS_C
	command += SystemCallBinDir::instance ().getSystemCallBinDir ();
	command += SLASH;
#endif
	command += "gunzip";
	command += " ";
	command += "\"" + filename + "\"";
	genSystem ( command, "", false );
	string suffix = genSuffixFromPath ( filename );
	const char* csuffix = suffix.c_str ();
	if ( !genStrcasecmp ( csuffix, "gz" ) ) return filename.substr ( 0, filename.length () - 3 );
	else if ( !genStrcasecmp ( csuffix, "z" ) ) return filename.substr ( 0, filename.length () - 2 );
	else if ( !genStrcasecmp ( csuffix, "tgz" ) ) return filename.substr ( 0, filename.length () - 4 ) + ".tar";
	else if ( !genStrcasecmp ( csuffix, "taz" ) ) return filename.substr ( 0, filename.length () - 4 ) + ".tar";
	else return filename.substr ( 0, filename.length () - 4 ); // bz2
}
string genPreprocessFile ( const string& filename )
{
	string retFilename;
	string suffix = genSuffixFromPath ( filename );
	string originalSuffix = suffix;
	const char* csuffix = suffix.c_str ();
	if ( !genStrcasecmp ( csuffix, "zip" ) ) {
		retFilename = genUnzip ( filename, true, true );	// Unzip the files into a subdirectory and delete the original archive
	}
	else if ( !genStrcasecmp ( csuffix, "7z" ) ) {
		retFilename = gen7za ( filename, true, true );		// 7za the files into a subdirectory and delete the original archive
	}
	else if ( !genStrcasecmp ( csuffix, "rar" ) ) {
		retFilename = genUnrar ( filename, true, true );	// unrar the files into a subdirectory and delete the original archive
	}
	else if ( !genStrcasecmp ( csuffix, "tar" ) ) {
		retFilename = gen7za ( filename, true, true, ".tar" );	// 7za the files into a subdirectory and delete the original archive
	}
	else if ( !genStrcasecmp ( csuffix, "gz" ) || !genStrcasecmp ( csuffix, "tgz" ) || !genStrcasecmp ( csuffix, "z" ) || !genStrcasecmp ( csuffix, "taz" ) || !genStrcasecmp ( csuffix, "bz2" ) || !genStrcasecmp ( csuffix, "cmn" ) ) {
		if ( !genStrcasecmp ( csuffix, "cmn" ) ) {
			retFilename = genCommon ( filename, true );		// common the file
		}
		else {
			retFilename = genGunzip ( filename );			// gunzip the file
		}
		suffix = genSuffixFromPath ( retFilename );
		if ( !genStrcasecmp ( suffix.c_str (), "tar" ) ) {
			retFilename = gen7za ( retFilename, true, true, ".tar" );	// 7za the files into a subdirectory and delete the original archive
		}
		else {			// This is to allow an upload that was just compressed to be moved to a directory
			const char* osuffix = originalSuffix.c_str ();
			if ( !genStrcasecmp ( osuffix, "gz" ) || !genStrcasecmp ( osuffix, "z" ) || !genStrcasecmp ( osuffix, "bz2" ) || !genStrcasecmp ( osuffix, "cmn" ) ) {
				genRename ( retFilename, retFilename + "." + originalSuffix );
				retFilename +=  "." + originalSuffix;
			}
		}
	}
	else {
		retFilename = filename;		// The file is not processed
	}
	return retFilename;
}
