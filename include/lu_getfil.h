/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_getfil.h                                                   *
*                                                                             *
*  Created    : July 20th 1996                                                *
*                                                                             *
*  Purpose    : Include file for:                                             *
*                                                                             *
*               lu_get_file.cpp                                               *
*                                                                             *
*               Utility function for reading parameter and database files.    *
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

#ifndef __lu_getfil_h
#define __lu_getfil_h

#include <string>
#include <lg_io.h>
#include <lgen_define.h>

std::string adjustPPOutputPath ( const std::string& path );

class Version {
	std::string version;
	static void getVersionParts ( const std::string& v, int& vMain, int& vSub, int& vInc );
	Version ();
public:
	static Version& instance ();
	std::string getVersion () const { return version; }
	bool isOlderVersion ( const std::string& v ) const;
	bool isNewerVersion ( const std::string& v ) const;
	static bool isOlderVersion ( const std::string& v1, const std::string& v2 );
	static bool isNewerVersion ( const std::string& v1, const std::string& v2 );
};

class MsparamsDir {
	std::string msparamsDir;
	std::string globalMSParamsDir;
	bool globalFlag;
	MsparamsDir ();
public:
	static MsparamsDir& instance ();
	std::string getParamPath ( const std::string& filename ) const;
};

class HTMLDir {
	std::string htmlDir;
	HTMLDir ();
public:
	static HTMLDir& instance ();
	std::string getVirtualHTMLDir () const;
	std::string getVirtualHTMLJavaDir () const;
	std::string getVirtualHTMLImagesDir () const;
	std::string getHTMLFullJavaDir () const;
};

class BinDir {
	std::string binDir;
	BinDir ();
public:
	static BinDir& instance ();
	std::string getBinDir () const { return binDir; }
	std::string getVirtualBinDir () const;
};

class SystemCallBinDir {
	std::string systemCallBinDir;
	SystemCallBinDir ();
public:
	static SystemCallBinDir& instance ();
	std::string getSystemCallBinDir () const { return systemCallBinDir; }
};

std::string getSystemCall ( const std::string& executable );

class SeqdbDir {
	std::string seqdbDir;
	SeqdbDir ();
public:
	static SeqdbDir& instance ();
	std::string getSeqdbDir () const { return seqdbDir; }
	StringVector getDatabaseList ( bool dbSearchFlag ) const;
	StringVector getUserDatabaseList ( bool dbSearchFlag ) const;
	std::string getDatabaseSuffix ( const std::string& databaseName ) const;
	std::string getDatabasePath ( const std::string& databaseName ) const;
	std::string getDatabasePathCreateOrAppend ( const std::string& databaseName ) const;
};

char* getParamsFileInfo ( const std::string& filename );
char* getParamsFileInfo ( const std::string& filename, int* numEntries );
char* getFileInfo ( const std::string& filename, char separator, int separatorsPerEntry, bool deleteComments );
char* getFileInfo ( const std::string& filename, char separator, int separatorsPerEntry, bool deleteComments, int* numEntries );
char* getFileAsCharPtr ( const std::string& filename );
std::string adjustedPPCurrentWorkingDirectory ( const std::string& path );
std::string ppBaseDir ();

class PPTempFile {
	std::string fullPathDir;
	std::string directory;
	std::string adjustedDir;
	std::string filename;
public:
	PPTempFile ( const std::string& directory, const std::string& suffix, bool dummy );
	PPTempFile ( const std::string& prefix, const std::string& suffix );
	std::string getAdjustedPath () const { return adjustedDir + filename; }
	std::string getFullPath () const { return fullPathDir + SLASH + filename; }
	std::string getFullPathDir () const { return fullPathDir; }
	std::string getURL () const;
};

class InfoParams : public GenNameValueStream {
	InfoParams ();
	std::string getPlatformValue ( const std::string& param ) const;
public:
	static InfoParams& instance ( bool reset = false );
	static time_t getLastModifiedTime ();
	std::string getSeqdbDir () const		{ return getPlatformValue ( "seqdb" ); }
	std::string getUploadTemp () const		{ return getPlatformValue ( "upload_temp" ); }
	std::string getRCommand () const		{ return getPlatformValue ( "r_command" ); }
	std::string getCentroidDir () const		{ return getPlatformValue ( "centroid_dir" ); }
	std::string getRawDir () const			{ return getPlatformValue ( "raw_dir" ); }
	std::string getUserRepository () const	{ return getPlatformValue ( "user_repository" ); }
};

void pidLogOutput ( const std::string& message, bool showTime = true );

#endif /* ! __lu_getfil_h */
