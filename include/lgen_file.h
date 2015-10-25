/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lgen_file.h                                                   *
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

#ifndef __lgen_file_h
#define __lgen_file_h

#include <string>
#include <ctime>
#include <lgen_define.h>

GENINT64 genFileSize ( const std::string& full_name );
GENINT64 genDirectorySize ( const std::string& name );
std::string stripFilenameChars ( const std::string& filename );
bool isFileType ( const std::string& filename, const std::string& type );
bool genFileExists ( const std::string& full_name );
time_t genLastModifyTime ( const std::string& full_name );
void genCreateDotExtensionFileList ( char*** name_list, int* num_names, const std::string& directory, const std::string& extension );
void genFreeNameList ( char** name_list, int num_names );
double genFreeDiskSpace ( const std::string& fileName );
double genVolumeSize ( const std::string& fileName );
bool isTextFile ( const std::string& fileName, int maxLen );
bool checkForCR ( const std::string& file_name );
bool checkForCR ( const std::string& file_name, int maxLen );
bool macTextFile ( const std::string& file_name );
bool macTextFile ( const std::string& file_name, int maxLen );
void dos2Unix ( const std::string& name1 );
void dos2Unix ( const std::string& name1, const std::string& name2 );
void mac2Unix ( const std::string& name1 );
void mac2Unix ( const std::string& name1, const std::string& name2 );
bool convertMacFiles ( const std::string& name );
void convertToUnixTextWithMessage ( const std::string& name );
void convertToUnixText ( const std::string& name );
void copyFileWithMessage ( const std::string& name1, const std::string& name2 );
void copyFile ( const std::string& name1, const std::string& name2 );
bool genIsFullPath ( const std::string& path );
std::string genSuffixFromPath ( const std::string& path );
std::string genFilenameFromPath ( const std::string& path );
std::string genShortFilenameFromPath ( const std::string& path );
std::string genDirectoryFromPath ( const std::string& path );
std::string genCurrentWorkingDirectory ();
void genChangeWorkingDirectory ( const std::string& workingDirectory );
#ifdef VIS_C
std::string genCurrentDOSDrive ();
#endif

class FileList {
	std::string directory;
	StringVector nameList;
#ifdef VIS_C
void addNames ( const std::string& filespec, const std::string& extensionToDelete = "" );
void addSubDirectoryNames ( const std::string& filespec );
#endif
public:
	FileList ( const std::string& directory, const std::string& prefix, const std::string& extension, bool deleteExtension );
	FileList ( const std::string& directory, const std::string& prefix = "" );
	FileList ( const std::string& directory, const StringVector& regExpList );
	typedef StringVector::size_type size_type;
	StringVector& getNameList () { return nameList; }
	void unlink () const;
	std::string& operator [] ( size_type i ) { return nameList [i]; }
	const std::string& operator [] ( size_type i ) const { return nameList [i]; }
	size_type size () const { return nameList.size (); }
};
std::string genNextFreeFileName ( const std::string& directory, const std::string& prefix, const std::string& suffix );
bool genCreateDirectory ( const std::string& dir );
bool genCreateNewDirectory ( const std::string& dir );
void genCreateDirectory ( const std::string& directory, bool errorIfPresent, bool deleteFiles );
std::string genMakeSubDirectory ( const std::string& dir, const std::string& subDir );
void genCreateDirectoryPath ( const std::string& ssName );
void genRmdir ( const std::string& dir );
void genUnlink ( const std::string& file );
void genUnlink ( const StringVector& files );
void genUnlinkDirectory ( const std::string& name );
bool genIsDirectory ( const std::string& filenameFullPath );
std::string genTranslateSlash ( const std::string& s );
void genCopyDirectoryContents ( const std::string& oldName, const std::string& newName, bool message = false );
int genRename ( const std::string& oldName, const std::string& newName, bool message = false );
bool genDifferentPartitions ( const std::string& file1, const std::string& file2 );
StringVector getRepositoryKeyList ( const std::string& d, int keyLen );

#endif /* ! __lgen_file_h */
