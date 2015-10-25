/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_fas_ind.h                                                  *
*                                                                             *
*  Created    : July 19th 2000                                                *
*                                                                             *
*  Purpose    : Functions to create the index file for fasta formatted        *
*               databases.                                                    *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2000-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_fas_ind_h
#define __lu_fas_ind_h

#include <string>
#include <vector>
#include <lgen_define.h>
#include <lgen_mmap.h>

class DatabaseIndicies {
	MMapFile <char>* databaseMap;
	MMapFile <GENINT64>* commentIndexMap;
	MMapFile <GENINT64>* proteinIndexMap;
	std::string fullDatabasePath;
	std::string fullIDCPath;
	std::string fullIDPPath;
	std::string fullIDIPath;

	unsigned int maxCommentLength;
	unsigned int maxProteinLength;
	unsigned int numEntries;

	void findOffsets ();
	void readIndexFile ();
	void writeIndexFile ();
	DatabaseIndicies& operator= ( DatabaseIndicies& rhs );
	DatabaseIndicies ( const DatabaseIndicies& rhs );
public:
	DatabaseIndicies ( const std::string& filename, bool createFlag = false );
	~DatabaseIndicies ();
	char* getCommentPointer ( unsigned int serialNumber, int* length );
	char* getProteinPointer ( unsigned int serialNumber, int* length );
	char* getEntryPointer ( unsigned int serialNumber, int* length );
	unsigned int getMaxCommentLength () const { return maxCommentLength; }
	unsigned int getMaxProteinLength () const { return maxProteinLength; }
	unsigned int getNumEntries () const { return numEntries; }
};
std::string getIndexFileLastModifiedTime ( const std::string& filename );	

#endif /* ! __lu_fas_ind_h */
