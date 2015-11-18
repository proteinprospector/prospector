/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_fasta.h                                                    *
*                                                                             *
*  Created    : July 20th 1996                                                *
*                                                                             *
*  Purpose    : Functions for extracting accession number, name, species, and *
*               sequence from various dialects of fasta formatted databases.  *
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

#ifndef __lu_fasta_h
#define __lu_fasta_h

#include <cstdio>
#include <string>
#include <vector>
#include <lgen_define.h>

class ParameterList;
class DatabaseEntry;

class DatabaseIndicies;
class CommentLine;
class SequenceReader;

class FastaServer {

	std::string filePath;
	std::string usedFileName;
	bool random;
	bool reverse;
	bool decoy;
	DatabaseIndicies* dbIndicies;
	CommentLine* commentLine;
	SequenceReader* sequenceReader;

	int numEntries;
	bool dna_database;

	CharPtrVector reading_frames;

	int cur;
	StringVectorSizeType line;
	void openUserDatabase ( const std::string& filename, bool createIndicies );
	void openPEFFDatabase ( const std::string& filename, bool createIndicies );
	void openDatabase ( const std::string& filename, bool createIndicies );
	std::string getUsedFilename ( const std::string& filename );
	CommentLine* getCommentLine ( const std::string& usedFileName, DatabaseIndicies* dbIndicies );
	SequenceReader* getSequenceReader ( const std::string& usedFileName, DatabaseIndicies* dbIndicies );
	void closeDatabase ();
public:
	FastaServer ( const std::string& file_name, bool createIndicies = false );
	~FastaServer ();
	void loadIntoMemory ();
	int getNumEntries () const { return numEntries; }
	std::string getFileName () const { return usedFileName; }
	bool getRandom () const { return random; }
	bool getReverse () const { return reverse; }
	bool getDecoy () const { return decoy; }
	std::string getFilePath () const { return filePath; }
	int getUnreadableSpeciesCount () const;
	bool getDNADatabase () const { return dna_database; }

	char* get_fasta_full_entry ( int serial_number, int* length );

	CharPtrVector& split_protein_to_reading_frames ( char* protein );

	char* get_fasta_protein ( int serialNumber, int dnaReadingFrame );
	char* getProtein (  const DatabaseEntry& databaseEntry );

	const char* getAccessionNumber ( int serialNumber );
	const char* getAccessionInfo ( int serialNumber );
	const char* getSpecies ( int serialNumber );
	const char* getName ( int serialNumber );
	const char* getUniprotID ( int serialNumber );

	void firstLine ( int serialNumber );
	void nextLine ();
	bool isDoneLine () const;

	std::string getLineAccessionNumber ();
	std::string getLineAccessionInfo ();
	std::string getLineSpecies ();
	std::string getLineName ();

	std::string getLineOrganismName ();			// OS
	std::string getLineGeneName ();				// GN
	std::string getLineProteinExistence ();		// PE
	std::string getLineSequenceVersion ();		// SV

	std::string getLineUniprotID ();

	time_t getDatabaseTime ();
	std::string getDatabasePath ();
	void setMaxNTermAA ( int m );
	int getMaxNTermAA () const;
};

typedef std::vector <FastaServer*> FastaServerPtrVector;

class PreloadedDatabases {
	FastaServerPtrVector fs;
public:
	PreloadedDatabases ();
	~PreloadedDatabases ();
};

void wrCommentLine ( FILE* fp, const std::string& filename, const std::string& accession_number, const std::string& name, const std::string& species );
void readProteinFromDNA ( int dnaReadingFrame, int maxNTermAA, char* fpointer, int length, char* protein, int& numUnknowns );

#endif /* ! __lu_fasta_h */
