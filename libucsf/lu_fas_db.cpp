/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_fas_db.cpp                                                 *
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
#ifdef VIS_C
#pragma warning( disable : 4503 )	// Disable warnings about decorated name length being too long
#endif
#ifndef VIS_C
#include <stdexcept>
#endif
#include <lg_new.h>
#include <lg_string.h>
#include <lgen_file.h>
#include <lgen_error.h>
#include <lu_fas_ind.h>
#include <lu_fasta.h>
#include <lu_check_db.h>
#include <lu_param_list.h>
#include <lu_html_form.h>
#include <lu_cgi_val.h>
#include <lu_db_entry.h>
using std::string;
using std::runtime_error;

#define MAX_NON_PRINT 31

class CommentLine {
protected:
	static const char* unreadable;
	DatabaseIndicies* dbIndicies;
	char* top_line_start;
	char* point;
	char* accessionNumber;
	char* accessionInfo;
	char* name;
	char* species;
	char* uniprotID;
	char* uniprotInfo;
	StringVector accessionNumberList;
	StringVector accessionInfoList;
	StringVector nameList;
	StringVector speciesList;
	StringVector uniprotIDList;
	StringVector uniprotInfoList;
	int unreadable_species_count;
	void genpeptLine ( bool speciesUndefined );
	void uniprotLine ();
	void spaceDashSpaceLine ( bool speciesUndefined );
	void get_default_top_line ( int serialNumber );
	void unreadableSpecies ();
	void copyTo ( char* to, int delimiter );
	void copyToStripTrailing ( char* to, int delimiter );
	void copyToUpper ( char* to, char* start, int len );
	void goPastNext ( int delimiter );
	void advance1 ();
	void advance4 ();
	void advance5 ();
	bool checkChar ( int c );
	void setNumericAccessionNumber ();
	void setAccessionNumber ( int delimiter );
	void setAccessionNumber ( int delimiter1, int delimiter2 );
	void setSpecies ( int delimiter );
	void setUniprotID ();
	void setName ( int delimiter );
	void setName ( int delimiter1, int delimiter2 );
	void copyToUpperStripTrailing ( char* to, int delimiter );
	void copyTo ( char* to, int delimiter1, int delimiter2 );
	void copyTo ( char* to, char* start, int len );
	void copyToPtr ( char* to, char* start, char* end );
	void copyNumericTo ( char* to );
	char* copyToCheckUnderscore ( char* to, int delimiter );
	friend class FastaServer;
public:
	CommentLine ( DatabaseIndicies* dbIndicies );
	virtual ~CommentLine ();
	virtual void getCommentLine ( int serialNumber, bool allLines = false ) = 0;
	virtual void getAccessionInfo ( int serialNumber, bool allLines = false );
	int getUnreadableSpeciesCount () const { return unreadable_species_count; }
	const char* getAccessionNumber () const { return accessionNumberList.empty () ? accessionNumber : accessionNumberList [0].c_str (); }
	const char* getAccessionInfo () const { return accessionInfoList.empty () ? accessionInfo : accessionInfoList [0].c_str (); }
	const char* getName () const { return nameList.empty () ? name : nameList [0].c_str (); }
	const char* getSpecies () const { return speciesList.empty () ? species : speciesList [0].c_str (); }
	const char* getUniprotID () const { return uniprotIDList.empty () ? uniprotID : uniprotIDList [0].c_str (); }
};
void CommentLine::getAccessionInfo ( int serialNumber, bool allLines )
{
	int length;
	if ( !accessionInfoList.empty () ) accessionInfoList.clear ();
	point = dbIndicies->getCommentPointer ( serialNumber, &length );
	do {
		point++;
		int i = 0;
		for ( ;  ; ) {
			char c = *point;
			if ( c == ' ' ) break;
			if ( c <= MAX_NON_PRINT ) break;
			accessionInfo [i++] = c;
			point++;
		}
		if ( accessionInfo [i-1] == '|' ) accessionInfo [i-1] = 0;
		else accessionInfo [i] = 0;
		for ( ; ; ) {						// Find the end of the comment line
			char c = *point;
			if ( c == '\n' || c == 1 ) break;
			point++;
		}
		if ( allLines ) accessionInfoList.push_back ( accessionInfo );
		else break;
	} while ( *point != '\n' );
	goPastNext ( '\n' );
}

class GenpeptCommentLine : public CommentLine {
public:
	GenpeptCommentLine ( DatabaseIndicies* dbIndicies );
	~GenpeptCommentLine ();
	void getCommentLine ( int serialNumber, bool allLines = false );
};

class SwissProtCommentLine : public CommentLine {
	bool checkStr2Advance ( const char* str );
public:
	SwissProtCommentLine ( DatabaseIndicies* dbIndicies );
	~SwissProtCommentLine ();
	void getCommentLine ( int serialNumber, bool allLines = false );
};

class LudwigCommentLine : public CommentLine {
public:
	LudwigCommentLine ( DatabaseIndicies* dbIndicies );
	~LudwigCommentLine ();
	void getCommentLine ( int serialNumber, bool allLines = false );
};

class OwlCommentLine : public CommentLine {
	bool setAccessionNumberCheckUnderscore ( int delimiter );
public:
	OwlCommentLine ( DatabaseIndicies* dbIndicies );
	~OwlCommentLine ();
	void getCommentLine ( int serialNumber, bool allLines = false );
};

class NCBINRCommentLine : public CommentLine {
	bool checkCharAdvance ( int c );
public:
	NCBINRCommentLine ( DatabaseIndicies* dbIndicies );
	~NCBINRCommentLine ();
	void getCommentLine ( int serialNumber, bool allLines = false );
};

class IPICommentLine : public CommentLine {
public:
	IPICommentLine ( DatabaseIndicies* dbIndicies );
	~IPICommentLine ();
	void getCommentLine ( int serialNumber, bool allLines = false );
};

class UniprotCommentLine : public CommentLine {
public:
	UniprotCommentLine ( DatabaseIndicies* dbIndicies );
	~UniprotCommentLine ();
	void getCommentLine ( int serialNumber, bool allLines = false );
};

class DBESTCommentLine : public CommentLine {
	char** dbest_species_list;
	char** uc_dbest_species_list;
	int num_dbest_species;
	void read_dbest_species_list_file ();
public:
	DBESTCommentLine ( DatabaseIndicies* dbIndicies );
	~DBESTCommentLine ();
	void getCommentLine ( int serialNumber, bool allLines = false );
};

class GenericCommentLine : public CommentLine {
	void skipSpace ();
	void skipChar ( int c );
	bool lineEnd ();
public:
	GenericCommentLine ( DatabaseIndicies* dbIndicies );
	~GenericCommentLine ();
	void getCommentLine ( int serialNumber, bool allLines = false );
	void getAccessionInfo ( int serialNumber, bool allLines = false )
	{
		accessionInfoList.clear ();
		accessionInfoList.push_back ( "" );
		if ( accessionInfoList.size () < accessionNumberList.size () ) accessionInfoList.resize ( accessionNumberList.size () );
	}
};

class DefaultCommentLine : public CommentLine {
public:
	DefaultCommentLine ( DatabaseIndicies* dbIndicies );
	~DefaultCommentLine ();
	void getCommentLine ( int serialNumber, bool allLines = false );
	void getAccessionInfo ( int serialNumber, bool allLines = false )
	{
		accessionInfoList.clear ();
		accessionInfoList.push_back ( "" );
		if ( accessionInfoList.size () < accessionNumberList.size () ) accessionInfoList.resize ( accessionNumberList.size () );
	}
};

class SequenceReader {
protected:
	DatabaseIndicies* dbIndicies;
	char* protein;
	int numUnknowns;
	int maxNTermAA;
public:
	SequenceReader ( DatabaseIndicies* dbIndicies );
	virtual ~SequenceReader ();
	virtual void readProtein ( int serial_number, int dnaReadingFrame ) = 0;
	char* getProtein () const { return protein; }
	int getNumUnknowns () const { return numUnknowns; }
	void setMaxNTermAA ( int m ) { maxNTermAA = m; }
	int getMaxNTermAA () const { return maxNTermAA; }
};

class ProteinSequenceReader : public SequenceReader {
public:
	ProteinSequenceReader ( DatabaseIndicies* dbIndicies );
	~ProteinSequenceReader ();
	void readProtein ( int serial_number, int dnaReadingFrame );
};

class DNASequenceReader : public SequenceReader {
public:
	DNASequenceReader ( DatabaseIndicies* dbIndicies );
	~DNASequenceReader ();
	void readProtein ( int serial_number, int dnaReadingFrame );
};

class DNAProteinSequenceReader : public SequenceReader {
public:
	DNAProteinSequenceReader ( DatabaseIndicies* dbIndicies );
	~DNAProteinSequenceReader ();
	void readProtein ( int serial_number, int dnaReadingFrame );
};

FastaServer::FastaServer ( const string& filename, bool createIndicies )
{
	if ( filename.find ( "UserProtein" ) != string::npos )
		openUserDatabase ( filename, createIndicies );
	else
		openDatabase ( filename, createIndicies );
	random	= ( usedFileName.find ( ".random" ) != string::npos );
	reverse	= ( usedFileName.find ( ".reverse" ) != string::npos );
	decoy = random || reverse;
}
FastaServer::~FastaServer ()
{
	closeDatabase ();
}
void FastaServer::openUserDatabase ( const string& filename, bool createIndicies )
{
	dna_database = false;
	usedFileName = filename.substr ( filename.find ( "UserProtein" ) );
	filePath = filename;

	dbIndicies = new DatabaseIndicies ( filename, createIndicies );
	numEntries = dbIndicies->getNumEntries ();

	commentLine = getCommentLine ( usedFileName, dbIndicies );

	sequenceReader = getSequenceReader ( usedFileName, dbIndicies );

	cur = -1;
}
void FastaServer::openDatabase ( const string& filename, bool createIndicies )
{
	if ( getDatabasePrefix ( filename ).empty () ) {
		ErrorHandler::genError ()->error ( "Invalid database filename prefix.\nSee the manual for information on naming databases." );
	}
	dna_database = is_dna_database ( filename );
	usedFileName = getUsedFilename ( filename );

	dbIndicies = new DatabaseIndicies ( usedFileName, createIndicies );
	numEntries = dbIndicies->getNumEntries ();

	commentLine = getCommentLine ( usedFileName, dbIndicies );

	sequenceReader = getSequenceReader ( usedFileName, dbIndicies );

	cur = -1;
}
string FastaServer::getUsedFilename ( const string& filename )
{
	if ( !genFileExists ( SeqdbDir::instance ().getDatabasePath ( filename ) ) ) {
		string substituteFileName = getBestSubstituteDatabase ( filename );
		if ( substituteFileName.empty () ) {
			ErrorHandler::genError ()->error ( "The database you have selected does not exist on the server.\nPerhaps the databases have been recently updated." );
			return "";
		}
		else
			return substituteFileName;
	}
	else
		return filename;
}
CommentLine* FastaServer::getCommentLine ( const string& usedFileName, DatabaseIndicies* dbIndicies )
{
	if		( is_genpept_database ( usedFileName ) )	return new GenpeptCommentLine ( dbIndicies );
	else if	( is_swissprot_database	( usedFileName ) )	return new SwissProtCommentLine ( dbIndicies );
	else if	( is_ludwig_database ( usedFileName ) )		return new LudwigCommentLine ( dbIndicies );
	else if	( is_owl_database ( usedFileName ) )		return new OwlCommentLine ( dbIndicies );
	else if	( is_ncbinr_database ( usedFileName ) )		return new NCBINRCommentLine ( dbIndicies );
	else if	( is_ipi_database ( usedFileName ) )		return new IPICommentLine ( dbIndicies );
	else if	( is_uniprot_database ( usedFileName ) )	return new UniprotCommentLine ( dbIndicies );
	else if	( is_est_database ( usedFileName ) )		return new DBESTCommentLine ( dbIndicies );
	else if ( is_generic_database ( usedFileName ) )	return new GenericCommentLine ( dbIndicies );
	else if ( is_default_database ( usedFileName ) )	return new DefaultCommentLine ( dbIndicies );
	else {
		ErrorHandler::genError ()->error ( "Invalid database filename prefix.\nSee the manual for information on naming databases." );
		return 0;
	}
}
SequenceReader* FastaServer::getSequenceReader ( const string& usedFileName, DatabaseIndicies* dbIndicies )
{
	if ( is_dna_format_database ( usedFileName ) )		return new DNASequenceReader ( dbIndicies );
	else if ( is_pdna_format_database ( usedFileName ) )return new DNAProteinSequenceReader ( dbIndicies );
	else												return new ProteinSequenceReader ( dbIndicies );
}
void FastaServer::closeDatabase ()
{
	delete dbIndicies;
	delete commentLine;
	delete sequenceReader;
}
void FastaServer::loadIntoMemory ()
{
	for ( int i = 1 ; i <= numEntries ; i++ ) {
		for ( firstLine ( i ) ; isDoneLine () ; nextLine () );
	}
}
int FastaServer::getUnreadableSpeciesCount () const
{
	return ( commentLine->getUnreadableSpeciesCount () );
}
char* FastaServer::get_fasta_full_entry ( int serialNumber, int* length )
{
	return ( dbIndicies->getEntryPointer ( serialNumber, length ) );
}
CharPtrVector& FastaServer::split_protein_to_reading_frames ( char* protein )
{
	reading_frames.resize ( 0 );
	if ( dna_database ) {
		reading_frames.push_back ( protein );
		bool dot = false;
		for ( ; ;  ) {
			char val = *protein;
			if ( val == 0 ) break;
			if ( val == '.' ) {
				*protein = 0;
				dot = true;
			}
			else {
				if ( dot ) {
					reading_frames.push_back ( protein );
					dot = false;
				}
			}
			protein++;
		}
	}
	else {
		reading_frames.push_back ( protein );
	}
	return ( reading_frames );
}

const char* FastaServer::getAccessionNumber ( int serialNumber )
{
	if ( cur != serialNumber ) {
		cur = serialNumber;
		commentLine->getCommentLine ( serialNumber );
	}
	return commentLine->getAccessionNumber ();
}
const char* FastaServer::getAccessionInfo ( int serialNumber )
{
	return commentLine->getAccessionInfo ();
}
const char* FastaServer::getSpecies ( int serialNumber )
{
	if ( cur != serialNumber ) {
		cur = serialNumber;
		commentLine->getCommentLine ( serialNumber );
	}
	return commentLine->getSpecies ();
}
const char* FastaServer::getUniprotID ( int serialNumber )
{
	if ( cur != serialNumber ) {
		cur = serialNumber;
		commentLine->getCommentLine ( serialNumber );
	}
	return commentLine->getUniprotID ();
}
const char* FastaServer::getName ( int serialNumber )
{
	if ( cur != serialNumber ) {
		cur = serialNumber;
		commentLine->getCommentLine ( serialNumber );
	}
	return commentLine->getName ();
}
void FastaServer::firstLine ( int serialNumber )
{
	line = 0;
	cur = serialNumber;
	commentLine->getCommentLine ( serialNumber, true );
	commentLine->getAccessionInfo ( serialNumber, true );
}
void FastaServer::nextLine ()
{
	line++;
}
bool FastaServer::isDoneLine () const
{
	return line < commentLine->accessionNumberList.size ();
}
string FastaServer::getLineAccessionNumber ()
{
	return commentLine->accessionNumberList [line];
}
string FastaServer::getLineAccessionInfo ()
{
	return commentLine->accessionInfoList [line];
}
string FastaServer::getLineSpecies ()
{
	return commentLine->speciesList [line];
}
string FastaServer::getLineName ()
{
	return commentLine->nameList [line];
}
string FastaServer::getLineUniprotID ()
{
	if ( commentLine->uniprotIDList.size () ) {
		return commentLine->uniprotIDList [line];
	}
	else return "";
}
/*
>db|UniqueIdentifier|EntryName ProteinName OS=OrganismName[ GN=GeneName]PE=ProteinExistence SV=SequenceVersion
>sp|Q4U9M9|104K_THEAN 104 kDa microneme/rhoptry antigen OS=Theileria annulata GN=TA08425 PE=3 SV=1
*/
string FastaServer::getLineOrganismName ()
{
	if ( commentLine->uniprotInfoList.size () ) {
		string s1 = commentLine->uniprotInfoList [line];
		size_t p1 = s1.find ( "OS=" );
		if ( p1 != string::npos ) {
			p1 += 3;
			size_t p2 = s1.find ( " GN=", p1 );
			if ( p2 != string::npos )		return s1.substr ( p1, p2-p1 );
			else {
				size_t p3 = s1.find ( " PE=", p1 );
				if ( p3 != string::npos )	return s1.substr ( p1, p3-p1 );
				else						return "";
			}
		}
		else return "";
	}
	else return "";
}
string FastaServer::getLineGeneName ()
{
	if ( commentLine->uniprotInfoList.size () ) {
		string s1 = commentLine->uniprotInfoList [line];
		size_t p1 = s1.find ( " GN=" );
		size_t p2 = s1.find ( " PE=" );
		if ( p1 != string::npos && p2 != string::npos ) {
			p1 += 4;
			return s1.substr ( p1, p2-p1 );
		}
		else return "";
	}
	else return "";
}
string FastaServer::getLineProteinExistence ()
{
	if ( commentLine->uniprotInfoList.size () ) {
		string s1 = commentLine->uniprotInfoList [line];
		size_t p1 = s1.find ( " PE=" );
		if ( p1 != string::npos ) {
			p1 += 4;
			size_t p2 = s1.find ( " SV=", p1 );
			if ( p2 != string::npos ) {
				string pExist = s1.substr ( p1, p2-p1 );
				if ( pExist == "1" )		return "Evidence at protein level";
				else if ( pExist == "2" )	return "Evidence at transcript level";
				else if ( pExist == "3" )	return "Inferred from homology";
				else if ( pExist == "4" )	return "Predicted";
				else if ( pExist == "5" )	return "Uncertain";
				else						return "";
			}
			else return "";
		}
		else return "";
	}
	else return "";
}
string FastaServer::getLineSequenceVersion ()
{
	if ( commentLine->uniprotInfoList.size () ) {
		string s1 = commentLine->uniprotInfoList [line];
		size_t p1 = s1.find ( " SV=" );
		if ( p1 != string::npos ) {
			p1 += 4;
			return s1.substr ( p1 );
		}
		else return "";
	}
	else return "";
}
char* FastaServer::get_fasta_protein ( int serialNumber, int dnaReadingFrame )
{
	sequenceReader->readProtein ( serialNumber, dnaReadingFrame );
	return ( sequenceReader->getProtein () );
}
char* FastaServer::getProtein ( const DatabaseEntry& databaseEntry )
{
	int serialNumber = databaseEntry.getIndex ();
	sequenceReader->readProtein ( serialNumber, databaseEntry.getDNAReadingFrame () );
	char* seq = sequenceReader->getProtein ();
	if ( dna_database ) {
		CharPtrVectorSizeType openReadingFrame = databaseEntry.getOpenReadingFrame ();
		CharPtrVector readingFrames = split_protein_to_reading_frames ( seq );
		if ( openReadingFrame >= readingFrames.size () ) {
			string err;
			err += "Invalid open reading frame ";
			err += gen_itoa ( openReadingFrame );
			err += ".\n";
			err += "Index Number ";
			err += gen_itoa ( serialNumber );
			err += " ";
			err += "(DNA reading frame ";
			err += gen_itoa ( databaseEntry.getDNAReadingFrame () );
			err += ") has ";
			err += gen_itoa ( readingFrames.size () );
			err += " open reading frames.";
			ErrorHandler::genError ()->error ( err );
		}
		return ( readingFrames [openReadingFrame] );
	}
	else
		return seq;
}
time_t FastaServer::getDatabaseTime ()
{
	return genLastModifyTime ( getDatabasePath () );
}
string FastaServer::getDatabasePath ()
{
	return SeqdbDir::instance ().getDatabasePath ( usedFileName );
}
void FastaServer::setMaxNTermAA ( int m )
{
	sequenceReader->setMaxNTermAA ( m );
}
int FastaServer::getMaxNTermAA () const
{
	return sequenceReader->getMaxNTermAA ();
}
SequenceReader::SequenceReader ( DatabaseIndicies* dbIndicies )
{
	maxNTermAA = 0;						// This is the default value
	this->dbIndicies = dbIndicies;
	unsigned int maxProteinLength = dbIndicies->getMaxProteinLength ();
	protein = new char [maxProteinLength];
}
SequenceReader::~SequenceReader ()
{
	delete [] protein;
}
ProteinSequenceReader::ProteinSequenceReader ( DatabaseIndicies* dbIndicies ) :
	SequenceReader ( dbIndicies ) {}

ProteinSequenceReader::~ProteinSequenceReader () {}

void ProteinSequenceReader::readProtein ( int serialNumber, int dnaReadingFrame )
{
	int length;
	char* fpointer = dbIndicies->getProteinPointer ( serialNumber, &length );
	char* fpointer_next = fpointer + length;
	char* cppointer = protein;
	for ( ; ; ) {
		if ( fpointer == fpointer_next ) break;
		char val = *fpointer++;

		if ( val > MAX_NON_PRINT ) {
//			if ( val == 'X' ) numUnknowns++;	Taken out for speed reasons
			*cppointer++ = val;
		}
	}
	if ( maxNTermAA )	*(protein+maxNTermAA) = 0;
	else				*cppointer = 0;
}
DNASequenceReader::DNASequenceReader ( DatabaseIndicies* dbIndicies ) :
	SequenceReader ( dbIndicies ) {}
DNASequenceReader::~DNASequenceReader () {}

void readProteinFromDNA ( int dnaReadingFrame, int maxNTermAA, char* fpointer, int length, char* protein, int& numUnknowns )
{
	char* sequence_end = fpointer + length;
	char* sequence_start = fpointer;
	char aa;
	char p1, p2, p3;
	numUnknowns = 0;
	char* cppointer = protein;
	if ( length == 0 ) {
		*cppointer = 0;
		return;
	}
	int direction = ( dnaReadingFrame > 3 );
	int reading_frame = ( direction ) ? dnaReadingFrame - 4 : dnaReadingFrame - 1;

	if ( direction ) {	/* Backwards */
		char* backwards_read_end = sequence_start - 1;
		fpointer = sequence_end - 1;
		while ( *fpointer <= MAX_NON_PRINT || *fpointer == '>' ) fpointer--;
		for ( int i = 0 ; i < reading_frame ; i++ ) {
			while ( *fpointer <= MAX_NON_PRINT ) fpointer--;
			fpointer--;
		}
		for ( ; ; ) {
			do {
				if ( fpointer == backwards_read_end ) goto label;
				p1 = *fpointer--;
			} while ( p1 <= MAX_NON_PRINT );
			do {
				if ( fpointer == backwards_read_end ) goto label;
				p2 = *fpointer--;
			} while ( p2 <= MAX_NON_PRINT );
			do {
				if ( fpointer == backwards_read_end ) goto label;
				p3 = *fpointer--;
			} while ( p3 <= MAX_NON_PRINT );
			switch ( p1 ) {
				case 'T':
					switch ( p2 ) {
						case 'T':
							switch ( p3 ) {
								case 'T':case 'C':
									aa = 'K';
									goto aa_set1;
								case 'A':case 'G':case 'R':
									aa = 'N';goto aa_set1;
								case 'Y':
									aa = 'K';goto aa_set1;
								default:
									aa = 'X';goto aa_set1;
							}
						case 'C':
							switch ( p3 ) {
								case 'T':case 'C':
									aa = 'R';goto aa_set1;
								case 'A':case 'G':case 'R':
									aa = 'S';goto aa_set1;
								case 'Y':
									aa = 'R';goto aa_set1;
								default:
									aa = 'X';goto aa_set1;
							}
						case 'A':
							switch ( p3 ) {
								case 'T':case 'A':case 'G':
									aa = 'I';goto aa_set1;
								case 'C':
									aa = 'M';goto aa_set1;
								case 'D':case 'K':case 'R':case 'W':
									aa = 'I';goto aa_set1;
								default:
									aa = 'X';goto aa_set1;
							}
						case 'G':
							aa = 'T';goto aa_set1;
						default:
							aa = 'X';goto aa_set1;
					}
				case 'C':
					switch ( p2 ) {
						case 'T':
							switch ( p3 ) {
								case 'T':case 'C':
									aa = 'E';goto aa_set1;
								case 'A':case 'G':case 'R':
									aa = 'D';goto aa_set1;
								case 'Y':
									aa = 'E';goto aa_set1;
								default:
									aa = 'X';goto aa_set1;
							}
						case 'C':
							aa = 'G';goto aa_set1;
						case 'A':
							aa = 'V';goto aa_set1;
						case 'G':
							aa = 'A';goto aa_set1;
						default:
							aa = 'X';goto aa_set1;
					}
				case 'A':
					switch ( p2 ) {
						case 'T':
							switch ( p3 ) {
								case 'T':case 'C':
									aa = '.';goto aa_set1;
								case 'A':case 'G':case 'R':
									aa = 'Y';goto aa_set1;
								case 'Y':
									aa = '.';goto aa_set1;
								default:
									aa = 'X';goto aa_set1;
							}
						case 'C':
							switch ( p3 ) {
								case 'T':
									aa = '.';goto aa_set1;
								case 'C':
									aa = 'W';goto aa_set1;
								case 'A':case 'G':case 'R':
									aa = 'C';goto aa_set1;
								default:
									aa = 'X';goto aa_set1;
							}
						case 'A':
							switch ( p3 ) {
								case 'T':case 'C':
									aa = 'L';goto aa_set1;
								case 'A':case 'G':case 'R':
									aa = 'F';goto aa_set1;
								case 'Y':
									aa = 'L';goto aa_set1;
								default:
									aa = 'X';goto aa_set1;
							}
						case 'G':
							aa = 'S';goto aa_set1;
						case 'Y':
							switch ( p3 ) {
								case 'T':
									aa = '.';goto aa_set1;
								default:
									aa = 'X';goto aa_set1;
							}
						default:
							aa = 'X';goto aa_set1;
					}
				case 'G':
					switch ( p2 ) {
						case 'T':
							switch ( p3 ) {
								case 'T':case 'C':
									aa = 'Q';goto aa_set1;
								case 'A':case 'G':case 'R':
									aa = 'H';goto aa_set1;
								case 'Y':
									aa = 'Q';goto aa_set1;
								default:
									aa = 'X';goto aa_set1;
							}
						case 'C':
							aa = 'R';goto aa_set1;
						case 'A':
							aa = 'L';goto aa_set1;
						case 'G':
							aa = 'P';goto aa_set1;
						default:
							aa = 'X';goto aa_set1;
					}
				case 'K':
					switch ( p2 ) {
						case 'C':
							switch ( p3 ) {
								case 'T':case 'C':case 'Y':
									aa = 'R';goto aa_set1;
								default:
									aa = 'X';goto aa_set1;
							}
						default:
							aa = 'X';goto aa_set1;
					}
				case 'R':
					switch ( p2 ) {
						case 'A':
							switch ( p3 ) {
								case 'T':case 'C':case 'Y':
									aa = 'L';goto aa_set1;
								default:
									aa = 'X';goto aa_set1;
							}
						default:
							aa = 'X';goto aa_set1;
					}
				default:
					aa = 'X';goto aa_set1;
			}
			aa_set1:;
			*cppointer++ = aa;
			if ( aa == 'X' ) numUnknowns++;
		}
	}
	else {
		fpointer = sequence_start;
		for ( int i = 0 ; i < reading_frame ; i++ ) {
			while ( *fpointer <= MAX_NON_PRINT ) fpointer++;
			fpointer++;
		}
		for ( ; ; ) {
			do {
				if ( fpointer == sequence_end ) goto label;
				p1 = *fpointer++;
			} while ( p1 <= MAX_NON_PRINT );
			do {
				if ( fpointer == sequence_end ) goto label;
				p2 = *fpointer++;
			} while ( p2 <= MAX_NON_PRINT );
			do {
				if ( fpointer == sequence_end ) goto label;
				p3 = *fpointer++;
			} while ( p3 <= MAX_NON_PRINT );
			switch ( p1 ) {
				case 'T':
					switch ( p2 ) {
						case 'T':
							switch ( p3 ) {
								case 'T':case 'C':
									aa = 'F';goto aa_set2;
								case 'A':case 'G':case 'R':
									aa = 'L';goto aa_set2;
								case 'Y':
									aa = 'F';goto aa_set2;
								default:
									aa = 'X';goto aa_set2;
							}
						case 'C':
							aa = 'S';goto aa_set2;
						case 'A':
							switch ( p3 ) {
								case 'T':case 'C':
									aa = 'Y';goto aa_set2;
								case 'A':case 'G':case 'R':
									aa = '.';goto aa_set2;
								case 'Y':
									aa = 'Y';goto aa_set2;
								default:
									aa = 'X';goto aa_set2;
							}
						case 'G':
							switch ( p3 ) {
								case 'T':case 'C':
									aa = 'C';goto aa_set2;
								case 'A':
									aa = '.';goto aa_set2;
								case 'G':
									aa = 'W';goto aa_set2;
								case 'Y':
									aa = 'C';goto aa_set2;
								default:
									aa = 'X';goto aa_set2;
							}
						case 'R':
							switch ( p3 ) {
								case 'A':
									aa = '.';goto aa_set2;
								default:
									aa = 'X';goto aa_set2;
							}
						default:
							aa = 'X';goto aa_set2;
					}
				case 'C':
					switch ( p2 ) {
						case 'T':
							aa = 'L';goto aa_set2;
						case 'C':
							aa = 'P';goto aa_set2;
						case 'A':
							switch ( p3 ) {
								case 'T':case 'C':
									aa = 'H';goto aa_set2;
								case 'A':case 'G':case 'R':
									aa = 'Q';goto aa_set2;
								case 'Y':
									aa = 'H';goto aa_set2;
								default:
									aa = 'X';goto aa_set2;
							}
						case 'G':
							aa = 'R';goto aa_set2;
						default:
							aa = 'X';goto aa_set2;
					}
				case 'A':
					switch ( p2 ) {
						case 'T':
							switch ( p3 ) {
								case 'T':case 'C':case 'A':
									aa = 'I';goto aa_set2;
								case 'G':
									aa = 'M';goto aa_set2;
								case 'H':case 'M':case 'W':case 'Y':
									aa = 'I';goto aa_set2;
								default:
									aa = 'X';goto aa_set2;
							}
						case 'C':
							aa = 'T';goto aa_set2;
						case 'A':
							switch ( p3 ) {
								case 'T':case 'C':
									aa = 'N';goto aa_set2;
								case 'A':case 'G':case 'R':
									aa = 'K';goto aa_set2;
								case 'Y':
									aa = 'N';goto aa_set2;
								default:
									aa = 'X';goto aa_set2;
							}
						case 'G':
							switch ( p3 ) {
								case 'T':case 'C':
									aa = 'S';goto aa_set2;
								case 'A':case 'G':case 'R':
									aa = 'R';goto aa_set2;
								case 'Y':
									aa = 'S';goto aa_set2;
								default:
									aa = 'X';goto aa_set2;
							}
						default:
							aa = 'X';goto aa_set2;
					}
				case 'G':
					switch ( p2 ) {
						case 'T':
							aa = 'V';goto aa_set2;
						case 'C':
							aa = 'A';goto aa_set2;
						case 'A':
							switch ( p3 ) {
								case 'T':case 'C':
									aa = 'D';goto aa_set2;
								case 'A':case 'G':case 'R':
									aa = 'E';goto aa_set2;
								case 'Y':
									aa = 'D';goto aa_set2;
								default:
									aa = 'X';goto aa_set2;
							}
						case 'G':
							aa = 'G';goto aa_set2;
						default:
							aa = 'X';goto aa_set2;
					}
				case 'M':
					switch ( p2 ) {
						case 'G':
							switch ( p3 ) {
								case 'A':case 'G':case 'R':
									aa = 'R';goto aa_set2;
								default:
									aa = 'X';goto aa_set2;
							}
						default:
							aa = 'X';goto aa_set2;
					}
				case 'Y':
					switch ( p2 ) {
						case 'T':
							switch ( p3 ) {
								case 'A':case 'G':case 'R':
									aa = 'L';goto aa_set2;
								default:
									aa = 'X';goto aa_set2;
							}
						default:
							aa = 'X';goto aa_set2;
					}
				default:
					aa = 'X';goto aa_set2;
			}
			aa_set2:;
			*cppointer++ = aa;
			if ( aa == 'X' ) numUnknowns++;
		}
	}
	label:;
	if ( maxNTermAA )	*(protein+maxNTermAA) = 0;
	else				*cppointer = 0;
}
void readProteinFromDNALowerCase ( int dnaReadingFrame, int maxNTermAA, char* fpointer, int length, char* protein, int& numUnknowns )
{
	char* sequence_end = fpointer + length;
	char* sequence_start = fpointer;
	char aa;
	char p1, p2, p3;
	numUnknowns = 0;
	char* cppointer = protein;
	if ( length == 0 ) {
		*cppointer = 0;
		return;
	}
	int direction = ( dnaReadingFrame > 3 );
	int reading_frame = ( direction ) ? dnaReadingFrame - 4 : dnaReadingFrame - 1;

	if ( direction ) {	/* Backwards */
		char* backwards_read_end = sequence_start - 1;
		fpointer = sequence_end - 1;
		while ( *fpointer <= MAX_NON_PRINT || *fpointer == '>' ) fpointer--;
		for ( int i = 0 ; i < reading_frame ; i++ ) {
			while ( *fpointer <= MAX_NON_PRINT ) fpointer--;
			fpointer--;
		}
		for ( ; ; ) {
			do {
				if ( fpointer == backwards_read_end ) goto label;
				p1 = *fpointer--;
			} while ( p1 <= MAX_NON_PRINT );
			do {
				if ( fpointer == backwards_read_end ) goto label;
				p2 = *fpointer--;
			} while ( p2 <= MAX_NON_PRINT );
			do {
				if ( fpointer == backwards_read_end ) goto label;
				p3 = *fpointer--;
			} while ( p3 <= MAX_NON_PRINT );
			switch ( p1 ) {
				case 't':
					switch ( p2 ) {
						case 't':
							switch ( p3 ) {
								case 't':case 'c':
									aa = 'K';
									goto aa_set1;
								case 'a':case 'g':case 'r':
									aa = 'N';goto aa_set1;
								case 'y':
									aa = 'K';goto aa_set1;
								default:
									aa = 'X';goto aa_set1;
							}
						case 'c':
							switch ( p3 ) {
								case 't':case 'c':
									aa = 'R';goto aa_set1;
								case 'a':case 'g':case 'r':
									aa = 'S';goto aa_set1;
								case 'y':
									aa = 'R';goto aa_set1;
								default:
									aa = 'X';goto aa_set1;
							}
						case 'a':
							switch ( p3 ) {
								case 't':case 'a':case 'g':
									aa = 'I';goto aa_set1;
								case 'c':
									aa = 'M';goto aa_set1;
								case 'd':case 'k':case 'r':case 'w':
									aa = 'I';goto aa_set1;
								default:
									aa = 'X';goto aa_set1;
							}
						case 'g':
							aa = 'T';goto aa_set1;
						default:
							aa = 'X';goto aa_set1;
					}
				case 'c':
					switch ( p2 ) {
						case 't':
							switch ( p3 ) {
								case 't':case 'c':
									aa = 'E';goto aa_set1;
								case 'a':case 'g':case 'r':
									aa = 'D';goto aa_set1;
								case 'y':
									aa = 'E';goto aa_set1;
								default:
									aa = 'X';goto aa_set1;
							}
						case 'c':
							aa = 'G';goto aa_set1;
						case 'a':
							aa = 'V';goto aa_set1;
						case 'g':
							aa = 'A';goto aa_set1;
						default:
							aa = 'X';goto aa_set1;
					}
				case 'a':
					switch ( p2 ) {
						case 't':
							switch ( p3 ) {
								case 't':case 'c':
									aa = '.';goto aa_set1;
								case 'a':case 'g':case 'r':
									aa = 'Y';goto aa_set1;
								case 'y':
									aa = '.';goto aa_set1;
								default:
									aa = 'X';goto aa_set1;
							}
						case 'c':
							switch ( p3 ) {
								case 't':
									aa = '.';goto aa_set1;
								case 'c':
									aa = 'W';goto aa_set1;
								case 'a':case 'g':case 'r':
									aa = 'C';goto aa_set1;
								default:
									aa = 'X';goto aa_set1;
							}
						case 'a':
							switch ( p3 ) {
								case 't':case 'c':
									aa = 'L';goto aa_set1;
								case 'a':case 'g':case 'r':
									aa = 'F';goto aa_set1;
								case 'y':
									aa = 'L';goto aa_set1;
								default:
									aa = 'X';goto aa_set1;
							}
						case 'g':
							aa = 'S';goto aa_set1;
						case 'y':
							switch ( p3 ) {
								case 't':
									aa = '.';goto aa_set1;
								default:
									aa = 'X';goto aa_set1;
							}
						default:
							aa = 'X';goto aa_set1;
					}
				case 'g':
					switch ( p2 ) {
						case 't':
							switch ( p3 ) {
								case 't':case 'c':
									aa = 'Q';goto aa_set1;
								case 'a':case 'g':case 'r':
									aa = 'H';goto aa_set1;
								case 'y':
									aa = 'Q';goto aa_set1;
								default:
									aa = 'X';goto aa_set1;
							}
						case 'c':
							aa = 'R';goto aa_set1;
						case 'a':
							aa = 'L';goto aa_set1;
						case 'g':
							aa = 'P';goto aa_set1;
						default:
							aa = 'X';goto aa_set1;
					}
				case 'k':
					switch ( p2 ) {
						case 'c':
							switch ( p3 ) {
								case 't':case 'c':case 'y':
									aa = 'R';goto aa_set1;
								default:
									aa = 'X';goto aa_set1;
							}
						default:
							aa = 'X';goto aa_set1;
					}
				case 'r':
					switch ( p2 ) {
						case 'a':
							switch ( p3 ) {
								case 't':case 'c':case 'y':
									aa = 'L';goto aa_set1;
								default:
									aa = 'X';goto aa_set1;
							}
						default:
							aa = 'X';goto aa_set1;
					}
				default:
					aa = 'X';goto aa_set1;
			}
			aa_set1:;
			*cppointer++ = aa;
			if ( aa == 'X' ) numUnknowns++;
		}
	}
	else {
		fpointer = sequence_start;
		for ( int i = 0 ; i < reading_frame ; i++ ) {
			while ( *fpointer <= MAX_NON_PRINT ) fpointer++;
			fpointer++;
		}
		for ( ; ; ) {
			do {
				if ( fpointer == sequence_end ) goto label;
				p1 = *fpointer++;
			} while ( p1 <= MAX_NON_PRINT );
			do {
				if ( fpointer == sequence_end ) goto label;
				p2 = *fpointer++;
			} while ( p2 <= MAX_NON_PRINT );
			do {
				if ( fpointer == sequence_end ) goto label;
				p3 = *fpointer++;
			} while ( p3 <= MAX_NON_PRINT );
			switch ( p1 ) {
				case 't':
					switch ( p2 ) {
						case 't':
							switch ( p3 ) {
								case 't':case 'c':
									aa = 'F';goto aa_set2;
								case 'a':case 'g':case 'r':
									aa = 'L';goto aa_set2;
								case 'y':
									aa = 'F';goto aa_set2;
								default:
									aa = 'X';goto aa_set2;
							}
						case 'c':
							aa = 'S';goto aa_set2;
						case 'a':
							switch ( p3 ) {
								case 't':case 'c':
									aa = 'Y';goto aa_set2;
								case 'a':case 'g':case 'r':
									aa = '.';goto aa_set2;
								case 'y':
									aa = 'Y';goto aa_set2;
								default:
									aa = 'X';goto aa_set2;
							}
						case 'g':
							switch ( p3 ) {
								case 't':case 'c':
									aa = 'C';goto aa_set2;
								case 'a':
									aa = '.';goto aa_set2;
								case 'g':
									aa = 'W';goto aa_set2;
								case 'y':
									aa = 'C';goto aa_set2;
								default:
									aa = 'X';goto aa_set2;
							}
						case 'r':
							switch ( p3 ) {
								case 'a':
									aa = '.';goto aa_set2;
								default:
									aa = 'X';goto aa_set2;
							}
						default:
							aa = 'X';goto aa_set2;
					}
				case 'c':
					switch ( p2 ) {
						case 't':
							aa = 'L';goto aa_set2;
						case 'c':
							aa = 'P';goto aa_set2;
						case 'a':
							switch ( p3 ) {
								case 't':case 'c':
									aa = 'H';goto aa_set2;
								case 'a':case 'g':case 'r':
									aa = 'Q';goto aa_set2;
								case 'y':
									aa = 'H';goto aa_set2;
								default:
									aa = 'X';goto aa_set2;
							}
						case 'g':
							aa = 'R';goto aa_set2;
						default:
							aa = 'X';goto aa_set2;
					}
				case 'a':
					switch ( p2 ) {
						case 't':
							switch ( p3 ) {
								case 't':case 'c':case 'a':
									aa = 'I';goto aa_set2;
								case 'g':
									aa = 'M';goto aa_set2;
								case 'h':case 'm':case 'w':case 'y':
									aa = 'I';goto aa_set2;
								default:
									aa = 'X';goto aa_set2;
							}
						case 'c':
							aa = 'T';goto aa_set2;
						case 'a':
							switch ( p3 ) {
								case 't':case 'c':
									aa = 'N';goto aa_set2;
								case 'a':case 'g':case 'r':
									aa = 'K';goto aa_set2;
								case 'y':
									aa = 'N';goto aa_set2;
								default:
									aa = 'X';goto aa_set2;
							}
						case 'g':
							switch ( p3 ) {
								case 't':case 'c':
									aa = 'S';goto aa_set2;
								case 'a':case 'g':case 'r':
									aa = 'R';goto aa_set2;
								case 'y':
									aa = 'S';goto aa_set2;
								default:
									aa = 'X';goto aa_set2;
							}
						default:
							aa = 'X';goto aa_set2;
					}
				case 'g':
					switch ( p2 ) {
						case 't':
							aa = 'V';goto aa_set2;
						case 'c':
							aa = 'A';goto aa_set2;
						case 'a':
							switch ( p3 ) {
								case 't':case 'c':
									aa = 'D';goto aa_set2;
								case 'a':case 'g':case 'r':
									aa = 'E';goto aa_set2;
								case 'y':
									aa = 'D';goto aa_set2;
								default:
									aa = 'X';goto aa_set2;
							}
						case 'g':
							aa = 'G';goto aa_set2;
						default:
							aa = 'X';goto aa_set2;
					}
				case 'm':
					switch ( p2 ) {
						case 'g':
							switch ( p3 ) {
								case 'a':case 'g':case 'r':
									aa = 'R';goto aa_set2;
								default:
									aa = 'X';goto aa_set2;
							}
						default:
							aa = 'X';goto aa_set2;
					}
				case 'y':
					switch ( p2 ) {
						case 't':
							switch ( p3 ) {
								case 'a':case 'g':case 'r':
									aa = 'L';goto aa_set2;
								default:
									aa = 'X';goto aa_set2;
							}
						default:
							aa = 'X';goto aa_set2;
					}
				default:
					aa = 'X';goto aa_set2;
			}
			aa_set2:;
			*cppointer++ = aa;
			if ( aa == 'X' ) numUnknowns++;
		}
	}
	label:;
	if ( maxNTermAA )	*(protein+maxNTermAA) = 0;
	else				*cppointer = 0;
}
void DNASequenceReader::readProtein ( int serial_number, int dnaReadingFrame )
{
	int length;
	char* fpointer = dbIndicies->getProteinPointer ( serial_number, &length );
	if ( length && islower ( *fpointer ) )
		readProteinFromDNALowerCase ( dnaReadingFrame, maxNTermAA, fpointer, length, protein, numUnknowns );
	else
		readProteinFromDNA ( dnaReadingFrame, maxNTermAA, fpointer, length, protein, numUnknowns );
}
DNAProteinSequenceReader::DNAProteinSequenceReader ( DatabaseIndicies* dbIndicies ) :
	SequenceReader ( dbIndicies ) {}
DNAProteinSequenceReader::~DNAProteinSequenceReader () {}

void DNAProteinSequenceReader::readProtein ( int serialNumber, int dnaReadingFrame )
{
	static char* fpointer;
	char* cppointer;
	char val;
	int i;
	static int saved_serial_number = -1;
	static int saved_frame = -1;

	if ( serialNumber != saved_serial_number || dnaReadingFrame != saved_frame + 1 ) {
		int length;
		fpointer = dbIndicies->getProteinPointer ( serialNumber, &length );
		for ( i = 1 ; i < dnaReadingFrame ; i++ ) while ( *fpointer++ != '\n' );
	}
	cppointer = protein;
	for ( ; ; ) {
		val = *fpointer++;
		if ( val <= MAX_NON_PRINT ) break;
		if ( val == 'X' ) numUnknowns++;
		*cppointer++ = val;
	}
	if ( maxNTermAA )	*(protein+maxNTermAA) = 0;
	else				*cppointer = 0;
	saved_serial_number = serialNumber;
	saved_frame = dnaReadingFrame;
}
const char* CommentLine::unreadable = "UNREADABLE";
CommentLine::CommentLine ( DatabaseIndicies* dbIndicies )
{
	this->dbIndicies = dbIndicies;
	unreadable_species_count = 0;
	unsigned int maxCommentLength = dbIndicies->getMaxCommentLength ();
	accessionNumber = new char [maxCommentLength];
	accessionInfo = new char [maxCommentLength];
	accessionInfo [0] = 0;
	name = new char [maxCommentLength];
	species = new char [maxCommentLength];
	uniprotID = new char [maxCommentLength];
	uniprotID [0] = 0;
	uniprotInfo = new char [maxCommentLength];
	uniprotInfo [0] = 0;
}
CommentLine::~CommentLine ()
{
	delete [] accessionNumber;
	delete [] accessionInfo;
	delete [] name;
	delete [] species;
	delete [] uniprotID;
	delete [] uniprotInfo;
}
void CommentLine::genpeptLine ( bool speciesUndefined )
{
	char* nameStart = point;
	char* commentEnd;
	bool set = false;

	for ( ; ; ) {						// Find the end of the comment line
		char c = *point;
		if ( c == '\n' || c == 1 ) {
			commentEnd = point;
			break;
		}
		point++;
	}
	if ( point [-1] == ']' ) {				// Possible correct format
		point--;
		int bracket = 0;
		for ( ; ; ) {
			point--;
			char c = *point;
			if ( c == '[' ) {
				if ( bracket == 0 ) {
					int nameLen = point - nameStart;
					if ( speciesUndefined ) {
						copyToUpper ( species, point + 1, commentEnd - point - 2 );
					}
					if ( point == nameStart )
						strcpy ( name, "No name." );
					else
						copyTo ( name, nameStart, nameLen );
					set = true;
					break;
				}
				bracket--;
			}
			if ( c == ']' ) bracket++;
			if ( point == nameStart ) break;
		}
	}
	if ( set == false ) {
		if ( speciesUndefined )	unreadableSpecies ();
		else					copyToPtr ( name, nameStart, commentEnd );
	}
	point = commentEnd;
}
void CommentLine::uniprotLine ()
{
	char* savePoint = point;
	advance1 ();
	bool old;
	bool newest;
	for ( ; ; ) {
		char c = *point++;
		if ( c == '_' ) {
			old = true;
			break;
		}
		if ( c == '|' ) {
			old = false;
			for ( ; ; ) {
				char c2 = *point++;
				if ( c2 == '|' ) {
					newest = true;
					break;
				}
				if ( c2 == ' ' ) {
					newest = false;
					break;
				}
			}
			break;
		}
	}
	point = savePoint;
/*
Newest format http://www.uniprot.org/help/fasta-headers
>db|UniqueIdentifier|EntryName ProteinName OS=OrganismName[ GN=GeneName]PE=ProteinExistence SV=SequenceVersion
*/
	if ( old ) {	// >104K_THEPA (P15711) 104 kDa microneme-rhoptry antigen
		setUniprotID ();
		goPastNext ( '_' );
		setSpecies ( ' ' );
		goPastNext ( '(' );
		setAccessionNumber ( ')' );
		goPastNext ( ' ' );
		genpeptLine ( false );
	}
	else {
		if ( newest ) {	// >sp|Q4U9M9|104K_THEAN 104 kDa microneme-rhoptry antigen precursor (p104) - Theileria annulata
			goPastNext ( '|' );
		}
		else {			// >Q4U9M9|104K_THEAN 104 kDa microneme-rhoptry antigen precursor (p104) - Theileria annulata
			advance1 ();
		}
		setAccessionNumber ( '|' );
		setUniprotID ();
		goPastNext ( '_' );
		setSpecies ( ' ' );
		goPastNext ( ' ' );
		spaceDashSpaceLine ( false );
	}
}
void CommentLine::spaceDashSpaceLine ( bool speciesUndefined )
{
	char* name_start = point;
	int name_len = 0;
	char* species_start = NULL;
	int species_len = 0;
	char* commentEnd;

	int phase = 1;
	for ( ; ; ) {
		char c = *point;
		if ( c == '\n' || c == 1 ) {
			species_len = point - species_start;
			commentEnd = point;
			break;
		}
		switch ( phase ) {
			case 1:
				if ( c == ' ' ) phase = 2;
				break;
			case 2:
				if ( c == '-' ) phase = 3;
				else phase = 1;
				break;
			case 3:
				if ( c == ' ' ) phase = 4;
				else phase = 1;
				break;
		}
		if ( phase == 4 ) {
			name_len = point - 2 - name_start;
			species_start = point + 1;
			phase = 1;
		}
		point++;
	}
	if ( speciesUndefined ) {
		if ( species_start == NULL || species_len <= 0 || name_len < 0 ) {
			unreadableSpecies ();
		}
		else if ( !strncmp ( point - 3, "...", 3 ) ) {	// Truncated species
			unreadableSpecies ();
		}
		else if ( point [-1] == '>' ) {
			unreadableSpecies ();
		}
		else {
			if ( point [-1] == '.' ) species_len--;	// Species with full stop at end
			copyToUpper ( species, species_start, species_len );
			if ( name_len == 0 )
				strcpy ( name, "No name." );
			else
				copyTo ( name, name_start, name_len );
		}
	}
	else {
		if ( name_len == 0 )
			copyToPtr ( name, name_start, point );
		else
			copyTo ( name, name_start, name_len );
	}
	char* osp = strstr ( name, " OS=" );
	if ( osp ) {
		strcpy ( uniprotInfo, osp+1 );
		*osp = 0;
	}
	point = commentEnd;
}
void CommentLine::get_default_top_line ( int serialNumber )
{
	sprintf ( accessionNumber, "%d", serialNumber );
	unreadableSpecies ();
}
void CommentLine::unreadableSpecies ()
{
	strcpy ( species, "UNREADABLE" );
	point = top_line_start + 1;
	copyTo ( name, '\n' );
	unreadable_species_count++;
}
void CommentLine::copyTo ( char* to, int delimiter )
{
	for ( ;  ; ) {
		char c = *point;
		if ( c == delimiter ) break;
		if ( c <= MAX_NON_PRINT ) break;
		*to++ = c;
		point++;
	}
	*to = 0;
}
void CommentLine::copyToStripTrailing ( char* to, int delimiter )
{
	for ( ;  ; ) {
		char c = *point;
		if ( c == delimiter ) break;
		if ( c <= MAX_NON_PRINT ) break;
		*to++ = c;
		point++;
	}
	*to = 0;
	--to;
	while ( *to == ' ' || *to == '\t' ) {
		*to-- = 0;
	}
}
void CommentLine::copyToUpper ( char* to, char* start, int len )
{
	char* end = start + len;
	for ( point = start ; point < end ; ) {
		*to++ = toupper ( *point++ );
	}
	*to = 0;
}
void CommentLine::goPastNext ( int delimiter )
{
	while ( *point++ != delimiter );
}
void CommentLine::advance1 ()
{
	point++;
}
void CommentLine::advance4 ()
{
	point += 4;
}
void CommentLine::advance5 ()
{
	point += 5;
}
bool CommentLine::checkChar ( int c )
{
	return ( *point == c );
}
void CommentLine::setNumericAccessionNumber ()
{
	copyNumericTo ( accessionNumber );
}
void CommentLine::setAccessionNumber ( int delimiter )
{
	copyTo ( accessionNumber, delimiter );
}
void CommentLine::setAccessionNumber ( int delimiter1, int delimiter2 )
{
	copyTo ( accessionNumber, delimiter1, delimiter2 );
}
void CommentLine::setSpecies ( int delimiter )
{
	copyTo ( species, delimiter );
}
void CommentLine::setUniprotID ()
{
	advance1 ();
	char* savePoint = point;
	copyTo ( uniprotID, ' ' );
	point = savePoint;
}
void CommentLine::setName ( int delimiter )
{
	copyTo ( name, delimiter );
}
void CommentLine::setName ( int delimiter1, int delimiter2 )
{
	copyTo ( name, delimiter1, delimiter2 );
}
void CommentLine::copyToUpperStripTrailing ( char* to, int delimiter )
{
	for ( ;  ; ) {
		char c = *point;
		if ( c == delimiter ) break;
		if ( c <= MAX_NON_PRINT ) break;
		*to++ = toupper ( c );
		point++;
	}
	*to = 0;
	--to;
	while ( *to == ' ' || *to == '\t' ) {
		*to-- = 0;
	}
}
void CommentLine::copyTo ( char* to, int delimiter1, int delimiter2 )
{
	for ( ;  ; ) {
		char c = *point;
		if ( c == delimiter1 ) break;
		if ( c == delimiter2 ) break;
		if ( c <= MAX_NON_PRINT ) break;
		*to++ = c;
		point++;
	}
	*to = 0;
}
void CommentLine::copyTo ( char* to, char* start, int len )
{
	char* end = start + len;
	for ( point = start ; point < end ; ) {
		*to++ = *point++;
	}
	*to = 0;
}
void CommentLine::copyToPtr ( char* to, char* start, char* end )
{
	for ( point = start ; point < end ; ) {
		*to++ = *point++;
	}
	*to = 0;
}
void CommentLine::copyNumericTo ( char* to )
{
	for ( ; ; ) {
		char c = *point;
		if ( isdigit ( c ) == false && c != '-' ) break;
		*to++ = c;
		point++;
	}
	*to = 0;
}
char* CommentLine::copyToCheckUnderscore ( char* to, int delimiter )
{
	char* underscore = NULL;
	for ( ;  ; ) {
		char c = *point;
		if ( c == '_' ) underscore = point + 1;
		if ( c == delimiter ) break;
		if ( c <= MAX_NON_PRINT ) break;
		*to++ = c;
		point++;
	}
	*to = 0;
	return underscore;
}
GenpeptCommentLine::GenpeptCommentLine ( DatabaseIndicies* dbIndicies ) :
	CommentLine ( dbIndicies ) {}

GenpeptCommentLine::~GenpeptCommentLine () {}
void GenpeptCommentLine::getCommentLine ( int serialNumber, bool allLines )
{
//	Sample entries:
//
//	>gi|216790 (D13314) arginine deiminase [Mycoplasma hominis]
//	>gi|261706|bbs|120303 (S50809) protein LG=immunoglobulin binding protein {immunoglobulin binding domains} [streptococcus, Peptide Recombinant, 455 aa]
//	>gi|404694 (L13970) lctR contains a helix-turn-helix domain shared by several DNA-binding proteins, suggesting that it may encode a regulator for other lct genes [Escherichia coli]
//	>gi|1387979 (L77099) 44% identity over 302 residues with hypothetical protein from Synechocystis sp, accession D64006_CD; expression induced by environmental stress; some similarity to glycosyl transferases; two potential membrane-spanning helices [Bacillus subtil>
//	>gi|2065214|gnl|PID|e314286 (Z95117) MLC1351.03c, unknown, len: 256 aa, similar to eg. YCIL_ECOLI P37765 hypothetical 32.7 kd protein in trpl- (291 aa), fasta clones, opt: 481 z-score: 570.9 E(): 8.5e-25, (42.4% identity in 229 aa overlap); contains PS01149 Hypothetical yciL/yejD>
//	>gi|1123088 (U42436) coded for by C. elegans cDNA yk56a1.3; coded for by C. elegans cDNA CEMSG41FB; coded for by C. elegans cDNA yk81f4.5; coded for by C. elegans cDNA yk56a1.5; coded for by C. elegans cDNA yk81f4.3;  similar to the S5P family of ribosomal proteins
//	>gi|2330745|gnl|PID|e334350 (Z98598) SPAC1B3.11c, ras-related protein, len:234aa, similar eg. to RB4B_RAT, P51146, ra
//	Zero name length - species OK.
//	>gi|1575686 (U70379)  [Synechococcus PCC7942]
//	Two bracketed fields
//	>gi|3928871 (AF093627) poly(ADP)-ribose polymerase [Zea mays]
//	Identical 2nd accession numbers.
//	>gi|3928875 (AF093611) putative chloroplast desaturase [Acetabularia acetabulum]
//	>gi|3928876 (AF093611) putative chloroplast desaturase [Acetabularia acetabulum]
//	No accession number.
//	>gi|3928883 unknown
//	Bracketed zone before line end.
//	>gi|3881286|gnl|PID|e1350785 (AL021507) [980325 dl] : Prediction spanned chimera, modified based on new 3' sequence information (o/l with F14D1); cDNA EST EMBL:D34402 comes from this gene; cDNA EST EMBL:D37454 comes from this gene; cDNA EST EMBL:D68054 comes from this gene; cDNA E>
//	No 0x01 characters as line extenders.
//
	int length;
	if ( !speciesList.empty () ) {
		accessionNumberList.clear ();
		speciesList.clear ();
		nameList.clear ();
	}
	point = dbIndicies->getCommentPointer ( serialNumber, &length );
	do {
		top_line_start = point;
		advance4 ();					// Skip past >gi|
		setNumericAccessionNumber ();
		goPastNext ( ' ' );
		genpeptLine ( true );
		if ( allLines ) {
			accessionNumberList.push_back ( accessionNumber );
			speciesList.push_back ( species );
			nameList.push_back ( name );
		}
		else break;
	} while ( *point != '\n' );
	goPastNext ( '\n' );
}
SwissProtCommentLine::SwissProtCommentLine ( DatabaseIndicies* dbIndicies ) :
	CommentLine ( dbIndicies ) {}

SwissProtCommentLine::~SwissProtCommentLine () {}
void SwissProtCommentLine::getCommentLine ( int serialNumber, bool allLines )
{
//	Sample entries:
//
//	>sp|P15394|REPA_AGRTU REPLICATING PROTEIN 
//
//	This is the output of the sp2fasta program
//
//	>gi|122068|sp|P16105|H32_BOVIN HISTONE H3 (H3.2)
//
//	This is what you get if you download the database from NCBI.
//
//	956 entries are of the following form:
//	>gi|400027|sp||HYEP_PSESP_2 [Segment 2 of 3] EPOXIDE HYDROLASE (EPOXIDE HYDRATASE)
//	In these cases the accession number is taken as being the first number ie. 400027
//	for the example shown.
//
//	The species is extracted from the code between the last vertical bar and the first space.
//	It appears after the first underscore and before the second underscore (if present).
//	The name is the rest of the line.
//
//  Uniprot style entry - September 2005
//
//  >104K_THEPA (P15711) 104 kDa microneme-rhoptry antigen
//
//  New Uniprot style entry - December 2006
//
//  >Q4U9M9|104K_THEAN 104 kDa microneme-rhoptry antigen precursor (p104) - Theileria annulata
//
//  Even newer Uniprot style entry - August 2008
//
//  >sp|Q4U9M9|104K_THEAN 104 kDa microneme/rhoptry antigen OS=Theileria annulata GN=TA08425 PE=3 SV=1
//
//  Newest line style from SwissProt
//
//  >gi|13878750|sp|Q9CDN0.1|RS18_LACLA RecName: Full=30S ribosomal protein S18
//
	int length;
	if ( !speciesList.empty () ) {
		accessionNumberList.clear ();
		speciesList.clear ();
		nameList.clear ();
		uniprotIDList.clear ();
		uniprotInfoList.clear ();
	}
	point = dbIndicies->getCommentPointer ( serialNumber, &length );
	do {
		top_line_start = point;
		if ( point [3] == '|' ) {	// Old style or newest style
			advance1 ();
			if ( checkStr2Advance ( "gi" ) ) {
				advance1 ();					// Skip past |
				setNumericAccessionNumber ();	// Deals with the empty SwissProt accession number bug.
				advance4 ();					// Skip past |sp|
			}
			else {								// Old style must be sp
				advance1 ();					// Skip past |
			}
			if ( checkChar ( '|' ) ) {			// Check for empty accession number
				setUniprotID ();
				goPastNext ( '_' );
				setSpecies ( '_' );
			}
			else {
				setAccessionNumber ( '|' );		// Overwrite numeric accession number
				char* dot = strchr ( accessionNumber, '.' );
				if ( dot ) *dot = 0;
				setUniprotID ();
				goPastNext ( '_' );
				setSpecies ( ' ' );
			}
			goPastNext ( ' ' );					// Skip past space
			if ( !strncmp ( point, "RecName: Full=", 14 ) ) {
				goPastNext ( '=' );
				setName ( ';', '\n' );
				for ( ;  ; ) {
					char c = *point;
					if ( c == '\n' ) break;
					if ( c <= MAX_NON_PRINT ) break;
					point++;
				}
			}
			else {
				setName ( '\n' );
			}
			char* osp = strstr ( name, " OS=" );
			if ( osp ) {
				strcpy ( uniprotInfo, osp+1 );
				*osp = 0;
			}
		}
		else
			uniprotLine ();
		if ( allLines ) {
			accessionNumberList.push_back ( accessionNumber );
			speciesList.push_back ( species );
			nameList.push_back ( name );
			uniprotIDList.push_back ( uniprotID );
			uniprotInfoList.push_back ( uniprotInfo );
		}
		else break;
	} while ( *point != '\n' );
	goPastNext ( '\n' );
}
bool SwissProtCommentLine::checkStr2Advance ( const char* str )
{
	if ( *point++ == *str++ ) {
		if ( *point++ == *str ) return true;
		else return false;
	}
	*point++;
	return false;
}
LudwigCommentLine::LudwigCommentLine ( DatabaseIndicies* dbIndicies ) :
	CommentLine ( dbIndicies ) {}

LudwigCommentLine::~LudwigCommentLine () {}
void LudwigCommentLine::getCommentLine ( int serialNumber, bool allLines )
{
// Info gathered April 29th 2008
//
// cat Aaegypti_nr.seq Agambiae_nr.seq Amellifera_nr.seq Btaurus_nr.seq Cbriggsae_nr.seq Celegans_nr.seq Cfamiliaris_nr.seq 
// Cintestinalis_nr.seq Cporcellus_nr.seq Csavignyi_nr.seq Dmelanogaster_nr.seq Dnovemcinctus_nr.seq Drerio_nr.seq Ecaballus_nr.seq 
// Eeuropaeus_nr.seq Etelfairi_nr.seq Fcatus_nr.seq Gaculeatus_nr.seq Ggallus_nr.seq Hsapiens_nr.seq Lafricana_nr.seq Mdomestica_nr.seq 
// Mlucifugus_nr.seq Mmulatta_nr.seq Mmurinus_nr.seq Mmusculus_nr.seq Oanatinus_nr.seq Ocuniculus_nr.seq Ogarnettii_nr.seq 
// Olatipes_nr.seq Oprinceps_nr.seq Pberghei_nr.seq Pchabaudi_nr.seq Pfalciparum_nr.seq Pknowlesi_nr.seq Ppygmaeus_nr.seq 
// Ptroglodytes_nr.seq Pvivax_nr.seq Pyoelii_nr.seq Rnorvegicus_nr.seq Saraneus_nr.seq Scerevisiae_nr.seq Stridecemlineatus_nr.seq
// Tbelangeri_nr.seq Tgondii_nr.seq Tnigroviridis_nr.seq Trubripes_nr.seq Xtropicalis_nr.seq sludge_aus_nr.seq sludge_us1_nr.seq 
// sludge_us2_nr.seq swiss_nr.seq swiss_varsplic_nr.seq trembl_nr.seq wormpep_nr.seq yeastpep_nr.seq > Ludwignr.04.29.2008
//
//
// >DBtag|AccessionNb|OtherId (GeneName)description[species]
//
// Aaegypti_nr.seq
//
// >ens|AAEL015636-PA pep:novel supercontig:AaegL1:supercont1.3990:5173:5592:1 gene:AAEL015636 transcript:AAEL015636-RA
//
//
// Agambiae_nr.seq
//
// >ens|AAEL015636-PA pep:novel supercontig:AaegL1:supercont1.3990:5173:5592:1 gene:AAEL015636 transcript:AAEL015636-RA
//
//
// Pberghei_nr.seq
//
// >plasmo|Plasmodium_berghei_strain_ANKA|PB_RP3411|PB107738.00.0|Annotation|Plasmodium_berghei_Sanger|(protein coding) conserved hypothetical protein
//
//
// Pfalciparum_nr.seq
//
// >plasmo|Plasmodium_falciparum_3D7|MAL8|MAL8P1.155|Annotation|Plasmodium_falciparum_Sanger_Stanford_TIGR|(protein coding) hypothetical protein
//
//
// sludge_aus_nr.seq
//
// >sludge|2000346600 Predicted ATPase (AAA+ superfamily) [Sludge/Australian, Phrap Assembly]
// >sludge|2000346610  [Sludge/Australian, Phrap Assembly]
//
//
// sludge_us1_nr.seq
//
// >sludge|2001000010  [Sludge/US, Jazz Assembly]
// >sludge|2001000030 Predicted metal-binding protein [Sludge/US, Jazz Assembly]
//
//
// swiss_nr.seq
//
// >DBtag|AccessionNb|OtherId (GeneName)description[species]
// >sp|Q4U9M9|104K_THEAN 104 kDa microneme/rhoptry antigen precursor (p104).[Theileria annulata]
//
//
// swiss_varsplic_nr.seq
//
// >sp_vs|P15455-2|12S1_ARATH-2 (CRA1)Isoform 2 of P15455.[Arabidopsis thaliana]
//
//
// trembl_nr.seq
//
// >tr|A0AQI8|A0AQI8_9ARCH (amoA)Putative ammonia monooxygenase (Fragment).[uncultured archaeon]
// >tr|A0B530|A0B530_METTP Quinolinate synthetase complex, A subunit.[Methanosaeta thermophila]
//
//
// wormpep_nr.seq
//
// >wp|CE03921|C04F6.1 WBGene00006929 locus:vit-5 status:Partially_confirmed SW:P06125 protein_id:AAA83587.1[C. elegans]
// >wp|CE15495|ZK1010.1 WBGene00006728 locus:ubq-2 UBQ-2 ubiquitin\; 60S Ribosomal protein L40 status:Confirmed[C. elegans]
// >wp|CE37074|ZC101.2f WBGene00006787 locus:unc-52 status:Partially_confirmed SW:Q06561 protein_id:CAH04744.1[C. elegans]
// >wp|CE39776|T14B4.4b WBGene00006636 locus:tsp-10 status:Confirmed[C. elegans]
// >wp|CE36390|K04H4.1b WBGene00001263 locus:emb-9 collagen status:Partially_confirmed SW:P17139 protein_id:CAE52901.2[C. elegans]
//
//
// yeastpep_nr.seq
//
// >yp|YIR031C DAL7 SGDID:S000001470, Chr IX from 414676-413012, reverse complement
// , Verified ORF, "Malate synthase, role in allantoin degradation unknown; express
// ion sensitive to nitrogen catabolite repression and induced by allophanate, an i
// ntermediate in allantoin degradation"[S. cerevisiae]
// >yp|YLR248W RCK2 SGDID:S000004238, Chr XII from 634254-636086, Verified ORF, "Pr
// otein kinase involved in the response to oxidative and osmotic stress; identifie
// d as suppressor of S. pombe cell cycle checkpoint mutations"[S. cerevisiae]
// >yp|YHR092C HXT4 SGDID:S000001134, Chr VIII from 288814-287132, reverse compleme
// nt, Verified ORF, "High-affinity glucose transporter of the major facilitator su
// perfamily, expression is induced by low levels of glucose and repressed by high
// levels of glucose"[S. cerevisiae]
//
// Older database format
//
// db|accno|ID|CRC Description[species]
//
// db - database
// CRC - 64bit cyclic redundancy check
//
// >gp|M84711|182775|000037AE195F7A9D v-fos transformation effector protein [Homo sapiens]
// >gp|AL391014|9716128|0006579AD1B1EEE8 putative DNA-binding protein [Streptomyces coelicolor A3(2)]
// >pir|A91719|GGIC1A|0027F62F6F36BA36 globin CTT-IA - midge (Chironomus thummi thummi)[Chironomus thummi thummi]
// >pir|JX0361|JX0361|00013E4475F84453 subtilisin-trypsin inhibitor, SIL10 - Streptomyces sp.[Streptomyces sp.]
// >pir|JC7193|PC7055|0154D83E82AA822B cell division protein FtsQ - Streptomyces collinus (fragment)[Streptomyces collinus]
// >pir|A29526|A29526|02AC2025766BCBC7 ubiquitin B processed pseudogene - human[Homo sapiens]
// >sp|P55820|SN25_RABIT|00014F740FEB29C5 (SNAP..)SYNAPTOSOMAL-ASSOCIATED PROTEIN 25 (SNAP-25) (SUPER PROTEIN) (SUP) (FRAGMENTS).[Oryctolagus cuniculus]
// >sp_vs|P16157-01|P16157|004EDB42F81EBDE8 ISOFORM 2.2 OF P16157[Homo sapiens]
// >tr|AF247519|AAF71733|0001F06BB33BD2E8 Gag protein (Fragment).[Human immunodeficiency virus type 1]
// >tr|U83613|O09751|0000148C132C06BD (POL)REVERSE TRANSCRIPTASE (FRAGMENT).[Human immunodeficiency virus type 1]
// >tr_vs|P70390-01|P70390|0172F8C6825A0023 ISOFORM OG12B/PRX3B OF P70390[Mus musculus]
// >wp|CE24847|C44C3.3|0205CAE438EE8B14 (ST.LOUIS) TR:P91157 protein_id:AAB37360.1[C. elegans]
// >yp|ORFP:YDR094W|0642CC1F954A58E2 YDR094W, Chr IV from 635833-636168[S. cerevisiae]
//
	int length;
	if ( !speciesList.empty () ) {
		accessionNumberList.clear ();
		speciesList.clear ();
		nameList.clear ();
	}
	point = dbIndicies->getCommentPointer ( serialNumber, &length );
	do {
		top_line_start = point;
		goPastNext ( '|' );
		setAccessionNumber ( '|' );
		goPastNext ( ' ' );					// Skip past space
		genpeptLine ( true );
		if ( allLines ) {
			accessionNumberList.push_back ( accessionNumber );
			speciesList.push_back ( species );
			nameList.push_back ( name );
		}
		else break;
	} while ( *point != '\n' );
	goPastNext ( '\n' );
}
OwlCommentLine::OwlCommentLine ( DatabaseIndicies* dbIndicies ) :
	CommentLine ( dbIndicies ) {}

OwlCommentLine::~OwlCommentLine () {}
void OwlCommentLine::getCommentLine ( int serialNumber, bool allLines )
{
//
//	The accession number is always between the second vertical bar and the first space.
//	The Dec 1998 version of Owl contained 279796 entries.
//
//	The first 74015 entries are from Swiss Prot and are of the following format:
//	>owl|Q62671|100K_RAT 100 KD PROTEIN (EC 6.3.2.-). - RATTUS NORVEGICUS (RAT).
//	The species is encoded after the underscore in the accession number.
//	The name is the rest of the line. However if the species is added at the end of
//	the line after space dash space then this is not included in the name.
//	The name is frequently trunctated in which case there are 3 dots at the end of the line.
//	>owl|P15455|12S1_ARATH 12S SEED STORAGE PROTEIN PRECURSOR. - ARABIDOPSIS THALIANA (MOUSE-EAR...
//	The species is sometimes truncated in this fashion. Also a second variant of the species is
//	often included in brackets.
//
//	The next 53488 entries are of the following format:
//	>owl|B40638|B40638 isocytochrome c2 - Rhodobacter sphaeroides
//
//	The next 150959 entries are of the following format:
//	>owl|Z31371|A7120FTSZ1 A7120FTSZ NID: g1100793 - Anabaena PCC7120.
//	Note the dot following the species name.
//
//	For the above two cases the species has to be extracted from the end of the line. If the
//	entry is truncated. Eg:
//	>owl|E69371|E69371 bile acid-inducible operon protein F (baiF-1) homolog - Archaeoglobus...
//	or
//	>owl|U70848|CELC43G24 CELC43G2 NID: g1572755 - Caenorhabditis elegans strain=Bristol...
//	The species has to be set to UNREADABLE. Note in this second case the first few letters
//	of the accession number probably denotes the sequence.
//
//
//	The last 1334 entries are of the following format:
//	>owl||NRL_1A00B hemoglobin beta chain mutant (V1M, W37Y) (deoxy), chain B - human
//	The only way to find these is to look for an empty field between the two vertical bars. The
//	NRL before the underscore can't be used to find these entries as some of the Swiss Prot
//	contain this string. Eg:
//	>owl|P32961|NRL1_ARATH NITRILASE 1 (EC 3.5.5.1). - ARABIDOPSIS THALIANA (MOUSE-EAR CRESS).
//	>owl|P54845|NRL_HUMAN NEURAL RETINA-SPECIFIC LEUCINE ZIPPER PROTEIN (NRL) (D14S46E)....
//	The species isn't always present. If it is present it is at the end of the line after a
//	space dash space. A second way of writing the species is sometimes included in brackets.
//
//	Old style entries:
//
//	>10KD_VIGUN 10 KD PROTEIN PRECURSOR (CLONE PSAS10). - VIGNA UNGUICULATA (COWPEA).
//	>AEOHFPA AEOHFPA NID: g141875 - A.hydrophila DNA, clone pPH4.
//	>pir|P15711|104K_THEPA 104 KD MICRONEME-RHOPTRY ANTIGEN. - THEILERIA PARVA.
//
	int length;
	if ( !speciesList.empty () ) {
		accessionNumberList.clear ();
		speciesList.clear ();
		nameList.clear ();
	}
	point = dbIndicies->getCommentPointer ( serialNumber, &length );
	do {
		top_line_start = point;
		goPastNext ( '|' );
		if ( checkChar ( '|' ) ) {			// Check for empty 2nd field.
			advance1 ();
			setAccessionNumber ( ' ' );
			advance1 ();
			spaceDashSpaceLine ( true );
		}
		else {
			goPastNext ( '|' );
			if ( setAccessionNumberCheckUnderscore ( ' ' ) ) {
				setSpecies ( ' ' );
				advance1 ();
				spaceDashSpaceLine ( false );
			}
			else {
				advance1 ();
				spaceDashSpaceLine ( true );
			}
		}
		if ( allLines ) {
			accessionNumberList.push_back ( accessionNumber );
			speciesList.push_back ( species );
			nameList.push_back ( name );
		}
		else break;
	} while ( *point != '\n' );
	goPastNext ( '\n' );
}
bool OwlCommentLine::setAccessionNumberCheckUnderscore ( int delimiter )
{
	char* underscore = copyToCheckUnderscore ( accessionNumber, delimiter );
	if ( underscore != NULL ) {
		point = underscore;
		return true;
	}
	return false;
}
NCBINRCommentLine::NCBINRCommentLine ( DatabaseIndicies* dbIndicies ) :
	CommentLine ( dbIndicies ) {}

NCBINRCommentLine::~NCBINRCommentLine () {}
void NCBINRCommentLine::getCommentLine ( int serialNumber, bool allLines )
{
/*
//	Sample entries:
//
// 356315/563276 entries
// 153509/0 entries
>gi|149575 (M76708) L(+)-lactate dehydrogenase [Lactobacillus casei]
>gi|45803 (X04609) gamma subunit (3'terminus); pid:g45803 [thermophilic bacterium PS3]
>gi|289135 (L10036) unknown [Anabaena PCC7120]
>gi|402254 (U01238) beta subunit of the molybdenum-iron nitrogenase [Frankia sp.]
>gi|414523 (U02284) beta-lactamase [Cloning vector pSP65]
>gi|439619 (L25848) [Salmonella typhimurium IS200 insertion sequence from SARA17, partial.], gene product [Salmonella typhimurium]
>gi|431128 (L15633) start [Transposon Tn916]
>gi|466378 (U07618) SSB [Unknown]
>gi|403947 (U01693) (M90060);  Homology to GenBank Accession numORF-X from STRATPASEA [Mycoplasma genitalium]
>gi|405516 (L22217) This ORF is homologous to nitroreductase from Enterobacter cloacae, Accession Number A38686, and Salmonella, Accession Number P15888. [Mycoplasma-like organism]
>gi|457139 (L29100) transposase [Insertion sequence IS150 homolog]
>gi|468279 (L31491) nreA [pTOM9]
>gi|413733 (L25424) orf 1 [Plasmid pCB2.4]
>gi|144453 (M94320) very similar to DNA polymerase of Bacillus subtilis bacteriophage SPO2; potential DNA polymerase; putative [Citrus greening disease-associated bacterium-like organism]
>gi|520517 (U10338) RNA polymerase II, largest subunit [Ilyanassa obsoleta ]
>gi|971400 (X88862) immunogenic polyprotein with 2A protease [Foot-and-mouth disease virus]
>gi|1008449 (L19624) envelope glycoprotein [Human immunodeficiency virus type 1]
>gi|1718307 (U75698) ORF 54; dUTPase homolog; EBV BLLF3 homolog [Kaposi's sarcoma-associated herpesvirus]
>gi|2271117 (AF008696) hemagglutinin [influenza A virus (A/South_Australia/68/92(H3N2))]
>gi|2429520 (AF025469) Similar to acetyl-CoA carboxylase; coded for by C. elegans cDNA yk16c3.3; coded for by C. elegans cDNA yk36b11.3; coded for by C. elegans cDNA yk43h8.3; coded for by C. elegans cDNA yk24d2.3; coded for by C. elegans cDNA yk24d2.5;...
>gi|2444119 (U88974) ORF40 [Streptococcus thermophilus temperate bacteriophage O1205]
>gi|2662546 (AF036688) No definition line found [Caenorhabditis elegans]
>gi|4206510 (AF066801) ribulose 1,5-bisphosphate carboxylase [Dictamnus sp. M.W.Chase-1820K]
// 1 entry (done)
>gi|3928883 unknown
// sp 77419/73385 entries
>gi|132349|sp|P15394|REPA_AGRTU REPLICATING PROTEIN
>gi|123494|sp|P22291|SULD_STRPN BIFUNCTIONAL FOLATE SYNTHESIS PROTEIN (DIHYDRONEOPTERIN ALDOLASE (DHNA) / 2-AMINO-4-HYDROXY-6-HYDROXYMETHYLDIHYDROPTERIDINE PYROPHOSPHOKINASE (7,8-DIHYDRO-6-HYDROXYMETHYLPTERIN PYROPHOSPHOKINASE) (HPPK) (6-HYDROXYMETHYL-7...
>gi|4033439|sp||LEC_VICVI_1 [Segment 1 of 4] LECTIN B4 (VVLB4)
// 61785/28 entries
>gi|216351|gnl|PID|d1003451 (D13793) ORF [Bacillus subtilis]
// 44139/105601 entries
>gi|282349|pir||A41961 chitinase (EC 3.2.1.14) D - Bacillus circulans
>gi|80297|pir||JN0146 hypothetical protein (div+ 3' region) - Bacillus subtilis (fragment)
>gi|77616|pir||A36125 branched-chain amino acid transport protein braC - Pseudomonas aeruginosa (strain PAO)
>gi|538696|pir||A40613 avirulence protein avrRpt2 - Pseudomonas syringae (strain DC3000, pv. tomato)
>gi|98505|pir||S21241 oligo-1,6-glucosidase (EC 3.2.1.10) - Bacillus "thermoamyloliquefaciens" (strain KP1071) (fragment)
>gi|320384|pir||A37388 probable DNA-binding protein 1A - Thermus aquaticus (strain HB8) insertion sequence IS1000
>gi|477498|pir||A49131 releasechannel homolog - fruit fly (Drosophila melanogaster) (fragment)
// 9875/0 entries (done)
>gi|3712669|bbs|85194 (S85224) vascular endothelial growth factor; VEGF 206 [Homo sapiens]
>gi|386065|bbs|133195 cytochrome c3 {N-terminal} [Desulfovibrio vulgaris, NCIMB 8303, Peptide Partial, 22 aa]
>gi|386067|bbs|133197 cytochrome c3
>gi|236142|bbs|57690 (S57688) EF-G=elongation factor G [Thermotoga maritima, Peptide, 682 aa] [Thermotoga maritima]
>gi|435743|bbs|139151 (S66567) alpha-atrial natriuretic factor/coat protein, alpha-ANF/coat protein=fusion polypeptide(coat protein, alpha-atrial natriuretic factor, alpha-ANF) [human, bacteriophage fr, expression vector pFAN15, Peptide PlasmidSynthetic...
>gi|913316|bbs|163145 (S76565) T-cell receptor beta chain VJ region {clone N4} [not specified, vesicular stomatitis virus-specific CTL, Peptide Partial, 15 aa] [unidentified]
>gi|833965|bbs|160632 (S75335) polyprotein(structural protein C, structural protein E, structural protein M, structural protein PreM, nonstructural protein NS1) [dengue type 1 D1 virus, Mochizuki, Peptide Partial, 50 aa, segment 2 of 2] [Dengue virus ty...
// 5116/7761 entries (done)
>gi|230242|pdb|1PFK|A Escherichia coli
>gi|4139942|pdb|1BC5|T Chain T, Chemotaxis Receptor Recognition By Protein Methyltransferase Cher
>gi|231004|pdb|4ER4|I synthetic construct
>gi|494001|pdb|1EGF|  Epidermal Growth Factor (Egf) (Nmr, 16 Structures)
>gi|493782|pdb|146L|  Lysozyme (E.C.3.2.1.17) Mutant With Cys 54 Replaced By Thr, Cys 97 Replaced By Ala, Leu 121 Replaced By Met, Ala 129 Replaced By Leu, Leu 133 Replaced By Met, Val 149 Replaced By Ile, Phe 153 Replaced By Trp (C54t,C97a,L121m,A129l,...
>gi|230275|pdb|1R1A|1 Human rhinovirus 1A
// 3545/3543 entries (done)
>gi|742246|prf||2009326A beta glucosidase [Cellvibrio gilvus]
>gi|225172|prf||1210227A amylase subtilisin inhibitor alpha [Hordeum vulgare var. distichum]
// 884/37698 entries (done)
>gi|2440229|dbj||AB006689_5 (AB006689) ORF13 [Agrobacterium rhizogenes]
>gi|1805521|dbj||D90852_18 (D90852) ORF_ID:o250#11; similar to [SwissProt Accession Number P19779]; start codon is not identified yet [Escherichia coli]
// 42/232258 entries (done)
>gi|1680564|gb||S58174_1 (S58174) putative RNA polymerase [Pelargonium leaf curl virus]
>gi|1683178|gb||S69825_2 (S69825) coat/capsid protein [Sweet potato feathery mottle virus (strain CH)]
>gi|1683615|gb||S81342_1 (S81342) unnamed protein product [Mus sp.]
//  0/68424 entries
>gi|6|emb|CAA42669.1| (X60065) beta-2-glycoprotein  I [Bos taurus]
>gi|6065756|emb|CAB58425.1| (AJ238324) Clostridium difficile binary toxin A
>gi|6018922|emb|CAB58111.1| (AL121806) /prediction=(method:""genefinder"", version:""084"", score:""32.36"")~/prediction=(method:""genscan"", version:""1.0"")~/match=(desc:""EUKARYOTIC TRANSLATION INITIATION FACTOR 4E (EIF-4E) (EIF4E) (MRNA CAP-BINDING PROTEIN) (EIF-4F 25 KD SUBU>
//  0/34578 entries
>gi|5713315|ref|NP_002060.1| guanine nucleotide binding protein (G protein), alpha inhibiting activity polypeptide 1
*/
	int length;
	if ( !speciesList.empty () ) {
		accessionNumberList.clear ();
		speciesList.clear ();
		nameList.clear ();
	}
	point = dbIndicies->getCommentPointer ( serialNumber, &length );
	do {
		top_line_start = point;
		advance4 ();					// Skip past >gi|
		setNumericAccessionNumber ();
		if ( checkCharAdvance ( ' ' ) ) {			// Space after accession number
			genpeptLine ( true );
		}
		else {										// Must be | (vertical bar)
			if ( !strncmp ( point, "sp", 2 ) ) {
				goPastNext ( '|' );
				if ( checkChar ( '|' ) ) {			// Check for empty accession number
					goPastNext ( '_' );
					setSpecies ( '_' );
				}
				else {
					goPastNext ( '_' );
					setSpecies ( ' ' );
				}
				goPastNext ( ' ' );					// Skip past space
				setName ( '\n' );
			}
			else if ( !strncmp ( point, "gnl", 3 ) ) {
				goPastNext ( ' ' );
				genpeptLine ( true );
			}
			else if ( !strncmp ( point, "pir", 3 ) ) {
				goPastNext ( ' ' );
				spaceDashSpaceLine ( true );
			}
			else if ( !strncmp ( point, "bbs", 3 ) ) {
				goPastNext ( ' ' );
				genpeptLine ( true );
			}
			else if ( !strncmp ( point, "pdb", 3 ) ) {
				unreadableSpecies ();
			}
			else if ( !strncmp ( point, "prf", 3 ) || !strncmp ( point, "dbj", 3 )  || !strncmp ( point, "gb", 2 ) ) {
				goPastNext ( ' ' );
				genpeptLine ( true );
			}
			else if ( !strncmp ( point, "emb", 3 ) ) {
				goPastNext ( ' ' );
				genpeptLine ( true );
			}
			else if ( !strncmp ( point, "ref", 3 ) ) {
				goPastNext ( ' ' );
				genpeptLine ( true );
			}
			else if ( !strncmp ( point, "tpg", 3 ) || !strncmp ( point, "tpe", 3 ) || !strncmp ( point, "tpd", 3 ) ) {
				goPastNext ( ' ' );
				genpeptLine ( true );
			}
			else {
				unreadableSpecies ();
			}
		}
		if ( allLines ) {
			accessionNumberList.push_back ( accessionNumber );
			speciesList.push_back ( species );
			nameList.push_back ( name );
		}
		else break;
	} while ( *point != '\n' );
	goPastNext ( '\n' );
	//gen_error ( "The database has an unknown comment line format.\n" );
}
bool NCBINRCommentLine::checkCharAdvance ( int c )
{
	return ( *point++ == c );
}
IPICommentLine::IPICommentLine ( DatabaseIndicies* dbIndicies ) :
	CommentLine ( dbIndicies ) {}

IPICommentLine::~IPICommentLine () {}
void IPICommentLine::getCommentLine ( int serialNumber, bool allLines )
{
/*
//	Sample entries:
//
>IPI:IPI00177321.1|REFSEQ_XP:XP_168060|ENSEMBL:ENSP00000343431 Tax_Id=9606 similar to NOD3 protein
>IPI:IPI00015171.1|UniProt/Swiss-Prot:O43931 Tax_Id=9606 AFG3-like protein 1
*/
	int length;
	if ( !speciesList.empty () ) {
		accessionNumberList.clear ();
		speciesList.clear ();
		nameList.clear ();
	}
	point = dbIndicies->getCommentPointer ( serialNumber, &length );
	do {
		advance5 ();					// Skip past >IPI:
		setAccessionNumber ( '|' );
		goPastNext ( ' ' );					// Skip past space
		goPastNext ( ' ' );					// Skip past space
		unreadable_species_count++;
		copyTo ( name, '\n' );
		if ( allLines ) {
			accessionNumberList.push_back ( accessionNumber );
			speciesList.push_back ( "UNREADABLE" );
			nameList.push_back ( name );
		}
		else break;
	} while ( *point != '\n' );
	goPastNext ( '\n' );
}
UniprotCommentLine::UniprotCommentLine ( DatabaseIndicies* dbIndicies ) :
	CommentLine ( dbIndicies ) {}

UniprotCommentLine::~UniprotCommentLine () {}
void UniprotCommentLine::getCommentLine ( int serialNumber, bool allLines )
{
/*
//	Sample entries:
//
>104K_THEPA (P15711) 104 kDa microneme-rhoptry antigen
>O05152_SULAC (O05152) Glycogen debranching enzyme
//
//  Even newer Uniprot style entry - August 2008
//
>sp|Q4U9M9|104K_THEAN 104 kDa microneme/rhoptry antigen OS=Theileria annulata GN=TA08425 PE=3 SV=1
>tr|A0AQI4|A0AQI4_9ARCH Putative ammonia monooxygenase (Fragment) OS=uncultured archaeon GN=amoA PE=4 SV=1
//
Newest format http://www.uniprot.org/help/fasta-headers

>db|UniqueIdentifier|EntryName ProteinName OS=OrganismName[ GN=GeneName]PE=ProteinExistence SV=SequenceVersion
>sp|P15711|104K_THEPA 104 kDa microneme/rhoptry antigen OS=Theileria parva GN=TP04_0437 PE=2 SV=1

db					is 'sp' for UniProtKB/Swiss-Prot and 'tr' for UniProtKB/TrEMBL.
UniqueIdentifier	is the primary accession number of the UniProtKB entry.
EntryName			is the entry name of the UniProtKB entry.
ProteinName			is the recommended name of the UniProtKB entry as annotated in the RecName field from release 14.0 on.
					For UniProtKB/TrEMBL entries without a RecName field, the SubName field is used. The 'precursor'
					attribute is excluded, 'Fragment' is included with the name if applicable.

OS					OrganismName is the scientific name of the organism of the UniProtKB entry.
GN					GeneName is the first gene name of the UniProtKB entry. If there is no gene name, OrderedLocusName
					or ORFname, the GN field is not listed.
PE					ProteinExistence is the numerical value describing the evidence for the existence of the protein.
						1. Evidence at protein level
						2. Evidence at transcript level
						3. Inferred from homology
						4. Predicted
						5. Uncertain
SV					SequenceVersion is the version number of the sequence.
*/
	int length;
	if ( !speciesList.empty () ) {
		accessionNumberList.clear ();
		speciesList.clear ();
		nameList.clear ();
		uniprotIDList.clear ();
		uniprotInfoList.clear ();
	}
	point = dbIndicies->getCommentPointer ( serialNumber, &length );
	do {
		uniprotLine ();
		if ( allLines ) {
			accessionNumberList.push_back ( accessionNumber );
			speciesList.push_back ( species );
			nameList.push_back ( name );
			uniprotIDList.push_back ( uniprotID );
			uniprotInfoList.push_back ( uniprotInfo );
		}
		else break;
	} while ( *point != '\n' );
	goPastNext ( '\n' );
}
DBESTCommentLine::DBESTCommentLine ( DatabaseIndicies* dbIndicies ) :
	CommentLine ( dbIndicies )
{
	read_dbest_species_list_file ();
}

DBESTCommentLine::~DBESTCommentLine () {}
void DBESTCommentLine::read_dbest_species_list_file ()
{
	char* fileInfo = getParamsFileInfo ( "dbEST.spl.txt", &num_dbest_species );
	dbest_species_list = new char* [num_dbest_species];
	uc_dbest_species_list = new char* [num_dbest_species];

	for ( int i = 0 ; i < num_dbest_species ; i++ ) {
		dbest_species_list [i] = ( i == 0 ) ? strtok ( fileInfo, "\n" ) : strtok ( NULL, "\n" );
		uc_dbest_species_list [i] = gen_new_string ( dbest_species_list [i] );
		gen_strupr ( uc_dbest_species_list [i] );
	}
}
void DBESTCommentLine::getCommentLine ( int serialNumber, bool allLines )
{
//	Sample entries:
//
//	>gi|1705383|gb|N20717|N20717 SMNHADA002044SK SmAW Schistosoma mansoni cDNA 5' 
//	>gi|3771232|gb|AI209290|AI209290 SWOvAFCAP09G09SK Onchocerca volvulus adult female cDNA (SAW98MLW-OvAF) Onchocerca volvulus cDNA clone SWOvAFCAP09G09 5', mRNA sequence [Onchocerca volvulus].gi|3789602|gb|AI216948|AI216948 SWOvAFCAP10G11SK Onchocerca volvulus adult female cDNA (SAW98MLW-OvAF) Onchocerca volvulus cDNA clone SWOvAFCAP10G115', mRNA sequence [Onchocerca volvulus]
	int length;
	if ( !speciesList.empty () ) {
		accessionNumberList.clear ();
		speciesList.clear ();
		nameList.clear ();
	}
	point = dbIndicies->getCommentPointer ( serialNumber, &length );
	do {
		top_line_start = point;
		advance4 ();					// Skip past >gi|
		setNumericAccessionNumber ();
		goPastNext ( ' ' );
		setName ( '\n' );
		species = const_cast <char*> (unreadable);
		char* n = name;
		bool found = false;
		for ( int i = 0 ; i < num_dbest_species ; i++ ) {
			if ( strstr ( n, dbest_species_list [i] ) ) {
				species = uc_dbest_species_list [i];
				found = true;
				break;
			}
		}
		if ( found == false ) unreadable_species_count++;
		if ( allLines ) {
			accessionNumberList.push_back ( accessionNumber );
			speciesList.push_back ( species );
			nameList.push_back ( name );
		}
		else break;
	} while ( *point != '\n' );
	goPastNext ( '\n' );
}
GenericCommentLine::GenericCommentLine ( DatabaseIndicies* dbIndicies ) :
	CommentLine ( dbIndicies ) {}

GenericCommentLine::~GenericCommentLine () {}
void GenericCommentLine::getCommentLine ( int serialNumber, bool allLines )
{
//	Sample entries:
//
//	the longest ">" line is 2986 characters long and starts out 
//	> 417909| zampli05rt1|Mouse|pancreas|
//	consisting of a number of "|" delimited fields.
//
	int length;
	if ( !speciesList.empty () ) {
		accessionNumberList.clear ();
		speciesList.clear ();
		nameList.clear ();
	}
	point = dbIndicies->getCommentPointer ( serialNumber, &length );
	do {
		top_line_start = point;
		advance1 ();					// Skip past >
		skipSpace ();					// Skip any leading space
		if ( lineEnd () ) {
			get_default_top_line ( serialNumber );	// Maybe should say no name - empty comment line.
		}
		else {
			setAccessionNumber ( '|', ' ' );
			skipSpace ();					// Skip any trailing space
			skipChar ( '|' );				// Skip '|' if there
			skipSpace ();					// Skip any space after a '|'
			if ( lineEnd () ) {
				unreadableSpecies ();
			}
			else {
				setName ( '|' );
				skipChar ( '|' );				// Skip '|' if there
				skipSpace ();					// Skip any space after a '|'
				if ( lineEnd () )
					unreadableSpecies ();
				else
					copyToUpperStripTrailing ( species, '|' );
			}
		}
		if ( allLines ) {
			accessionNumberList.push_back ( accessionNumber );
			speciesList.push_back ( species );
			nameList.push_back ( name );
		}
		else break;
	} while ( *point != '\n' );
	goPastNext ( '\n' );
}
void GenericCommentLine::skipSpace ()
{
	while ( *point == ' ' || *point == '\t' ) point++;
}
void GenericCommentLine::skipChar ( int c )
{
	while ( *point == c ) point++;
}
bool GenericCommentLine::lineEnd ()
{
	return ( *point <= MAX_NON_PRINT );
}
DefaultCommentLine::DefaultCommentLine ( DatabaseIndicies* dbIndicies ) :
	CommentLine ( dbIndicies ) {}

DefaultCommentLine::~DefaultCommentLine () {}
void DefaultCommentLine::getCommentLine ( int serialNumber, bool allLines )
{
	int length;
	if ( !speciesList.empty () ) {
		accessionNumberList.clear ();
		speciesList.clear ();
		nameList.clear ();
	}
	point = dbIndicies->getCommentPointer ( serialNumber, &length );
	do {
		top_line_start = point;
		get_default_top_line ( serialNumber );
		if ( allLines ) {
			accessionNumberList.push_back ( accessionNumber );
			speciesList.push_back ( species );
			nameList.push_back ( name );
		}
		else break;
	} while ( *point != '\n' );
	goPastNext ( '\n' );
}
PreloadedDatabases::PreloadedDatabases ()
{
	StringVector preloadDatabases = InfoParams::instance ().getStringVectorValue ( "preload_database" );
	for ( int i = 0 ; i < preloadDatabases.size () ; i++ ) {
		fs.push_back ( new FastaServer ( preloadDatabases [i] ) );
		fs.back ()->loadIntoMemory ();
	}
}
PreloadedDatabases::~PreloadedDatabases ()
{
	for ( int i = 0 ; i < fs.size () ; i++ ) {
		delete fs [i];
	}
}

void wrCommentLine ( FILE* fp, const string& filename, const string& accession_number, const string& name, const string& species )
{
	const char* a = accession_number.c_str ();
	const char* n = name.c_str ();
	const char* s = species.c_str ();
	int numeric_accession_number = atoi ( a );

	if ( is_genpept_database ( filename ) )
		fprintf ( fp, ">gi|%s %s [%s]\n", a, n, s );
	else if ( is_swissprot_database ( filename ) ) {
		if ( strchr ( s, '_' ) == NULL )
			fprintf ( fp, ">sp|%s|_%s %s\n", a, s, n );
		else
			fprintf ( fp, ">sp|%s|%s %s\n", a, s, n );
	}
	else if ( is_ludwig_database ( filename ) )
		fprintf ( fp, ">gp|%s|%s %s [%s]\n", a, a, n, s );
	else if ( is_owl_database ( filename ) )
		fprintf ( fp, ">owl||%s %s - %s.\n", a, n, s );
	else if ( is_ncbinr_database ( filename ) )
		fprintf ( fp, ">gi|%s %s [%s]\n", a, n, s );
	else if ( is_ipi_database ( filename ) )
		fprintf ( fp, ">IPI:%s| Tax_Id=0 %s\n", a, n );
	else if ( is_uniprot_database ( filename ) )
		fprintf ( fp, ">%s_%s (%s) %s\n", a, s, a, n );
	else if ( is_dbest_database ( filename ) ) {
		if ( numeric_accession_number <= 0 ) ErrorHandler::genError ()->error ( "The accession number entered must be numeric and greater than zero.\n" );
		fprintf ( fp, ">gi|%u|gb|1|%s %s\n", numeric_accession_number, n, s );
	}
	else if ( is_generic_numeric_acc_number_database ( filename ) ) {
		if ( numeric_accession_number <= 0 ) ErrorHandler::genError ()->error ( "The accession number entered must be numeric and greater than zero.\n" );
		fprintf ( fp, ">%u|%s|%s\n", numeric_accession_number, n, s );
	}
	else if ( is_generic_string_acc_number_database ( filename ) )
		fprintf ( fp, ">%s|%s|%s\n", a, n, s );
	else
		throw runtime_error ( "Invalid database file name.\n" );
}
