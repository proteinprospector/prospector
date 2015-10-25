/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_faindex.cpp                                                *
*                                                                             *
*  Created    : October 27th 1997                                             *
*                                                                             *
*  Purpose    : Functions for creating the database index files.              *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1997-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifndef VIS_C
#include <stdexcept>
#endif
#include <algorithm>
#include <lg_new.h>
#include <lg_stdio.h>
#include <lgen_file.h>
#include <lu_fasta.h>
#include <lu_acc_num.h>
#include <lu_mass_seq.h>
#include <lu_pi.h>
#include <lu_getfil.h>
#include <lu_mass.h>
#include <lu_html.h>
#include <lu_check_db.h>
#include <lg_string.h>
#include <lu_species.h>
using std::string;
using std::vector;
using std::ios_base;
using std::endl;
using std::sort;
using std::cout;
using std::reverse;
using std::random_shuffle;
using std::unique;
using std::runtime_error;

static const int REPORTING_FREQUENCY = 100000;

class UniprotAndNCBISpeciesMap {
	UniProtSpeciesMap upsm;
	NCBISpeciesMap nsm;
	bool uniProtDefault;
public:
	UniprotAndNCBISpeciesMap ( bool uniProtDefault );
	string getActualSpecies ( const string& sp );
	int getNode ( const string& sp );
	int getDefaultNode ( const string& sp ) { return uniProtDefault ? upsm.getNode ( sp ) : nsm.getNode ( sp ); }
	int getNonDefaultNode ( const string& sp ) { return uniProtDefault ? nsm.getNode ( sp ) : upsm.getNode ( sp ); }
};
UniprotAndNCBISpeciesMap::UniprotAndNCBISpeciesMap ( bool uniProtDefault ) :
	 uniProtDefault ( uniProtDefault ) 
{
}

/*
	"PEPTIDE PARTIAL"							// 5591 entries
	"PEPTIDE"									// 1952 entries
	"PEPTIDE CHLOROPLAST"						// 9 entries
	"PEPTIDE CHLOROPLAST PARTIAL"				// 84 entries
	"PEPTIDE INSERTION PARTIAL"					// 2 entries
	"PEPTIDE KINETOPLAST"						// 1 entry
	"PEPTIDE KINETOPLAST PARTIAL"				// 1 entry
	"PEPTIDE MITOCHONDRIAL"						// 41 entries
	"PEPTIDE MITOCHONDRIAL KINETOPLAST PARTIAL"	// 1 entry
	"PEPTIDE MITOCHONDRIAL MUTANT"				// 2 entries
	"PEPTIDE MITOCHONDRIAL PARTIAL"				// 91 entries
	"PEPTIDE MITOCHONDRIAL PARTIAL MUTANT"		// 4 entries
	"PEPTIDE MUTAGENESIS"						// 4 entries
	"PEPTIDE MUTANT"							// 56 entries
	"PEPTIDE PARTIAL MUTAGENESIS"				// 9 entries
	"PEPTIDE PARTIAL MUTANT"					// 246 entries
	"PEPTIDE PARTIALMUTANT"						// 4 entries
	"PEPTIDE PLASMID"							// 12 entries
	"PEPTIDE PLASMID INSERTION MUTANT"			// 1 entry
	"PEPTIDE PLASMID MITOCHONDRIAL PARTIAL"		// 1 entry
	"PEPTIDE PLASMID PARTIAL"					// 11 entries
	"PEPTIDE PLASMID RECOMBINANT PARTIAL"		// 3 entries
	"PEPTIDE RECOMBINANT"						// 9 entries
	"PEPTIDE RECOMBINANT PARTIAL"				// 61 entries
	"PEPTIDE RECOMBINANT PARTIAL MUTANT"		// 1 entry
	"PEPTIDE SYNTHETIC"							// 5 entries
	"PEPTIDE TRANSPOSON"						// 7 entries
	"PEPTIDE TRANSPOSON PARTIAL"				// 32 entries
												// total - 8241 entries - 11 Dec 2007
*/
/*
Function that looks in the taxonomy node map for a matching species or species grouping. If it can't
find one it returns an empty string.
*/
string UniprotAndNCBISpeciesMap::getActualSpecies ( const string& sp )
{
	int node = getDefaultNode ( sp );
	if ( node == 0 ) node = getNonDefaultNode ( sp );
	if ( node == 0 ) {
		if ( sp.find ( ", PEPTIDE" ) != string::npos ) {
			int x = sp.find ( ',' );
			if ( x != string::npos ) {
				string sp2 = sp.substr ( 0, x );
				int node = getDefaultNode ( sp2 );
				if ( node == 0 ) node = getNonDefaultNode ( sp2 );
				if ( node == 0 ) {
					int y = sp2.find ( '=' );
					if ( y != string::npos ) {
						string sp3 = sp2.substr ( 0, y );
						int node = getDefaultNode ( sp3 );
						if ( node == 0 ) node = getNonDefaultNode ( sp3 );
						if ( node == 0 ) {
							return "";		// Unknown species	
						}
						else return sp3;	// The taxonomy species
					}
					else return "";		// Unknown species	
				}
				else return sp2;
			}
			else return "";
		}
		else if ( isSuffix ( sp, " (FRAGMENT)" ) ||											// 1071 entries
			isSuffix ( sp, " (FRAGMENTS)" ) ) {												// 201 entries
																							// total = 1272 entries
			string sp2 = gen_strtrim ( sp.substr ( 0, sp.length () - 11 ) );
			int node = getDefaultNode ( sp2 );
			if ( node == 0 ) node = getNonDefaultNode ( sp );
			if ( node == 0 ) {
				return "";
			}
			else return sp2;
		}
		else return "";
	}
	else return sp;
}
int UniprotAndNCBISpeciesMap::getNode ( const string& sp )
{
	int node = getDefaultNode ( sp );
	if ( node == 0 ) {
		node = getNonDefaultNode ( sp );
		if ( node == 0 ) {
			if ( sp.find ( ", PEPTIDE" ) != string::npos ) {
				int x = sp.find ( ',' );
				if ( x != string::npos ) {
					string sp2 = sp.substr ( 0, x );
					int node = getDefaultNode ( sp2 );
					if ( node == 0 ) {
						node = getNonDefaultNode ( sp2 );
						if ( node == 0 ) {
							int y = sp2.find ( '=' );
							if ( y != string::npos ) {
								string sp3 = sp2.substr ( 0, y );
								int node = getDefaultNode ( sp3 );
								if ( node == 0 ) { 
									node = getNonDefaultNode ( sp3 );
									if ( node == 0 ) return 0;		// Unknown species	
									else return node;				// The taxonomy species
								}
								else return node;
							}
							else return 0;		// Node not found	
						}
						else return node;
					}
					else return node;
				}
				else return 0;	// Node not found
			}
			else if ( isSuffix ( sp, " (FRAGMENT)" ) ||											// 1071 entries
				isSuffix ( sp, " (FRAGMENTS)" ) ) {												// 201 entries
																								// total = 1272 entries
				string sp2 = gen_strtrim ( sp.substr ( 0, sp.length () - 11 ) );
				int node = getDefaultNode ( sp2 );
				if ( node == 0 ) {
					node = getNonDefaultNode ( sp2 );
					if ( node == 0 ) return 0;	// Node not found
					else return node;
				}
				else return node;
			}
			else return 0;	// Node not found
		}
		else return node;
	}
	else return node;
}
static void faindexSerial ( FastaServer* fs, const string& fileName );
static void faindexParallel ( FastaServer* fs, const string& fileName );
static void processSpecies ( FastaServer* fs, const string& fileName );
static void processAccessionNumbers ( FastaServer* fs, const string& fileName );
static void processAlphaNumericAccessionNumbers ( FastaServer* fs, const string& fileName );
static void processNumericAccessionNumbers ( FastaServer* fs, const string& fileName );
static void processMW ( FastaServer* fs, const string& fileName );
static void processPI ( FastaServer* fs, const string& fileName );
static string convertDNADatabaseToProteinDatabase ( const string& database );
static void deleteDatabaseFiles ( const string& database );
static void deleteDatabaseIndexFiles ( const string& database );
static void deleteCRCharactersFromDatabaseFile ( const string& fileName );
static void writeTaxAndTL ( const string& fileName, VectorIndexedInt& species );
static void writeUnk ( const string& fileName, const SetString& unkSpecies );
static void writeACN ( const string& fileName, vector <IndexedInt>& numericAccessionNumber );
static void writeACC ( const string& fileName, VectorIndexedConstCString& accessionNumber );
static void writeMW ( const string& fileName, IndexedDouble* mole_wt, int numEntries );
static void writePI ( const string& fileName, IndexedDouble* pi, int numEntries );

void create_new_database_from_indicies_list ( FastaServer* fs, const string& database, const string& subDatabaseID, const IntVector& indicies )
{
	if ( indicies.size () == 0 ) ErrorHandler::genError ()->error ( "Can't create the new database as the number of indicies is zero.\n" );

	string databaseSuffix = SeqdbDir::instance ().getDatabaseSuffix ( database );		// Does the original database have a suffix
	GenOFStream ost ( SeqdbDir::instance ().getSeqdbDir () + database + subDatabaseID + databaseSuffix, ios_base::binary );

	for ( IntVectorSizeType i = 0 ; i < indicies.size () ; i++ ) {
		int length;
		char* entry = fs->get_fasta_full_entry ( indicies [i], &length );
 		ost.write ( (char*) entry, length * sizeof (char) );
	}
}
void create_database_files ( const string& fileName, bool parallel )
{
	deleteDatabaseIndexFiles ( fileName );
	
	deleteCRCharactersFromDatabaseFile ( fileName );

	ErrorHandler::genError ()->message ( "Counting entries in database and creating index files (.idc, .idp and .idi).\n" );
	FastaServer* fs = new FastaServer ( fileName, true );

	ErrorHandler::genError ()->message ( gen_itoa ( fs->getNumEntries () ) + " entries in the database.\n" );

	initialise_amino_acid_weights ( MapStringConstModPtr (), ElementalFormulaVector (), false );

	if ( parallel )	faindexParallel ( fs, fileName );
	else			faindexSerial ( fs, fileName );

	delete fs;
}
static void faindexSerial ( FastaServer* fs, const string& fileName )
{
	bool dna_database = is_dna_database ( fileName );
	if ( !dna_database ) {
		processMW ( fs, fileName );
		processPI ( fs, fileName );
	}
	processAccessionNumbers ( fs, fileName );
	processSpecies ( fs, fileName );
}
static void faindexParallel ( FastaServer* fs, const string& fileName )
{
	bool spFlag = is_swissprot_database ( fileName ) || is_uniprot_database ( fileName );
	bool dna_database = is_dna_database ( fileName );
	int numEntries = fs->getNumEntries ();
	bool numAccessionNumber = is_numeric_acc_number_database ( fileName );

	VectorIndexedInt species;
	species.reserve ( numEntries );

	GenOFStream ost1 ( SeqdbDir::instance ().getSeqdbDir () + fileName + ".unr" );	// Unreadable entries
	UpdatingJavascriptMessage ujm;
	UniprotAndNCBISpeciesMap unspm ( spFlag );
	SetString unkSpecies;
	bool unr = false;
	vector <IndexedInt> numericAccessionNumber;
	VectorIndexedConstCString accessionNumber;
	if ( numAccessionNumber )
		numericAccessionNumber.reserve ( numEntries );
	else
		accessionNumber.reserve ( numEntries );
	IndexedDouble* mole_wt;
	IndexedDouble* pi;
	if ( !dna_database ) {
		mole_wt = new IndexedDouble [numEntries];
		pi = new IndexedDouble [numEntries];
	}
	for ( int n = 0 ; n < numEntries ; n++ ) {
		int n_plus_1 = n+1;
		if ( n_plus_1 % REPORTING_FREQUENCY == 0 ) ujm.writeMessage ( cout, n_plus_1 );
		for ( fs->firstLine ( n_plus_1 ) ; fs->isDoneLine () ; fs->nextLine () ) {
			string lineSpecies = fs->getLineSpecies ();
			IndexedInt speciesEntry;
			if ( lineSpecies == "UNREADABLE" ) {
				speciesEntry.number = -1;
				ost1 << fs->getLineName () << '\n';
				unr = true;
			}
			else {
				speciesEntry.number = unspm.getNode ( lineSpecies );
				if ( speciesEntry.number == 0 ) {									// This species doesn't appear in the taxonomy file
					unkSpecies.insert ( lineSpecies );
				}
			}
			speciesEntry.index = n_plus_1;
			species.push_back ( speciesEntry );
			string tmpAccessionNum = fs->getLineAccessionNumber(); 
			const char* curAccessionNumber = tmpAccessionNum.c_str(); 
			if ( numAccessionNumber ) {
				IndexedInt nan;
				nan.number = atoi ( curAccessionNumber );
				nan.index = n_plus_1;
				numericAccessionNumber.push_back ( nan );
			}
			else {
				IndexedConstCString an;
				an.name = gen_new_string ( curAccessionNumber );
				an.index = n_plus_1;
				accessionNumber.push_back ( an );
			}
		}
		if ( !dna_database ) {
			char* frame = fs->get_fasta_protein ( n_plus_1, 1 );
			mole_wt [n].index = n_plus_1;
			ProteinMW pmw ( frame );
			mole_wt [n].number = pmw.getMass ();

			pi [n].index = n_plus_1;
			ProteinPI ppi ( frame );
			pi [n].number = ppi.getProteinPI ();
		}
	}
	ost1.close ();
	ujm.deletePreviousMessage ( cout );

	if ( unr ) {
		ErrorHandler::genError ()->message ( gen_itoa ( fs->getUnreadableSpeciesCount () ) + " comment lines contained unreadable species. These are listed in the .unr file.\n" );
	}
	else genUnlink ( SeqdbDir::instance ().getSeqdbDir () + fileName + ".unr" );


	if ( !dna_database ) {
		writeMW ( fileName, mole_wt, numEntries );
		writePI ( fileName, pi, numEntries );
	}
	writeUnk ( fileName, unkSpecies );

	writeTaxAndTL ( fileName, species );

	if ( numAccessionNumber )
		writeACN ( fileName, numericAccessionNumber );
	else
		writeACC ( fileName, accessionNumber );
}
static void processSpecies ( FastaServer* fs, const string& fileName )
{
	int numEntries = fs->getNumEntries ();

	ErrorHandler::genError ()->message ( "Processing entries to find taxonomy.\n" );
	VectorIndexedInt species;
	species.reserve ( numEntries );

	GenOFStream ost ( SeqdbDir::instance ().getSeqdbDir () + fileName + ".unr" );	// Unreadable entries
	UpdatingJavascriptMessage ujm;
	bool spFlag = is_swissprot_database ( fileName ) || is_uniprot_database ( fileName );
	UniprotAndNCBISpeciesMap unspm ( spFlag );
	SetString unkSpecies;
	bool unr = false;
	for ( int n = 0 ; n < numEntries ; n++ ) {
		int n_plus_1 = n+1;
		if ( n_plus_1 % REPORTING_FREQUENCY == 0 ) ujm.writeMessage ( cout, n_plus_1 );
		for ( fs->firstLine ( n_plus_1 ) ; fs->isDoneLine () ; fs->nextLine () ) {
			string lineSpecies = fs->getLineSpecies ();
			IndexedInt speciesEntry;
			if ( lineSpecies == "UNREADABLE" ) {
				speciesEntry.number = -1;
				ost << fs->getLineName () << '\n';
				unr = true;
			}
			else {
				speciesEntry.number = unspm.getNode ( lineSpecies );
				if ( speciesEntry.number == 0 ) {									// This species doesn't appear in the taxonomy file
					unkSpecies.insert ( lineSpecies );
				}
			}
			speciesEntry.index = n_plus_1;
			species.push_back ( speciesEntry );
		}
	}
	ujm.deletePreviousMessage ( cout );
	ost.close ();
	if ( unr ) {
		ErrorHandler::genError ()->message ( gen_itoa ( fs->getUnreadableSpeciesCount () ) + " comment lines contained unreadable species. These are listed in the .unr file.\n" );
	}
	else genUnlink ( SeqdbDir::instance ().getSeqdbDir () + fileName + ".unr" );

	writeUnk ( fileName, unkSpecies );

	writeTaxAndTL ( fileName, species );
}
static void processAccessionNumbers ( FastaServer* fs, const string& fileName )
{
	ErrorHandler::genError ()->message ( "Processing entries to find accession numbers.\n" );

	if ( is_numeric_acc_number_database ( fileName ) ) 
		processNumericAccessionNumbers ( fs, fileName );
	else
		processAlphaNumericAccessionNumbers ( fs, fileName );
}
static void processNumericAccessionNumbers ( FastaServer* fs, const string& fileName )
{
	int numEntries = fs->getNumEntries ();
	vector <IndexedInt> numericAccessionNumber;
	numericAccessionNumber.reserve ( numEntries );

	UpdatingJavascriptMessage ujm;
	for ( int j = 0 ; j < numEntries ; j++ ) {
		int j_plus_1 = j+1;
		if ( j_plus_1 % REPORTING_FREQUENCY == 0 ) ujm.writeMessage ( cout, j_plus_1 );
		for ( fs->firstLine ( j_plus_1 ) ; fs->isDoneLine () ; fs->nextLine () ) {
			string tmpAccessionNum = fs->getLineAccessionNumber(); 
			const char* curAccessionNumber = tmpAccessionNum.c_str(); 
			IndexedInt nan;
			nan.number = atoi ( curAccessionNumber );
			nan.index = j_plus_1;
			numericAccessionNumber.push_back ( nan );
		}
	}
	ujm.deletePreviousMessage ( cout );
	writeACN ( fileName, numericAccessionNumber );
}
static void processAlphaNumericAccessionNumbers ( FastaServer* fs, const string& fileName )
{
	int numEntries = fs->getNumEntries ();
	VectorIndexedConstCString accessionNumber;

	accessionNumber.reserve ( numEntries );

	UpdatingJavascriptMessage ujm;
	for ( int j = 0 ; j < numEntries ; j++ ) {
		int j_plus_1 = j+1;
		if ( j_plus_1 % REPORTING_FREQUENCY == 0 ) ujm.writeMessage ( cout, j_plus_1 );
		for ( fs->firstLine ( j_plus_1 ) ; fs->isDoneLine () ; fs->nextLine () ) {
			string tmpAccessionNum = fs->getLineAccessionNumber(); 
			const char* curAccessionNumber = tmpAccessionNum.c_str(); 
			IndexedConstCString an;
			an.name = gen_new_string ( curAccessionNumber );
			an.index = j_plus_1;
			accessionNumber.push_back ( an );
		}
	}
	ujm.deletePreviousMessage ( cout );
	writeACC ( fileName, accessionNumber );
}
static void processMW ( FastaServer* fs, const string& fileName )
{
	ErrorHandler::genError ()->message ( "Calculating Protein MWs.\n" );
	int numEntries = fs->getNumEntries ();
	IndexedDouble* mole_wt = new IndexedDouble [numEntries];

	UpdatingJavascriptMessage ujm;
	for ( int k = 0 ; k < numEntries ; k++ ) {
		int k_plus_1 = k+1;
		if ( k_plus_1 % REPORTING_FREQUENCY == 0 ) ujm.writeMessage ( cout, k_plus_1 );
		mole_wt [k].index = k_plus_1;
		char* frame = fs->get_fasta_protein ( k_plus_1, 1 );
		ProteinMW pmw ( frame );
		mole_wt [k].number = pmw.getMass ();
	}
	ujm.deletePreviousMessage ( cout );
	writeMW ( fileName, mole_wt, numEntries );
}
static void processPI ( FastaServer* fs, const string& fileName )
{
	ErrorHandler::genError ()->message ( "Calculating Protein pIs.\n" );
	int numEntries = fs->getNumEntries ();
	IndexedDouble* pi = new IndexedDouble [numEntries];

	UpdatingJavascriptMessage ujm;
	for ( int m = 0 ; m < numEntries ; m++ ) {
		int m_plus_1 = m+1;
		if ( m_plus_1 % REPORTING_FREQUENCY == 0 ) ujm.writeMessage ( cout, m_plus_1 );
		pi [m].index = m_plus_1;
		char* frame = fs->get_fasta_protein ( m_plus_1, 1 );
		ProteinPI ppi ( frame );
		pi [m].number = ppi.getProteinPI ();
	}
	ujm.deletePreviousMessage ( cout );
	writePI ( fileName, pi, numEntries );
}
static void deleteCRCharactersFromDatabaseFile ( const string& fileName )
{
	string name2 = SeqdbDir::instance ().getDatabasePath ( fileName );

	if ( !genFileExists ( name2 ) ) {
		ErrorHandler::genError ()->error ( "\nDatabase file does not exist: " + name2 + ". Function: deleteCRCharactersFromDatabaseFile\n" );
	}
	bool dosFlag = checkForCR ( name2 );
	bool macFlag = dosFlag && macTextFile ( name2, 5000 );
	if ( dosFlag ) {
		double size = (double) genFileSize ( name2 );
		double freeDiskSpace = genFreeDiskSpace ( name2 );

		if ( size < freeDiskSpace ) {
			ErrorHandler::genError ()->message ( "\nRemoving CR characters from the main database file.\n" );
			ErrorHandler::genError ()->message ( "\nThe original database file will be renamed with a .old suffix and will be deleted if the operation is successful.\n" );

			string name1 = name2 + ".old";

			genRename ( name2, name1 );

			if ( macFlag )	mac2Unix ( name1, name2 );
			else			dos2Unix ( name1, name2 );

			genUnlink ( name1 );
		}
		else ErrorHandler::genError ()->error ( "\nInsufficient disk space to proceed.\n" );
	}
}
void add_single_database_entry ( const string& filename, const string& protein, const string& name, const string& species, const string& accessionNumber )
{
	StringSizeType len = protein.length ();
	if ( strspn ( protein.c_str (), "ABCDEFGHIKLMNPQRSTUVWXYZ" ) < len ) {
		ErrorHandler::genError ()->error ( "The protein contains invalid amino acids.\n" );
	}
	AccessionNumberMap* am = getAccessionNumberMap ( filename );

	if ( am->accessionNumberUnique ( accessionNumber.c_str () ) == false ) {
		ErrorHandler::genError ()->error ( "The accession number you have used already exists in the database\n" );
	}
	delete am;

	FILE* fp = gen_fopen_binary ( SeqdbDir::instance ().getDatabasePathCreateOrAppend ( filename ), "a", "Function: add_single_database_entry" );

	try {
		wrCommentLine ( fp, filename, accessionNumber, name, species );
	}
	catch ( runtime_error e ) {		// Catch database login problems
		ErrorHandler::genError ()->error ( e );
	}

	for ( int i = 1 ; i <= len ; i++ ) {
		putc ( protein [i-1], fp );
		if ( i % 70 == 0 || i == len ) putc ( '\n', fp );
	}
	gen_fclose ( fp, "Function: add_single_database_entry" );
}
string dna_database_to_protein_database ( const string& database, bool delete_dna_database )
{
	if ( is_dna_format_database ( database ) == false ) {
		ErrorHandler::genError ()->error ( "The original database is not a DNA database.\n" );
	}
	string dnaFileName = SeqdbDir::instance ().getDatabasePath ( database );
	double size = (double) genFileSize ( dnaFileName );
	double freeDiskSpace = genFreeDiskSpace ( dnaFileName );

	if ( size * 2.0 > freeDiskSpace ) ErrorHandler::genError ()->error ( "\nInsufficient disk space to proceed.\n" );

	string pdatabase = convertDNADatabaseToProteinDatabase ( database );

	if ( delete_dna_database ) deleteDatabaseFiles ( database );	// This function must be outside convertDNADatabaseToProteinDatabase as the files need to close before they can be deleted.

	return pdatabase;
}
static string convertDNADatabaseToProteinDatabase ( const string& database )
{
	FastaServer fs ( database );
	int numEntries = fs.getNumEntries ();

	string databaseSuffix = SeqdbDir::instance ().getDatabaseSuffix ( database );
	string pdatabase = "p" + database;
	string protFileName = SeqdbDir::instance ().getSeqdbDir () + pdatabase + databaseSuffix;

	ErrorHandler::genError ()->message ( "\nCreating Protein Database " + pdatabase + ".\n" );

	FILE* fp = gen_fopen_binary ( protFileName, "w", "convertDNADatabaseToProteinDatabase" );

	UpdatingJavascriptMessage ujm;
	for ( int i = 0 ; i < numEntries ; i++ ) {
		int i_plus_1 = i+1;
		if ( i_plus_1 % REPORTING_FREQUENCY == 0 ) ujm.writeMessage ( cout, i_plus_1 );
		int length;
		char* entry = fs.get_fasta_full_entry ( i+1, &length );
		do {
			fputc ( *entry, fp );
		} while ( *entry++ != '\n' );
		for ( int frame = 1 ; frame <= 6 ; frame++ ) {
			char* protein = fs.get_fasta_protein ( i+1, frame );
			gen_fwrite ( protein, strlen ( protein ) * sizeof (char), 1, fp, "convertDNADatabaseToProteinDatabase" );
			fputc ( '\n', fp );
		}
	}
	ujm.deletePreviousMessage ( cout );
	gen_fclose ( fp, "convertDNADatabaseToProteinDatabase" );
	return pdatabase;
}
string convertDatabaseToRandomDatabase ( const string& database, bool reverseFlag, bool concatFlag )
{
	if ( is_dna_format_database ( database ) ) ErrorHandler::genError ()->error ( "The original database must be a protein database.\n" );
	string fileName = SeqdbDir::instance ().getDatabasePath ( database );
	double size = (double) genFileSize ( fileName );
	double freeDiskSpace = genFreeDiskSpace ( fileName );

	double freeSpaceFactor = ( concatFlag ) ? 2.4 : 1.2;
	if ( size * freeSpaceFactor > freeDiskSpace ) ErrorHandler::genError ()->error ( "\nInsufficient disk space to proceed.\n" );

	FastaServer fs ( database );
	int numEntries = fs.getNumEntries ();

	string databaseSuffix = SeqdbDir::instance ().getDatabaseSuffix ( database );
	string rdatabase = database;
	if ( reverseFlag ) rdatabase += ".reverse";
	else {
		rdatabase += ".random";
		srand ( static_cast<unsigned int> ( time ( 0 ) ) );
	}
	if ( concatFlag ) rdatabase += ".concat";
	string protFileName = SeqdbDir::instance ().getSeqdbDir () + rdatabase + databaseSuffix;

	ErrorHandler::genError ()->message ( "\nCreating Protein Database " + rdatabase + ".\n" );

	FILE* fp = gen_fopen_binary ( protFileName, "w", "convertDatabaseToRandomDatabase" );

	UpdatingJavascriptMessage ujm;
	if ( concatFlag ) {	// Duplicate the normal database
		for ( int i = 0 ; i < numEntries ; i++ ) {
			int i_plus_1 = i+1;
			if ( i_plus_1 % REPORTING_FREQUENCY == 0 ) ujm.writeMessage ( cout, i_plus_1 );
			int length;
			char* entry = fs.get_fasta_full_entry ( i+1, &length );
			do {
				fputc ( *entry, fp );
			} while ( *entry++ != '\n' );
			char* protein = fs.get_fasta_protein ( i+1, 1 );
			int len = strlen ( protein );
			if ( len ) {	// Account for the empty protein case
				gen_fwrite ( protein, len * sizeof (char), 1, fp, "convertDatabaseToRandomDatabase" );
				fputc ( '\n', fp );
			}
		}
	}
	for ( int i = 0 ; i < numEntries ; i++ ) {
		int i_plus_1 = i+1;
		if ( i_plus_1 % REPORTING_FREQUENCY == 0 ) ujm.writeMessage ( cout, i_plus_1 );
		int length;
		if ( concatFlag ) {
			try {
				wrCommentLine ( fp, database, string ( "-" ) + fs.getAccessionNumber ( i+1 ), fs.getName ( i+1 ), fs.getSpecies ( i+1 ) );
			}
			catch ( runtime_error e ) {		// Error so delete the database created so far
				gen_fclose ( fp, "convertDatabaseToRandomDatabase" );
				genUnlink ( protFileName );
				ErrorHandler::genError ()->error ( e );
			}
		}
		else {
			char* entry = fs.get_fasta_full_entry ( i+1, &length );
			do {
				fputc ( *entry, fp );
			} while ( *entry++ != '\n' );
		}
		char* protein = fs.get_fasta_protein ( i+1, 1 );
		int len = strlen ( protein );
		if ( reverseFlag ) reverse ( protein, protein + len );
		else random_shuffle ( protein, protein + len );
		if ( len ) {	// Account for the empty protein case
			gen_fwrite ( protein, len * sizeof (char), 1, fp, "convertDatabaseToRandomDatabase" );
			fputc ( '\n', fp );
		}
	}
	ujm.deletePreviousMessage ( cout );
	gen_fclose ( fp, "convertDatabaseToRandomDatabase" );
	return rdatabase;
}
static void deleteDatabaseFiles ( const string& database )
{
	genUnlink ( SeqdbDir::instance ().getDatabasePath ( database ) );
	deleteDatabaseIndexFiles ( database );
}
static void deleteDatabaseIndexFiles ( const string& database )
{
	static string database_extension [] = {
		"acc", "acn", "idx", "idc", "idp", "idi", "mw", "pi", "sl", "sp", "tax", "tl", "unk", "unr", "usp", "END"
	};
	for ( int i = 0 ; database_extension [i] != "END" ; i++ ) {
		string fileName = SeqdbDir::instance ().getSeqdbDir () + database + string ( "." ) + database_extension [i];
		genUnlink ( fileName );
	}
}
static void writeTaxAndTL ( const string& fileName, VectorIndexedInt& species )
{
	ErrorHandler::genError ()->message ( "Sorting species.\n" );
	sort ( species.begin (), species.end (), SortIndexedIntAscending () );
	species.erase ( unique ( species.begin (), species.end () ), species.end () );		// This is the list of sorted nodes
	ErrorHandler::genError ()->message ( "Creating taxonomy (.tax) and taxonomy list (.tl) files.\n" );
	GenOFStream ost1 ( SeqdbDir::instance ().getSeqdbDir () + fileName + ".tax" );
	GenOFStream ost2 ( SeqdbDir::instance ().getSeqdbDir () + fileName + ".tl" );
	int jj;

	for ( VectorIndexedIntSizeType ii = 0 ; ii < species.size () ; ) {
		int numRepeats = 1;
		for ( jj = 1 ; ii + jj < species.size () && ( species [ii].number == species [ii+jj].number ) ; jj++ ) {
			numRepeats++;
		}
		ost1 << species [ii].number << ' ' << numRepeats << '\n';
		ost2 << species [ii].number << ' ' << numRepeats << '\n';
		for ( jj = 0 ; jj < numRepeats ; jj++ ) {
			ost1 << species [ii++].index << '\n';
		}
	}
}
static void writeUnk ( const string& fileName, const SetString& unkSpecies )
{
	if ( !unkSpecies.empty () ) {
		GenOFStream ost ( SeqdbDir::instance ().getSeqdbDir () + fileName + ".unk" );
		for ( SetStringConstIterator xx = unkSpecies.begin () ; xx != unkSpecies.end () ; xx++ ) {
			ost << *xx << '\n';
		}
		ost.close ();
		ErrorHandler::genError ()->message ( gen_itoa ( unkSpecies.size () ) + " species names are not present in the taxonomy files. These are listed in the .unk file." );
	}
}
static void writeACN ( const string& fileName, vector <IndexedInt>& numericAccessionNumber )
{
	ErrorHandler::genError ()->message ( "Sorting accession numbers.\n" );
	sort ( numericAccessionNumber.begin (), numericAccessionNumber.end (), SortIndexedIntAscending () );	// sort causes dbEST to hang
	ErrorHandler::genError ()->message ( "Creating accession number file (.acn).\n" );
	GenOFStream ost ( SeqdbDir::instance ().getSeqdbDir () + fileName + ".acn", ios_base::binary );
	if ( !numericAccessionNumber.empty () ) {
 		ost.write ( (char*) &numericAccessionNumber [0], numericAccessionNumber.size () * sizeof (IndexedInt) );
	}
}
static void writeACC ( const string& fileName, VectorIndexedConstCString& accessionNumber )
{
	ErrorHandler::genError ()->message ( "Sorting accession numbers.\n" );
	sort ( accessionNumber.begin (), accessionNumber.end (), SortIndexedStrcasecmpAscending () );
	ErrorHandler::genError ()->message ( "Creating accession number file (.acc).\n" );
	GenOFStream ost ( SeqdbDir::instance ().getSeqdbDir () + fileName + ".acc" );
	for ( VectorIndexedConstCStringSizeType i = 0 ; i < accessionNumber.size () ; i++ ) {
		ost << accessionNumber [i].name << ' ' << accessionNumber [i].index << endl;
	}
	ost.close ();
	for ( VectorIndexedConstCStringSizeType ii = 0 ; ii < accessionNumber.size () ; ii++ ) {
		delete [] const_cast <char*> (accessionNumber [ii].name);
	}
}
static void writeMW ( const string& fileName, IndexedDouble* mole_wt, int numEntries )
{
	ErrorHandler::genError ()->message ( "Sorting protein molecular weights.\n" );
	sort ( mole_wt, mole_wt + numEntries, SortIndexedDoubleAscending () );
	ErrorHandler::genError ()->message ( "Creating protein molecular weight file (.mw).\n" );
	GenOFStream ost ( SeqdbDir::instance ().getSeqdbDir () + fileName + ".mw", ios_base::binary );
 	ost.write ( (char*) mole_wt, numEntries * sizeof (IndexedDouble) );
	ost.close ();
	delete [] mole_wt;
}
static void writePI ( const string& fileName, IndexedDouble* pi, int numEntries )
{
	ErrorHandler::genError ()->message ( "Sorting protein pIs.\n" );
	sort ( pi, pi + numEntries, SortIndexedDoubleAscending () );
	ErrorHandler::genError ()->message ( "Creating protein pI file (.pi).\n" );
	GenOFStream ost ( SeqdbDir::instance ().getSeqdbDir () + fileName + ".pi", ios_base::binary );
 	ost.write ( (char*) pi, numEntries * sizeof (IndexedDouble) );
	ost.close ();
	delete [] pi;
}
