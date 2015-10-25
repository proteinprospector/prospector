/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_check_db.cpp                                               *
*                                                                             *
*  Created    : November 26th 2003                                            *
*                                                                             *
*  Purpose    : Functions for dealing with multiply charged data.             *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1997-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lg_string.h>
#include <lgen_file.h>
#include <lgen_error.h>
#include <lu_getfil.h>
#include <lu_check_db.h>
using std::make_pair;
using std::string;
static const string GENPEPT1	= "gen";
static const string GENPEPT2	= "Genpept";
static const string SWISSPROT1	= "swp";
static const string SWISSPROT2	= "SwissProt";
static const string UNIPROT		= "UniProt";
static const string LUDWIGNR	= "Ludwignr";
static const string OWL1		= "owl";
static const string OWL2		= "Owl";
static const string NCBINR1		= "NCBInr";
static const string NCBINR2		= "nr";
static const string IPI1		= "ipi";
static const string IPI2		= "IPI";
static const string DBEST1		= "dbest";
static const string DBEST2		= "dbEST";
static const string DA			= "DA";
static const string DN			= "DN";
static const string PA			= "PA";
static const string PN			= "PN";
static const string PDA			= "pDA";
static const string PDN			= "pDN";
static const string PDBEST1		= "pdbest";
static const string PDBEST2		= "pdbEST";
static const string DDEFAULT	= "Ddefault";
static const string PDDEFAULT	= "pDdefault";
static const string PDEFAULT	= "Pdefault";
static const string USERPROTEIN	= "UserProtein";

string getDatabasePrefix ( const string& filename )
{
	string prefix;
	if ( isPrefix ( filename, GENPEPT1 ) )		prefix = GENPEPT1;
	if ( isPrefix ( filename, GENPEPT2 ) )		prefix = GENPEPT2;
	if ( isPrefix ( filename, SWISSPROT1 ) )	prefix = SWISSPROT1;
	if ( isPrefix ( filename, SWISSPROT2 ) )	prefix = SWISSPROT2;
	if ( isPrefix ( filename, UNIPROT ) )		prefix = UNIPROT;
	if ( isPrefix ( filename, LUDWIGNR ) )		prefix = LUDWIGNR;
	if ( isPrefix ( filename, OWL1 ) )			prefix = OWL1;
	if ( isPrefix ( filename, OWL2 ) )			prefix = OWL2;
	if ( isPrefix ( filename, NCBINR1 ) )		prefix = NCBINR1;
	if ( isPrefix ( filename, NCBINR2 ) )		prefix = NCBINR2;
	if ( isPrefix ( filename, DBEST1 ) )		prefix = DBEST1;
	if ( isPrefix ( filename, DBEST2 ) )		prefix = DBEST2;
	if ( isPrefix ( filename, DA ) )			prefix = DA;
	if ( isPrefix ( filename, DN ) )			prefix = DN;
	if ( isPrefix ( filename, PA ) )			prefix = PA;
	if ( isPrefix ( filename, PN ) )			prefix = PN;
	if ( isPrefix ( filename, PDA ) )			prefix = PDA;
	if ( isPrefix ( filename, PDN ) )			prefix = PDN;
	if ( isPrefix ( filename, IPI1 ) )			prefix = IPI1;
	if ( isPrefix ( filename, IPI2 ) )			prefix = IPI2;
	if ( isPrefix ( filename, PDBEST1 ) )		prefix = PDBEST1;
	if ( isPrefix ( filename, PDBEST2 ) )		prefix = PDBEST2;
	if ( isPrefix ( filename, DDEFAULT ) )		prefix = DDEFAULT;
	if ( isPrefix ( filename, PDDEFAULT ) )		prefix = PDDEFAULT;
	if ( isPrefix ( filename, PDEFAULT ) )		prefix = PDEFAULT;
	return prefix;
}
string getBestSubstituteDatabase ( const string& filename )
{
	string prefix = getDatabasePrefix ( filename );
	bool ranFilename = isSuffix ( filename, ".random" );
	bool revFilename = isSuffix ( filename, ".reverse" );
	if ( prefix.empty () ) {
		ErrorHandler::genError ()->error ( "Invalid database filename prefix.\nSee the manual for information on naming databases." );
	}
	string bestSubstituteFilename;
	FileList fList ( SeqdbDir::instance ().getSeqdbDir (), prefix, ".idc", true );
	GENINT64 maxFileSize = 0;
	for ( FileList::size_type i = 0 ; i < fList.size () ; i++ ) {
		string f = fList [i];
		bool ranF = isSuffix ( f, ".random" );
		bool revF = isSuffix ( f, ".reverse" );
		bool conF = isSuffix ( f, ".concat" );
		if ( !ranF && !revF && !conF ) {
			GENINT64 size = genFileSize ( SeqdbDir::instance ().getDatabasePath ( f ) );
			if ( size > maxFileSize ) {
				maxFileSize = size;
				bestSubstituteFilename = f;
			}
		}
	}
	return bestSubstituteFilename;
}
bool is_genpept_database ( const string& fileName )
{
	return ( isPrefix ( fileName, GENPEPT1 ) || isPrefix ( fileName, GENPEPT2 ) );
}
bool is_swissprot_database ( const string& fileName )
{
	return ( isPrefix ( fileName, SWISSPROT1 ) || isPrefix ( fileName, SWISSPROT2 ) );
}
bool is_ludwig_database ( const string& fileName )
{
	return ( isPrefix ( fileName, LUDWIGNR ) );
}
bool is_owl_database	( const string& fileName )
{
	return ( isPrefix ( fileName, OWL1 ) || isPrefix ( fileName, OWL2 ) );
}
bool is_ncbinr_database ( const string& fileName )
{
	return ( isPrefix ( fileName, NCBINR1 ) || isPrefix ( fileName, NCBINR2 ) );
}
bool is_ipi_database ( const string& fileName )
{
	return ( isPrefix ( fileName, IPI1 ) || isPrefix ( fileName, IPI2 ) );
}
bool is_uniprot_database ( const string& fileName )
{
	return ( isPrefix ( fileName, UNIPROT ) );
}
bool is_user_protein_database ( const string& fileName )
{
	return ( isPrefix ( fileName, USERPROTEIN ) );
}
bool is_dbest_database ( const string& fileName )
{
	return ( isPrefix ( fileName, DBEST1 ) || isPrefix ( fileName, DBEST2 ) );
}

static bool is_da_database ( const string& fileName )
{
	return ( isPrefix ( fileName, DA ) );
}
static bool is_dn_database ( const string& fileName )
{
	return ( isPrefix ( fileName, DN ) );
}
static bool is_pa_database ( const string& fileName )
{
	return ( isPrefix ( fileName, PA ) );
}
static bool is_pn_database ( const string& fileName )
{
	return ( isPrefix ( fileName, PN ) );
}

static bool is_pda_database	( const string& fileName )
{
	return ( isPrefix ( fileName, PDA ) );
}
static bool is_pdn_database	( const string& fileName )
{
	return ( isPrefix ( fileName, PDN ) );
}
static bool is_pdbest_database ( const string& fileName )
{
	return ( isPrefix ( fileName, PDBEST1 ) || isPrefix ( fileName, PDBEST2 ) );
}

bool is_est_database ( const string& fileName )
{
	return ( is_dbest_database ( fileName ) || is_pdbest_database ( fileName ) );
}

static bool is_generic_dna_database	( const string& fileName )
{
	return ( is_da_database ( fileName ) || is_dn_database ( fileName ) );
}
static bool is_generic_pdna_database ( const string& fileName )
{
	return ( is_pda_database ( fileName ) || is_pdn_database ( fileName ) );
}
static bool is_generic_protein_database ( const string& fileName )
{
	return ( is_pa_database ( fileName ) || is_pn_database ( fileName ) );
}
bool is_generic_database ( const string& fileName )
{
	return ( is_generic_dna_database ( fileName ) || is_generic_pdna_database ( fileName ) || is_generic_protein_database ( fileName ) );
}

static bool is_default_dna_database	( const string& fileName )
{
	return ( isPrefix ( fileName, DDEFAULT ) );
}
static bool is_default_pdna_database ( const string& fileName )
{
	return ( isPrefix ( fileName, PDDEFAULT ) );
}
static bool is_default_protein_database ( const string& fileName )
{
	return ( isPrefix ( fileName, PDEFAULT ) );
}
bool is_default_database ( const string& fileName )
{
	return ( is_default_dna_database ( fileName ) || is_default_pdna_database ( fileName ) || is_default_protein_database ( fileName ) || is_user_protein_database ( fileName ) );
}

bool is_generic_string_acc_number_database ( const string& fileName )
{
	return ( is_da_database ( fileName ) || is_pa_database ( fileName ) );
}
bool is_generic_numeric_acc_number_database ( const string& fileName )
{
	return ( is_dn_database ( fileName ) || is_pn_database ( fileName ) );
}

bool is_dna_format_database	( const string& fileName )
{
	return ( is_dbest_database ( fileName ) || is_generic_dna_database ( fileName ) || is_default_dna_database ( fileName ) );
}
bool is_pdna_format_database	( const string& fileName )
{
	return ( is_pdbest_database ( fileName ) || is_generic_pdna_database ( fileName ) || is_default_pdna_database ( fileName ) );
}

bool is_numeric_acc_number_database ( const string& fileName )
{
	return ( is_est_database ( fileName ) || is_ncbinr_database ( fileName ) || is_generic_numeric_acc_number_database ( fileName ) || is_default_database ( fileName ) || is_genpept_database ( fileName ) );
}
bool is_dna_database ( const string& fileName )
{
	return ( is_dna_format_database ( fileName ) || is_pdna_format_database ( fileName ) );
}
bool isFullyDecoyDatabase ( const string& database )
{
	bool concat = isSuffix ( database, ".concat" );
	bool random = database.find ( ".random" ) != string::npos;
	bool reverse = database.find ( ".reverse" ) != string::npos;
	if ( ( !concat && ( random || reverse ) ) || database == "User Protein Random" || database == "User Protein Reverse" )
		return true;
	else
		return false;
}
PairIntInt getDNAFrameTranslationPair ( const string& db, int dft )
{
	int start = 1;
	int end = 1;
	if ( is_dna_database ( db ) ) {
		start	= ( dft > 0 ) ? 1 : 4;
		end		= ( dft > 0 ) ? dft : ( dft == -1 ) ? 4 : 6;
	}
	return make_pair ( start, end );
}
bool getConcatDBPair ( const string& db, PairStringString& pss )
{
	if ( isSuffix ( db, ".concat" ) && !genFileExists ( SeqdbDir::instance ().getDatabasePath ( db ) ) ) {
		string d1 = db.substr ( 0, db.length () - 7 );
		int len = isSuffix ( d1, ".random" ) ? 7 : 8; 
		string d2 = d1.substr ( 0, d1.length () - len );
		pss = make_pair ( d2, d1 );
		return true;
	}
	return false;
}

DBSearchFlags::DBSearchFlags ( const StringVector& dbList ) :
	concatFlag ( false ),
	randomFlag ( false ),
	reverseFlag ( false ),
	userFlag ( false )
{
	for ( StringVectorSizeType i = 0 ; i < dbList.size () ; i++ ) {
		const string& d = dbList [i];
		if ( isSuffix ( d, ".concat" ) )			concatFlag = true;
		if ( d.find ( ".random" ) != string::npos )	randomFlag = true;
		if ( d.find ( ".reverse" ) != string::npos )reverseFlag = true;
		if ( d == "User Protein" )					userFlag = true;
	}
}
