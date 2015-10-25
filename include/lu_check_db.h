/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_check_db.h                                                 *
*                                                                             *
*  Created    : December 3rd 2002                                             *
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
*  Copyright (2002-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_check_db_h
#define __lu_check_db_h

#include <string>
#include <lgen_define.h>

std::string getDatabasePrefix ( const std::string& filename );
std::string getBestSubstituteDatabase ( const std::string& filename );
bool is_genpept_database ( const std::string& fileName );
bool is_swissprot_database ( const std::string& fileName );
bool is_ludwig_database ( const std::string& fileName );
bool is_owl_database ( const std::string& fileName );
bool is_ncbinr_database ( const std::string& fileName );
bool is_ipi_database ( const std::string& fileName );
bool is_uniprot_database ( const std::string& fileName );
bool is_user_protein_database ( const std::string& fileName );
bool is_est_database ( const std::string& fileName );
bool is_generic_database ( const std::string& fileName );
bool is_default_database ( const std::string& fileName );
bool is_pdna_format_database ( const std::string& fileName );
bool is_dbest_database ( const std::string& fileName );
bool is_generic_string_acc_number_database ( const std::string& fileName );
bool is_generic_numeric_acc_number_database ( const std::string& fileName );
bool is_dna_format_database	( const std::string& fileName );
bool is_numeric_acc_number_database ( const std::string& fileName );
bool is_dna_database ( const std::string& fileName );
bool isFullyDecoyDatabase ( const std::string& database );
PairIntInt getDNAFrameTranslationPair ( const std::string& db, int dft );
bool getConcatDBPair ( const std::string& db, PairStringString& pss );

class DBSearchFlags {
	bool concatFlag;
	bool randomFlag;
	bool reverseFlag;
	bool userFlag;
public:
	DBSearchFlags ( const StringVector& dbList );
	bool getConcatFlag () const { return concatFlag; }
	bool getRandomFlag () const { return randomFlag; }
	bool getReverseFlag () const { return reverseFlag; }
	bool getUserFlag () const { return userFlag; } 
	bool getRandomAndReverseFlag () const { return randomFlag && reverseFlag; } 
};

#endif /* ! __lu_check_db_h */
