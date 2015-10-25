/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_hom_par.cpp                                                *
*                                                                             *
*  Created    : September 3rd 2001                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lgen_error.h>
#include <lu_pq_vector.h>
#include <lu_hom_par.h>
#include <lu_param_list.h>
using std::ostream;
using std::string;
using std::runtime_error;

MSHomologyParameters::MSHomologyParameters ( const ParameterList* params ) :

	MSSearchParameters ( params ),

	minMatches		( params->getIntValue	( "min_matches", 1 ) ),
	scoreMatrix		( params->getStringValue	( "score_matrix", "BLOSUM62MS" ) ),
	maxPeptideHits	( params->getIntValue	( "max_hits", 2000000 ) ),

	productMassTolerance( "fragment_masses", params ),
	previousAA			( params->getIntValue	( "previous_aa", 1 ) ),
	nextAA				( params->getIntValue	( "next_aa", 1 ) )
{
	const char* value;
	if ( params->getValue ( "possible_sequences", value ) ) {
		try {
			StringVector temp1, temp2, temp3;
			getPostQuery3Vectors ( value, temp1, temp2, temp3 );	// Should throw exception if more than 3 fields entered

			bool flag = getPostQuery3Vectors ( value, sequenceSet, possibleSequences, maxSeqErrors );
			if ( !flag ) {
				flag = getPostQuery2Vectors ( value, possibleSequences, maxSeqErrors );		// Assume sequence set not entered
				if ( flag ) {
					sequenceSet.assign ( possibleSequences.size (), 1 );					// Assume sequence set is 1.
				}
				else {
					getPostQueryVector ( value, possibleSequences, '\n' );	// Assume sequence set and possible errors not entered
					sequenceSet.assign ( possibleSequences.size (), 1 );	// Assume sequence set is 1.
					maxSeqErrors.assign ( possibleSequences.size (), 0 );	// Assume no sequence errors.
				}
			}
		}
		catch ( runtime_error ) {
			ErrorHandler::genError ()->error ( "Possible sequences field is incorrectly formatted.\n" );
		}
		if ( possibleSequences.empty () ) {
			ErrorHandler::genError ()->error ( "Possible sequences field contains no correctly formatted sequences.\n" );
		}
		for ( StringVectorSizeType i = 0 ; i < possibleSequences.size () ; i++ ) {
			string s = possibleSequences [i];
			static string legalCharacters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789.[]{}|";
			string::size_type ind = s.find_first_not_of ( legalCharacters );	// Check for illegal characters in sequence
			if ( ind != string::npos ) {
				ErrorHandler::genError ()->error (
					"Possible sequences field contains the illegal character: '" + string ( 1, s[ind] ) + "'" );
			}
			if ( !genAnyUpper ( s ) ) possibleSequences [i] = genToUpper ( s );	// Try converting to upper case if all lower case
		}
	}
}
void MSHomologyParameters::printHTML ( ostream& os ) const
{
	MSSearchParameters::printHTML ( os );
	ParameterList::printHTML ( os, "Min matches", minMatches );
	ParameterList::printHTML ( os, "Score matrix", scoreMatrix );
	if ( possibleSequences.size () ) {
		os << "List of Sequences: ";
		os << "<b>";
		for ( StringVectorSizeType i = 0 ; i < possibleSequences.size () ; i++ ) {
			os << sequenceSet [i] << " " << possibleSequences [i] << " " << maxSeqErrors [i] << "<br />";
		}
		os << "</b>";
	}
	productMassTolerance.getTolerance ()->printHTML ( os, "Mass Tolerance" );
}
