/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_hom_srch.cpp                                               *
*                                                                             *
*  Created    : Septmber 5th 2001                                             *
*                                                                             *
*  Purpose    : MS-Homology search and hit classes.                           *
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
#include <lp_frame.h>
#include <lu_mat_score.h>
#include <lu_hom_srch.h>
#include <lu_hom_par.h>
#include <lu_unk_rexp.h>
#include <lu_param_list.h>
#include <lu_table.h>
#include <lu_delim.h>
using std::string;
using std::ostream;
using std::endl;
using std::sort;

int HomologyHit::numPreviousAA = 1;
int HomologyHit::numNextAA = 1;

HomologyHit::HomologyHit ( FastaServer* fs, int index, int dnaReadingFrame, int openReadingFrame, const char* protein, const string& measuredSequence, int sequenceSet, int numSubstitutions ) :
	ProteinHit ( fs, index, dnaReadingFrame, openReadingFrame ),
	sequence ( get_regular_expression_match () ),
	previousAA ( getPreviousAA ( protein ) ),
	nextAA ( getNextAA () ),
	startAA ( get_regular_expression_loc1 () - protein + 1 ),
	measuredSequence ( measuredSequence ),
	sequenceSet ( sequenceSet ),
	numSubstitutions ( numSubstitutions ),
	scoreSet ( false )
{
}
double HomologyHit::getScore () const
{
	if ( scoreSet == false ) {
		score = matrix_score ( measuredSequence.c_str (), sequence.c_str () );
		scoreSet = true;
	}
	return score;
}
void HomologyHit::printHTMLHeader ( ostream& os ) const
{
	tableRowStart ( os );
		tableHeader ( os, "Protein Score" );
		tableHeader ( os, "Peptide Score" );
		tableHeader ( os, "Peptide Sequence" );
		tableHeader ( os, "Matching Sequence" );
		tableHeader ( os, "Start AA" );
		ProteinHit::printHTMLHeader ( os );
	tableRowEnd ( os );
}
void HomologyHit::printDelimitedHeader ( ostream& os ) const
{
	delimitedRowStart ( os );
		delimitedHeader ( os, "Protein Score" );
		delimitedHeader ( os, "Peptide Score" );
		delimitedHeader ( os, "Peptide Sequence" );
		delimitedHeader ( os, "Previous AA" );
		delimitedHeader ( os, "Matching Sequence" );
		delimitedHeader ( os, "Next AA" );
		delimitedHeader ( os, "Start AA" );
		ProteinHit::printDelimitedHeader ( os );
	delimitedRowEnd ( os );
}
void HomologyHit::printHTML ( ostream& os, bool newProtein ) const
{
	tableRowStart ( os );
		if ( newProtein ) {
			tableHeaderStart ( os );
				genPrintSigFig ( os, proteinScore, 0 );
				os << endl;
			tableHeaderEnd ( os );
		}
		else tableEmptyCell ( os );
		tableHeaderStart ( os );
			genPrintSigFig ( os, score, 0 );
			os << endl;
		tableHeaderEnd ( os );
		tableCell ( os, measuredSequence );
		tableDataStart ( os, "", "", true );
			os << "(" << previousAA << ")";
			int len = sequence.length ();
			for ( int j = 0 ; j < len ; j++ ) {
				if ( sequence [j] == measuredSequence [j] )
					os << sequence [j];
				else
					os << "<font size=\"+1\"><b>" << sequence [j] << "</b></font>";
			}
			os << "(" << nextAA << ")";
		tableDataEnd ( os );
		tableDataStart ( os, "", "center" );
			os << startAA;
		tableDataEnd ( os );
		if ( newProtein )
			printHTMLHit2 ( os );
	tableRowEnd ( os );
}
void HomologyHit::printDelimited ( ostream& os, bool newProtein ) const
{
	delimitedRowStart ( os );
		if ( newProtein ) {
			delimitedCellSigFig ( os, proteinScore, 0 );
		}
		else delimitedEmptyCell ( os );
		delimitedCellSigFig ( os, score, 0 );
		delimitedCell ( os, measuredSequence );
		delimitedCell ( os, previousAA );
		delimitedCell ( os, sequence );
		delimitedCell ( os, nextAA );
		delimitedCell ( os, startAA );
		if ( newProtein )
			printDelimitedHit2 ( os );
	delimitedRowEnd ( os );
}
void HomologyHit::printXML ( ostream& os, bool newProtein ) const
{
	if ( newProtein ) {
		ParameterList::printXML ( os, "protein_score", proteinScore );
		printXMLHit ( os );
	}
	os << "<peptide>" << endl;
		ParameterList::printXML ( os, "peptide_score", score );
		ParameterList::printXML ( os, "measured_sequence", measuredSequence );
		ParameterList::printXML ( os, "previous_aa", previousAA );
		ParameterList::printXML ( os, "database_sequence", sequence );
		ParameterList::printXML ( os, "next_aa", nextAA );
		ParameterList::printXML ( os, "start_aa", startAA );
	os << "</peptide>" << endl;
}
string HomologyHit::getPreviousAA ( const char* protein ) const
{
	const char* loc1 = get_regular_expression_loc1 ();
	int startAA = loc1 - protein;
	int numDisp = genMin ( numPreviousAA, startAA );
	string ret;
	if ( numDisp < numPreviousAA ) ret += string ( "-" );
	for ( int i = numDisp ; i-- ; ) {
		ret += *(loc1-i-1);
	}
	return ret;
}
string HomologyHit::getNextAA () const
{
	const char* loc2 = get_regular_expression_loc2 ();
	string ret;
	int i = 0;
	for ( ; i < numNextAA ; i++ ) {
		char c = loc2 [i];
		if ( c == 0 ) break;
		ret += c;
	}
	if ( i < numNextAA ) ret += '-';
	return ret;
}
void HomologyHit::setNumPreviousAA ( int num )
{
	numPreviousAA = num;
}
void HomologyHit::setNumNextAA ( int num )
{
	numNextAA = num;
}
void HomologyHits::calculateScores ( int minMatches, const string& scoreMatrix )
{
	sort ( peptideHits.begin (), peptideHits.end (), sort_multi_sequence_list_hits () );
	int numUniqueHits;
	int start_i;
	int i, j, k;
	int numPeptideHits = peptideHits.size ();

	set_score_matrix ( scoreMatrix );
	for ( i = 0, k = 0 ; i < numPeptideHits ; i++ ) {
		if ( i == 0 || peptideHits [i].getIndex () != peptideHits [i-1].getIndex () ) {
			numUniqueHits = 1;
			start_i = i;
		}
		else {
			if ( peptideHits [i].getSequenceSet () != peptideHits [i-1].getSequenceSet () ) numUniqueHits++;
		}
		if ( i == numPeptideHits - 1 || peptideHits [i].getIndex () != peptideHits [i+1].getIndex () ) {
			if ( numUniqueHits >= minMatches ) {
				int end_i = i;
				double proteinScore = 0.0;
				double maxScore = 0.0;
				for ( j = start_i ; j <= end_i ; j++ ) {

					maxScore = genMax ( peptideHits [j].getScore (), maxScore );
					if ( j == end_i || peptideHits [j].getSequenceSet () != peptideHits [j+1].getSequenceSet () ) {
						proteinScore += maxScore;
						maxScore = 0.0;
					}
				}
				for ( j = start_i ; j <= end_i ; j++ ) {
					peptideHits [j].proteinScore = proteinScore;
					peptideHits [k++] = peptideHits [j];
				}
				HomologyHit& h = peptideHits [i];
				proteinHits.push_back ( ProteinHomologyHit ( h.getDBIndex (), h.getIndex (), h.getAccessionNumber (), h.getDNAReadingFrame (), h.getOpenReadingFrame (), proteinScore ) );
			}
		}
	}
	peptideHits.erase ( peptideHits.begin () + k, peptideHits.end () );
	sort ( peptideHits.begin (), peptideHits.end (), sort_hits_by_score () );
}
HomologySearch::HomologySearch ( const MSHomologyParameters& params ) :
	DatabaseSearch ( params ),
	homologyParams ( params ),
	maxPeptideHits ( params.getMaxPeptideHits () )
{
	HomologyHit::setNumPreviousAA ( params.getPreviousAA () );
	HomologyHit::setNumNextAA ( params.getNextAA () );
	homologyHits = new HomologyHits ();
	initialise_regular_expression_hit_list ( params.getPossibleSequences (), params.getSequenceSet (), params.getMaxSeqErrors (), params.getProductMassTolerance (), !params.getNoEnzyme () );
	sequenceSearch ();
	homologyHits->calculateScores ( params.getMinMatches (), params.getScoreMatrix () );

	numHits = homologyHits->size ();
	databaseHits = homologyHits;
}
void HomologySearch::sequenceSearch ()
{
	for ( int i = 0 ; i < fs.size () ; i++ ) {
		doSearch ( fs [i], i );
		FrameIterator::resetElapsedTime ( 1 );
	}
}
void HomologySearch::doSearch ( FastaServer* fsPtr, int i )
{
	bool noEnzyme = params.getNoEnzyme ();
	ProteinHit::addFS ( fsPtr, i );
	FrameIterator fi ( fsPtr, params.getIndicies ( i ), dnaFrameTranslationPairVector [i], params.getTempOverride () );

	if ( !noEnzyme ) init_fasta_enzyme_function ( params.getEnzyme () );
	char* reading_frame;
	while ( ( reading_frame = fi.getNextFrame () ) != NULL ) {
		if ( noEnzyme ) {
			while ( multi_unknome_hits_present ( reading_frame ) ) {
				add_hit ( fsPtr, fi.getEntry (), fi.getFrameTranslation (), fi.getFrame (), reading_frame );
			}
		}
		else {
			const IntVector& cleavageIndex = enzyme_fragmenter ( reading_frame );

			for ( IntVectorSizeType m = 0 ; m < cleavageIndex.size () ; m++ ) {
				int offset = ( m == 0 ) ? 0 : cleavageIndex [m-1] + 1;
				int n = m;
				if ( multi_unknome_hits_present ( reading_frame + offset ) ) {
					const char* loc = get_regular_expression_loc2 ();
					while ( loc > reading_frame + cleavageIndex [n] + 1 ) n++;
					if ( loc == reading_frame + cleavageIndex [n] + 1 ) 
						add_hit ( fsPtr, fi.getEntry (), fi.getFrameTranslation (), fi.getFrame (), reading_frame );
				}
			}
		}
	}
}
void HomologySearch::add_hit ( FastaServer* fs, int index, int dnaReadingFrame, int openReadingFrame, char* protein )
{
	char* measuredSequence;
	int sequenceSet;
	int numSubstitutions;
	bool reset = true;

	while ( get_unknome_rexp_results ( &measuredSequence, &sequenceSet, &numSubstitutions, reset ) ) {
		reset = false;
		homologyHits->addHit ( HomologyHit ( fs, index, dnaReadingFrame, openReadingFrame, protein, measuredSequence, sequenceSet, numSubstitutions ) );
	}
	if ( homologyHits->peptideHitsSize () > maxPeptideHits ) {
		ErrorHandler::genError ()->error ( "The maximum number of hits has been exceeded.\n" );
	}
}
void HomologySearch::printHTMLHitsTable ( ostream& os ) const
{
	int numPeptideHits = homologyHits->peptideHitsSize ();
	int start_i;

	os << "Number of Peptide Hits: <b>" << numPeptideHits << "</b>" << endl;

	os << "<table cellspacing=\"3\">" << endl;
	for ( int i = 0, j = 0 ; i < numPeptideHits ; i++ ) {
		if ( i == 0 ) {
			homologyHits->printHTMLHeader ( os, i );
		}
		if ( i == 0 || homologyHits->differentProteins ( i, i-1 ) ) {
			start_i = i;
		}
		if ( i == numPeptideHits - 1 || homologyHits->differentProteins ( i, i+1 ) ) {
			for ( int k = start_i ; k <= i ; k++ ) {
				homologyHits->printHTMLHit2 ( os, k == start_i, k );
			}
			j++;
			if ( j >= numHits ) break;
		}
	}
	os << "</table>" << endl;
}
void HomologySearch::printDelimitedHits ( ostream& os ) const
{
	int numPeptideHits = homologyHits->peptideHitsSize ();
	int start_i;
	for ( int i = 0, j = 0 ; i < numPeptideHits ; i++ ) {
		if ( i == 0 ) {
			homologyHits->printDelimitedHeader ( os, i );
		}
		if ( i == 0 || homologyHits->differentProteins ( i, i-1 ) ) {
			start_i = i;
		}
		if ( i == numPeptideHits - 1 || homologyHits->differentProteins ( i, i+1 ) ) {
			for ( int k = start_i ; k <= i ; k++ ) {
				homologyHits->printDelimitedHit2 ( os, k == start_i, k );
			}
			j++;
			if ( j >= numHits ) break;
		}
	}
}
void HomologySearch::printXMLHits ( ostream& os ) const
{
	int numPeptideHits = homologyHits->peptideHitsSize ();
	int start_i;
	os << "<mshomology_results>" << endl;
	for ( int i = 0, j = 0 ; i < numPeptideHits ; i++ ) {
		if ( i == 0 || homologyHits->differentProteins ( i, i-1 ) ) {
			start_i = i;
		}
		if ( i == numPeptideHits - 1 || homologyHits->differentProteins ( i, i+1 ) ) {
			os << "<hit>" << endl;
			for ( int k = start_i ; k <= i ; k++ ) {
				homologyHits->printXMLHit ( os, k == start_i, k );
			}
			os << "</hit>" << endl;
			j++;
			if ( j >= numHits ) break;
		}
	}
	os << "</mshomology_results>" << endl;
}
void HomologySearch::printParamsBodyHTML ( ostream& os ) const
{
	homologyParams.printHTML ( os );
}
