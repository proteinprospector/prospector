/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_patt_srch.cpp                                              *
*                                                                             *
*  Created    : Septmber 5th 2001                                             *
*                                                                             *
*  Purpose    : MS-Pattern search and hit classes.                            *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lgen_error.h>
#include <lgen_reg_exp.h>
#include <lp_frame.h>
#include <lu_patt_par.h>
#include <lu_patt_srch.h>
#include <lu_param_list.h>
#include <lu_inst.h>
#include <lu_table.h>
using std::ostream;
using std::endl;
using std::string;

PatternHit::PatternHit ( FastaServer* fs, int index, int dnaReadingFrame, int openReadingFrame, RegularExpression* regExp ) :
	ProteinHit ( fs, index, dnaReadingFrame, openReadingFrame ),
	regExp ( regExp ),
	sequence ( regExp->getMatch () ),
	previousAA ( regExp->getPrevious () ),
	nextAA ( regExp->getNext () ),
	startAA ( regExp->getStartAA () ),
	matchedPepMW ( peptide_formula_to_molecular_weight ( sequence.c_str () ) ),
	numSubstitutions ( regExp->getNErr () )
{
}
void PatternHit::printHTMLHeader ( ostream& os, bool reportSubstitutions ) const
{
	tableRowStart ( os );
		if ( reportSubstitutions ) tableHeader ( os, "Number of Substitutions" );
		tableHeader ( os, "Matching Sequence" );
		tableHeader ( os, "Start AA" );
		tableHeader ( os, "Peptide M+H" );
		ProteinHit::printHTMLHeader ( os );
	tableRowEnd ( os );
}
void PatternHit::printHTML ( ostream& os, bool reportSubstitutions ) const
{
	tableRowStart ( os );
		if ( reportSubstitutions ) tableCell ( os, numSubstitutions, true );

		tableDataStart ( os, "", "", true );
			os << "(" << previousAA << ")";
			if ( numSubstitutions )
				printSequenceWithErrors ( os, sequence, (RegularExpressionWithErrors*)regExp );
			else
				os << sequence;

			os << "(" << nextAA << ")";
			os << endl;
		tableDataEnd ( os );
			
		tableDataStart ( os, "", "center" );
			os << startAA;
		tableDataEnd ( os );

		tableDataStart ( os, "", "right" );
			genPrint ( os, matchedPepMW, instInf->getParentPeakPrecision ().getMassDecimalPlaces () );
			os << endl;
		tableDataEnd ( os );
		ProteinHit::printHTMLHit2 ( os );
	tableRowEnd ( os );
}
void PatternHit::printXML ( ostream& os ) const
{
	os << "<hit>" << endl;
	ParameterList::printXML ( os, "num_substitutions", numSubstitutions );
	ParameterList::printXML ( os, "previous_aa", previousAA );
	ParameterList::printXML ( os, "next_aa", nextAA );
	ParameterList::printXML ( os, "start_aa", startAA );
	printXMLHit ( os );
	os << "</hit>" << endl;
}
void PatternHit::printSequenceWithErrors ( ostream& os, const string& sequenceHit, RegularExpressionWithErrors* regExp )
{
	const char* p = regExp->errorLocation ( sequenceHit.c_str (), true );
	for ( StringSizeType i = 0 ; i < sequenceHit.length () ; i++ ) {
		if ( &(sequenceHit[i]) != p )
			os << sequenceHit[i];
		else {
			os << "<font size=\"+1\"><b>" << sequenceHit[i] << "</b></font>";
			p = regExp->errorLocation ( sequenceHit.c_str (), false );
		}
	}
}
PatternSearch::PatternSearch ( const MSPatternParameters& params ) :
	DatabaseSearch ( params ),
	patternParams ( params ),
	maxPeptideHits ( params.getMaxPeptideHits () )
{
	PatternHits* patternHits = new PatternHits ( params.getMaxAASubstitutions () != 0 );
	for ( int i = 0 ; i < fs.size () ; i++ ) {
		doSearch ( patternHits, fs [i], i );
		FrameIterator::resetElapsedTime ( 1 );
	}
	patternHits->sortHits ();

	numHits = patternHits->size ();
	databaseHits = patternHits;
}
void PatternSearch::doSearch ( PatternHits* patternHits, FastaServer* fsPtr, int i )
{
	ProteinHit::addFS ( fsPtr, i );
	FrameIterator fi ( fsPtr, params.getIndicies ( i ), dnaFrameTranslationPairVector [i], params.getTempOverride () );

	int maxErrs = patternParams.getMaxAASubstitutions ();
	string regularExp = patternParams.getRegularExpression ();
	bool noEnzyme = params.getNoEnzyme ();
	if ( !noEnzyme ) regularExp = "^" + regularExp;

	RegularExpression* regexp;
	if ( maxErrs )	
		regexp = new RegularExpressionWithErrors ( regularExp, maxErrs );
	else
		regexp = new RegularExpression ( regularExp );

	if ( !noEnzyme ) init_fasta_enzyme_function ( params.getEnzyme () );

	char* readingFrame;
	while ( ( readingFrame = fi.getNextFrame () ) != NULL ) {
		if ( noEnzyme ) {
			while ( regexp->isPresentMultiOverlap ( readingFrame ) ) {
				patternHits->addHit ( PatternHit ( fsPtr, fi.getEntry (), fi.getFrameTranslation (), fi.getFrame (), regexp ) );
			}
		}
		else {
			const IntVector& cleavage_index = enzyme_fragmenter ( readingFrame );

			for ( IntVectorSizeType m = 0 ; m < cleavage_index.size () ; m++ ) {
				int offset = ( m == 0 ) ? 0 : cleavage_index [m-1] + 1;
				int n = m;
				if ( regexp->isPresent ( readingFrame + offset ) ) {
					const char* loc = regexp->getMatchEndPointer ();
					while ( loc > readingFrame + cleavage_index [n] + 1 ) n++;
					if ( loc == readingFrame + cleavage_index [n] + 1 )
						patternHits->addHit ( PatternHit ( fsPtr, fi.getEntry (), fi.getFrameTranslation (), fi.getFrame (), regexp ) );
				}
			}
		}
		if ( patternHits->size () > maxPeptideHits ) {
			ErrorHandler::genError ()->error ( "The maximum number of hits has been exceeded.\n" );
		}
	}
}
void PatternSearch::printParamsBodyHTML ( ostream& os ) const
{
	patternParams.printHTML ( os );
}
