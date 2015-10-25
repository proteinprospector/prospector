/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_spep_srch.cpp                                              *
*                                                                             *
*  Created    : February 12th 2002                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2002-2010) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifndef VIS_C
#include <stdexcept>
#endif
#include <lgen_error.h>
#include <lu_mat_score.h>
#include <lu_spep_srch.h>
#include <lu_param_list.h>
#include <lu_table.h>
using std::string;
using std::ostream;
using std::endl;
using std::sort;
using std::runtime_error;

ModificationHit::ModificationHit ( const EnzymeFragment& enzFrag, const PeakMatch& peakMatch, const PeptideSequence& mutatedSequence ) :
	enzFrag ( enzFrag ),
	peakMatch ( peakMatch ),
	mutatedSequence ( mutatedSequence )
{
}
bool ModificationHit::operator!= ( const ModificationHit& rhs ) const
{
	if ( mutatedSequence.getSequence () != rhs.mutatedSequence.getSequence () || enzFrag.getFragment () != rhs.enzFrag.getFragment () )
		return true;
	else
		return false;
}
bool ModificationHit::isHomologyMatch ( const ModificationHit& rhs ) const
{
	const string& frag = enzFrag.getFragment ();
	const string& frag2 = rhs.enzFrag.getFragment ();

	if ( frag.length () == frag2.length () ) {
		if ( (double)matrix_score ( frag.c_str (), frag2.c_str () ) / (double)frag.length () > 4.0 ) {
			return true;
		}
	}
	return false;
}
void ModificationHit::printHTMLHeader ( ostream& os, const PeakMatchContext& peakMatchContext, bool reportModifications )
{
	tableRowStart ( os );
		peakMatch.printHeaderHTML ( os, peakMatchContext );
		EnzymeFragment::printHeaderHTML ( os );
	tableRowEnd ( os );
}
void ModificationHit::printHTML ( ostream& os, const PeakMatchContext& peakMatchContext, const MSProductLink& productLink, bool reportModifications, bool hideLinks ) const
{
	tableRowStart ( os );
		peakMatch.printHTML ( os, peakMatchContext );
		int maxCharge = peakMatch.getCharge ();
		if ( mutatedSequence.getSequence ().empty () ) {
			enzFrag.printHTML ( os, hideLinks, productLink, maxCharge, 0 );
		}
		else {
			enzFrag.printHTML ( os, true, productLink, maxCharge, 0 );
			productLink.write2 ( os, enzFrag.getFragment (), mutatedSequence, hideLinks, maxCharge );
		}
	tableRowEnd ( os );
}
void ModificationHit::printXML ( ostream& os, const PeakMatchContext& peakMatchContext ) const
{
	os << "<peptide_match>" << endl;
	peakMatch.printXML ( os, peakMatchContext );
	ParameterList::printXML ( os, "mutated_sequence", mutatedSequence.getSequence () );
	enzFrag.printXML ( os );
	os << "</peptide_match>" << endl;
}

SinglePeptideSearch::~SinglePeptideSearch () {}

void SinglePeptideSearch::calculateAACovered ( int index, int protLen )
{
	coverageMap [index].resize ( protLen );
	for ( ModificationHitVectorSizeType i = 0 ; i < hits.size () ; i++ ) {
		const EnzymeFragment& ef = hits [i].getEnzymeFragment ();
		try {
			coverageMap [index].setCoverage ( ef.getStartAA (), ef.getEndAA () );
		}
		catch ( runtime_error e ) {
			ErrorHandler::genError ()->message ( e );
		}
	}
}
bool SinglePeptideSearch::areHits () const
{
	return !hits.empty ();
}
void SinglePeptideSearch::calculateNumUnique ()
{
	numUnique = 1;
	for ( ModificationHitVectorSizeType i = 1 ; i < hits.size () ; i++ ) {
		if ( hits [i] != hits [i-1] ) {
			numUnique++;
		}
	}
}
void SinglePeptideSearch::printHTMLBody ( ostream& os, bool hideLinks )
{
	for ( ModificationHitVectorSizeType i = 0 ; i < hits.size () ; i++ ) {
		if ( i == 0 ) hits [i].printHTMLHeader ( os, peakMatchContext, getReportModifications () );
		hits [i].printHTML ( os, peakMatchContext, productLink, getReportModifications (), hideLinks );
	}
}
void SinglePeptideSearch::printXMLBody ( ostream& os )
{
	for ( ModificationHitVectorSizeType i = 0 ; i < hits.size () ; i++ ) {
		hits [i].printXML ( os, peakMatchContext );
	}
}
bool SinglePeptideSearch::getReportModifications () const
{
	if ( reportModificationsSet == false ) {
		reportModifications = false;
		for ( ModificationHitVectorSizeType i = 0 ; i < hits.size () ; i++ ) {
			if ( hits [i].getMutated () ) {
				reportModifications = true;
				break;
			}
		}
	}
	return reportModifications;
}
ModificationSearch::ModificationSearch ( const EnzymeFragmentContainer& enzFrags, const ModificationTable& modificationTable, const PeakContainer& peaks ) :
	SinglePeptideSearch ( peaks )
{
	for ( PeakContainerSizeType i = 0 ; i < peaks.size () ; i++ ) {
		Peak* peak = peaks [i];
		bool matched = false;

		if ( peak->getAverageMass () ) {
			modified_mass_convert ( false );
			peak->setMassAverage ();
		}
		for ( EnzymeFragmentIterator efi ( enzFrags ) ; efi.more () ; efi.advance () ) {
			const EnzymeFragment& ef = efi.getEnzFrag ();
			double unmodifiedMass = ef.getMass ();

           	if ( peak->isMatch ( unmodifiedMass ) ) {
				matched = true;
				hits.push_back ( ModificationHit ( ef, PeakMatch ( peak, unmodifiedMass ), PeptideSequence ( "" ) ) );
			}
			else if ( peak->isInRange ( unmodifiedMass, modificationTable.getMaxParentError () ) ) {
				double diff = peak->getMass () - unmodifiedMass;
				const bool firstFragment = ef.getFirstFragment ();
				const bool lastFragment = ef.getLastFragment ();

				PeptideSequenceVector mutatedSequences;
				modificationTable.getMutatedSequences ( diff - peak->getTolerance (), diff + peak->getTolerance (), ef.getFragment (), firstFragment, lastFragment, mutatedSequences );
				for ( PeptideSequenceVectorSizeType j = 0 ; j < mutatedSequences.size () ; j++ ) {
					matched = true;
					const PeptideSequence& ps = mutatedSequences [j];
					PeakMatch pm ( peak, ps.getMW () );
					hits.push_back ( ModificationHit ( ef, pm, ps ) );
					errors.push_back ( pm.getError ( peakMatchContext ) );
				}
//if ( get_enzyme_terminal_specificity () == 'N' && ef.previousAA != '-' && (*e)->getTerminalSpecificity () != 'C' ) start_f++;
//if ( get_enzyme_terminal_specificity () == 'C' && ef.nextAA != '-' && (*e)->getTerminalSpecificity () != 'N' ) end_f--;
			}
		}
		peakUsed [i] = matched;
		if ( peak->getAverageMass () ) modified_mass_convert ( true );
	}
	sort ( hits.begin (), hits.end (), sortModificationHits () );
	coverageMap.resize ( 1 );
	calculateAACovered ( 0, enzFrags.getProtLen () );
	calculateNumUnique ();
	calculateTIC ( peaks );
	calculateStats ();
}
ModificationSearch::~ModificationSearch () {}
int ModificationSearch::calculateNumHomologyMatches ( const ModificationSearch* ls ) const
{
	int numMatches = 0;
	if ( areHits () && ls->areHits () ) {
		if ( hits.size () > 500 && ls->hits.size () > 500 ) {	// Make sure the algorithm doesn't get stuck here if someone chose stupid parameters
			return 0;
		}
		ModificationHit lastHit = hits [0];
		for ( ModificationHitVectorSizeType i = 0 ; i < hits.size () ; i++ ) {
			const ModificationHit& currentHit = hits [i];
			if ( i == 0 || currentHit != lastHit ) {
				ModificationHit lastHomHit = ls->hits [0];
				for ( ModificationHitVectorSizeType j = 0 ; j < ls->hits.size () ; j++ ) {
					const ModificationHit& homHit = ls->hits [j];
					if ( j == 0 || homHit != lastHomHit ) {
						if ( currentHit.isHomologyMatch ( homHit ) ) {
							numMatches++;
							break;
						}
					}
					lastHomHit = homHit;
				}
			}
			lastHit = currentHit;
		}
	}
	return numMatches;
}
NonSpecificSearch::NonSpecificSearch ( const char* frame, int frameLen, const PeakContainer& peaks, int maxNonSpecificHits ) :
	SinglePeptideSearch ( peaks )
{
	double* startProtein = get_protein_double_mass_array ( frame );
	double* endProtein = startProtein + frameLen;
	for ( PeakContainerSizeType i = 0 ; i < peaks.size () ; i++ ) {
		Peak* peak = peaks [i];
		bool matched = false;

		if ( peak->getAverageMass () ) {
			modified_mass_convert ( false );
			peak->setMassAverage ();
		}
		double startMass= peak->getMassMinusTol ();
		double endMass	= peak->getMassPlusTol ();
		double* start = startProtein;
		double* end = startProtein;
		double sum = terminal_wt;
		if ( sum < endMass ) {	// This is only a problem if someone has entered a stupid mass
			for ( ; end <= endProtein || sum >= startMass ; ) {
				if ( sum < startMass ) {
					sum += *end++;
				}
				else if ( sum > endMass ) {
					sum -= *start++;
					while ( sum > startMass && end - start > 1 ) {
						sum -= *--end;
					}
				}
				else {
					if ( end - start > 0 ) {
						matched = true;
						EnzymeFragment ef ( frame, start-startProtein, end-startProtein-1, 0 );
						hits.push_back ( ModificationHit ( ef, PeakMatch ( peak, sum ), PeptideSequence ( "" ) ) );
						if ( hits.size () > maxNonSpecificHits ) {
							throw runtime_error ( "The maximum number of hits has been exceeded." );
						}
						if ( end == endProtein ) {
							sum -= *start++;
							while ( sum > startMass && end - start > 1 ) {
								sum -= *--end;
							}
						}
						else {
							sum += *end++;
						}
					}
					else {
						if ( end != endProtein ) sum += *end++;
						else break;
					}
				}
			}
		}
		peakUsed [i] = matched;
		if ( peak->getAverageMass () ) modified_mass_convert ( true );
	}
	sort ( hits.begin (), hits.end (), sortNSpecHits () );
	coverageMap.resize ( 1 );
	calculateAACovered ( 0, frameLen );
	calculateNumUnique ();
	calculateTIC ( peaks );
	calculateStats ();
}
NonSpecificSearch::~NonSpecificSearch () {}
