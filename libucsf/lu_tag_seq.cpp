/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_tag_seq.cpp                                                *
*                                                                             *
*  Created    : July 21st 1996                                                *
*                                                                             *
*  Purpose    : Contains most of the search algorithms for MS-Seq.            *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1996-2009) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lgen_error.h>
#include <lu_tag_frag.h>
#include <lu_tag_par.h>
using std::string;

MSSeqSearch::MSSeqSearch ( MSMSDataPoint* file, MSTagParameters& params ) :
	MSMSSearch ( file, params )
{
	numPeaks = peaks.size ();
	if ( !peaks.getSpectrumRetained () ) {
		ErrorHandler::genError ()->error ( "The spectrum has been rejected because the precursor m/z is too low.\nNote that the precursor ion must be specified before the fragment ions.\n" );
	}
	if ( numPeaks < 2 ) {
		ErrorHandler::genError ()->error ( "This program requires at least two fragment peaks after peak filtering.\n" );
	}
	if ( peaks.getPrecursorLessThanFragments ( parentMass ) ) {
		ErrorHandler::genError ()->warning ( "The precursor ion mass is less than the lowest fragment ion mass.\nNote that the precursor ion must be specified before the fragment ions.\n" );
	}
	numMasses = numPeaks - 1;

	string ionSeries = params.getIonSeries ();
	if ( ionSeries == "a" ) {
		tagMass1 = peaks [0]->getMass () + c_o;
		tagMass3 = parentMass - peaks [numPeaks - 1]->getMass () - h1 - c_o;
	}
	if ( ionSeries == "b" ) {
		tagMass1 = peaks [0]->getMass ();
		tagMass3 = parentMass - peaks [numPeaks - 1]->getMass () - h1;
	}
	if ( ionSeries == "c" ) {
		tagMass1 = peaks [0]->getMass () - n_h3;
		tagMass3 = parentMass - peaks [numPeaks - 1]->getMass () - h1 + n_h3;
	}
	if ( ionSeries == "a" || ionSeries == "b" || ionSeries == "c" ) {
		tagMass1Tolerance = peaks [0]->getTolerance ();
		for ( int i = 0 ; i < numMasses ; i++ ) {
			seqTag.push_back ( peaks [i+1]->getMass () - peaks [i]->getMass () );
			seqTagTolerance .push_back ( peaks [i+1]->getTolerance () + peaks [i]->getTolerance () );
		}
		tagMass3Tolerance = parentTolerance + peaks [numPeaks - 1]->getTolerance ();
	}
	if ( ionSeries == "y" ) {
		seqTag.resize ( numPeaks );
		seqTagTolerance.resize ( numPeaks );
		tagMass1 = parentMass - peaks [numPeaks - 1]->getMass () + h1; 
		tagMass1Tolerance = parentTolerance + peaks [numPeaks - 1]->getTolerance ();
		for ( int i = 0 ; i < numMasses ; i++ ) {
			seqTag [numMasses - 1 - i] = peaks [i+1]->getMass () - peaks [i]->getMass ();
			seqTagTolerance [numMasses - 1 - i] = peaks [i+1]->getTolerance () + peaks [i]->getTolerance ();
		}
		tagMass3 = peaks [0]->getMass () - h2;
		tagMass3Tolerance = peaks [0]->getTolerance ();
	}
}
MSSeqSearch::~MSSeqSearch ()
{
}
bool MSSeqSearch::doMatch ( const string& peptide, bool nTermPeptide, bool cTermPeptide, double mol_wt, TagMatchVector& tagMatch, const ScoreType& minScore )
{
	if ( compositionSearch ) {
		if ( compositionSearch->doCompositionSearch ( peptide ) == false ) return ( 0 );
	}
	return doSeqMatch ( peptide, mol_wt, tagMatch );
}
bool MSSeqSearch::doSeqMatch ( const string& peptide, double mol_wt, TagMatchVector& tagMatch )
{
	int i, j;

	double ion = n_terminus_wt;
	int len = peptide.length ();
	for ( i = 0 ; i < len ; i++ ) {
		ion += amino_acid_wt [peptide [i]];
		if ( genAbsDiff ( ion, tagMass1 ) < tagMass1Tolerance ) {	/* Tag 1 matches */
			double current_aa_sum = 0.0;
			for ( i++, j = 0 ; i < len && j < numMasses ; i++ ) {
				current_aa_sum += amino_acid_wt [peptide [i]];
				if ( current_aa_sum < seqTag [j] - seqTagTolerance [j] ) continue;
				if ( genAbsDiff ( current_aa_sum, seqTag[j] ) >= seqTagTolerance [j] ) {	/* Sequence filter */
					return false;
				}
				current_aa_sum = 0.0;
				j++;
			}
			break;	/* Sequence matches */
		}
		if ( tagMass1 < ion ) return false;
	}
	tagMatch.clear ();
	tagMatch.push_back ( TagMatch ( &parentPeak, 0, 0, PeptideSequence ( peptide ) ) );
	return true;
}
MSSeqMass1Mass3Search::~MSSeqMass1Mass3Search ()
{
}
bool MSSeqMass1Mass3Search::doSeqMatch ( const string& peptide, double mol_wt, TagMatchVector& tagMatch )
{
	double ion = n_terminus_wt;
	int len = peptide.length ();
	for ( int i = 0 ; i < len ; i++ ) {
		ion += amino_acid_wt [peptide[i]];
		if ( genAbsDiff ( ion, tagMass1 ) < tagMass1Tolerance ) {				/* Tag 1 matches */
			ion = c_terminus_wt;
			for ( int j = len ; j-- ; ) {
				if ( j <= i ) return false;
				ion += amino_acid_wt [peptide [j]];
				if ( genAbsDiff ( ion, tagMass3 ) < tagMass3Tolerance ) {		/* Tag 3 matches */
					tagMatch.clear ();
					tagMatch.push_back ( TagMatch ( &parentPeak, 0, 0, PeptideSequence ( peptide ) ) );
					return true;
				}
				if ( tagMass3 < ion ) return false;
			}
			return false;
		}
		if ( tagMass1 < ion ) return false;
	}
	return false;
}
MSSeqMass1SeqSearch::~MSSeqMass1SeqSearch ()
{
}
bool MSSeqMass1SeqSearch::doSeqMatch ( const string& peptide, double mol_wt, TagMatchVector& tagMatch )
{
	double ion = n_terminus_wt;
	int i, j;

	int len = peptide.length ();
	for ( i = 0 ; i < len ; i++ ) {
		ion += amino_acid_wt [peptide [i]];
		if ( genAbsDiff ( ion, tagMass1 ) < tagMass1Tolerance ) {
			double current_aa_sum = 0.0;
			for ( j = 0 ; i < len && j < numMasses ; i++ ) {
				current_aa_sum += amino_acid_wt [peptide[i+j]];
				if ( current_aa_sum < seqTag [j] - seqTagTolerance [j] ) continue;
				if ( genAbsDiff ( current_aa_sum, seqTag[j] ) >= seqTagTolerance [j] ) {	/* Sequence filter */
					return false;
				}
				current_aa_sum = 0.0;
				j++;
			}
			if ( j == numMasses ) {
				tagMatch.clear ();
				tagMatch.push_back ( TagMatch ( &parentPeak, 0, 0, PeptideSequence ( peptide ) ) );
				return true;
			}
			else return false;
		}
		if ( tagMass1 < ion ) return false;
	}
	return false;
}
MSSeqMass3SeqSearch::~MSSeqMass3SeqSearch ()
{
}
bool MSSeqMass3SeqSearch::doSeqMatch ( const string& peptide, double mol_wt, TagMatchVector& tagMatch )
{
	double ion = c_terminus_wt;
	int i, j;

	int len = peptide.length ();
	for ( i = len ; i-- ; ) {
		ion += amino_acid_wt [peptide [i]];
		if ( genAbsDiff ( ion, tagMass3 ) < tagMass3Tolerance ) {
			double current_aa_sum = 0.0;
			for ( j = numMasses - 1 ; j >= 0 ; ) {
				if ( i-- == 0 ) break;
				current_aa_sum += amino_acid_wt [peptide [i]];
				if ( current_aa_sum < seqTag [j] - seqTagTolerance [j] ) continue;
				if ( genAbsDiff ( current_aa_sum, seqTag [j] ) >= seqTagTolerance [j] ) {	/* Sequence filter */
					return false;
				}
				current_aa_sum = 0.0;
				j--;
			}
			if ( j == -1 ) {
				tagMatch.clear ();
				tagMatch.push_back ( TagMatch ( &parentPeak, 0, 0, PeptideSequence ( peptide ) ) );
				return true;
			}
			else return false;
		}
		if ( tagMass3 < ion ) return false;
	}
	return false;
}
