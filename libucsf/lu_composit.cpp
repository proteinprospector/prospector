/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_composit.cpp                                               *
*                                                                             *
*  Created    : January 7th 1997                                              *
*                                                                             *
*  Purpose    : Composition searches based on bit maps.                       *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1997-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lgen_error.h>
#include <lu_immonium.h>
#include <lu_mass.h>
#include <lu_composit.h>
#include <lu_param_list.h>
#include <lu_charge.h>
using std::ostream;
using std::string;

CompositionSearchParameters::CompositionSearchParameters ( const ParameterList* params ) :
	compIons			( params->getStringVectorValue ( "comp_ion" ) ),
	compositionExclude	( params->getStringValue ( "composition_exclude" , "" ) ),
	excludeFlag			( params->getBoolValue ( "exclude_flag", false ) ),
	searchFlag			( params->getBoolValue ( "composition_search", false ) )
{
}
void CompositionSearchParameters::printHTML ( ostream& os ) const
{
	ParameterList::printHTMLContainer ( os, "Composition Include", compIons );
	ParameterList::printHTML ( os, "Exclude Flag", excludeFlag );
	ParameterList::printHTML ( os, "Composition Exclude", compositionExclude );
	ParameterList::printHTML ( os, "Composition Search Flag", searchFlag );
}
CompositionSearch::CompositionSearch ( const CompositionSearchParameters& compSearchParams, PeakContainer& peaks )
{
	if ( compSearchParams.getSearchFlag () ) {
		if ( compSearchParams.getExcludeFlag () ) {
			excludedComposition = ImmoniumInfo::getExcludeListFromPeaks ( peaks, compSearchParams.getCompositionExclude () );
		}
		else {
			excludedComposition = compSearchParams.getCompositionExclude ();
		}
		bool removeDuplicates = !peaks.empty ();
		StringVector includeList = getIncludeList ( peaks, compSearchParams.getCompIons (), removeDuplicates );

		compositionSearchFlag = ( includeList.size () != 0 || excludedComposition.length () != 0 );

		knownCompositionSequence = calculateKnownCompositionSequence ( includeList );
		includeMask = calculateIncludeMask ( includeList );
		eitherOrMasks = calculateEitherOrMasks ( includeList, includeMask );

		specifiedComposition = calculateCompositionSequence ( includeList );

		excludeMask = string_to_mask ( excludedComposition );
	}
	else compositionSearchFlag = false;
}
bool CompositionSearch::doCompositionSearch ( const string& peptide ) const
{
	if ( compositionSearchFlag == false ) return true;

	unsigned int compositionMask = string_to_mask ( peptide );

	if ( ( includeMask & compositionMask ) != includeMask ) return false;

	if ( excludeMask & compositionMask ) return false;

	for ( UIntVectorSizeType i = 0 ; i < eitherOrMasks.size () ; i++ ) {
		if ( ( eitherOrMasks [i] & compositionMask ) == 0 ) return false;
	}
	return true;
}
double CompositionSearch::subtractKnownComposition ( double mass ) const
{
	for ( StringSizeType i = 0 ; i < knownCompositionSequence.length () ; i++ ) {
		mass -= amino_acid_wt [knownCompositionSequence[i]];
	}
	return ( mass );
}
void CompositionSearch::printHTML ( ostream& os ) const
{
	if ( specifiedComposition != "" ) {
		ParameterList::printHTML ( os, "Composition Ions Present", specifiedComposition );
	}
	if ( excludedComposition != "" ) {
		ParameterList::printHTML ( os, "Composition Ions Excluded", excludedComposition );
	}
}
void CompositionSearch::printXML ( ostream& os ) const
{
	if ( specifiedComposition != "" ) {
		ParameterList::printXML ( os, "composition_ions_present", specifiedComposition );
	}
	if ( excludedComposition != "" ) {
		ParameterList::printXML ( os, "composition_ions_excluded", excludedComposition );
	}
}
StringVector CompositionSearch::getIncludeList ( PeakContainer& peaks, const StringVector& compAA, bool removeDuplicates ) const
{
	StringVector includeList;
	for ( StringVectorSizeType f = 0 ; f < compAA.size () ; f++ ) includeList.push_back ( compAA [f] );

	const StringVector& includeCompositions = ImmoniumInfo::getMajorIonCompositions ();
	double immoniumTolerance = ImmoniumInfo::getTolerance ();
	const DoubleVector& includeMasses = ImmoniumInfo::getMajorIonMasses ();
	int includeNumEntries = includeMasses.size ();
	double minimumFragmentMass = ImmoniumInfo::getMinFragmentMass ();
	double maxImmoniumMass = includeMasses [includeNumEntries-1] + immoniumTolerance;

	PeakContainerSizeType i;
	int j;
	for ( i = 0, j = 0 ; i < peaks.size () ; i++ ) {
		if ( peaks [i]->getMass () > maxImmoniumMass && peaks [i]->getMass () > minimumFragmentMass ) {
			peaks [j++] = peaks [i];
		}
		else {
			if ( peaks [i]->getCharge () != 1 ) {
				ErrorHandler::genError ()->error ( "Ions in the immonium region of the spectrum must be singly charged.\n" );
			}
			string inc;
			for ( int k = 0 ; k < includeNumEntries ; k++ ) {

				if ( genAbsDiff ( includeMasses [k], peaks [i]->getMass () ) <= immoniumTolerance ) {
					inc += includeCompositions [k];
				}
			}
			if ( inc.length () ) {
				if ( inc.length () == 1 ) includeList.push_back ( inc );
				else includeList.push_back ( '[' + inc + ']' );
			}
			else {
				if ( peaks [i]->getMass () > minimumFragmentMass ) {
					peaks [j++] = peaks [i];
				}
			}
		}
	}
	peaks.truncate ( j );

	if ( removeDuplicates ) {
		int z = 0;
		for ( StringVectorSizeType x = 0 ; x < includeList.size () ; x++ ) {			/* Cut out repeat composition info */
			bool flag = false;
			for ( StringVectorSizeType y = x + 1 ; y < includeList.size () ; y++ ) {
				if ( includeList [x] == includeList [y] ) flag = true;
			}
			if ( flag == false ) includeList [z++] = includeList [x];
		}
		includeList.resize ( z );
	}

	int d = 0;
	for ( StringVectorSizeType a = 0 ; a < includeList.size () ; a++ ) {	/* If all the aa's of an either/or combination are present as single aa's delete the either/or combination */
		int len = includeList [a].length ();
		if ( len > 1 ) {
			int count = 0;
			for ( int b = 1 ; b < len - 1 ; b++ ) {
				for ( StringVectorSizeType c = 0 ; c < includeList.size () ; c++ ) {
					if ( includeList [c].length () == 1 ) {
						if ( includeList [c][0] == includeList [a][b] ) count++;
					}
				}
			}
			if ( count < len - 2 ) includeList [d++] = includeList [a];
		}
		else includeList [d++] = includeList [a];
	}
	includeList.resize ( d );
	return includeList;
}
string CompositionSearch::calculateKnownCompositionSequence ( const StringVector& includeList ) const
{
	string sequence;
	for ( StringVectorSizeType i = 0 ; i < includeList.size () ; i++ ) {
		if ( includeList [i].length () == 1 ) {
			sequence += includeList [i];
		}
	}
	return sequence;
}
string CompositionSearch::calculateCompositionSequence ( const StringVector& compList ) const
{
	string sequence;
	for ( StringVectorSizeType i = 0 ; i < compList.size () ; i++ ) {
		sequence += compList [i];
	}
	return sequence;
}
unsigned int CompositionSearch::calculateIncludeMask ( const StringVector& includeList )
{
	unsigned int mask = 0;
	for ( StringVectorSizeType i = 0 ; i < includeList.size () ; i++ ) {
		if ( includeList [i].length () == 1 ) {
			mask |= aa_composition_mask [includeList [i][0]];
		}
	}
	return mask;
}
UIntVector CompositionSearch::calculateEitherOrMasks ( const StringVector& includeList, unsigned int includeMask ) const
{
	UIntVector maskList;
	for ( StringVectorSizeType i = 0 ; i < includeList.size () ; i++ ) {
		if ( includeList [i].length () > 1 ) {
			unsigned int mask = 0;
			for ( StringSizeType j = 0 ; j < includeList [i].length () ; j++ ) {
				mask |= aa_composition_mask [includeList [i][j]];
			}
			if ( ( mask & includeMask ) == 0 ) maskList.push_back ( mask );
		}
	}
	return maskList;
}
