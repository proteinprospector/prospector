/******************************************************************************
*                                                                             *
*  Program    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_dig_srch.cpp                                               *
*                                                                             *
*  Created    : September 27th 2001                                           *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lu_html.h>
#include <lu_dig_srch.h>
#include <lu_dig_par.h>
using std::ostream;
using std::string;

MSSingleSearch::MSSingleSearch ( const MSDigestParameters& digParams ) :
	MSProgram ( digParams ),
	digParams ( digParams ),
	se ( getSingleEntry ( digParams.getSingleEntryParameters () ) )
{
	for ( SingleEntryPtrVectorSizeType i = 0 ; i < se.size () ; i++ ) {
		enzFrags.push_back ( EnzymeFragmentContainer ( se [i]->getProtein (), digParams.getEnzymeParameters () ) );
	}
}
void MSSingleSearch::printProteinHTML ( ostream& os, int searchIndex, int proteinIndex ) const
{
	string protein;
	if ( digParams.getSeparateProteinsFlag () )
		protein = se [searchIndex]->getProtein ();
	else
		protein = se [proteinIndex]->getProtein ();
	if ( digParams.getHideProteinSequence () == false ) {
		const IntVector& cleavageIndicies = enzyme_fragmenter ( protein );
		htmlPrintProteinSequence ( os, protein, cleavageIndicies, digParams.getCoverageMap (), true );
		printProteinCoverage ( os, digParams.getCoverageMap () );
	}
}
void MSSingleSearch::printProteinCoverage ( ostream& os, const CharVector& aaCovered ) const
{
}
void MSSingleSearch::printBodyHTML ( ostream& os )
{
	if ( digParams.getSeparateProteinsFlag () ) {
		for ( SingleEntryPtrVectorSizeType i = 0 ; i < se.size () ; i++ ) {
			se [i]->printHTML ( os );
			printProteinHTML ( os, i, 0 );
			printResultsHTML ( os, i );
		}
	}
	else {
		for ( SingleEntryPtrVectorSizeType i = 0 ; i < se.size () ; i++ ) {
			se [i]->printHTML ( os );
			printProteinHTML ( os, 0, i );
		}
		printResultsHTML ( os, 0 );
	}
}
void MSSingleSearch::printBodyTabDelimitedText ( ostream& os )
{
	printHeaderDelimited ( os );
	if ( digParams.getSeparateProteinsFlag () ) {
		for ( SingleEntryPtrVectorSizeType i = 0 ; i < se.size () ; i++ ) {
			printResultsDelimited ( os, i );
		}
	}
	else
		printResultsDelimited ( os, 0 );
}
void MSSingleSearch::printBodyXML ( ostream& os )
{
	if ( digParams.getSeparateProteinsFlag () ) {
		for ( SingleEntryPtrVectorSizeType i = 0 ; i < se.size () ; i++ ) {
			se [i]->printXML ( os );
			printResultsXML ( os, i );
		}
	}
	else {
		for ( SingleEntryPtrVectorSizeType i = 0 ; i < se.size () ; i++ ) {
			se [i]->printXML ( os );
		}
		printResultsXML ( os, 0 );
	}
}

MSDigestSearch::MSDigestSearch ( const MSDigestParameters& digParams ) :
	MSSingleSearch ( digParams )
{
	int maxHits = digParams.getMaxHits ();
	if ( digParams.getSeparateProteinsFlag () ) {
		for ( EnzymeFragmentContainerVectorSizeType i = 0 ; i < enzFrags.size () ; i++ ) {
			potentialMSFragments.push_back ( PotentialMSFragmentContainer ( enzFrags [i], digParams.getDigestFragmentParameters (), digParams.getReportMultCharge () ? 0 : 1, maxHits, i ) );
		}
	}
	else {
		potentialMSFragments.push_back ( PotentialMSFragmentContainer ( enzFrags, digParams.getDigestFragmentParameters (), digParams.getReportMultCharge () ? 0 : 1, maxHits ) );
	}
}
void MSDigestSearch::printResultsHTML ( ostream& os, int i )
{
	potentialMSFragments [i].printHTML ( os, digParams.getHideHTMLLinks () );
}
void MSDigestSearch::printResultsDelimited ( ostream& os, int i )
{
	potentialMSFragments [i].printDelimited ( os );
}
void MSDigestSearch::printResultsXML ( ostream& os, int i )
{
	potentialMSFragments [i].printXML ( os );
}
void MSDigestSearch::printProteinCoverage ( ostream& os, const CharVector& aaCovered ) const
{
	int protLen = aaCovered.size ();
	if ( protLen ) {
		int coverCount = 0;
		for ( int m = 0; m < protLen ; m++ ) {
			if ( aaCovered [m] > 0 ) coverCount++;
		}
		os << "The&nbsp;matched&nbsp;peptides&nbsp;cover&nbsp;<b>" << 100 * coverCount / protLen << "%</b>";
		os << "&nbsp;";
		os << "(" << coverCount << "/" << protLen << "&nbsp;AA's)&nbsp;of&nbsp;the&nbsp;protein.";
		os << "<p />";
	}
}
void MSDigestSearch::printParamsBodyHTML ( ostream& os ) const
{
	digParams.printHTML ( os );
}
void MSDigestSearch::printHeaderDelimited ( ostream& os )
{
	PotentialMSFragment::printHeaderDelimited ( os, true );
}
