/******************************************************************************
*                                                                             *
*  Program    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_bdg_srch.cpp                                               *
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
*  Copyright (2001-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lu_html.h>
#include <lu_bdg_srch.h>
#include <lu_bdg_par.h>
#include <lu_aa_calc.h>
using std::string;
using std::ostream;

MSSingleSearch* getMSBridgeSearch ( MSBridgeParameters& bridgeParams )
{
	if ( bridgeParams.getDataSetInfo ()->getNoDataFlag () ) {
		return new MSSingleSearch ( bridgeParams );
	}
	else {
		return new MSBridgeSearch ( bridgeParams );
	}
}
MSBridgeSearch::MSBridgeSearch ( MSBridgeParameters& bridgeParams ) :
	MSSingleSearch ( bridgeParams ),
	bridgeParams ( bridgeParams ),
	parentPeaks ( getParentPeaks ( bridgeParams ) ),
	linkInfo ( bridgeParams.getLinkInfo () )
{
	PotentialLinkFragment::setAccessionNumbers ( bridgeParams.getAccessionNumbers () );
	AACalculator aaCalc ( bridgeParams.getPeakContainerInfo ().getMonoisotopicFlag (), bridgeParams.getAAInitInfo ().getConstMods () );
	if ( digParams.getSeparateProteinsFlag () ) {
		for ( EnzymeFragmentContainerVectorSizeType i = 0 ; i < enzFrags.size () ; i++ ) {
			linksSearch.push_back ( new LinksSearch ( enzFrags [i], parentPeaks, linkInfo, aaCalc, bridgeParams.getKnownSequences () ) );
		}
	}
	else {
		linksSearch.push_back ( new LinksSearch ( enzFrags, parentPeaks, linkInfo, aaCalc, bridgeParams.getKnownSequences () ) );
	}
}
void MSBridgeSearch::printParamsBodyHTML ( ostream& os ) const
{
	bridgeParams.printHTML ( os );
}
PeakContainer MSBridgeSearch::getParentPeaks ( MSBridgeParameters& bridgeParams )
{
	return PeakContainer ( bridgeParams.getDataSetInfo (), bridgeParams.getMSPeakFilterOptions (), bridgeParams.getPeakContainerInfo () );
}
void MSBridgeSearch::printProteinHTML ( ostream& os, int searchIndex, int proteinIndex ) const
{
	string protein;
	if ( bridgeParams.getSeparateProteinsFlag () )
		protein = se [searchIndex]->getProtein ();
	else
		protein = se [proteinIndex]->getProtein ();
	if ( bridgeParams.getHideProteinSequence () == false ) {
		CharVector coverageMap = linksSearch [searchIndex]->getAACovered ( proteinIndex );
		const IntVector& cleavageIndicies = enzyme_fragmenter ( protein );
		htmlPrintProteinSequence ( os, protein, cleavageIndicies, coverageMap, true );
	}
	if ( bridgeParams.getSeparateProteinsFlag () )
		linksSearch [searchIndex]->printCoverCountHTML ( os, 0 );
	else
		linksSearch [0]->printCoverCountHTML ( os, proteinIndex );
}
void MSBridgeSearch::printResultsHTML ( ostream& os, int i )
{
	linksSearch [i]->printHTML ( os, bridgeParams.getHideHTMLLinks () );
}
void MSBridgeSearch::printResultsXML ( ostream& os, int i )
{
	linksSearch [i]->printXML ( os );
}
