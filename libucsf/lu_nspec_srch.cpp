/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_nspec_srch.cpp                                             *
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
*  Copyright (2002-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifndef VIS_C
#include <stdexcept>
#endif
#include <lu_html.h>
#include <lu_nspec_par.h>
#include <lu_nspec_srch.h>
using std::string;
using std::ostream;
using std::runtime_error;

MSNonSpecificSearch::MSNonSpecificSearch ( MSNonSpecificParameters& nonSpecificParams ) :
	MSSingleSearch ( nonSpecificParams ),
	nonSpecificParams ( nonSpecificParams ),
	parentPeaks ( nonSpecificParams.getDataSetInfo ()->getDataSet ( 0 ), nonSpecificParams.getMSPeakFilterOptions (), nonSpecificParams.getPeakContainerInfo () )
{
	try {
		for ( SingleEntryPtrVectorSizeType i = 0 ; i < se.size () ; i++ ) {
			nonSpecificSearch.push_back ( new NonSpecificSearch ( se [i]->getProtein (), se [i]->getProteinLength (), parentPeaks, nonSpecificParams.getMaxHits () ) );
		}
	}
	catch ( runtime_error e ) {
		ErrorHandler::genError ()->error ( e );
	}
}
MSNonSpecificSearch::~MSNonSpecificSearch ()
{
	for ( SingleEntryPtrVectorSizeType i = 0 ; i < se.size () ; i++ ) {
		delete nonSpecificSearch [i];
	}
}
void MSNonSpecificSearch::printProteinHTML ( ostream& os, int searchIndex, int proteinIndex ) const
{
	string protein;
	if ( nonSpecificParams.getSeparateProteinsFlag () )
		protein = se [searchIndex]->getProtein ();
	else
		protein = se [proteinIndex]->getProtein ();
	if ( nonSpecificParams.getHideProteinSequence () == false ) {
		CharVector coverageMap (protein.length ());
		//for ( int i = 0 ; i < protein.length () ; i++ ) {
		//	coverageMap [i] = (protein [i] == linkAA) ? 1 : 0;
		//}
		const IntVector& cleavageIndicies = enzyme_fragmenter ( protein );
		htmlPrintProteinSequence ( os, protein, cleavageIndicies, coverageMap, true );
	}
}
void MSNonSpecificSearch::printResultsHTML ( ostream& os, int i )
{
	nonSpecificSearch [i]->printHTML ( os, nonSpecificParams.getHideHTMLLinks () );
}
void MSNonSpecificSearch::printResultsXML ( ostream& os, int i )
{
	nonSpecificSearch [i]->printXML ( os );
}
void MSNonSpecificSearch::printParamsBodyHTML ( ostream& os ) const
{
	nonSpecificParams.printHTML ( os );
}
void MSNonSpecificSearch::printBodyHTML ( ostream& os )
{
	if ( nonSpecificParams.getSeparateProteinsFlag () ) {
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
		for ( SingleEntryPtrVectorSizeType j = 0 ; j < se.size () ; j++ ) {
			printResultsHTML ( os, j );
		}
	}
}
void MSNonSpecificSearch::printBodyXML ( ostream& os )
{
	if ( nonSpecificParams.getSeparateProteinsFlag () ) {
		for ( SingleEntryPtrVectorSizeType i = 0 ; i < se.size () ; i++ ) {
			se [i]->printXML ( os );
			printResultsXML ( os, i );
		}
	}
	else {
		for ( SingleEntryPtrVectorSizeType i = 0 ; i < se.size () ; i++ ) {
			se [i]->printXML ( os );
		}
		for ( SingleEntryPtrVectorSizeType j = 0 ; j < se.size () ; j++ ) {
			printResultsXML ( os, j );
		}
	}
}
