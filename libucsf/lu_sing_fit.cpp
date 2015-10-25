/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_sing_fit.cpp                                               *
*                                                                             *
*  Created    : October 29th 2001                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2010) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lu_html.h>
#include <lu_sing_fit.h>
#include <lu_app_gr.h>
#include <lu_param_list.h>
#include <lu_table.h>
using std::ostream;
using std::endl;
using std::showpos;
using std::noshowpos;
using std::stringstream;
using std::string;

SingleFitSearch::SingleFitSearch ( const PeakContainer& peaks ) :
	peaks ( peaks ),
	productLink ( "msdigest" ),
	isotopeLink (),
	peakMatchContext ( peaks.getTolerance (), instInf->getParentPeakPrecision (), peaks.getNonUnitChargeData () ),
	peakUsed ( peaks.size () )
{
}
SingleFitSearch::~SingleFitSearch () {}
void SingleFitSearch::printHTML ( ostream& os, bool hideLinks )
{
	if ( areHits () ) {
		startJavascript ( os );
		productLink.printHTML ( os );
		isotopeLink.printHTML ( os );
		endJavascript ( os );

		os << "<table>" << endl;
		printHTMLBody ( os, hideLinks );
		os << "</table>" << endl;
		os << "<br />" << endl;

		StringVector categories;
		categories.push_back ( "Unmatched" );
		categories.push_back ( "Matched" );
		BoolDeque showCategory ( 2, false );
		showCategory [1] = true;
		LabelledCatagorizedGraphData graphData ( categories, showCategory );
		PeakPrecision pp = peakMatchContext.getPeakPrecision ();
		for ( PeakContainerSizeType i = 0 ; i < peaks.size () ; i++ ) {
			int color = peakUsed [i] ? 0 : 1;
			string category = peakUsed [i] ? "1|" : "0|";
			int charge = peaks [i]->getCharge ();
			stringstream label;
			genPrint ( label, peaks [i]->getMOverZ (), pp.getMassDecimalPlaces () );
			if ( charge != 1 ) label << "<sup>" << showpos << charge << noshowpos << "</sup>";
			graphData.add ( peaks [i]->getMOverZ (), peaks [i]->getIntensity (), label.str (), color, category );
		}
		if ( SpectrumGraph::getDrawGraph () ) {
			SpectrumGraph s ( "fit_graph.par.txt" );
			s.drawGraph ( os, graphData, false );
			os << "<br />" << endl;
		}
	}
	else {
		os << "No hits." << endl;
		os << "<br />" << endl;
	}
}
void SingleFitSearch::printXML ( ostream& os )
{
	printXMLBody ( os );
	printStatsXML ( os );
}
void SingleFitSearch::calculateTIC ( const PeakContainer& peaks )
{
	double matchedTic = 0.0;
	double totalTic = 0.0;
	for ( PeakContainerSizeType i = 0 ; i < peaks.size () ; i++ ) {
		if ( peakUsed [i] ) matchedTic += peaks [i]->getIntensity ();
		totalTic += peaks [i]->getIntensity ();
	}
	percentTic = matchedTic * 100.0 / totalTic;
}
void SingleFitSearch::calculateStats ()
{
	if ( errors.size () > 1 ) {
		double adev, svar, skew, curt;

		try {
			moment ( &errors[0]-1, errors.size (), &aveError, &adev, &sdevError, &svar, &skew, &curt );
		}
		catch ( lNrecMomentZeroVariance ) {}	// Not bothered
	}
}
void SingleFitSearch::printCoverageHTML ( ostream& os ) const
{
	for ( BoolDequeSizeType i = 0 ; i < peakUsed.size () ; i++ ) {
		if ( peakUsed [i] ) tableCell ( os, "x", false, true );
		else tableCell ( os, ".", false, true );
	}
}
void SingleFitSearch::printCoverCountHTML ( ostream& os, int i ) const
{
	coverageMap [i].printCoverCountHTML ( os );
}
void SingleFitSearch::printStatsHeaderHTML ( ostream& os ) const
{
	tableHeader ( os, "%<br />Cov" );
	tableHeader ( os, "%<br />TIC" );
	string toleranceUnits = peaks.getTolerance ()->getUnitsString ();
	tableHeaderStart ( os );
		os << "Mean<br />Err<br />" << toleranceUnits << endl;
	tableHeaderEnd ( os );
	tableHeaderStart ( os );
		os << "Data<br />Tol<br />" << toleranceUnits << endl;
	tableHeaderEnd ( os );
}
void SingleFitSearch::printStatsHTML ( ostream& os, int i ) const
{
	tableHeaderStart ( os );
		genPrint ( os, coverageMap [i].getPercentCoverage (), 1 );
		os << endl;
	tableHeaderEnd ( os );
	tableHeaderStart ( os );
		genPrint ( os, percentTic, 1 );
		os << endl;
	tableHeaderEnd ( os );
	if ( errors.size () > 1 ) {
		tableHeaderStart ( os, "", "", true );
			genPrintSigFig ( os, aveError, 3 );
			os << endl;
		tableHeaderEnd ( os );
		tableHeaderStart ( os );
			genPrintSigFig ( os, 2 * sdevError, 3 );
			os << endl;
		tableHeaderEnd ( os );
	}
	else tableEmptyNCells ( os, 2 );
}
void SingleFitSearch::printStatsXML ( ostream& os ) const
{
	os << "<statistics>" << endl;
	ParameterList::printDoubleXMLFixed ( os, "percent_coverage", coverageMap [0].getPercentCoverage (), 1 );
	ParameterList::printDoubleXMLFixed ( os, "percent_tic", percentTic, 1 );
	if ( errors.size () > 1 ) {
		string toleranceUnits = peaks.getTolerance ()->getUnitsString ();
		ParameterList::printDoubleXMLSigFig ( os, "mean_error", aveError, 3 );
		ParameterList::printDoubleXMLSigFig ( os, "data_tolerance", 2 * sdevError, 3 );
	}
	os << "</statistics>" << endl;
}
