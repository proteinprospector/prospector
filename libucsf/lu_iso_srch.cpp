/******************************************************************************
*                                                                             *
*  Program    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_iso_srch.cpp                                               *
*                                                                             *
*  Created    : October 22nd 2001                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2011) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lu_app_gr.h>
#include <lu_delim.h>
#include <lu_iso_par.h>
#include <lu_iso_srch.h>
#include <lu_param_list.h>
#include <lu_iso.h>
#include <lu_table.h>
using std::vector;
using std::ostream;
using std::endl;

MSIsotopeSearch::MSIsotopeSearch ( const MSIsotopeParameters& params ) :
	MSProgram ( params ),
	isoParams ( params )
{
	graphData = 0;
	if ( isoParams.getDisplayGraph () ) setGraphData ();
}
MSIsotopeSearch::~MSIsotopeSearch ()
{
	delete graphData;
}
void MSIsotopeSearch::setGraphData ()
{
	IsotopeProfile* ip = isoParams.getIsotopeProfile ();
	graphData = new GraphData ( *ip );
	if ( isoParams.getResolution () >= 1000000 ) graphData->setXPrecision ( 5 );
	if ( isoParams.getResolution () >= 10000000 ) graphData->setXPrecision ( 6 );
	delete ip;
}
void MSIsotopeSearch::printBodyHTML ( ostream& os )
{
	vector <IsotopePeakStats*> ips = isoParams.getIsotopePeakStats ();

	if ( ips.size () > 1 && graphData ) {
		SpectrumGraph s1 ( "sp_graph.par.txt" );
		s1.drawGraph ( os, *graphData );
		os << "<p />" << endl;
	}
	for ( int i = 0 ; i < ips.size () ; i++ ) {
		os << "<hr /><p />" << endl;
		ElementalFormula ef = ips [i]->getElementalFormula ();
		ParameterList::printHTMLNC ( os, "Elemental Composition", ef );
		ParameterList::printDoubleHTMLFixed ( os, "Monoisotopic M/Z", ips [i]->getIdealMonoisotopicMZ (), 5 );

		if ( ips [i]->getNoDistribution () ) {
			os << "No isotope distribution, monoisotopic peak is 100%.<p />" << endl;
		}
		else {
			if ( ips [i]->getMonoisotopicPeakTooSmall () ) os << "Monoisotopic peak less than minimum probability.<p />" << endl;
			else {
				os << "Total Abundance: <b>";
				genPrint ( os, ips [i]->getTotalAbundance () * 100.0, 2 );
				os << "%</b><br />" << endl;
				os << "<hr /><p />" << endl;
				if ( ips.size () == 1 ) {
					tableStart ( os );
						tableRowStart ( os );
							os << "<td valign=\"top\" align=\"left\">" << endl << endl << endl;
								printSummaryReportHTML ( os, ips [i] );
							tableDataEnd ( os );
							if ( graphData ) {
								os << "<td valign=\"bottom\" align=\"right\">" << endl;
									SpectrumGraph s1 ( "sp_graph.par.txt" );
									s1.drawGraph ( os, *graphData );
								tableDataEnd ( os );
							}
						tableRowEnd ( os );
					tableEnd ( os );
				}
				else {
					printSummaryReportHTML ( os, ips [i] );
				}
				if ( isoParams.getDetailedReport () ) {
					os << "<hr />" << endl;
					os << "<p />" << endl;
					printDetailedReportHTML ( os, ips [i] );
				}
			}
		}
	}
}
void MSIsotopeSearch::printBodyXML ( ostream& os )
{
	vector <IsotopePeakStats*> vips = isoParams.getIsotopePeakStats ();
	os << "<results>" << endl;
		if ( graphData ) graphData->writeXML ( os );
		for ( int i = 0 ; i < vips.size () ; i++ ) {
			os << "<distribution>" << endl;
				ParameterList::printXML ( os, "elemental_composition", vips [i]->getElementalFormula ().getFormula () );
				ParameterList::printDoubleXMLFixed ( os, "monoisotopic_m_over_z", vips [i]->getIdealMonoisotopicMZ (), 5 );

				if ( vips [i]->getNoDistribution () ) {
					ParameterList::printXML ( os, "error", "No isotope distribution, monoisotopic peak is 100%." );
				}
				else {
					if ( vips [i]->getMonoisotopicPeakTooSmall () ) {
						ParameterList::printXML ( os, "error", "Monoisotopic peak less than minimum probability." );
					}
					else {
						ParameterList::printDoubleXMLFixed ( os, "total_abundance", vips [i]->getTotalAbundance () * 100.0, 2 );
						printSummaryReportXML ( os, vips [i] );
						if ( isoParams.getDetailedReport () ) printDetailedReportXML ( os, vips [i] );
					}
				}
			os << "</distribution>" << endl;
		}
	os << "</results>" << endl;
}
void MSIsotopeSearch::printBodyTabDelimitedText ( ostream& os )
{
	vector <IsotopePeakStats*> vips = isoParams.getIsotopePeakStats ();

	if ( graphData ) graphData->writeTabDelimitedText ( os );
	os << endl << endl << endl;
	for ( int i = 0 ; i < vips.size () ; i++ ) {
		delimitedRowStart ( os );
			delimitedCell ( os, "Elemental Composition:" );
			delimitedCell ( os, vips [i]->getElementalFormula ().getFormula () );
		delimitedRowEnd ( os );

		delimitedRowStart ( os );
			delimitedCell ( os, "Monoisotopic M/Z:" );
			delimitedCell ( os, vips [i]->getIdealMonoisotopicMZ (), 5 );
		delimitedRowEnd ( os );

		if ( vips [i]->getNoDistribution () ) {
			delimitedRowStart ( os );
				delimitedCell ( os, "No isotope distribution, monoisotopic peak is 100%." );
			delimitedRowEnd ( os );
		}
		else {
			if ( vips [i]->getMonoisotopicPeakTooSmall () ) {
				delimitedRowStart ( os );
					delimitedCell ( os, "Monoisotopic peak less than minimum probability." );
				delimitedRowEnd ( os );
			}
			else {
				delimitedRowStart ( os );
					delimitedCell ( os, "Total Abundance:" );
					delimitedCell ( os, vips [i]->getTotalAbundance () * 100.0, 2 );
				delimitedRowEnd ( os );
				os << endl;
				printSummaryReportTabDelimitedText ( os, vips [i] );
				if ( isoParams.getDetailedReport () ) {
					os << endl;
					printDetailedReportTabDelimitedText ( os, vips [i] );
				}
			}
		}
		if ( i != vips.size () - 1 ) os << endl << endl << endl;
	}
}
void MSIsotopeSearch::printSummaryReportHTML ( ostream& os, const IsotopePeakStats* ips )
{
	double maximumProbability = ips->getGroupMaximumProbability ();

	tableStart ( os );
		tableRowStart ( os );
			tableHeader ( os, "Isotope<br />Number", "", "center" );
			tableHeader ( os, "m/z", "", "center" );
			tableHeader ( os, "Percent<br />Total", "", "center" );
			tableHeader ( os, "Percent<br />Maximum", "", "center" );
		tableRowEnd ( os );

		IsotopePeakStatsConstIterator i ( ips );
		for ( int j = 0 ; i.more () ; i.advance (), j++ ) {
			tableRowStart ( os );
				tableDataStart ( os, "", "left" );
					os << j << endl;
				tableDataEnd ( os );
				tableDataStart ( os, "", "right" );
					genPrint ( os, i.averageMass (), 5 );
					os << endl;
				tableDataEnd ( os );
				tableDataStart ( os, "", "right" );
					genPrint ( os, i.totalProbability () * 100.0, 2 );
					os << endl;
				tableDataEnd ( os );
				double y = i.totalProbability () * 100.0 / maximumProbability;
				tableDataStart ( os, "", "right" );
					genPrint ( os, y, 2 );
					os << endl;
				tableDataEnd ( os );
			tableRowEnd ( os );
		}
	tableEnd ( os );
}
void MSIsotopeSearch::printSummaryReportXML ( ostream& os, const IsotopePeakStats* ips )
{
	os << "<summary_report>" << endl;
		double maximumProbability = ips->getGroupMaximumProbability ();
		IsotopePeakStatsConstIterator i ( ips );
		for ( int j = 0 ; i.more () ; i.advance (), j++ ) {
			os << "<isotope";
			os << " ";
			os << "number=\"" << j << "\""; 
			os << " ";
			os << "m_over_z=\""; 
			genPrint ( os, i.averageMass (), 5 );
			os << "\" ";
			os << "percent_total=\""; 
			genPrint ( os, i.totalProbability () * 100.0, 2 );
			os << "\" ";
			os << "percent_maximum=\""; 
			genPrint ( os, i.totalProbability () * 100.0 / maximumProbability, 2 );
			os << "\" />";
			os << endl;
		}
	os << "</summary_report>" << endl;
}
void MSIsotopeSearch::printSummaryReportTabDelimitedText ( ostream& os, const IsotopePeakStats* ips )
{
	double maximumProbability = ips->getGroupMaximumProbability ();
	delimitedRowStart ( os );
		delimitedHeader ( os, "Isotope Number" );
		delimitedHeader ( os, "m/z" );
		delimitedHeader ( os, "Percent Total" );
		delimitedHeader ( os, "Percent Maximum" );
	delimitedRowEnd ( os );
	IsotopePeakStatsConstIterator i ( ips );
	for ( int j = 0 ; i.more () ; i.advance (), j++ ) {
		delimitedRowStart ( os );
			delimitedCell ( os, j );
			delimitedCell ( os, i.averageMass (), 5 );
			delimitedCell ( os, i.totalProbability () * 100.0, 2 );
			delimitedCell ( os, i.totalProbability () * 100.0 / maximumProbability, 2 );
		delimitedRowEnd ( os );
	}
}
void MSIsotopeSearch::printDetailedReportHTML ( ostream& os, const IsotopePeakStats* ips )
{
	os << "<table>" << endl;

	tableRowStart ( os );
		tableHeader ( os, "Isotope<br />Content" );
		tableHeader ( os, "m/z" );
		tableHeader ( os, "Percent<br />Total" );
		tableHeader ( os, "Percent<br />Group Max" );
	tableRowEnd ( os );

	IsotopicDistributionConstIterator id ( ips );
	IsotopePeakStatsConstIterator ipsci ( ips );
	for ( int i = 0 ; id.more () ; id.advance (), i++ ) {
		if ( i == 0 ) {
			tableRowStart ( os );
			tableHeader ( os, "Monoisotopic", "", "left" );
		}
		else {
			if ( i == ipsci.numComponents () ) {
				tableEmptyRow ( os );
				ipsci.advance ();
				i = 0;
			}
			tableRowStart ( os );
			tableHeaderStart ( os, "", "left" );
				os << id.formula () << endl;
			tableHeaderEnd ( os );
		}
		tableDataStart ( os, "", "right" );
			genPrint ( os, id.massOffset (), 5 );
			os << endl;
		tableDataEnd ( os );

		tableDataStart ( os, "", "right" );
			genPrint ( os, id.probability () * 100.0, 2 );
			os << endl;
		tableDataEnd ( os );

		tableDataStart ( os, "", "right" );
			genPrint ( os, id.probability () * 100.0 / ipsci.maximumProbability (), 2 );
			os << endl;
		tableDataEnd ( os );

		tableRowEnd ( os );
	}
	os << "</table>" << endl;
}
void MSIsotopeSearch::printDetailedReportXML ( ostream& os, const IsotopePeakStats* ips )
{
	os << "<detailed_report>" << endl;
		IsotopicDistributionConstIterator id ( ips );
		IsotopePeakStatsConstIterator ipsci ( ips );
		int group = 0;
		for ( int i = 0 ; id.more () ; id.advance (), i++ ) {
			os << "<detailed_isotope";
			os << " ";
			os << "content=\""; 
			if ( i == 0 ) {
				os << "Monoisotopic\"";
			}
			else {
				if ( i == ipsci.numComponents () ) {
					ipsci.advance ();
					group++;
					i = 0;
				}
				os << id.formula ();
			}
			os << "\" ";
			os << "group=\"" << group << "\""; 
			os << " ";
			os << "m_over_z=\""; 
			genPrint ( os, id.massOffset (), 5 );

			os << "\" ";
			os << "percent_total=\""; 
			genPrint ( os, id.probability () * 100.0, 2 );

			os << "\" ";
			os << "percent_group_max=\""; 
			genPrint ( os, id.probability () * 100.0 / ipsci.maximumProbability (), 2 );

			os << "\" />";
			os << endl;
		}
	os << "</detailed_report>" << endl;
}
void MSIsotopeSearch::printDetailedReportTabDelimitedText ( ostream& os, const IsotopePeakStats* ips )
{
	delimitedRowStart ( os );
		delimitedHeader ( os, "Isotope Content" );
		delimitedHeader ( os, "Group" );
		delimitedHeader ( os, "m/z" );
		delimitedHeader ( os, "Percent Total" );
		delimitedHeader ( os, "Percent Group Max" );
	delimitedRowEnd ( os );

	IsotopicDistributionConstIterator id ( ips );
	IsotopePeakStatsConstIterator ipsci ( ips );
	int group = 0;
	for ( int i = 0 ; id.more () ; id.advance (), i++ ) {
		if ( i == 0 ) {
			delimitedRowStart ( os );
			delimitedCell ( os, "Monoisotopic" );
		}
		else {
			if ( i == ipsci.numComponents () ) {
				ipsci.advance ();
				i = 0;
				group++;
			}
			delimitedRowStart ( os );
			delimitedCell ( os, id.formula () );
		}
		delimitedCell ( os, group );
		delimitedCell ( os, id.massOffset (), 5 );
		delimitedCell ( os, id.probability () * 100.0, 2 );
		delimitedCell ( os, id.probability () * 100.0 / ipsci.maximumProbability (), 2 );

		delimitedRowEnd ( os );
	}
}
void MSIsotopeSearch::printParamsBodyHTML ( ostream& os ) const
{
	isoParams.printHTML ( os );
}
