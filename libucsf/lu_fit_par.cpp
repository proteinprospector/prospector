/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_fit_par.cpp                                                *
*                                                                             *
*  Created    : July 12th 2001                                                *
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
#include <lg_string.h>
#include <lu_fit_par.h>
#include <lu_cgi_val.h>
#include <lu_getfil.h>
#include <lu_html.h>
#include <lu_html_form.h>
#include <lu_param_list.h>
using std::ostream;
using std::count;
using std::string;
using std::endl;

MowseInfo::MowseInfo ( const ParameterList* params ) :
	mowseOn		( params->getBoolValue	( "mowse_on", false ) ),
	mowsePFactor( params->getDoubleValue( "mowse_pfactor", 0.4 ) )
{
}
void MowseInfo::printHTML ( ostream& os ) const
{
	ParameterList::printHTML ( os, "MOWSE On", mowseOn );
	ParameterList::printHTML ( os, "MOWSE P Factor", mowsePFactor );
}
void MowseInfo::copyToCGI ( ostream& os, const ParameterList* params )
{
	if ( params->copyToCGI ( os, "mowse_on" ) ) {
		params->copyToCGI ( os, "mowse_pfactor" );
	}
}
void MowseInfo::copyToHiddenFormEntry ( ostream& os, const ParameterList* params )
{
	if ( params->copyToHiddenFormEntry ( os, "mowse_on" ) ) {
		params->copyToHiddenFormEntry ( os, "mowse_pfactor" );
	}
}
void MowseInfo::copyToHiddenFormJavascriptEntry ( ostream& os, const ParameterList* params )
{
	if ( params->copyToHiddenFormJavascriptEntry ( os, "mowse_on" ) ) {
		params->copyToHiddenFormJavascriptEntry ( os, "mowse_pfactor" );
	}
}
int MSFitParameters::maxReportedHitsLimit = InfoParams::instance ().getIntValue ( "msfit_max_reported_hits_limit", 500 );
MSFitParameters::MSFitParameters ( const ParameterList* params ) :
	MSSearchParameters	( params, "ms_" ),

	msPeakFilterOptions	( params ),

	dataSetInfo			( params ),
	peakContainerInfo	( params ),
	specID				( params ),

	minMatches			( params->getIntValue	( "min_matches", 5 ) ),
	mowseInfo			( params ),
	sortType			( params->getStringValue( "sort_type", "Score Sort" ) ),
	multipleModification( params ),
	multipleModification2( params, true ),
	minParentIonMatches	( params->getIntValue	( "min_parent_ion_matches", 1 ) ),
	reportHomologousProteins	( params->getStringValue ( "ms_report_homologous_proteins", "Interesting" ) ),

	modificationParameters		( params, "ms_" )
{
	if ( getMaxReportedHits () > maxReportedHitsLimit ) {
		ErrorHandler::genError ()->error ( "The maximum number of reported hits cannot be greater than " + gen_itoa ( maxReportedHitsLimit ) + ".\n" );
	}
	init_fasta_enzyme_function ( enzymeParameters.getEnzyme () );
	multipleModification2.setUserMods ( aaInitInfo );
}
void MSFitParameters::printHTML ( ostream& os ) const
{
	MSSearchParameters::printHTML ( os );
	ParameterList::printHTML ( os, "Minimum Matches", minMatches );
	ParameterList::printHTML ( os, "Sort Type", sortType );
	multipleModification.printHTML ( os );
	ParameterList::printHTML ( os, "Min Parent Ion Matches", minParentIonMatches );
	mowseInfo.printHTML ( os );
	ParameterList::printHTML ( os, "Report Homologous Proteins", reportHomologousProteins );
}
MSFitLink::MSFitLink ()
{
}
void MSFitLink::write ( ostream& os, int minMatches, const PeakContainer& peaks, const BoolDeque& peakUsed, const string& str )
{
	static int idx = 1;
	string variableName = "msfit" + gen_itoa ( idx ) + "form";
	idx++;
	int numUnmatchedMasses = count ( peakUsed.begin (), peakUsed.end (), false );

	ProgramLink::putForm ( os, variableName, "mssearch.cgi" );
	startJavascript ( os );
	os << "printMSFitLinkItems ();" << endl;
	os << "printUnmatchedPeaks" << str << " ();" << endl;
	endJavascript ( os );
	printHTMLFORMHidden ( os, "min_matches", 1 + ( minMatches * numUnmatchedMasses / peaks.size () ) );

	os << "</form>" << endl;
	ProgramLink::putFormLink ( os, variableName, "Search for another component." );
}
void MSFitLink::putHidden ( ostream& os ) const
{
	const ParameterList* params = ProgramLink::getParams ();
	os << "function printMSFitLinkItems () {" << endl;
		printHTMLFORMJavascriptHidden ( os, "search_name", "msfit" );
		printHTMLFORMJavascriptHidden ( os, "output_type", "HTML" );
		printHTMLFORMJavascriptHidden ( os, "report_title", "MS-Fit" );
		printHTMLFORMJavascriptHidden ( os, "version", Version::instance ().getVersion () );
		MSSearchParameters searchParameters ( params, "ms_" );
		searchParameters.putHiddenFormJavascriptEntry ( os );

		MultipleModification::copyToHiddenFormJavascriptEntry ( os, params );

		ModificationParameters::copyToHiddenFormJavascriptEntry ( os, params, "ms_" );
		params->copyToHiddenFormJavascriptEntry ( os, "min_parent_ion_matches" );

		MassInfo::copyToHiddenFormJavascriptEntry ( os, params );

		MowseInfo::copyToHiddenFormJavascriptEntry ( os, params );

		params->copyToHiddenFormJavascriptEntry ( os, "detailed_report" );
		params->copyToHiddenFormJavascriptEntry ( os, "display_graph" );
	os << "}" << endl;
}
