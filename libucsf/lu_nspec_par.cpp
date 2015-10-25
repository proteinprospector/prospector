/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_nspec_par.cpp                                              *
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
*  Copyright (2002-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lu_nspec_par.h>
#include <lu_check_db.h>
#include <lu_getfil.h>
#include <lu_cgi_val.h>
#include <lu_html_form.h>
#include <lu_html.h>
#include <lu_param_list.h>
using std::string;
using std::ostream;
using std::endl;

MSNonSpecificParameters::MSNonSpecificParameters ( const ParameterList* params ) :

	MSDigestParameters ( params ),

	msPeakFilterOptions	( params ),
	dataSetInfo			( params ),
	peakContainerInfo	( params ),

	maxHits	( params->getIntValue	( "max_hits", 50000 ) )
{
}
void MSNonSpecificParameters::printHTML ( ostream& os ) const
{
	MSDigestParameters::printHTML ( os );
}
int MSNonSpecificLink::num = 0;
StringVector MSNonSpecificLink::database;
MSNonSpecificLink::MSNonSpecificLink ( const string& db )
{
	num++;
	ind = num;
	database.push_back ( db );
}
void MSNonSpecificLink::write ( ostream& os, int indexNumber, int dnaReadingFrame, int openReadingFrame, const string& str, int index )
{
	static int idx = 1;
	string variableName = "msnonspecific" + gen_itoa ( idx ) + "form";
	idx++;
	ProgramLink::putForm ( os, variableName, "mssearch.cgi" );
	startJavascript ( os );
	os << "printMSNonSpecificLinkItems" << index+1 << " ();" << endl;
	os << "printUnmatchedPeaks" << str << " ();" << endl;
	endJavascript ( os );

	bool dnaFlag = is_dna_database ( database [index] );
	SingleEntryParameters::putIndexHiddenFormEntry ( os, dnaFlag, indexNumber, dnaReadingFrame, openReadingFrame );

	os << "</form>" << endl;
	ProgramLink::putFormLink ( os, variableName, "Do a non-specific cleavage search." );
}
void MSNonSpecificLink::putHidden ( ostream& os ) const
{
	const ParameterList* params = ProgramLink::getParams ();
	string db = database [ind-1];
	os << "function printMSNonSpecificLinkItems" << ind << " () {" << endl;
		printHTMLFORMJavascriptHidden ( os, "search_name", "msnonspecific" );
		printHTMLFORMJavascriptHidden ( os, "output_type", "HTML" );
		printHTMLFORMJavascriptHidden ( os, "report_title", "MS-NonSpecific" );
		printHTMLFORMJavascriptHidden ( os, "version", Version::instance ().getVersion () );
		if ( !db.empty () )	printHTMLFORMJavascriptHidden2 ( os, "database", db );
		else				params->copyToHiddenFormJavascriptEntry ( os, "database", ind-1 );
		params->copyToHiddenFormJavascriptEntry ( os, "n_term_aa_limit" );
		AAInitInfo::copyToHiddenFormJavascriptEntry ( os, params );
		EnzymeParameters::copyToHiddenFormJavascriptEntry ( os, params, true );
		MassInfo::copyToHiddenFormJavascriptEntry ( os, params );
	os << "}" << endl;
}
