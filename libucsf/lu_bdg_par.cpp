/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_bdg_par.cpp                                                *
*                                                                             *
*  Created    : June 14th 2001                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lgen_math.h>
#include <lu_bdg_par.h>
#include <lu_getfil.h>
#include <lu_check_db.h>
#include <lu_cgi_val.h>
#include <lu_html_form.h>
#include <lu_html.h>
#include <lu_param_list.h>
using std::string;
using std::ostream;
using std::endl;
using std::istringstream;
using std::getline;

MSBridgeParameters::MSBridgeParameters ( const ParameterList* params ) :

	MSDigestParameters ( params ),

	msPeakFilterOptions	( params ),
	knownSequences		( initKnownSequences ( params ) ),
	dataSetInfo			( params ),
	peakContainerInfo	( params ),
	specID				( params ),
	linkInfo ( new LinkInfo ( params, true ) )
{
	linkInfo->setUserMods ( aaInitInfo );
}
void MSBridgeParameters::printHTML ( ostream& os ) const
{
	MSDigestParameters::printHTML ( os );
	linkInfo->printHTML ( os );
}
StringVector MSBridgeParameters::initKnownSequences ( const ParameterList* params )
{
	StringVector sv;
	string data = params->getFileStringValue ( "data", "" );
	istringstream isstr ( data );
	int c;
	string line;
	string newData;
	while ( ( c = isstr.peek () ) != EOF ) {
		if ( !getline ( isstr, line ) ) break;
		line = gen_strtrim ( line );	// Trim the line
		if ( line.length () != 0 ) {
			if ( line [line.length () - 1] == '\r' ) line = line.substr ( 0, line.length () - 1 );
			if ( line [0] == '>' || line [0] == 'E' ) continue;
			if ( isdigit ( line [0] ) ) {
				string::size_type idx1 = line.find ( "<" );
				string::size_type idx2 = line.rfind ( ">" );
				if ( idx1 == string::npos || idx2 == string::npos ) {
					sv.push_back ( "" );
					newData.append ( line + "\n" );
				}
				else {
					sv.push_back ( line.substr ( idx1, idx2-idx1+1 ) ); //10,39+1-10=32
					newData.append ( gen_strtrim ( line.substr ( 0, idx1 ) ) + "\n" );
				}
			}
		}
	}
	const_cast <ParameterList*> (params)->setValue ( "data", newData );
	return sv;
}
int MSBridgeLink::num = 0;
StringVector MSBridgeLink::database;
MSBridgeLink::MSBridgeLink ( const string& db )
{
	num++;
	ind = num;
	database.push_back ( db );
}
void MSBridgeLink::write ( ostream& os, int indexNumber, int dnaReadingFrame, int openReadingFrame, const MultipleModification2& mm2, const string& str, int index )
{
	static int idx = 1;
	if ( idx == 1 ) {		// These are items which are passed in but unchanged between calls
		startJavascript ( os );
		os << "function printMSBridgeLinkItems () {" << endl;
			mm2.putHiddenFormJavascriptEntry ( os );
		os << "}" << endl;
		endJavascript ( os );
	}
	string variableName = "msbridge" + gen_itoa ( idx ) + "form";
	idx++;
	ProgramLink::putForm ( os, variableName, "mssearch.cgi" );
	startJavascript ( os );
	os << "printMSBridgeLinkItems ();" << endl;
	os << "printMSBridgeLinkItems" << index+1 << " ();" << endl;
	os << "printUnmatchedPeaks" << str << " ();" << endl;
	endJavascript ( os );

	bool dnaFlag = is_dna_database ( database [index] );
	SingleEntryParameters::putIndexHiddenFormEntry ( os, dnaFlag, indexNumber, dnaReadingFrame, openReadingFrame );

	os << "</form>" << endl;
	ProgramLink::putFormLink ( os, variableName, "Search for disulfide linked peptides." );
}
void MSBridgeLink::putHidden ( ostream& os ) const
{
	const ParameterList* params = ProgramLink::getParams ();
	string db = database [ind-1];
	os << "function printMSBridgeLinkItems" << ind << " () {" << endl;
		printHTMLFORMJavascriptHidden ( os, "search_name", "msbridge" );
		printHTMLFORMJavascriptHidden ( os, "output_type", "HTML" );
		printHTMLFORMJavascriptHidden ( os, "report_title", "MS-Bridge" );
		printHTMLFORMJavascriptHidden ( os, "version", Version::instance ().getVersion () );
		printHTMLFORMJavascriptHidden ( os, "link_search_type", "Disulfide (C)" );
/* Parameters got from html page */
		if ( !db.empty () )	printHTMLFORMJavascriptHidden2 ( os, "database", db );
		else				params->copyToHiddenFormJavascriptEntry ( os, "database", ind-1 );
		params->copyToHiddenFormJavascriptEntry ( os, "n_term_aa_limit" );
		AAInitInfo::copyToHiddenFormJavascriptEntry ( os, params, true );
		EnzymeParameters::copyToHiddenFormJavascriptEntry ( os, params, true );
		MassInfo::copyToHiddenFormJavascriptEntry ( os, params );
	os << "}" << endl;
}
