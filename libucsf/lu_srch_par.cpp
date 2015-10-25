/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_srch_par.cpp                                               *
*                                                                             *
*  Created    : June 16th 2001                                                *
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
#include <lu_check_db.h>
#include <lu_srch_par.h>
#include <lu_cgi_val.h>
#include <lu_html_form.h>
#include <lu_param_list.h>
#include <lu_const_mod.h>
using std::ostream;
using std::string;

MSSearchParameters::MSSearchParameters ( const ParameterList* params, const string& prefix, const StringVector& db ) :
	MSProgramParameters	( params ),

	prefix			( prefix ),

	aaInitInfo		( params ),
	database		( db.empty () ? params->getStringVectorValue ( "database" ) : db ),
	enzymeParameters( params ),
	tempOverride	( params->getBoolValue	( "temp_override", false ) ),
	maxReportedHits	( params->getIntValue	( prefix + "max_reported_hits", 50 ) ),
	maxHits			( params->getIntValue	( "max_hits", 200 ) ),
	comment			( params->getStringValue	( "comment", "" ) ),

	preSearchInfo	( params, prefix ),
	preSearchOnly	( params->getBoolValue	( "pre_search_only", false ) ),

	detailedReport	( params->getBoolValue	( "detailed_report", false ) ),
	dnaFrameTranslation	( params->getIntValue( "dna_frame_translation", 3 ) ),
	maxNTermAA ( params->getIntValue ( "n_term_aa_limit", 0 ) ),

	singleEntryParameters ( params ),
	preSearchDone ( false )
{
	if ( preSearchOnly ) PreSearchInfo::setReportTaxonomy ();
}
MSSearchParameters::~MSSearchParameters ()
{
}
void MSSearchParameters::putCGI ( ostream& os, bool noEnzyme ) const
{
	aaInitInfo.putCGI ( os );
	for ( int i = 0 ; i < database.size () ; i++ ) {
		printCGIString ( os, "database", database [i] );
	}
	printCGI ( os, "n_term_aa_limit", maxNTermAA );
	if ( noEnzyme )
		enzymeParameters.putNoEnzymeCGI ( os );
	else
		enzymeParameters.putCGI ( os );

	printCGI ( os, "dna_frame_translation", dnaFrameTranslation );
	printCGI ( os, prefix + "max_reported_hits", maxReportedHits );
	printCGI ( os, "max_hits", maxHits );
	printCGIString ( os, "comment", comment );
	preSearchInfo.putCGI ( os );
}
void MSSearchParameters::putHiddenFormJavascriptEntry ( ostream& os, bool noEnzyme ) const
{
	aaInitInfo.putHiddenFormEntryJavascript ( os );
	printHTMLFORMJavascriptHiddenContainer2 ( os, "database", database );
	printHTMLFORMJavascriptHidden ( os, "n_term_aa_limit", maxNTermAA );
	if ( noEnzyme )
		enzymeParameters.putNoEnzymeHiddenFormEntryJavascript ( os );
	else
		enzymeParameters.putHiddenFormEntryJavascript ( os );

	printHTMLFORMJavascriptHidden ( os, "dna_frame_translation", dnaFrameTranslation );
	printHTMLFORMJavascriptHidden ( os, prefix + "max_reported_hits", maxReportedHits );
	printHTMLFORMJavascriptHidden ( os, "max_hits", maxHits );
	printHTMLFORMJavascriptHidden ( os, "comment", comment );
	preSearchInfo.putHiddenFormEntryJavascript ( os );
}
void MSSearchParameters::printHTML ( ostream& os ) const
{
	ParameterList::printHTMLContainer ( os, "Database searched", database );
	enzymeParameters.printHTML ( os );
	aaInitInfo.printHTML ( os );
	for ( StringVectorSizeType i = 0 ; i < database.size () ; i++ ) {
		if ( is_dna_database ( database [i] ) ) {
			PairIntInt pii = getDNAFrameTranslationPair ( database [i], dnaFrameTranslation );
			ParameterList::printHTML ( os, "DNA Frame Translation Start", pii.first );
			ParameterList::printHTML ( os, "DNA Frame Translation End", pii.second );
			break;
		}
	}
	if ( comment != "" ) ParameterList::printHTML ( os, "Sample ID (comment)", comment );
}
void MSSearchParameters::doSearch ( const FastaServerPtrVector fs ) const
{
	if ( !preSearchDone ) {
		preSearchInfo.doSearch ( fs );
		preSearchDone = true;
	}
}
