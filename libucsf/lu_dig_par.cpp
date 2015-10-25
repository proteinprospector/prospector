/******************************************************************************
*                                                                             *
*  Program    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_dig_par.cpp                                                *
*                                                                             *
*  Created    : June 13th 2001                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lgen_error.h>
#include <lu_pq_vector.h>
#include <lu_dig_par.h>
#include <lu_getfil.h>
#include <lu_check_db.h>
#include <lu_param_list.h>
#include <lu_cgi_val.h>
using std::ostream;
using std::string;

#ifdef VIS_C
#pragma warning( disable : 4503 )	// Disable warnings about decorated name length being too long
#endif

MSDigestParameters::MSDigestParameters ( const ParameterList* params ) :

	MSProgramParameters		( params ),

	aaInitInfo				( params ),
	singleEntryParameters	( params ),

	enzymeParameters		( params, true ),

	hideProteinSequence		( params->getBoolValue	( "hide_protein_sequence", false ) ),

	reportMultCharge		( params->getBoolValue	( "report_mult_charge", false ) ),
	hideHtmlLinks			( params->getBoolValue	( "hide_html_links", false ) ),

	digestFragmentParameters( params ),

	instrumentName			( params->getStringValue	( "instrument_name", "" ) ),

#ifdef CHEM_SCORE
	chemScoreParams			( params ),
#endif

	multipleModification	( params ),

	separateProteinsFlag	( params->getBoolValue	( "separate_proteins", false ) ),

	maxHits	( params->getIntValue	( "max_hits", 50000 ) )
{
	init_fasta_enzyme_function ( enzymeParameters.getEnzyme () );
	multipleModification.setUserMods ( aaInitInfo );
	const char* value;

	if ( params->getValue ( "coverage_map", value ) ) {
		IntVector tempArray;

		getPostQueryVector ( value, tempArray, ' ' );

		int val = tempArray [0];
		coverageMap.clear ();
		for ( IntVectorSizeType i = 1 ; i < tempArray.size () ; i++ ) {
			for ( int j = 0 ; j < tempArray [i] ; j++ ) {
				coverageMap.push_back ( val );
			}
			val = !val;
		}
	}
	bullBreeseFlag = false;
	if ( params->getBoolValue ( "bull_breese", false ) ) {
		bullBreeseFlag = true;
	}
	PotentialMSFragment::setBullBreese ( bullBreeseFlag );
	HPLCFlag = false;
	if ( params->getBoolValue ( "hplc_index", false ) ) {
		HPLCFlag = true;
	}
	PotentialMSFragment::setHPLC ( HPLCFlag );
}
void MSDigestParameters::printHTML ( ostream& os ) const
{
	singleEntryParameters.printHTML ( os );
	multipleModification.printHTML ( os );
	enzymeParameters.printHTML ( os );
	aaInitInfo.printHTML ( os );
	digestFragmentParameters.printHTML ( os );
#ifdef CHEM_SCORE
	chemScoreParams.printHTML ( os );
#endif
}
bool MSDigestLink::init = false;
std::string MSDigestLink::urlParams;
std::string MSDigestLink::linkName = "msdigestLink";
int MSDigestLink::num = 0;
StringVector MSDigestLink::database;
MSDigestLink::MSDigestLink ( const string& programName, const string& db )
{
	num++;
	ind = num;
	database.push_back ( db );
	if ( init == false ) {
		initLinks ( programName );
		init = true;
	}
}
void MSDigestLink::printHTML ( ostream& os ) const
{
	os << linkName << ind << "=\"";
	os << ProgramLink::getURLStart ( "mssearch" );
	os << "?";
	if ( !urlParams.empty () ) os << urlParams << "&";	// User defined parameters
	const ParameterList* params = ProgramLink::getParams ();
	string db = database [ind-1];
	printCGIString ( os, "search_name", "msdigest" );
	printCGIString ( os, "output_type", "HTML" );
	printCGIString ( os, "report_title", "MS-Digest" );
	printCGIString ( os, "version", Version::instance ().getVersion () );
	if ( !db.empty () )	printCGIString ( os, "database", db );
	else				params->copyToCGI ( os, "database", ind-1 );
	params->copyToCGI ( os, "n_term_aa_limit" );
	AAInitInfo::copyToCGI ( os, params );
	EnzymeParameters::copyToCGI ( os, params, true );
#ifdef CHEM_SCORE
	ChemScoreParameters::copyToCGI ( os, params );
#endif
	MultipleModification2::copyToCGI ( os, params );
	os << "\";\n";
}
void MSDigestLink::write ( ostream& os, int indexNumber, int dnaReadingFrame, int openReadingFrame, int num )
{
	const ParameterList* params = ProgramLink::getParams ();
	bool dnaFlag = is_dna_database ( params->getStringValue ( "database", "" ) );
	ProgramLink::openLink ( os, linkName, num+1 );

	SingleEntryParameters::putIndexCGI ( os, dnaFlag, indexNumber, dnaReadingFrame, openReadingFrame );
	os << "\\\">" << indexNumber;
	ProgramLink::closeLink ( os );
}
void MSDigestLink::write ( ostream& os, int indexNumber, int dnaReadingFrame, int openReadingFrame, const CoverageMap& coverageMap, int num )
{
	const ParameterList* params = ProgramLink::getParams ();
	bool dnaFlag = is_dna_database ( params->getStringValue ( "database", "" ) );
	ProgramLink::openLink ( os, linkName, num+1 );

	coverageMap.putCGI ( os );

	SingleEntryParameters::putIndexCGI ( os, dnaFlag, indexNumber, dnaReadingFrame, openReadingFrame );
	os << "\\\">" << indexNumber;
	ProgramLink::closeLink ( os );
}
void MSDigestLink::initLinks ( const string& programName )
{
	GenCommentedIFStream ifs ( MsparamsDir::instance ().getParamPath ( "idxlinks.txt" ) );
	string line;
	for ( ; ; ) {
		if ( !ifs.getUncommentedLine ( line ) ) break;
		string fileProg = line;
		if ( !ifs.getUncommentedLine ( line ) ) break;
		int end = line.find ( "?" );
		urlParams = line.substr ( end + 1, line.length () - ( end + 1 ) );
		if ( fileProg == programName ) return;
	}
	ErrorHandler::genError ()->error ( "The program " + programName + " does not occur in the parameter file idxlinks.txt.\n" );
}
MSDigestLinkNameValueStream::MSDigestLinkNameValueStream ( bool process ) :
	GenCommentedIFStream ( MsparamsDir::instance ().getParamPath ( "idxlinks.txt" ) )
{
	string line;
	bool nameFlag = true;
	string name;
	string value;
	while ( getUncommentedLine ( line ) ) {
		if ( nameFlag ) {
			name = line;
			nameFlag = false;
		}
		else {
			if ( process ) {
				int end = line.find ( "?" );
				value = line.substr ( end + 1, line.length () - ( end + 1 ) );
			}
			else {
				value = line;
			}
			params [name].push_back ( value );
			nameFlag = true;
		}
	}
}
string MSDigestLinkNameValueStream::getStringValue ( const string& name, const string& defaultVal ) const
{
	MapStringToStringVectorConstIterator cur = params.find ( name );
	if ( cur != params.end () )
		return (*cur).second [0];
	else
		return defaultVal;
}
StringVector MSDigestLinkNameValueStream::getNameList () const
{
	StringVector sv;
	for ( MapStringToStringVectorConstIterator i = params.begin () ; i != params.end () ; i++ ) {
		sv.push_back ( (*i).first );
	}
	return sv;
}
