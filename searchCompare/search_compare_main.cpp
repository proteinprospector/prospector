/******************************************************************************
*                                                                             *
*  Program    : searchCompare                                                 *
*                                                                             *
*  Filename   : search_compare_main.cpp                                       *
*                                                                             *
*  Created    : March 13th 2003                                               *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2003-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <memory>
#ifdef VIS_C
#pragma warning( disable : 4503 )	// Disable warnings about decorated name length being too long
#endif
#ifndef VIS_C
#include <stdexcept>
#endif
#include <lu_getfil.h>
#include <lu_html.h>
#include <lu_prog_par.h>
#include <lu_proj_file.h>
#include <lu_version.h>
#include <lu_param_list.h>
#include <ld_init.h>
#include <lu_mzidentml.h>
#include <sc_sres_rep.h>
using std::vector;
using std::string;
using std::ostream;
using std::cout;
using std::endl;
using std::sort;
using std::set_intersection;
using std::back_inserter;
using std::runtime_error;
using std::auto_ptr;

static void getSearchResults ( const SearchCompareParams& params, vector <SearchResults*>& searchResults );
static void writeReport ( const SearchCompareParams& params, const vector <SearchResults*>& searchResults );
static void doReport ( ostream& os, const SearchCompareParams& params, const StringVector& idUsedList, const vector <SearchResultsReport*>& srr );
static void writeParamsXML ( ParameterList* paramList );

class sortIDs {
	public:
		bool operator () ( const string& lhs, const string& rhs ) const
		{
			int end = lhs.find ( '-' );
			if ( end != string::npos ) {
				int lhsFrac = atoi ( lhs.substr ( 0, end ).c_str () ); 
				int lhsSpot = atoi ( lhs.substr ( end+1 ).c_str () );
				end = rhs.find ( '-' );
				int rhsFrac = atoi ( rhs.substr ( 0, end ).c_str () );
				int rhsSpot = atoi ( rhs.substr ( end+1 ).c_str () );
				if ( lhsFrac == rhsFrac ) return lhsSpot < rhsSpot;
				else return lhsFrac < rhsFrac;
			}
			else {
				return atoi ( lhs.c_str () ) < atoi ( rhs.c_str () );
			}
		}
};

int main ( int argc, char** argv )
{
	initialiseProspector ();
	ParameterList pList ( argc, argv );
	try {
		if ( pList.getBoolValue ( "save_params" ) ) {
			init_html ( cout, "Saving Search Compare Parameters" );
			pList.removeName ( "save_params" );
			ParameterList cookieParamList ( "", false, false );		// Create empty parameter set
			cookieParamList.appendParameters ( pList );
			cookieParamList.removeName ( "accession_nums" );
			cookieParamList.removeName ( "remove" );
			cookieParamList.removeName ( "id_filter_list" );
			cookieParamList.removeName ( "data" );
			cookieParamList.removeName ( "version" );
			bool ret = cookieParamList.copyToCookie ( cout, "search_compare_params" );
			if ( !ret ) {
				ErrorHandler::genError ()->error ( "Could not save the parameters as their length exceeds the maximum cookie length.\n" );
			}
			cout << "<p>Settings saved</p>" << endl;
			cout << "<input type=\"button\" value=\"Search Form\" onclick=\"history.go(-1)\">" << endl;
			exit ( 0 );
		}
		if ( pList.empty () ) {
			ErrorHandler::genError ()->error ( "No parameters passed to Search Compare.\n" );
		}
		ProgramLink::setParams ( &pList );
		init_html_premature_stop ( "Search Compare", true );
		SearchCompareParams params ( &pList );
		sresFPR = params.getReportType () == "False Positive Rate";
		sresTime = params.getReportType () == "Time" || sresFPR;
		sresXLinks = params.getReportType () == "Crosslinked Peptides";
		sresProt = params.getReportType () == "Protein";
		sresKeepReplicates = params.getPeptideFilter () == "Keep Replicates" || sresTime;
		sresKeepCharges = params.getPeptideFilter () == "Best Per Charge";
		sresKeepTimeReplicates = params.getPeptideFilter () == "Keep Replicates";
		sresViewer = params.getSaveFormat () == 'V';
		ujm = new UpdatingJavascriptMessage;
		vector <SearchResults*> searchResults;
		getSearchResults ( params, searchResults );
		writeReport ( params, searchResults );
		delete ujm;
		printProgramInformationHTML ( cout, "Search Compare" );
		ProteinInfo::deleteTempDirs ();
	}
	catch ( runtime_error e ) {
		pList.writeLogError ( e.what () );
	}
	return 0;
}
static void getSearchResults ( const SearchCompareParams& params, vector <SearchResults*>& searchResults )
{
	sresSingleProject = true;
	sresMainAndSupplementary = ( params.getSaveFormat () != 'P'  && params.getSaveFormat () != 'B' );
	sresMergedFlag = params.getMergeOption () == "Merged";
	string proj;
	StringVector filenames = params.getFilenames ();
#ifdef MYSQL_DATABASE
	for ( StringVectorSizeType i = 0 ; i < filenames.size () ; i++ ) {
		BatchJobItem* bji;
		string searchEndTime;
		string searchTime;
		try {
			bji = MySQLPPSDDBase::instance ().getBatchJobByKey ( filenames [i] );
			searchEndTime = MySQLPPSDDBase::instance ().getSearchEndTimeByKey ( filenames [i] );
			searchTime = MySQLPPSDDBase::instance ().getSearchTimeByKey ( filenames [i] );
		}
		catch ( runtime_error e ) {		// Catch database login problems
			ErrorHandler::genError ()->error ( e );
		}
		if ( !bji ) {
			ErrorHandler::genError ()->error ( "Unknown search key.\n" );
		}
		string temp = bji->getUserID () + '/' + getUncalibratedProject ( bji->getProjectName () );
		if ( i != 0 && temp != proj ) sresSingleProject = false;
		proj = temp;
		MySQLPPSDDBase::instance ().updateProjectRecordUpdated ( bji->getProjectID () );
		searchResults.push_back ( new SearchResults ( bji->getProjectName (), bji->getResultsName (), bji->getResultsFullPath (), params, i, searchEndTime, searchTime ) );
	}
#endif
	if ( filenames.empty () ) {
		int numSearches = params.getNumCommandLineSearches ();
		StringVector projectNames = params.getCommandLineProjectNames ();
		StringVector resultsNames = params.getCommandLineResultsNames ();
		StringVector resultsFullPaths = params.getCommandLineResultsFullPaths ();
		for ( StringVectorSizeType i = 0 ; i < numSearches ; i++ ) {
			string temp = getUncalibratedProject ( projectNames [i] );
			if ( i != 0 && temp != proj ) sresSingleProject = false;
			proj = temp;
			searchResults.push_back ( new SearchResults ( projectNames [i], resultsNames [i], resultsFullPaths [i], params, i, "", "" ) );
		}
		if ( searchResults.empty () ) ErrorHandler::genError ()->error ( "No results to display.\n" );
	}
	if ( sresTime && !sresMergedFlag && !sresSingleProject ) {
		ErrorHandler::genError ()->error ( "Time reports using the separated option can only compare results files from the same project.\n" );
	}
	SearchResults::mergeResults ( params );
}
static void writeReport ( const SearchCompareParams& params, const vector <SearchResults*>& searchResults )
{
	bool remove = params.getRemove ();
	bool multisample = params.getMultiSample ();
	string reportHomologousProteins = params.getReportHomologousProteins ();
	string reportType = params.getReportType ();
	char format = params.getSaveFormat ();
	bool quanProteinFlag = PPProteinHitQuanInfo::getQuan ();

	StringVector idList = searchResults [0]->getIDList ();
	StringVector idFilterList = params.getIDFilterList ();
	sort ( idFilterList.begin (), idFilterList.end () );
	StringVector idUsedList;
	if ( !idFilterList.empty () && multisample )
		set_intersection ( idList.begin (), idList.end (), idFilterList.begin (), idFilterList.end (), back_inserter ( idUsedList ) );
	else
		idUsedList = idList;
	sort ( idUsedList.begin (), idUsedList.end (), sortIDs () );

	if ( reportType == "Calibration" ) {
		for ( StringVectorSizeType i = 0 ; i < idUsedList.size () ; i++ ) {
			SearchResultsPeptideReport srpr ( searchResults, idUsedList [i] );
			printAbortFunctions ( cout );
			srpr.printHTMLCalibration ( cout );
		}
	}
	else {
		vector <SearchResultsReport*> srr;
		for ( StringVectorSizeType i = 0 ; i < idUsedList.size () ; i++ ) {
			if ( reportType == "Protein" ) {
				if ( quanProteinFlag ) {
					srr.push_back ( new SearchResultsPeptideReport ( searchResults, remove, params.getAccessionNumbers (), params.getSortType (), params.getSortType2 (), params.getReportHitsType (), reportHomologousProteins, idUsedList [i] ) );
				}
				else {
					srr.push_back ( new SearchResultsProteinReport ( searchResults, remove, params.getAccessionNumbers (), params.getReportHitsType (), reportHomologousProteins, idUsedList [i] ) );
				}
			}
			else if ( reportType == "Peptide" || sresXLinks ) {
				srr.push_back ( new SearchResultsPeptideReport ( searchResults, remove, params.getAccessionNumbers (), params.getSortType (), params.getSortType2 (), params.getReportHitsType (), reportHomologousProteins, idUsedList [i] ) );
			}
			else if ( sresTime ) {
				srr.push_back ( new SearchResultsPeptideReport ( searchResults, remove, params.getAccessionNumbers (), params.getSortType (), params.getSortType2 (), params.getUnmatchedSpectra (), params.getReportHitsType (), reportHomologousProteins, idUsedList [i], format != 'P' && format != 'F' && format != 'V' ) );
			}
			else ErrorHandler::genError ()->error ( "Invalid report type.\n" );
		}
		printAbortFunctions ( cout );
		ujm->deletePreviousMessage ( cout );
		doReport ( cout, params, idUsedList, srr );
		for ( std::vector <SearchResultsReport*>::size_type j = 0 ; j < srr.size () ; j++ ) delete srr [j];
	}
}
static void doReport ( ostream& os, const SearchCompareParams& params, const StringVector& idUsedList, const vector <SearchResultsReport*>& srr )
{
	char format = params.getSaveFormat ();
#ifdef MYSQL_DATABASE
	bool checkboxes = params.getCheckboxes ();
#endif
	StringVectorSizeType numSections = idUsedList.size ();
	bool delimHeaderPrinted = false;
	for ( StringVectorSizeType i = 0 ; i < numSections ; i++ ) {
		if ( format == 'H' ) {
			if ( i == 0 ) {
				srr [i]->printHistogramHTML ( os );
#ifdef MYSQL_DATABASE
				if ( checkboxes ) srr [i]->printReportHeader ( os );
#endif
			}
			srr [i]->printHTML ( os );
			if ( i == numSections-1 ) {
#ifdef MYSQL_DATABASE
				if ( checkboxes ) srr [i]->printReportFooter ( os );
#endif
			}
		}
		else if ( format == 'P' ) {
			if ( sresSingleProject ) srr [i]->printPepXML ( os, params.getOutputDirectory (), params.getOutputFilename () );
			else
				ErrorHandler::genError ()->error ( "Pep XML reports can only be generated for results files from a single project.\n" );
		}
		else if ( format == 'M' ) {
			if ( sresSingleProject ) srr [i]->printMZIdentML ( os, params.getOutputDirectory (), params.getOutputFilename () );
			else
				ErrorHandler::genError ()->error ( "mzIdentML reports can only be generated for results files from a single project.\n" );
		}
		else if ( format == 'B' ) {
			srr [i]->writeBiblioSpec ( os, params.getOutputDirectory (), params.getOutputFilename (), idUsedList [i], params.getNormalizedBiblioSpec () );
		}
		else if ( format == 'F' || format == 'V' ) {
			if ( params.getReportType () == "Peptide" || params.getReportType () == "Time" || params.getReportType () == "Crosslinked Peptides" ) {
				srr [i]->printMGF ( os, params.getOutputDirectory () );
				if ( format == 'V' ) {
					srr [i]->printDelimitedReport ( os, i, i == ( numSections - 1 ), params.getOutputDirectory (), params.getOutputFilename (), delimHeaderPrinted );
				}
			}
			else
				ErrorHandler::genError ()->error ( "Peak list files can only be generated from Peptide and Time reports.\n" );
		}
		else
			srr [i]->printDelimitedReport ( os, i, i == ( numSections - 1 ), params.getOutputDirectory (), params.getOutputFilename (), delimHeaderPrinted );
	}
}
static void writeParamsXML ( ParameterList* paramList )
{
	string filename = paramList->getStringValue ( "settings_filename", "" );
	if ( !filename.empty () ) {
		paramList->removeName ( "save_settings" );
		paramList->removeName ( "settings_filename" );
		//GenOFStream ost ( ResultsDir::instance ().getResultsDir ( "search_compare" ) + SLASH + filename + ".xml" );
		//printXMLHeader ( ost );
		//printXMLVersion ( ost );
		//paramList->XMLParameters ( ost );
	}
	else {
		ErrorHandler::genError ()->error ( "You must specify a file name when saving the settings.\n" );
	}
}
