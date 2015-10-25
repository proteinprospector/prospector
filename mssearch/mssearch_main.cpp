/******************************************************************************
*                                                                             *
*  Program    : mssearch                                                      *
*                                                                             *
*  Filename   : mssearch_main.cpp                                             *
*                                                                             *
*  Created    : September 11th 2001                                           *
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
#ifndef VIS_C
#include <stdexcept>
#endif
#ifdef PP_MULTI
#include <lp_btag_mpi.h>
#else
#include <lp_frame.h>
#endif
#include <lu_comp_srch.h>
#include <lu_comp_par.h>
#include <lu_dbst_srch.h>
#include <lu_hom_srch.h>
#include <lu_hom_par.h>
#include <lu_patt_srch.h>
#include <lu_patt_par.h>
#include <lu_prot_srch.h>
#include <lu_tag_par.h>
#include <lu_tag_srch.h>
#include <lu_fit_srch.h>
#include <lu_filter_par.h>
#include <lu_filter_srch.h>
#include <lu_viewer_par.h>
#include <lu_viewer_srch.h>
#include <lu_iso_srch.h>
#include <lu_bdg_srch.h>
#include <lu_nspec_srch.h>
#include <lu_prod_srch.h>
#include <lu_script.h>
#include <lu_expec_par.h>
#include <lu_param_list.h>
#include <lu_getfil.h>
#include <lu_html.h>
#ifdef BATCHTAG
#include <lu_btag_run.h>
#include <ld_init.h>
#endif
using std::runtime_error;
using std::string;
using std::cout;
//using std::ios_base;

namespace {
#ifdef BATCHTAG
void runBTag ( ParameterList* paramList, int searchStage, const string& searchJobID, int startSerial )
{
	init_html_premature_stop ( paramList->getStringValue ( "report_title" ), true );
#ifdef MYSQL_DATABASE
	if ( !searchJobID.empty () ) MySQLPPSDDBase::instance ().updateSearchStage ( searchJobID, searchStage );
#endif
	runBatchTag ( paramList, 1, 0, searchJobID, startSerial );
	joinResultsFiles ( paramList, 1 );
}
#endif
void runProspectorProgram ( ParameterList* paramList, const string& searchName )
{
	string outputType = paramList->getStringValue ( "output_type" );
	string reportTitle = paramList->getStringValue ( "report_title" );
	bool resultsToFile = paramList->getBoolValue ( "results_to_file" );
	if ( outputType == "HTML" || outputType == "Tab delimited text" || outputType == "mgf" || ( resultsToFile && outputType == "XML" ) ) {
		if ( searchName == "msproduct" || searchName == "msisotope" || searchName == "msviewer" ) {
			if ( resultsToFile )
				init_html ( cout, reportTitle );
			else
				init_html ( cout, "" );
		}
		else
			init_html_premature_stop ( reportTitle, true );
	}
	MSProgram* ds;
	if ( searchName == "msfit" ) {
		MSFitParameters params ( paramList );

		ds = new FitSearch ( params );
		ds->outputResults ();
	}
	else if ( searchName == "mstag" || searchName == "msseq" ) {
		MSTagParameters params ( paramList );
		ds = getTagSearch ( params );
		ds->outputResults ();
	}
	else if ( searchName == "mspattern" ) {
		MSPatternParameters params ( paramList );
		if ( params.getPreSearchOnly () ) 
			ds = new ProteinSearch ( params );
		else
			ds = new PatternSearch ( params );

		ds->outputResults ();
	}
	else if ( searchName == "mshomology" ) {
		MSHomologyParameters params ( paramList );

		ds = new HomologySearch ( params );
		ds->outputResults ();
	}
	else if ( searchName == "dbstat" ) {
		DBStatParameters params ( paramList );

		ds = new DBStatSearch ( params );
		ds->outputResults ();
	}
	else if ( searchName == "mscomp" ) {
		MSCompParameters params ( paramList );
		ds = new MSCompSearch ( params );
		ds->outputResults ();
	}
	else if ( searchName == "msisotope" ) {
		MSIsotopeParameters params ( paramList );
		ds = new MSIsotopeSearch ( params );
		ds->outputResults ();
	}
	else if ( searchName == "msdigest" ) {
		MSDigestParameters params ( paramList );
		ds = new MSDigestSearch ( params );
		ds->outputResults ();
	}
	else if ( searchName == "msbridge" ) {
		MSBridgeParameters params ( paramList );
		ds = new MSBridgeSearch ( params );
		ds->outputResults ();
	}
	else if ( searchName == "msnonspecific" ) {
		MSNonSpecificParameters params ( paramList );
		ds = new MSNonSpecificSearch ( params );
		ds->outputResults ();
	}
	else if ( searchName == "msviewer" ) {
		MSViewerParameters params ( paramList );
		ds = new MSViewerSearch ( params );
		ds->outputResults ( false );
	}
	else if ( searchName == "msfilter" ) {
		MSFilterParameters params ( paramList );
		ds = new MSFilterSearch ( params );
		ds->outputResults ( false );
	}
	else if ( searchName == "msproduct" ) {
		MSProductParameters params ( paramList );
		ds = new MSProductSearch ( params );
		ds->outputResults ( false );
	}
	else if ( searchName == "msmultisample" ) {
		MSFitParameters fitParams ( paramList );
		ds = new FitSearch ( fitParams );
		ds->outputResults ();
	}
	else {
		ErrorHandler::genError ()->error ( "Invalid search_name parameter.\n" );
	}
}
}

int main ( int argc, char** argv )
{
#ifdef PP_MULTI
	if ( isSuffix ( argv [0], "mssearchmpi.cgi" ) ) {
		BatchTagMPI btmpi ( argc, argv );
		btmpi.run ();
	}
	else {
#else
#ifdef BATCHTAG
		if ( argc != 1 && string ( argv [1] ) == "-k" )			initialiseProspectorBTag ( argv [2] );	// Single Processor Batch-Tag
		else if ( argc != 1 && string ( argv [1] ) == "-c" )	initialiseProspectorBTag ();	// Single Processor Batch-Tag
		else
#endif
#endif
			initialiseProspector ();
		ParameterList paramList ( argc, argv );
		try {
			if ( paramList.empty () ) {
				ErrorHandler::genError ()->error ( "No parameters passed to Prospector program.\n" );
			}
			//string ver = paramList.getStringValue ( "version" );
			//if ( !ver.empty () && ver != Version::instance ().getVersion () ) {
			//	ErrorHandler::genError ()->error ( "Form version mismatch. Try reloading the form.\n" );
			//}
			ProgramLink::setParams ( &paramList );
			MSProgram::setParams ( &paramList );

			if ( paramList.getBoolValue ( "create_script", false ) ) {
				writeScript ( &paramList );
			}
			else if ( paramList.getBoolValue ( "create_params", false ) ) { 
				writeParamsXML ( &paramList, paramList.getStringValue ( "search_name", "" ) );
			}
			else {
				string searchName = paramList.getStringValue ( "search_name", "" );
#ifdef BATCHTAG
				if ( searchName == "batchtag" ) {
					string searchKey = paramList.getStringValue ( "search_key" );
					string searchJobID;
					int startSerial = 1;
					bool expectationSearchFirst = InfoParams::instance ().getBoolValue ( "expectation_search_first" );
					bool expectationSearchDone = false;
#ifdef MYSQL_DATABASE
					if ( !searchKey.empty () ) {
						JobItem* jobItem = MySQLPPSDDBase::instance ().getSearchJobByKey ( searchKey );
						searchJobID = jobItem->getSearchJobID ();
#ifndef PP_MULTI
						FrameIterator::setSearchJobID ( searchJobID );
#endif
						int searchStage = jobItem->getSearchStage ();
						startSerial = jobItem->getStartSerial ();
						expectationSearchDone = expectationSearchFirst && ( searchStage == 2 );
					}
#endif
					ParameterList* expParamList = 0;
					if ( paramList.getStringValue ( "expect_calc_method" ) != "None" && !expectationSearchDone ) {
						string outputFilename;
						expParamList = getExpectationParams ( &paramList, outputFilename, startSerial != 1 );
						paramList.addOrReplaceName ( "expect_coeff_file", outputFilename );
					}
					if ( !expectationSearchFirst )	runBTag ( &paramList, 2, searchJobID, startSerial );
					if ( expParamList )				runBTag ( expParamList, 1, searchJobID, startSerial );
					if ( expectationSearchFirst )	runBTag ( &paramList, 2, searchJobID, startSerial );
#ifdef MYSQL_DATABASE
					if ( !searchKey.empty () ) MySQLPPSDDBase::instance ().setJobDone ( searchJobID );
#endif
				}
				else
#endif
					runProspectorProgram ( &paramList, searchName );
			}
		}
		catch ( runtime_error e ) {
			paramList.writeLogError ( e.what () );
		}
		//catch ( ios_base::failure& e2 ) {
		//	paramList.writeLogError ( e2.what () );
		//}
#ifdef PP_MULTI
	}
#endif
	DBSearch::deleteTempDirs ();
	return 0;
}
