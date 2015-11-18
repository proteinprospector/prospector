/******************************************************************************
*                                                                             *
*  Program    : searchCompare                                                 *
*                                                                             *
*  Filename   : sc_params.cpp                                                 *
*                                                                             *
*  Created    : April 9th 2003                                                *
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
#ifdef VIS_C
#pragma warning( disable : 4503 )	// Disable warnings about decorated name length being too long
#endif
#include <lgen_error.h>
#ifdef MYSQL_DATABASE
#include <ld_init.h>				// Only required during pepXML development
#endif
#include <lu_getfil.h>
#include <lu_usermod.h>
#include <lu_pre_srch.h>
#include <lu_quan_ratio.h>
#include <lu_param_list.h>
#include <sc_params.h>
#include <sc_sres_rep.h>
#include <sc_xlink.h>
#include <sc_quan.h>
using std::string;
using std::endl;
using std::cout;

SearchCompareParams::SearchCompareParams ( const ParameterList* params ) :
	saveFormat				( params->getStringValue		( "save_format", "HTML" ) ),
	outputDirectory			( params->getStringValue		( "output_directory", "" ) ),
	outputFilename			( params->getStringValue		( "output_filename", "" ) ),
	reportType				( params->getStringValue		( "report_type", "Protein" ) ),
	reportHomologousProteins( params->getStringValue		( "report_homologous_proteins", "Interesting" ) ),
	databaseType			( params->getStringVectorValue	( "database_type" ) ),
	multiSample				( params->getBoolValue			( "multi_sample", false ) ),
	bestDiscOnly			( params->getBoolValue			( "best_disc_only" ) ),
	discScoreGraph			( params->getBoolValue			( "disc_score_graph" ) && saveFormat == "HTML" ),
	modificationScoreThreshold( params->getDoubleValue		( "slip_threshold", 0.0 ) ),
	maxProteinEValue		( params->getDoubleValue		( "max_protein_evalue", 0.01 ) ),
	maxPeptideEValue		( params->getDoubleValue		( "max_peptide_evalue", 0.05 ) ),
	minPPProteinScore		( params->getDoubleValue		( "min_protein_score", 0.0 ) ),
	minPPPeptideScore		( params->getDoubleValue		( "min_peptide_score", 0.0 ) ),
	xlMinLowScore			( params->getDoubleValue		( "xl_min_low_score", 0.0 ) ),
	xlMinScoreDiff			( params->getDoubleValue		( "xl_min_score_diff", 0.0 ) ),
	xlMaxLowExpectation		( params->getDoubleValue		( "xl_max_low_expectation", 1.0 ) ),
	remove					( params->getBoolValue			( "remove" ) ),
	peptideFilter			( params->getStringValue		( "peptide_filter" ) ),
	reportHitsType			( params->getStringValue		( "report_hits_type", "" ) ),
	mergeOption				( params->getStringValue		( "merge_option", "Separated" ) ),
	sortType				( params->getStringValue		( "sort_type", "" ) ),
	sortType2				( params->getStringValue		( "sort_type_2", "" ) ),
	unmatchedSpectra		( params->getBoolValue			( "unmatched_spectra" ) ),
	resolution				( params->getDoubleValue		( "resolution", 10000.0 ) ),
	reportMMods				( params->getBoolValue			( "report_mass_mod" ) ),
	reportCoverage			( params->getBoolValue			( "report_coverage" ) ),
	checkboxes				( reportType == "Peptide" ? params->getBoolValue ( "report_checkboxes", false ) : false ),
	aaii					( params ),
	isotopePurity			( params ),
	purityCorrection		( params ),
	reporterIonWindow		( params->getDoubleValue ( "reporter_ion_window", 0.4 ) )
{
	if ( isPrefix ( saveFormat, "pepXML" ) ) {
#ifdef MYSQL_DATABASE
		//if ( MySQLPPSDDBase::instance ().getUserName ( ParameterList::getCookieValue ( "key" ) ) != "pepxml" ) {
		//	ErrorHandler::genError ()->error ( "pepXML output is currently under development and is only available to certain users.\n" );
		//}
#endif
		reportType = "Time";
		sortType = "Time";
		sortType2 = "";
		unmatchedSpectra = false;
		mergeOption = "Merged";
		if ( saveFormat == "pepXML Unfiltered" ) {
			peptideFilter = "Keep Replicates";
			bestDiscOnly = false;
			maxProteinEValue = std::numeric_limits<double>::max();
			maxPeptideEValue = std::numeric_limits<double>::max();
			minPPProteinScore = -std::numeric_limits<double>::max();
			minPPPeptideScore = -std::numeric_limits<double>::max();
		}
	}
	if ( isPrefix ( saveFormat, "mzIdentML" ) || isPrefix ( saveFormat, "BiblioSpec" ) ) {
		reportType = "Time";
		sortType = "Time";
		sortType2 = "";
		unmatchedSpectra = false;
		mergeOption = "Merged";
		discScoreGraph = false;
	}
	if ( isPrefix ( saveFormat, "BiblioSpec" ) ) {
		reportCoverage = false;
	}
	Usermod::initialiseAllUsermodAAInfo ();
	string rawType = params->getStringValue ( "raw_type", "" );
	string quanType = params->getStringValue ( "quan_type", "" );
	QuantitationMulti::setPurityCorrection ( &purityCorrection );
	QuantitationMultiMassWindow::setReporterIonWindow ( reporterIonWindow );
	PeptidePositionQuan::initialiseAACalculator ( aaii.getConstMods (), rawType, quanType, resolution );
	PeptidePosition::initialiseAACalculator ( aaii.getConstMods () );
	if ( reportType == "Crosslinked Peptides" ) {
		SearchResultsCrosslinkPeptideHit::initialiseAACalculator ( aaii.getConstMods () );
	}
	PeptidePosition::initialiseComposition ( params->getStringVectorValue ( "comp_ion" ), params->getPQStringVectorValue ( "mass_comp_list" ), params->getStringValue ( "comp_mask_type", "AND" ) );
	filenames = params->getPQStringVectorValue ( "data" );
	idFilterList = params->getPQStringVectorValue ( "id_filter_list" );
	for ( StringVectorSizeType i = 0 ; i < idFilterList.size () ; i++ ) {
		idFilterSet.insert ( atoi ( idFilterList [i].c_str () ) );
	}
	setReportItems ( params );
	PreSearchInfo::setAccessionNumbers ( accessionNumbers, params->getPQStringVectorValue ( "accession_nums" ) );
	PeptidePosition::setRunMSProductFlag ( params->getBoolValue ( "run_msproduct" ) );
	PeptidePosition::initialiseDataSetInfo ( params->getDoubleValue ( "rt_int_start", 0.0 ), params->getDoubleValue ( "rt_int_end", 0.0 ) );
	StringVector instList = DiscriminantScoreInstrumentList::instance ().getNames ();
	for ( StringVectorSizeType j = 0 ; j < instList.size () ; j++ ) {
		if ( saveFormat == "pepXML Unfiltered" )
			minBestDiscScore [instList [j]] = -std::numeric_limits<double>::max();
		else
			minBestDiscScore [instList [j]] = params->getDoubleValue ( "min_best_disc_score_" + genTranslateCharacter(instList [j],'-','_'), 0.0 );
	}
	initPeptidesToDelete ( params->getStringVectorValue ( "cb" ) );

	projectNames	= params->getPQStringVectorValue ( "project_names" );
	resultsNames	= params->getPQStringVectorValue ( "results_names" );
	resultsFullPaths= params->getPQStringVectorValue ( "results_full_paths" );
	if ( reportType == "Modifications" ) SiteScores::init ( params->getStringVectorValue ( "comp_ion" ), true );
}
void SearchCompareParams::setReportItems ( const ParameterList* params )
{
	SearchResultsProteinInfo::setParams ( params );

	PPPeptideHitInfo::setReportUnmatched( params->getBoolValue ( "report_unmatched" ) );
	PPPeptideHitInfo::setReportNumPks	( params->getBoolValue ( "report_num_pks" ) );
	PPPeptideHitInfo::setReportRank		( params->getBoolValue ( "report_rank" ) );
	PPPeptideHitInfo::setReportScore	( params->getBoolValue ( "report_score" ) );
	PPPeptideHitInfo::setReportScoreDiff( params->getBoolValue ( "report_score_diff" ) );
	PPPeptideHitInfo::setReportExpectation( params->getBoolValue ( "report_expectation" ) );
	PPPeptideHitInfo::setReportPValue	( params->getBoolValue ( "report_p_value" ) );
	PPPeptideHitInfo::setReportMValue	( params->getBoolValue ( "report_nlog_p_value" ) );
	PPPeptideHitInfo::setReportNumPrecursor( params->getBoolValue ( "report_num_precursor" ) );
	PPPeptideHitInfo::setReportGradient	( params->getBoolValue ( "report_gradient" ) );
	PPPeptideHitInfo::setReportOffset	( params->getBoolValue ( "report_offset" ) );
	PPPeptideHitInfo::setReportDiscScore( params->getBoolValue ( "report_disc_score" ) );
	PPPeptideHitInfo::setReportRepeats	( params->getBoolValue ( "report_repeats" ) );

	PeptidePosition::setReportCheckboxes( checkboxes );
	PeptidePosition::setReportSearchNumber	( params->getBoolValue ( "report_search_number" ) );
	PeptidePosition::setReportMPlusH		( params->getBoolValue ( "report_m_plus_h" ) );
	PeptidePosition::setReportMOverZ		( params->getBoolValue ( "report_m_over_z" ) );
	PeptidePosition::setReportCharge		( params->getBoolValue ( "report_charge" ) );
	PeptidePosition::setReportMPlusHCalc	( params->getBoolValue ( "report_m_plus_h_calc" ) );
	PeptidePosition::setReportMOverZCalc	( params->getBoolValue ( "report_m_over_z_calc" ) );
	PeptidePosition::setReportIntensity		( params->getBoolValue ( "report_intensity" ) );
	PeptidePosition::setReportError			( params->getBoolValue ( "report_error" ) );
	PeptidePosition::setReportDBPeptide		( params->getBoolValue ( "report_db_peptide" ) );
	PeptidePosition::setPeptideModType		( params->getStringValue ( "peptide_mod_type" ) );
	PeptidePosition::setReportProteinMods	( params->getBoolValue ( "report_protein_mod" ) );
	PeptidePosition::setReportTime			( params->getBoolValue ( "report_time" ) );
	PeptidePosition::setReportMSMSInfo		( params->getBoolValue ( "report_msms_info" ) );
	PeptidePosition::setReportStartAA		( params->getBoolValue ( "report_start_aa" ) );
	PeptidePosition::setReportEndAA			( params->getBoolValue ( "report_end_aa" ) );
	PeptidePosition::setReportPreviousAA	( params->getIntValue ( "report_previous_aa" ) );
	PeptidePosition::setReportNextAA		( params->getIntValue ( "report_next_aa" ) );
	PeptidePosition::setReportElemComp		( params->getBoolValue ( "report_elem_comp" ) );
	PeptidePosition::setReportMissedCleavages( params->getBoolValue ( "report_missed_cleavages" ) );
	PeptidePosition::setReportMModValue		( params->getBoolValue ( "report_mass_mod" ) );
	PeptidePosition::setReportLength		( params->getBoolValue ( "report_length" ) );
	PeptidePosition::setReportComposition	( params->getBoolValue ( "report_composition" ) );
	PeptidePosition::setReportLinks			( params->getBoolValue ( "report_links" ) );

	PeakFitData::setReportPeakIntensity	( params->getBoolValue ( "rep_intensity" ) );
	PeakFitData::setReportPeakSNR		( params->getBoolValue ( "rep_snr" ) );
	PeakFitData::setReportPeakResolution( params->getBoolValue ( "rep_resolution" ) );
	PeakFitData::setReportPeakCSIntensity( params->getBoolValue ( "rep_cs_intensity" ) );
	PeakFitData::setReportPeakFWHM		( params->getBoolValue ( "rep_fwhm" ) );
	PeakFitData::setReportPeakArea		( params->getBoolValue ( "rep_area" ) );
	PeakFitData::setReportPeakCSArea	( params->getBoolValue ( "rep_cs_area" ) );
	PeakFitData::setReportNoiseMean		( params->getBoolValue ( "rep_n_mean" ) );
	PeakFitData::setReportStdDev		( params->getBoolValue ( "rep_n_stdev" ) );
	PeakFitData::setSNRThreshold		( params->getDoubleValue ( "snr_threshold" ) );
	PeakFitData::setReportFormulaString	( false );
	double areaThreshold = params->getDoubleValue ( "area_threshold" );
	PeakFitData::setAreaThreshold		( areaThreshold );
	double intensityThreshold = params->getDoubleValue ( "intensity_threshold" );
	PeakFitData::setIntensityThreshold	( intensityThreshold );
	if ( intensityThreshold != 0.0 || areaThreshold != 0.0 ) {
		if ( !isQuanMSMS ( params->getStringValue ( "quan_type", "" ) ) ) {
			ErrorHandler::genError ()->error ( "Intensity/Area thresholds not implemented for this quantitation type.\nUse a SNR threshold instead.\n" );
		}
	}
	QuantitationData::setReportActualLightHeavyIntensityRatio	( params->getBoolValue ( "rep_a_lh_int" ) );
	QuantitationData::setReportActualLightHeavyAreaRatio		( params->getBoolValue ( "rep_a_lh_area" ) );

	SearchResultsProteinLine::setReportNumber( params->getBoolValue ( "report_number" ) );
	SearchResultsProteinLine::setReportLinks ( params->getBoolValue ( "report_links" ) );

	SearchResultsPeptideLine::setReportNumber( params->getBoolValue ( "report_number" ) );
	SearchResultsPeptideLine::setReportLinks ( params->getBoolValue ( "report_links" ) );

	ProteinInfo::setReportUniprotID	( params->getBoolValue ( "report_uniprot_id" ) );
	ProteinInfo::setReportGeneName	( params->getBoolValue ( "report_gene_name" ) );
	ProteinInfo::setReportAccession	( params->getBoolValue ( "report_accession" ) );
	ProteinInfo::setReportVersion	( params->getBoolValue ( "report_version" ) );
	ProteinInfo::setReportIndex		( params->getBoolValue ( "report_index" ) );
	ProteinInfo::setReportLength	( params->getBoolValue ( "report_prot_len" ) );
	ProteinInfo::setReportMW		( params->getBoolValue ( "report_mw" ) );
	ProteinInfo::setReportPI		( params->getBoolValue ( "report_pi" ) );
	ProteinInfo::setReportSpecies	( params->getBoolValue ( "report_species" ) );
	ProteinInfo::setReportName		( params->getBoolValue ( "report_name" ) );
	ProteinInfo::setReportLinks		( params->getBoolValue ( "report_links" ) );
	ProteinInfo::setTaxonomyMatch	( params->getPQStringVectorValue ( "preferred_species" ) );

	if ( QuantitationRatio::getQuanReport () ) {
		PPProteinHitQuanInfo::setQuanParams ( params );
	}
}
void SearchCompareParams::initPeptidesToDelete ( const StringVector& s )
{
	if ( !s.empty () ) {
		peptidesToDelete.resize ( filenames.size () );
		for ( int i = 0 ; i < s.size () ; i++ ) {
			string str = s[i];
			string::size_type start = 0;
			string::size_type end;
			int fileNum = genNextInt ( str, "\t", start, end );
			string id = genNextString ( str, "\t", start, end );
			string specId = genNextString ( str, "\t", start, end );
			string pep = str.substr ( start );
			peptidesToDelete [fileNum][id].push_back ( PairStringString ( specId, pep ) );
		}
	}
}
VectorPairStringString SearchCompareParams::getRemovePeptides ( int fileNum, const std::string& id ) const
{
	std::map <string, VectorPairStringString>::const_iterator cur = peptidesToDelete [fileNum].find ( id );
	if ( cur != peptidesToDelete [fileNum].end () ) {
		return (*cur).second;
	}
	return VectorPairStringString ();
}
