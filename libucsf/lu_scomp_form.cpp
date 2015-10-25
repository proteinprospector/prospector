/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_scomp_form.cpp                                             *
*                                                                             *
*  Created    : October 29th 2004                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2004-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef MYSQL_DATABASE
#include <lgen_file.h>
#include <lu_param_list.h>
#include <lu_html.h>
#include <lu_quan_ratio.h>
#include <lu_quan_multi.h>
#include <lu_scomp_form.h>
#include <lu_table.h>
#include <lu_getfil.h>
#include <lu_check_db.h>
#include <lu_disc_sc.h>
#include <lu_inst.h>
#include <lu_usermod.h>
#include <lu_sim_ent.h>
using std::ostream;
using std::string;
using std::vector;
using std::endl;
using std::back_inserter;
using std::copy;
using std::find;

const char* FormItemSCFormat::options1 [] = { "HTML", "Tab delimited text", 0 };
const char* FormItemSCFormat::options2 [] = {	"HTML",
												"Tab delimited text",
												"Filtered Peak Lists",
												"MS-Viewer Files",
//												"mzIdentML",
												"pepXML Filtered",
												"pepXML Unfiltered",
												"BiblioSpec",
												"BiblioSpec (Normalized RT)", 0 };
const char* FormItemSCReportType::optionsCal []		= { "Calibration", 0 };
const char* FormItemSCReportType::optionsNormal []	= { "Protein", "Peptide", "Time", "False Positive Rate", 0 };
const char* FormItemSCReportType::optionsXLink []	= { "Protein", "Peptide", "Time", "Crosslinked Peptides", "False Positive Rate", 0 };

FormItemSCReportType::FormItemSCReportType ( bool cal, bool xlink, const string& type ) :
	FormItemSelect ( "Report Type", "", getName (), getOptions ( cal, xlink ), getSelect ( cal, xlink, type ) )
{
}
const char** FormItemSCReportType::getOptions ( bool cal, bool xlink )
{
	if ( cal )	return optionsCal;
	else {
		if ( xlink )	return optionsXLink;
		else			return optionsNormal;
	}
}
string FormItemSCReportType::getSelect ( bool cal, bool xlink, const string& type )
{
	if ( cal )	return optionsCal [0];
	else {
		if ( xlink )	return !type.empty () ? type : optionsXLink [0];
		else			return !type.empty () ? type : optionsNormal [0];
	}
}
static string PP_DIRECTORY				= "pp_directory";
static string SC_DIRECTORY				= "sc_directory";
static string MULTI_SAMPLE				= "multi_sample";
static string REMOVE					= "remove";
static string ACCESSION_NUMS			= "accession_nums";
static string ID_FILTER_LIST			= "id_filter_list";
static string PREFERRED_SPECIES			= "preferred_species";
static string MIN_BEST_DISC_SCORE		= "min_best_disc_score_";
static string MIN_PROTEIN_SCORE			= "min_protein_score";
static string MIN_PEPTIDE_SCORE			= "min_peptide_score";
static string MAX_PROTEIN_EVALUE		= "max_protein_evalue";
static string MAX_PEPTIDE_EVALUE		= "max_peptide_evalue";
static string PEPTIDE_FILTER			= "peptide_filter";
static string BEST_DISC_ONLY			= "best_disc_only";
static string DISC_SCORE_GRAPH			= "disc_score_graph";
static string SLIP_THRESHOLD			= "slip_threshold";
static string MOD_COMP_LIST				= "mod_comp_list";
static string MASS_COMP_LIST			= "mass_comp_list";
static string REPORT_SEARCH_NUMBER		= "report_search_number";
static string REPORT_M_PLUS_H			= "report_m_plus_h";
static string REPORT_M_OVER_Z			= "report_m_over_z";
static string REPORT_CHARGE				= "report_charge";
static string REPORT_M_PLUS_H_CALC		= "report_m_plus_h_calc";
static string REPORT_M_OVER_Z_CALC		= "report_m_over_z_calc";
static string REPORT_INTENSITY			= "report_intensity";
static string REPORT_ERROR				= "report_error";
static string REPORT_UNMATCHED			= "report_unmatched";
static string REPORT_NUM_PKS			= "report_num_pks";
static string REPORT_RANK				= "report_rank";
static string REPORT_SCORE				= "report_score";
static string REPORT_SCORE_DIFF			= "report_score_diff";
static string REPORT_EXPECTATION		= "report_expectation";
static string REPORT_P_VALUE			= "report_p_value";
static string REPORT_NLOG_P_VALUE		= "report_nlog_p_value";
static string REPORT_NUM_PRECURSOR		= "report_num_precursor";
static string REPORT_GRADIENT			= "report_gradient";
static string REPORT_OFFSET				= "report_offset";
static string REPORT_DISC_SCORE			= "report_disc_score";
static string REPORT_REPEATS			= "report_repeats";
static string REPORT_PROT_SCORE			= "report_prot_score";
static string REPORT_NUM_UNIQUE			= "report_num_unique";
static string REPORT_PEPTIDE_COUNT		= "report_peptide_count";
static string REPORT_BEST_SCORE			= "report_best_score";
static string REPORT_BEST_EXPECT		= "report_best_expect";
static string REPORT_COVERAGE			= "report_coverage";
static string REPORT_BEST_DISC_SCORE	= "report_best_disc_score";
static string REPORT_DB_PEPTIDE			= "report_db_peptide";
static string PEPTIDE_MOD_TYPE			= "peptide_mod_type";
static string REPORT_PROTEIN_MOD		= "report_protein_mod";
static string REPORT_MSMS_INFO			= "report_msms_info";
static string REPORT_LENGTH				= "report_length";
static string REPORT_COMPOSITION		= "report_composition";
static string REPORT_START_AA			= "report_start_aa";
static string REPORT_END_AA				= "report_end_aa";
static string REPORT_PREVIOUS_AA		= "report_previous_aa";
static string REPORT_NEXT_AA			= "report_next_aa";
static string REPORT_MISSED_CLEAVAGES	= "report_missed_cleavages";
static string REPORT_MASS_MOD			= "report_mass_mod";
static string REPORT_TIME				= "report_time";
static string REPORT_NUMBER				= "report_number";
static string REPORT_ACCESSION			= "report_accession";
static string REPORT_UNIPROT_ID			= "report_uniprot_id";
static string REPORT_GENE_NAME			= "report_gene_name";
static string REPORT_PROT_LEN			= "report_prot_len";
static string REPORT_MW					= "report_mw";
static string REPORT_PI					= "report_pi";
static string REPORT_SPECIES			= "report_species";
static string REPORT_NAME				= "report_name";
static string REPORT_LINKS				= "report_links";

static string QUAN_TYPE					= "quan_type";

static string REP_INTENSITY				= "rep_intensity";
static string INTENSITY_THRESHOLD		= "intensity_threshold";
static string REP_RESOLUTION			= "rep_resolution";
static string REP_CS_INTENSITY			= "rep_cs_intensity";
static string REP_A_LH_INT				= "rep_a_lh_int";
static string REP_AREA					= "rep_area";
static string REP_CS_AREA				= "rep_cs_area";
static string AREA_THRESHOLD			= "area_threshold";
static string REP_A_LH_AREA				= "rep_a_lh_area";

static string REP_SNR					= "rep_snr";
static string SNR_THRESHOLD				= "snr_threshold";
static string REP_N_MEAN				= "rep_n_mean";
static string REP_N_STDEV				= "rep_n_stdev";
static string RT_INT_START				= "rt_int_start";
static string RT_INT_END				= "rt_int_end";
static string RESOLUTION				= "resolution";

static string PURITY_CORRECTION		= "purity_correction";
static string REPORTER_ION_WINDOW	= "reporter_ion_window";

static string REP_Q_MEDIAN	= "rep_q_median";
static string REP_Q_IQR		= "rep_q_iqr";
static string REP_Q_MEAN	= "rep_q_mean";
static string REP_Q_N_SDV	= "rep_q_n_sdv";
static string REP_Q_STDEV	= "rep_q_stdev";
static string REP_Q_NUM		= "rep_q_num";

static string REPORT_XL_PEPTIDE		= "report_xl_peptide";
static string REPORT_XL_SCORE		= "report_xl_score";
static string REPORT_XL_EXPECTATION	= "report_xl_expectation";
static string REPORT_XL_NLOG_P		= "report_xl_nlog_p";
static string REPORT_XL_RANK		= "report_xl_rank";
static string REPORT_XL_LOW_SCORE	= "report_xl_low_score";
static string REPORT_XL_LOW_EXPECTATION	= "report_xl_low_expectation";
static string REPORT_XL_LOW_NLOG_P		= "report_xl_low_nlog_p";
static string REPORT_XL_AA			= "report_xl_aa";
static string XL_MIN_LOW_SCORE		= "xl_min_low_score";
static string XL_MAX_LOW_EXPECTATION= "xl_max_low_expectation";
static string XL_MIN_SCORE_DIFF		= "xl_min_score_diff";

static string DATA						= "data";
static string REPORT_HOMOLOGOUS_PROTEINS= "report_homologous_proteins";
static string REPORT_HITS_TYPE			= "report_hits_type";
static string MERGE_OPTION				= "merge_option";
static string SORT_TYPE					= "sort_type";
static string SORT_TYPE_2				= "sort_type_2";
static string UNMATCHED_SPECTRA			= "unmatched_spectra";
static string PARENT_MASS_CONVERT		= "parent_mass_convert";

static string RUN_MSPRODUCT = "run_msproduct";

SearchCompareForm::SearchCompareForm ( const VectorConstParameterListPtr& params, const string& reportDefault ) :
	fvj (),
	compIonForm ( params, "ACDEFGHIKLMNPQRSTUVWY" ),
	crosslinkingForm ( params, &fvj, true ),
	calForm ( reportDefault == "Calibration" ),
	numSearches ( params.size () )
{
	expectationFlag = true;
	xLinkFlag = false;
	for ( VectorConstParameterListPtrSizeType i = 0 ; i < params.size () ; i++ ) {
		const ParameterList* p = params [i];
		if ( p->getStringValue ( "expect_calc_method", "None" ) == "None" ) {
			expectationFlag = false;
			break;
		}
	}
	for ( VectorConstParameterListPtrSizeType j = 0 ; j < params.size () ; j++ ) {
		const ParameterList* p = params [j];
		string linkSearchType = p->getStringValue ( "link_search_type" );
		if ( !linkSearchType.empty () && linkSearchType != "No Link"  ) {
			xLinkFlag = true;
			break;
		}
	}
	create ( params );
}
void SearchCompareForm::createItems ()
{
	formItemMap [FormItemVersion::getName ()]	= new FormItemVersion ();
	formItemMap [FormItemSCFormat::getName ()]		= new FormItemSCFormat ( true );

	formItemMap [ACCESSION_NUMS]	= new FormItemTextArea ( "Accession<br />Numbers", "", ACCESSION_NUMS, 3, 20, StringVector () );
	formItemMap [REMOVE]			= new FormItemCheckbox ( "Remove", "", REMOVE, false );

	formItemMap [MULTI_SAMPLE]		= new FormItemCheckbox ( "Multi Sample", "", MULTI_SAMPLE, false );
	formItemMap [ID_FILTER_LIST]	= new FormItemTextArea ( "Spot/Fraction", "", ID_FILTER_LIST, 2, 20, StringVector () );

	formItemMap [PREFERRED_SPECIES] = new FormItemTextArea ( "Preferred<br />Species", "", PREFERRED_SPECIES, 4, 40, StringVector () );

	StringVector instList = DiscriminantScoreInstrumentList::instance ().getNames ();
	for ( StringVectorSizeType i = 0 ; i < instList.size () ; i++ ) {
		string inst = instList [i];
		if ( i == 0 ) initialiseInstrumentName ( inst );
		string in = genTranslateCharacter(inst,'-','_');
		formItemMap [MIN_BEST_DISC_SCORE+in] = new FormItemText ( inst, "", MIN_BEST_DISC_SCORE+in, 4, 6, "0.0", fvj.addSignedFloatingPointValidator ( MIN_BEST_DISC_SCORE+in, "Min Best Discriminant Score " + inst ) );
	}

	formItemMap [MIN_PROTEIN_SCORE]	= new FormItemText ( "Min Score Protein", "", MIN_PROTEIN_SCORE, 4, 10, "22.0", fvj.addSignedFloatingPointValidator ( MIN_PROTEIN_SCORE, "Min Prot Score" ) );
	formItemMap [MIN_PEPTIDE_SCORE]	= new FormItemText ( "Peptide", "", MIN_PEPTIDE_SCORE, 4, 10, "15.0", fvj.addSignedFloatingPointValidator ( MIN_PEPTIDE_SCORE, "Min Pep Score" ) );
	formItemMap [MAX_PROTEIN_EVALUE]= new FormItemText ( "Max E Value Protein", "", MAX_PROTEIN_EVALUE, 4, 10, "0.01", fvj.addPositiveFloatingPointOrExponentValidator ( MAX_PROTEIN_EVALUE, "Max Prot Expect Val" ) );
	formItemMap [MAX_PEPTIDE_EVALUE]= new FormItemText ( "Peptide", "", MAX_PEPTIDE_EVALUE, 4, 10, "0.01", fvj.addPositiveFloatingPointOrExponentValidator ( MAX_PEPTIDE_EVALUE, "Max Expect Val" ) );

	StringVector peptideFilterOptions;
	peptideFilterOptions.push_back ( "Best Peptide Only" );
	peptideFilterOptions.push_back ( "Best Per Charge" );
	peptideFilterOptions.push_back ( "Keep Replicates" );
	string peptideFilterOpt = calForm ? peptideFilterOptions [2] : peptideFilterOptions [0];
	formItemMap [PEPTIDE_FILTER]	= new FormItemSelect ( "", "", PEPTIDE_FILTER, peptideFilterOptions, peptideFilterOpt );
	formItemMap [BEST_DISC_ONLY]	= new FormItemCheckbox ( "Best Discr Only", "", BEST_DISC_ONLY, true );
	formItemMap [DISC_SCORE_GRAPH]	= new FormItemCheckbox ( "Discr Score Graph", "", DISC_SCORE_GRAPH, true );

	formItemMap [FormItemCompMaskType::getName ()]	= new FormItemCompMaskType ( "OR" );
	formItemMap [MOD_COMP_LIST] = new FormItemSelectMultiple ( "", "", "comp_ion", Usermod::getSCompNames (), StringVector (), 4 );
	formItemMap [MASS_COMP_LIST]= new FormItemTextArea ( "Mass<br />Modifications", "", MASS_COMP_LIST, 3, 15, StringVector () );

	formItemMap [REPORT_M_PLUS_H]		= new FormItemCheckbox ( "M+H", "", REPORT_M_PLUS_H, false );
	formItemMap [REPORT_M_OVER_Z]		= new FormItemCheckbox ( "M/Z", "", REPORT_M_OVER_Z, true );
	formItemMap [REPORT_CHARGE]			= new FormItemCheckbox ( "Charge", "", REPORT_CHARGE, true );
	formItemMap [REPORT_M_PLUS_H_CALC]	= new FormItemCheckbox ( "M+H Calc", "", REPORT_M_PLUS_H_CALC, false );
	formItemMap [REPORT_M_OVER_Z_CALC]	= new FormItemCheckbox ( "M/Z Calc", "", REPORT_M_OVER_Z_CALC, false );
	formItemMap [REPORT_INTENSITY]		= new FormItemCheckbox ( "Intensity", "", REPORT_INTENSITY, false );
	formItemMap [REPORT_ERROR]			= new FormItemCheckbox ( "Error", "", REPORT_ERROR, true );
	formItemMap [REPORT_UNMATCHED]		= new FormItemCheckbox ( "Unmatched", "", REPORT_UNMATCHED, false );
	formItemMap [REPORT_NUM_PKS]		= new FormItemCheckbox ( "Num Pks", "", REPORT_NUM_PKS, false );
	formItemMap [REPORT_RANK]			= new FormItemCheckbox ( "Rank", "", REPORT_RANK, false );
	formItemMap [REPORT_SEARCH_NUMBER]	= new FormItemCheckbox ( "Search #", "", REPORT_SEARCH_NUMBER, false );

	formItemMap [REPORT_SCORE]		= new FormItemCheckbox ( "Score", "", REPORT_SCORE, true );
	formItemMap [REPORT_SCORE_DIFF]	= new FormItemCheckbox ( "Score Difference", "", REPORT_SCORE_DIFF, !expectationFlag );
	formItemMap [REPORT_EXPECTATION]= new FormItemCheckbox ( "E Val", "", REPORT_EXPECTATION, expectationFlag );
	formItemMap [REPORT_P_VALUE]	= new FormItemCheckbox ( "P Val", "", REPORT_P_VALUE, false );
	formItemMap [REPORT_NLOG_P_VALUE]= new FormItemCheckbox ( "-10logP", "", REPORT_NLOG_P_VALUE, false );
	formItemMap [REPORT_NUM_PRECURSOR]= new FormItemCheckbox ( "# Precursor", "", REPORT_NUM_PRECURSOR, false );
	formItemMap [REPORT_GRADIENT]	= new FormItemCheckbox ( "Gradient", "", REPORT_GRADIENT, false );
	formItemMap [REPORT_OFFSET]		= new FormItemCheckbox ( "Offset", "", REPORT_OFFSET, false );
	formItemMap [REPORT_DISC_SCORE]	= new FormItemCheckbox ( "Discriminant Score", "", REPORT_DISC_SCORE, !expectationFlag );
	formItemMap [REPORT_REPEATS]	= new FormItemCheckbox ( "# in DB", "", REPORT_REPEATS, true );

	formItemMap [REPORT_PROT_SCORE]		= new FormItemCheckbox ( "Protein Score", "", REPORT_PROT_SCORE, false );
	formItemMap [REPORT_NUM_UNIQUE]		= new FormItemCheckbox ( "Num Unique", "", REPORT_NUM_UNIQUE, true );
	formItemMap [REPORT_PEPTIDE_COUNT]	= new FormItemCheckbox ( "Peptide Count", "", REPORT_PEPTIDE_COUNT, false );
	formItemMap [REPORT_BEST_SCORE]		= new FormItemCheckbox ( "Best Peptide Score", "", REPORT_BEST_SCORE, false );
	formItemMap [REPORT_BEST_EXPECT]	= new FormItemCheckbox ( "Best Expectation Value", "", REPORT_BEST_EXPECT, expectationFlag );
	formItemMap [REPORT_COVERAGE]		= new FormItemCheckbox ( "Coverage", "", REPORT_COVERAGE, true );
	formItemMap [REPORT_BEST_DISC_SCORE]= new FormItemCheckbox ( "Best Discriminant Score", "", REPORT_BEST_DISC_SCORE, true );

	formItemMap [REPORT_DB_PEPTIDE]		= new FormItemCheckbox ( "DB Peptide", "", REPORT_DB_PEPTIDE, true );
	StringVector peptideModTypeOptions;
	peptideModTypeOptions.push_back ( "Off" );
	peptideModTypeOptions.push_back ( "Mods In Peptide" );
	peptideModTypeOptions.push_back ( "Variable Mods Only" );
	//peptideModTypeOptions.push_back ( "Constant Mods Only" );
	peptideModTypeOptions.push_back ( "All Mods (1 column)" );
	peptideModTypeOptions.push_back ( "All Mods (2 columns)" );
	string peptideModTypeOpt = peptideModTypeOptions [1];
	formItemMap [PEPTIDE_MOD_TYPE]		= new FormItemSelect ( "Mod Reporting", "", PEPTIDE_MOD_TYPE, peptideModTypeOptions, peptideModTypeOpt );
	formItemMap [REPORT_PROTEIN_MOD]	= new FormItemCheckbox ( "Protein Mods", "", REPORT_PROTEIN_MOD, false );
	formItemMap [REPORT_MASS_MOD]		= new FormItemCheckbox ( "Mass Mods", "", REPORT_MASS_MOD, false );
	formItemMap [SLIP_THRESHOLD]		= new FormItemText ( "SLIP Threshold", "", SLIP_THRESHOLD, 3, 3, "0", fvj.addPositiveFloatingPointValidator ( SLIP_THRESHOLD, "SLIP Threshold" ) );

	formItemMap [REPORT_MISSED_CLEAVAGES]= new FormItemCheckbox ( "Missed Cleavages", "", REPORT_MISSED_CLEAVAGES, false );
	formItemMap [REPORT_TIME]			= new FormItemCheckbox ( "Time", "", REPORT_TIME, true );
	formItemMap [REPORT_MSMS_INFO]		= new FormItemCheckbox ( "MSMS Info", "", REPORT_MSMS_INFO, false );
	formItemMap [REPORT_LENGTH]			= new FormItemCheckbox ( "Length", "", REPORT_LENGTH, false );
	formItemMap [REPORT_COMPOSITION]	= new FormItemCheckbox ( "Composition", "", REPORT_COMPOSITION, false );
	formItemMap [REPORT_START_AA]		= new FormItemCheckbox ( "Start AA", "", REPORT_START_AA, false );
	formItemMap [REPORT_END_AA]			= new FormItemCheckbox ( "End AA", "", REPORT_END_AA, false );
	formItemMap [REPORT_PREVIOUS_AA]	= new FormItemText ( "Previous AA", "", REPORT_PREVIOUS_AA, 1, 1, "0", fvj.addPositiveIntegerValidator ( REPORT_PREVIOUS_AA, "Previous AA" ) );
	formItemMap [REPORT_NEXT_AA]		= new FormItemText ( "Next AA", "", REPORT_NEXT_AA, 1, 1, "0", fvj.addPositiveIntegerValidator ( REPORT_NEXT_AA, "Next AA" ) );

	formItemMap [REPORT_NUMBER]		= new FormItemCheckbox ( "Number", "", REPORT_NUMBER, true );
	formItemMap [REPORT_ACCESSION]	= new FormItemCheckbox ( "Accession", "", REPORT_ACCESSION, true );
	formItemMap [REPORT_UNIPROT_ID]	= new FormItemCheckbox ( "Uniprot ID", "", REPORT_UNIPROT_ID, false );
	formItemMap [REPORT_GENE_NAME]	= new FormItemCheckbox ( "Gene Name", "", REPORT_GENE_NAME, false );
	formItemMap [REPORT_PROT_LEN]	= new FormItemCheckbox ( "Protein Length", "", REPORT_PROT_LEN, false );
	formItemMap [REPORT_MW]			= new FormItemCheckbox ( "MW", "", REPORT_MW, true );
	formItemMap [REPORT_PI]			= new FormItemCheckbox ( "pI", "", REPORT_PI, false );
	formItemMap [REPORT_SPECIES]	= new FormItemCheckbox ( "Species", "", REPORT_SPECIES, true );
	formItemMap [REPORT_NAME]		= new FormItemCheckbox ( "Name", "", REPORT_NAME, true );
	formItemMap [REPORT_LINKS]		= new FormItemCheckbox ( "Links", "", REPORT_LINKS, true );
	formItemMap [FormItemSCCheckboxes::getName ()] = new FormItemSCCheckboxes ();

	formItemMap [FormItemRawType::getName ()] = new FormItemRawType ();
	StringVector quanTypeOptions = QuantitationInfo::instance ().getNames ();
	StringVector quanMSMSNames = getQuanMSMSNames ();
	copy ( quanMSMSNames.begin (), quanMSMSNames.end (), back_inserter ( quanTypeOptions ) );
	formItemMap [QUAN_TYPE] = new FormItemSelect ( "Quantitation", "", QUAN_TYPE, quanTypeOptions, quanTypeOptions [0] );

	formItemMap [REP_INTENSITY]	= new FormItemCheckbox ( "Intensity", "", REP_INTENSITY, false );
	formItemMap [INTENSITY_THRESHOLD]= new FormItemText ( "Threshold", "", INTENSITY_THRESHOLD, 5, 10, "0", fvj.addPositiveFloatingPointOrExponentValidator ( INTENSITY_THRESHOLD, "Intensity Threshold" ) );
	formItemMap [REP_RESOLUTION]= new FormItemCheckbox ( "Resolution", "", REP_RESOLUTION, false );
	formItemMap [REP_CS_INTENSITY]= new FormItemCheckbox ( "CS", "", REP_CS_INTENSITY, false );
	formItemMap [REP_A_LH_INT]	= new FormItemCheckbox ( "L/H Int", "", REP_A_LH_INT, false );
	formItemMap [REP_AREA]		= new FormItemCheckbox ( "Area", "", REP_AREA, false );
	formItemMap [REP_CS_AREA]	= new FormItemCheckbox ( "CS", "", REP_CS_AREA, false );
	formItemMap [AREA_THRESHOLD]= new FormItemText ( "Threshold", "", AREA_THRESHOLD, 5, 10, "0", fvj.addPositiveFloatingPointOrExponentValidator ( AREA_THRESHOLD, "Area Threshold" ) );
	formItemMap [REP_A_LH_AREA]	= new FormItemCheckbox ( "L/H Area", "", REP_A_LH_AREA, false );

	formItemMap [REP_SNR]		= new FormItemCheckbox ( "SNR", "", REP_SNR, false );
	formItemMap [SNR_THRESHOLD]	= new FormItemText ( "Threshold", "", SNR_THRESHOLD, 5, 10, "10.0", fvj.addPositiveFloatingPointValidator ( SNR_THRESHOLD, "SNR Threshold" ) );
	formItemMap [REP_N_MEAN]	= new FormItemCheckbox ( "Noise Mean", "", REP_N_MEAN, false );
	formItemMap [REP_N_STDEV]	= new FormItemCheckbox ( "Noise SD", "", REP_N_STDEV, false );
	formItemMap [RT_INT_START]	= new FormItemText ( "RT Int (sec)", "", RT_INT_START, 5, 10, "0.0", fvj.addSignedFloatingPointValidator ( RT_INT_START, "RT Int Start" ) );
	formItemMap [RT_INT_END]	= new FormItemText ( "to", "", RT_INT_END, 5, 10, "0.0", fvj.addSignedFloatingPointValidator ( RT_INT_END, "RT Int End" ) );
	formItemMap [RESOLUTION]	= new FormItemText ( "Resolution", "", RESOLUTION, 8, 10, "10000.0", fvj.addPositiveFloatingPointOrExponentValidator ( RESOLUTION, "Resolution" ) );

	formItemMap [REP_Q_MEDIAN]	= new FormItemCheckbox ( "Median", "", REP_Q_MEDIAN, false );
	formItemMap [REP_Q_IQR]		= new FormItemCheckbox ( "IQR", "", REP_Q_IQR, false );
	formItemMap [REP_Q_MEAN]	= new FormItemCheckbox ( "Mean", "", REP_Q_MEAN, false );
	formItemMap [REP_Q_N_SDV]	= new FormItemText ( "", "", REP_Q_N_SDV, 3, 3, "2.0", fvj.addPositiveFloatingPointValidator ( REP_Q_N_SDV, "Number of Std Dev" ) );
	formItemMap [REP_Q_STDEV]	= new FormItemCheckbox ( "Std Dev", "", REP_Q_STDEV, false );
	formItemMap [REP_Q_NUM]		= new FormItemCheckbox ( "Num", "", REP_Q_NUM, false );

	if ( xLinkFlag ) {
		formItemMap [REPORT_XL_PEPTIDE]		= new FormItemCheckbox ( "Crosslinked Peptide", "", REPORT_XL_PEPTIDE, true );

		formItemMap [REPORT_XL_SCORE]		= new FormItemCheckbox ( "MS-Tag Score", "", REPORT_XL_SCORE, true );
		formItemMap [REPORT_XL_EXPECTATION]	= new FormItemCheckbox ( "EVal", "", REPORT_XL_EXPECTATION, expectationFlag );
		formItemMap [REPORT_XL_NLOG_P]		= new FormItemCheckbox ( "-10logP", "", REPORT_XL_NLOG_P, expectationFlag );
		formItemMap [REPORT_XL_RANK]		= new FormItemCheckbox ( "Rank", "", REPORT_XL_RANK, true );

		formItemMap [REPORT_XL_LOW_SCORE]		= new FormItemCheckbox ( "MS-Tag Low Score", "", REPORT_XL_LOW_SCORE, true );
		formItemMap [REPORT_XL_LOW_EXPECTATION]	= new FormItemCheckbox ( "EVal", "", REPORT_XL_LOW_EXPECTATION, expectationFlag );
		formItemMap [REPORT_XL_LOW_NLOG_P]		= new FormItemCheckbox ( "-10logP", "", REPORT_XL_LOW_NLOG_P, expectationFlag );

		formItemMap [REPORT_XL_AA]			= new FormItemCheckbox ( "XL AA", "", REPORT_XL_AA, true );

		formItemMap [XL_MIN_LOW_SCORE]		= new FormItemText ( "Min MS-Tag Low Score", "", XL_MIN_LOW_SCORE, 4, 10, "10.0", fvj.addSignedFloatingPointValidator ( XL_MIN_LOW_SCORE, "Min MS-Tag Low Score" ) );
		formItemMap [XL_MAX_LOW_EXPECTATION]= new FormItemText ( "Max MS-Tag Low Expectation", "", XL_MAX_LOW_EXPECTATION, 4, 10, "1", fvj.addSignedFloatingPointValidator ( XL_MAX_LOW_EXPECTATION, "Max MS-Tag Low Expectation" ) );
		formItemMap [XL_MIN_SCORE_DIFF]		= new FormItemText ( "Min Score Diff", "", XL_MIN_SCORE_DIFF, 4, 10, "0.0", fvj.addSignedFloatingPointValidator ( XL_MIN_SCORE_DIFF, "Min Score Diff" ) );
	}

	formItemMap [FormItemIsotopePurity::getName ( "C", 13 )] = new FormItemIsotopePurity ( "C", 13, &fvj );
	formItemMap [FormItemIsotopePurity::getName ( "N", 15 )] = new FormItemIsotopePurity ( "N", 15, &fvj );
	formItemMap [FormItemIsotopePurity::getName ( "O", 18 )] = new FormItemIsotopePurity ( "O", 18, &fvj );
	StringVector purityOptions = getPurityNames ();
	formItemMap [PURITY_CORRECTION] = new FormItemSelect ( "Purity Corr", "", PURITY_CORRECTION, purityOptions, purityOptions [0] );
	formItemMap [REPORTER_ION_WINDOW] = new FormItemReporterIonWindow ( &fvj, "0.4" );

	formItemMap [DATA]	= new FormItemTextArea ( "", "", DATA, 6, 80, StringVector () );

	formItemMap [FormItemSCReportType::getName ()] = new FormItemSCReportType ( calForm, xLinkFlag );

	StringVector reportHomologousProteinsOptions;
	reportHomologousProteinsOptions.push_back ( "All" );
	reportHomologousProteinsOptions.push_back ( "Interesting" );
	reportHomologousProteinsOptions.push_back ( "None" );
	formItemMap [REPORT_HOMOLOGOUS_PROTEINS] = new FormItemSelect ( "Report Homologous Proteins", "", REPORT_HOMOLOGOUS_PROTEINS, reportHomologousProteinsOptions, "Interesting" );

	StringVector reportHitsTypeOptions;
	reportHitsTypeOptions.push_back ( "Union" );
	reportHitsTypeOptions.push_back ( "Protein Intersection" );
	reportHitsTypeOptions.push_back ( "Peptide Intersection" );
	reportHitsTypeOptions.push_back ( "Protein Difference" );
	reportHitsTypeOptions.push_back ( "Peptide Difference" );
	formItemMap [REPORT_HITS_TYPE] = new FormItemSelect ( "Report Hits Type", "", REPORT_HITS_TYPE, reportHitsTypeOptions, "Union" );
	StringVector mergeOptions;
	mergeOptions.push_back ( "Separated" );
	mergeOptions.push_back ( "Merged" );
	formItemMap [MERGE_OPTION] = new FormItemSelect ( "", "", MERGE_OPTION, mergeOptions, "Separated" );
	StringVector sortTypeOptions;
	sortTypeOptions.push_back ( "Discriminant Score" );
	sortTypeOptions.push_back ( "Peptide Score" );
	if ( expectationFlag ) sortTypeOptions.push_back ( "Expectation Value" );
	sortTypeOptions.push_back ( "Start Residue" );
	sortTypeOptions.push_back ( "End Residue" );
	sortTypeOptions.push_back ( "Crosslink AA" );
	sortTypeOptions.push_back ( "Fraction/RT" );
	sortTypeOptions.push_back ( "RT" );
	sortTypeOptions.push_back ( "m/z" );
	sortTypeOptions.push_back ( "M+H" );
	sortTypeOptions.push_back ( "Intensity" );
	sortTypeOptions.push_back ( "Error" );
	sortTypeOptions.push_back ( "Time" );
	sortTypeOptions.push_back ( "Charge/M+H" );
	sortTypeOptions.push_back ( "Mass Mod" );
	string sortDefault = expectationFlag ? "Expectation Value" : "Discriminant Score";
	StringVector sortTypeOptions2;
	sortTypeOptions2.push_back ( "" );
	copy ( sortTypeOptions.begin (), sortTypeOptions.end (), back_inserter ( sortTypeOptions2 ) );
	string sortDefault2 = "";
	formItemMap [SORT_TYPE]		= new FormItemSelect ( "Sort Type", "", SORT_TYPE, sortTypeOptions, sortDefault );
	formItemMap [SORT_TYPE_2]	= new FormItemSelect ( "", "", SORT_TYPE_2, sortTypeOptions2, sortDefault2 );
	formItemMap [UNMATCHED_SPECTRA] = new FormItemCheckbox ( "Unmatched Spectra", "", UNMATCHED_SPECTRA, false );

	formItemMap [FormItemSaveParameters::getName ()] = new FormItemSaveParameters ();

	StringVector parentMassConvertOptions;
	parentMassConvertOptions.push_back ( "monoisotopic" );
	formItemMap [PARENT_MASS_CONVERT]		= new FormItemSelect ( "", "", PARENT_MASS_CONVERT, parentMassConvertOptions, "monoisotopic" );
	formItemMap [FormItemMSMSPkFilter::getName ()] = new FormItemMSMSPkFilter ( true );
	formItemMap [FormItemMaxMSMSPeaks::getName ()] = new FormItemMaxMSMSPeaks ( &fvj, true );
	formItemMap [FormItemMaxReportedHits::getName("msms_")]	= new FormItemMaxReportedHits ("msms_", "", &fvj, true);

	formItemMap [RUN_MSPRODUCT] = new FormItemCheckbox ( "", "", RUN_MSPRODUCT, false );
}
void SearchCompareForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];
		formItemMap [MULTI_SAMPLE]->setValue ( p );

		StringVector data_value;
		for ( VectorConstParameterListPtrSizeType j = 0 ; j < params.size () ; j++ ) {
			data_value.push_back ( params [j]->getStringValue ( "search_key" ) );
			instrumentList.push_back ( params [j]->getStringValue ( "instrument_name" ) );
		}
		formItemMap [DATA]->setValue ( data_value );

		formItemMap [PARENT_MASS_CONVERT]->setValue ( p );

// Potential default parameters
		formItemMap [FormItemSCFormat::getName ()]->setValue ( p );
		formItemMap [REMOVE]->setValue ( p );
		formItemMap [ID_FILTER_LIST]->setValue ( p );
		formItemMap [PREFERRED_SPECIES]->setValue ( p );
		StringVector instList = DiscriminantScoreInstrumentList::instance ().getNames ();
		for ( StringVectorSizeType k = 0 ; k < instList.size () ; k++ ) {
			formItemMap [MIN_BEST_DISC_SCORE+genTranslateCharacter(instList [k],'-','_')]->setValue ( p );
		}
		formItemMap [MIN_PROTEIN_SCORE]->setValue ( p );
		formItemMap [MIN_PEPTIDE_SCORE]->setValue ( p );
		if ( expectationFlag ) {
			formItemMap [MAX_PROTEIN_EVALUE]->setValue ( p );
			formItemMap [MAX_PEPTIDE_EVALUE]->setValue ( p );
		}
		formItemMap [PEPTIDE_FILTER]->setValue ( p );
		formItemMap [BEST_DISC_ONLY]->setValue ( p );
		formItemMap [DISC_SCORE_GRAPH]->setValue ( p );
		formItemMap [FormItemCompMaskType::getName ()]->setValue ( p );
		formItemMap [MOD_COMP_LIST]->setValue ( p );
		formItemMap [MASS_COMP_LIST]->setValue ( p );
		formItemMap [REPORT_M_PLUS_H]->setValue ( p );
		formItemMap [REPORT_M_OVER_Z]->setValue ( p );
		formItemMap [REPORT_CHARGE]->setValue ( p );
		formItemMap [REPORT_M_PLUS_H_CALC]->setValue ( p );
		formItemMap [REPORT_M_OVER_Z_CALC]->setValue ( p );
		formItemMap [REPORT_INTENSITY]->setValue ( p );
		formItemMap [REPORT_ERROR]->setValue ( p );
		formItemMap [REPORT_UNMATCHED]->setValue ( p );
		formItemMap [REPORT_NUM_PKS]->setValue ( p );
		formItemMap [REPORT_RANK]->setValue ( p );
		formItemMap [REPORT_SEARCH_NUMBER]->setValue ( p );
		formItemMap [REPORT_SCORE]->setValue ( p );
		formItemMap [REPORT_SCORE_DIFF]->setValue ( p );
		if ( expectationFlag ) {
			formItemMap [REPORT_EXPECTATION]->setValue ( p );
			formItemMap [REPORT_P_VALUE]->setValue ( p );
			formItemMap [REPORT_NLOG_P_VALUE]->setValue ( p );
		}
		formItemMap [REPORT_NUM_PRECURSOR]->setValue ( p );
		if ( expectationFlag ) {
			formItemMap [REPORT_GRADIENT]->setValue ( p );
			formItemMap [REPORT_OFFSET]->setValue ( p );
		}
		formItemMap [REPORT_DISC_SCORE]->setValue ( p );
		formItemMap [REPORT_REPEATS]->setValue ( p );
		formItemMap [REPORT_PROT_SCORE]->setValue ( p );
		formItemMap [REPORT_NUM_UNIQUE]->setValue ( p );
		formItemMap [REPORT_PEPTIDE_COUNT]->setValue ( p );
		formItemMap [REPORT_BEST_SCORE]->setValue ( p );
		if ( expectationFlag ) {
			formItemMap [REPORT_BEST_EXPECT]->setValue ( p );
		}
		formItemMap [REPORT_COVERAGE]->setValue ( p );
		formItemMap [REPORT_BEST_DISC_SCORE]->setValue ( p );
		formItemMap [REPORT_DB_PEPTIDE]->setValue ( p );
		formItemMap [PEPTIDE_MOD_TYPE]->setValue ( p );
		formItemMap [REPORT_PROTEIN_MOD]->setValue ( p );
		formItemMap [REPORT_MASS_MOD]->setValue ( p );
		formItemMap [SLIP_THRESHOLD]->setValue ( p );

		formItemMap [REPORT_MISSED_CLEAVAGES]->setValue ( p );
		formItemMap [REPORT_TIME]->setValue ( p );
		formItemMap [REPORT_MSMS_INFO]->setValue ( p );
		formItemMap [REPORT_LENGTH]->setValue ( p );
		formItemMap [REPORT_COMPOSITION]->setValue ( p );
		formItemMap [REPORT_START_AA]->setValue ( p );
		formItemMap [REPORT_END_AA]->setValue ( p );
		formItemMap [REPORT_PREVIOUS_AA]->setValue ( p );
		formItemMap [REPORT_NEXT_AA]->setValue ( p );

		formItemMap [REPORT_NUMBER]->setValue ( p );
		formItemMap [REPORT_ACCESSION]->setValue ( p );
		formItemMap [REPORT_UNIPROT_ID]->setValue ( p );
		formItemMap [REPORT_GENE_NAME]->setValue ( p );
		formItemMap [REPORT_PROT_LEN]->setValue ( p );
		formItemMap [REPORT_MW]->setValue ( p );
		formItemMap [REPORT_PI]->setValue ( p );
		formItemMap [REPORT_SPECIES]->setValue ( p );
		formItemMap [REPORT_NAME]->setValue ( p );
		formItemMap [REPORT_LINKS]->setValue ( p );
		formItemMap [FormItemSCCheckboxes::getName ()]->setValue ( p );
		formItemMap [FormItemRawType::getName ()]->setValue ( p );
		formItemMap [QUAN_TYPE]->setValue ( p );
		formItemMap [REP_INTENSITY]->setValue ( p );
		formItemMap [INTENSITY_THRESHOLD]->setValue ( p );
		formItemMap [REP_RESOLUTION]->setValue ( p );
		formItemMap [REP_CS_INTENSITY]->setValue ( p );
		formItemMap [REP_A_LH_INT]->setValue ( p );
		formItemMap [REP_AREA]->setValue ( p );
		formItemMap [REP_CS_AREA]->setValue ( p );
		formItemMap [AREA_THRESHOLD]->setValue ( p );
		formItemMap [REP_A_LH_AREA]->setValue ( p );
		formItemMap [REP_SNR]->setValue ( p );
		formItemMap [SNR_THRESHOLD]->setValue ( p );
		formItemMap [REP_N_MEAN]->setValue ( p );
		formItemMap [REP_N_STDEV]->setValue ( p );
		formItemMap [RT_INT_START]->setValue ( p );
		formItemMap [RT_INT_END]->setValue ( p );
		formItemMap [RESOLUTION]->setValue ( p );
		formItemMap [REP_Q_MEDIAN]->setValue ( p );
		formItemMap [REP_Q_IQR]->setValue ( p );
		formItemMap [REP_Q_MEAN]->setValue ( p );
		formItemMap [REP_Q_N_SDV]->setValue ( p );
		formItemMap [REP_Q_STDEV]->setValue ( p );
		formItemMap [REP_Q_NUM]->setValue ( p );
		formItemMap [FormItemIsotopePurity::getName ( "C", 13 )]->setValue ( p );
		formItemMap [FormItemIsotopePurity::getName ( "N", 15 )]->setValue ( p );
		formItemMap [FormItemIsotopePurity::getName ( "O", 18 )]->setValue ( p );
		formItemMap [PURITY_CORRECTION]->setValue ( p );
		formItemMap [REPORTER_ION_WINDOW]->setValue ( p );
		if ( !calForm ) {
			formItemMap [FormItemSCReportType::getName ()]->setValue ( p );
		}
		if ( xLinkFlag ) {
			formItemMap [REPORT_XL_PEPTIDE]->setValue ( p );
			formItemMap [REPORT_XL_SCORE]->setValue ( p );
			if ( expectationFlag ) {
				formItemMap [REPORT_XL_EXPECTATION]->setValue ( p );
				formItemMap [REPORT_XL_NLOG_P]->setValue ( p );
			}
			formItemMap [REPORT_XL_RANK]->setValue ( p );
			formItemMap [REPORT_XL_LOW_SCORE]->setValue ( p );
			if ( expectationFlag ) {
				formItemMap [REPORT_XL_LOW_EXPECTATION]->setValue ( p );
				formItemMap [REPORT_XL_LOW_NLOG_P]->setValue ( p );
			}
			formItemMap [REPORT_XL_AA]->setValue ( p );
			formItemMap [XL_MIN_LOW_SCORE]->setValue ( p );
			if ( expectationFlag ) {
				formItemMap [XL_MAX_LOW_EXPECTATION]->setValue ( p );
			}
			formItemMap [XL_MIN_SCORE_DIFF]->setValue ( p );
		}
		formItemMap [REPORT_HOMOLOGOUS_PROTEINS]->setValue ( p );
		formItemMap [REPORT_HITS_TYPE]->setValue ( p );
		formItemMap [MERGE_OPTION]->setValue ( p );
		if ( expectationFlag ) {
			formItemMap [SORT_TYPE]->setValue ( p );
			formItemMap [SORT_TYPE_2]->setValue ( p );
		}
		else {
			if ( p->getStringValue ( SORT_TYPE ) != "Expectation Value" )	formItemMap [SORT_TYPE]->setValue ( p );
			if ( p->getStringValue ( SORT_TYPE_2 ) != "Expectation Value" )	formItemMap [SORT_TYPE_2]->setValue ( p );
		}
		formItemMap [UNMATCHED_SPECTRA]->setValue ( p );
		formItemMap [FormItemMSMSPkFilter::getName ()]->setValue ( p );
		formItemMap [FormItemMaxMSMSPeaks::getName ()]->setValue ( p );
		formItemMap [FormItemMaxReportedHits::getName("msms_")]->setValue ( p );
		formItemMap [RUN_MSPRODUCT]->setValue ( p );
	}
}
void SearchCompareForm::printRawDataValidationJavascript ( ostream& os ) const
{
	startJavascript ( os );
		os << "function validateCheckRawData( form ) {" << endl;
		os << "\t" << "var rd = false;" << endl;
		os << "\t" << "if(form.rep_q_median.checked){" << endl;
		os << "\t\t" << "rd = true;" << endl;
		os << "\t" << "}" << endl;
		os << "\t" << "if(form.rep_q_iqr.checked){" << endl;
		os << "\t\t" << "rd = true;" << endl;
		os << "\t" << "}" << endl;
		os << "\t" << "if(form.rep_q_mean.checked){" << endl;
		os << "\t\t" << "rd = true;" << endl;
		os << "\t" << "}" << endl;
		os << "\t" << "if(form.rep_q_stdev.checked){" << endl;
		os << "\t\t" << "rd = true;" << endl;
		os << "\t" << "}" << endl;
		os << "\t" << "if(form.rep_q_num.checked){" << endl;
		os << "\t\t" << "rd = true;" << endl;
		os << "\t" << "}" << endl;
		os << "\t" << "if(form.rep_intensity.checked){" << endl;
		os << "\t\t" << "rd = true;" << endl;
		os << "\t" << "}" << endl;
		os << "\t" << "if(form.rep_resolution.checked){" << endl;
		os << "\t\t" << "rd = true;" << endl;
		os << "\t" << "}" << endl;
		os << "\t" << "if(form.rep_cs_intensity.checked){" << endl;
		os << "\t\t" << "rd = true;" << endl;
		os << "\t" << "}" << endl;
		os << "\t" << "if(form.rep_a_lh_int.checked){" << endl;
		os << "\t\t" << "rd = true;" << endl;
		os << "\t" << "}" << endl;
		os << "\t" << "if(form.rep_area.checked){" << endl;
		os << "\t\t" << "rd = true;" << endl;
		os << "\t" << "}" << endl;
		os << "\t" << "if(form.rep_cs_area.checked){" << endl;
		os << "\t\t" << "rd = true;" << endl;
		os << "\t" << "}" << endl;
		os << "\t" << "if(form.rep_a_lh_area.checked){" << endl;
		os << "\t\t" << "rd = true;" << endl;
		os << "\t" << "}" << endl;
		os << "\t" << "if(form.rep_snr.checked){" << endl;
		os << "\t\t" << "rd = true;" << endl;
		os << "\t" << "}" << endl;
		os << "\t" << "if(form.rep_n_mean.checked){" << endl;
		os << "\t\t" << "rd = true;" << endl;
		os << "\t" << "}" << endl;
		os << "\t" << "if(form.rep_n_stdev.checked){" << endl;
		os << "\t\t" << "rd = true;" << endl;
		os << "\t" << "}" << endl;
		os << "\t" << "if(rd){" << endl;
		os << "\t\t" << "form.action = \"" << ProgramLink::getURLStart ( "searchCompareRawData" ) << "\";" << endl;
		os << "\t" << "}" << endl;
		os << "\t" << "else{" << endl;
		os << "\t\t" << "form.action = \"" << ProgramLink::getURLStart ( "searchCompare" ) << "\";" << endl;
		os << "\t" << "}" << endl;
		os << "\t" << "validateForms(form);" << endl;
		os << "}" << endl;
	endJavascript ( os );
}
void SearchCompareForm::printHTML ( ostream& os )
{
	init_html ( os, "Search Compare", "batchtagman.htm#search_compare" );
	fvj.print ( os );
	if ( InfoParams::instance ().getBoolValue ( "raw_data_forwarding" ) ) {
		printRawDataValidationJavascript ( os );
		os << "<form";
		os << " ";
		os << "method=\"post\"";
		os << " ";
		os << "onsubmit=\"return validateCheckRawData(this)\"";
		os << ">" << endl;
	}
	else printHTMLFORMStart ( os, "post", "searchCompare", false, true );
	tableStart ( os, true );

		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true );
				formItemMap [FormItemSCFormat::getName ()]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );

		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true );
				tableStart ( os, false );
					tableRowStart ( os );
						tableHeaderStart ( os, "", "left", true );
							formItemMap [ACCESSION_NUMS]->printHTML ( os );
							formItemMap [REMOVE]->printHTML ( os );
						tableHeaderEnd ( os );
						tableHeaderStart ( os, "", "left", true );
							formItemMap [MULTI_SAMPLE]->printHTML ( os );
							os << "<br />" << endl;
							formItemMap [ID_FILTER_LIST]->printHTML ( os );
						tableHeaderEnd ( os );
						tableHeaderStart ( os, "", "left", true );
							formItemMap [PREFERRED_SPECIES]->printHTML ( os );
						tableHeaderEnd ( os );
					tableRowEnd ( os );
				tableEnd ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );

		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true );
				printHTMLDiscriminentItems ( os );
				formItemMap [MIN_PROTEIN_SCORE]->printHTML ( os );
				formItemMap [MIN_PEPTIDE_SCORE]->printHTML ( os );
				if ( expectationFlag ) {
					formItemMap [MAX_PROTEIN_EVALUE]->printHTML ( os );
					formItemMap [MAX_PEPTIDE_EVALUE]->printHTML ( os );
				}
				formItemMap [PEPTIDE_FILTER]->printHTML ( os );
				formItemMap [BEST_DISC_ONLY]->printHTML ( os );
				formItemMap [DISC_SCORE_GRAPH]->printHTML ( os );
				os << "<br />" << endl;
				ExpandableJavascriptBlock ejb2 ( "Peptide Composition" );
				ejb2.printHeader ( os );
				tableStart ( os, true );
					tableRowStart ( os );
						tableHeaderStart ( os );
							formItemMap [FormItemCompMaskType::getName ()]->printHTML ( os );
						tableHeaderEnd ( os );
						tableHeaderStart ( os );
							compIonForm.printHTML2Rows ( os );
						tableHeaderEnd ( os );
						tableHeaderStart ( os );
							formItemMap [MOD_COMP_LIST]->printHTML ( os );
						tableHeaderEnd ( os );
						tableHeaderStart ( os );
							formItemMap [MASS_COMP_LIST]->printHTML ( os );
						tableHeaderEnd ( os );
					tableRowEnd ( os );
				tableEnd ( os );
				ejb2.printFooter ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );

		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true );
				formItemMap [REPORT_M_PLUS_H]->printHTML ( os );
				formItemMap [REPORT_M_OVER_Z]->printHTML ( os );
				formItemMap [REPORT_CHARGE]->printHTML ( os );
				formItemMap [REPORT_M_PLUS_H_CALC]->printHTML ( os );
				formItemMap [REPORT_M_OVER_Z_CALC]->printHTML ( os );
				formItemMap [REPORT_INTENSITY]->printHTML ( os );
				formItemMap [REPORT_ERROR]->printHTML ( os );
				formItemMap [REPORT_UNMATCHED]->printHTML ( os );
				formItemMap [REPORT_NUM_PKS]->printHTML ( os );
				formItemMap [REPORT_RANK]->printHTML ( os );
				formItemMap [REPORT_SEARCH_NUMBER]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [REPORT_SCORE]->printHTML ( os );
				formItemMap [REPORT_SCORE_DIFF]->printHTML ( os );
				if ( expectationFlag ) {
					formItemMap [REPORT_EXPECTATION]->printHTML ( os );
					formItemMap [REPORT_P_VALUE]->printHTML ( os );
					formItemMap [REPORT_NLOG_P_VALUE]->printHTML ( os );
				}
				formItemMap [REPORT_NUM_PRECURSOR]->printHTML ( os );
				if ( expectationFlag ) {
					formItemMap [REPORT_GRADIENT]->printHTML ( os );
					formItemMap [REPORT_OFFSET]->printHTML ( os );
				}
				formItemMap [REPORT_DISC_SCORE]->printHTML ( os );
				formItemMap [REPORT_REPEATS]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [REPORT_PROT_SCORE]->printHTML ( os );
				formItemMap [REPORT_NUM_UNIQUE]->printHTML ( os );
				formItemMap [REPORT_PEPTIDE_COUNT]->printHTML ( os );
				formItemMap [REPORT_BEST_SCORE]->printHTML ( os );
				if ( expectationFlag ) formItemMap [REPORT_BEST_EXPECT]->printHTML ( os );
				formItemMap [REPORT_COVERAGE]->printHTML ( os );
				formItemMap [REPORT_BEST_DISC_SCORE]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [REPORT_DB_PEPTIDE]->printHTML ( os );
				formItemMap [PEPTIDE_MOD_TYPE]->printHTML ( os );
				formItemMap [REPORT_PROTEIN_MOD]->printHTML ( os );
				formItemMap [SLIP_THRESHOLD]->printHTML ( os );
				formItemMap [REPORT_MASS_MOD]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [REPORT_MISSED_CLEAVAGES]->printHTML ( os );
				formItemMap [REPORT_TIME]->printHTML ( os );
				formItemMap [REPORT_MSMS_INFO]->printHTML ( os );
				formItemMap [REPORT_LENGTH]->printHTML ( os );
				formItemMap [REPORT_COMPOSITION]->printHTML ( os );
				formItemMap [REPORT_START_AA]->printHTML ( os );
				formItemMap [REPORT_END_AA]->printHTML ( os );
				formItemMap [REPORT_PREVIOUS_AA]->printHTML ( os );
				formItemMap [REPORT_NEXT_AA]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [REPORT_NUMBER]->printHTML ( os );
				formItemMap [REPORT_ACCESSION]->printHTML ( os );
				formItemMap [REPORT_UNIPROT_ID]->printHTML ( os );
				formItemMap [REPORT_GENE_NAME]->printHTML ( os );
				formItemMap [REPORT_PROT_LEN]->printHTML ( os );
				formItemMap [REPORT_MW]->printHTML ( os );
				formItemMap [REPORT_PI]->printHTML ( os );
				formItemMap [REPORT_SPECIES]->printHTML ( os );
				formItemMap [REPORT_NAME]->printHTML ( os );
				formItemMap [REPORT_LINKS]->printHTML ( os );
				if ( numSearches == 1 ) {
					formItemMap [FormItemSCCheckboxes::getName ()]->printHTML ( os );
				}
			tableHeaderEnd ( os );
		tableRowEnd ( os );

		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true );
				ExpandableJavascriptBlock ejb3 ( "Raw Data/Quantitation" );
				ejb3.printHeader ( os );
				formItemMap [FormItemRawType::getName ()]->printHTML ( os );
				formItemMap [QUAN_TYPE]->printHTML ( os );
				formItemMap [REP_Q_MEDIAN]->printHTML ( os );
				formItemMap [REP_Q_IQR]->printHTML ( os );
				formItemMap [REP_Q_MEAN]->printHTML ( os );
				formItemMap [REP_Q_N_SDV]->printHTML ( os );
				formItemMap [REP_Q_STDEV]->printHTML ( os );
				formItemMap [REP_Q_NUM]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [REP_INTENSITY]->printHTML ( os );
				formItemMap [INTENSITY_THRESHOLD]->printHTML ( os );
				formItemMap [REP_RESOLUTION]->printHTML ( os );
				formItemMap [REP_CS_INTENSITY]->printHTML ( os );
				formItemMap [REP_A_LH_INT]->printHTML ( os );
				formItemMap [REP_AREA]->printHTML ( os );
				formItemMap [REP_CS_AREA]->printHTML ( os );
				formItemMap [AREA_THRESHOLD]->printHTML ( os );
				formItemMap [REP_A_LH_AREA]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [REP_SNR]->printHTML ( os );
				formItemMap [SNR_THRESHOLD]->printHTML ( os );
				formItemMap [REP_N_MEAN]->printHTML ( os );
				formItemMap [REP_N_STDEV]->printHTML ( os );
				formItemMap [RT_INT_START]->printHTML ( os );
				formItemMap [RT_INT_END]->printHTML ( os );
				formItemMap [RESOLUTION]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [FormItemIsotopePurity::getName ( "C", 13 )]->printHTML ( os );
				formItemMap [FormItemIsotopePurity::getName ( "N", 15 )]->printHTML ( os );
				formItemMap [FormItemIsotopePurity::getName ( "O", 18 )]->printHTML ( os );
				formItemMap [PURITY_CORRECTION]->printHTML ( os );
				formItemMap [REPORTER_ION_WINDOW]->printHTML ( os );
				ejb3.printFooter ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );

		if ( xLinkFlag ) {
			tableRowStart ( os );
				tableHeaderStart ( os, "", "left", true );
					ExpandableJavascriptBlock ejb4 ( "Crosslinking" );
					ejb4.printHeader ( os );
					formItemMap [REPORT_XL_PEPTIDE]->printHTML ( os );
					os << "<br />" << endl;
					formItemMap [REPORT_XL_SCORE]->printHTML ( os );
					if ( expectationFlag ) {
						formItemMap [REPORT_XL_EXPECTATION]->printHTML ( os );
						formItemMap [REPORT_XL_NLOG_P]->printHTML ( os );
					}
					formItemMap [REPORT_XL_RANK]->printHTML ( os );
					os << "<br />" << endl;
					formItemMap [REPORT_XL_LOW_SCORE]->printHTML ( os );
					if ( expectationFlag ) {
						formItemMap [REPORT_XL_LOW_EXPECTATION]->printHTML ( os );
						formItemMap [REPORT_XL_LOW_NLOG_P]->printHTML ( os );
					}
					os << "<br />" << endl;
					formItemMap [REPORT_XL_AA]->printHTML ( os );
					os << "<br />" << endl;
					formItemMap [XL_MIN_LOW_SCORE]->printHTML ( os );
					if ( expectationFlag ) {
						formItemMap [XL_MAX_LOW_EXPECTATION]->printHTML ( os );
					}
					formItemMap [XL_MIN_SCORE_DIFF]->printHTML ( os );
					os << "<br />" << endl;
					ejb4.printFooter ( os );
				tableHeaderEnd ( os );
			tableRowEnd ( os );
		}

		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true );
				tableStart ( os, false );
					tableRowStart ( os );
						tableHeaderStart ( os );
							formItemMap [FormItemSCReportType::getName ()]->printHTML ( os );
						tableHeaderEnd ( os );
						tableHeaderStart ( os );
							formItemMap [REPORT_HOMOLOGOUS_PROTEINS]->printHTML ( os );
						tableHeaderEnd ( os );
					tableRowEnd ( os );
				tableEnd ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );

		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true );
				if ( numSearches > 1 ) {
					formItemMap [REPORT_HITS_TYPE]->printHTML ( os );
					formItemMap [MERGE_OPTION]->printHTML ( os );
				}
				formItemMap [SORT_TYPE]->printHTML ( os );
				formItemMap [SORT_TYPE_2]->printHTML ( os );
				formItemMap [UNMATCHED_SPECTRA]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );

		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				tableStart ( os );
					tableRowStart ( os );
						tableHeaderStart ( os, "", "center", true );
							printHTMLFORMSubmit ( os, "Compare Searches" );
						tableHeaderEnd ( os );
						tableHeaderStart ( os, "", "center", true );
							formItemMap [FormItemSaveParameters::getName ()]->printHTML ( os );
						tableHeaderEnd ( os );
					tableRowEnd ( os );
				tableEnd ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );

		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true );
				formItemMap [FormItemMSMSPkFilter::getName ()]->printHTML ( os );
				formItemMap [FormItemMaxMSMSPeaks::getName ()]->printHTML ( os );
				formItemMap [FormItemMaxReportedHits::getName ("msms_")]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
	tableEnd ( os );
	showHiddenItems ( os );
	printHTMLFORMEnd ( os );
	printHTMLFooter ( os, "2002" );
}
void SearchCompareForm::printHTMLDiscriminentItems ( ostream& os )
{
	os << "Min Best Discriminant Score" << endl;
	StringVector instList = DiscriminantScoreInstrumentList::instance ().getNames ();
	for ( StringVectorSizeType i = 0, j = 0 ; i < instList.size () ; i++ ) {
		if ( find ( instrumentList.begin (), instrumentList.end (), instList [i] ) != instrumentList.end () ) {
			string in = genTranslateCharacter(instList [i],'-','_');
			formItemMap [MIN_BEST_DISC_SCORE+in]->printHTML ( os );
			if ( (++j) % 5 == 0 ) os << "<br />" << endl;
		}
	}
	os << "<br />" << endl;
}
SearchCompareFormCalibrationForm::SearchCompareFormCalibrationForm ( const VectorConstParameterListPtr& params ) :
	SearchCompareForm ( params, "Calibration" )
{
}
void SearchCompareFormCalibrationForm::printHTML ( ostream& os )
{
	init_html ( os, "Search Compare", "batchtagman.htm#search_compare" );
	fvj.print ( os );
	printHTMLFORMStart ( os, "post", "searchCompare", false, true );
	tableStart ( os, true );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true );
				printHTMLDiscriminentItems ( os );
				formItemMap [MIN_PROTEIN_SCORE]->printHTML ( os );
				formItemMap [MIN_PEPTIDE_SCORE]->printHTML ( os );
				if ( expectationFlag ) {
					formItemMap [MAX_PROTEIN_EVALUE]->printHTML ( os );
					formItemMap [MAX_PEPTIDE_EVALUE]->printHTML ( os );
				}
				os << "<br />" << endl;
			tableHeaderEnd ( os );
		tableRowEnd ( os );

		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				printHTMLFORMSubmit ( os, "Calibrate" );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
	tableEnd ( os );
	showHiddenItems ( os );
	printHTMLFORMEnd ( os );
	printHTMLFooter ( os, "2002" );
}
#endif
