/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_tag_frag.cpp                                               *
*                                                                             *
*  Created    : July 21st 1996                                                *
*                                                                             *
*  Purpose    : Contains most of the search algorithms for MS-Tag.            *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1996-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <climits>
#include <lg_string.h>
#include <lu_get_link.h>
#include <lu_mass_conv.h>
#include <lu_tag_frag.h>
#include <lu_tag_par.h>
#include <lu_immonium.h>
using std::string;
using std::fill;

MSMSSearch* getMSTagSearch ( MSMSDataPoint* file, MSTagParameters& params )
{
	StringVector searchTypes = params.getModificationParameters ().getHomologyTypes ();

	if ( params.getModificationParameters ().getSimpleTagSearch () ) {
		return new MSTagSearch ( file, params );
	}
	if ( !searchTypes.empty () ) {
		if ( searchTypes [0] == "low mass error" ) {
			string ionSeries = params.getIonSeries ();
			if ( ionSeries == "a" || ionSeries == "b" || ionSeries == "c" )
				return new MSSeqMass3SeqSearch ( file, params );
			if ( ionSeries == "y" )
				return new MSSeqMass1SeqSearch ( file, params );
		}
		if ( searchTypes [0] == "high mass error" ) {
			string ionSeries = params.getIonSeries ();
			if ( ionSeries == "a" || ionSeries == "b" || ionSeries == "c" )
				return new MSSeqMass1SeqSearch ( file, params );
			if ( ionSeries == "y" )
				return new MSSeqMass3SeqSearch ( file, params );
		}
		if ( searchTypes [0] == "no errors" )
			return new MSSeqSearch ( file, params );
		if ( searchTypes [0] == "middle masses error" )
			return new MSSeqMass1Mass3Search ( file, params );
		if ( searchTypes [0] == "parent mass" )
			return new MSMSSearch ( file, params );
	}
	return new MSTagSearchAllowErrors ( file, params );	// default case
}
MSMSSearch::MSMSSearch ( MSMSDataPoint* file, MSTagParameters& params ) :
	parentTolerance ( file->getPrecursorTolerance () ),
	parentPeak ( file->getPrecursorMZ (), parentTolerance, file->getPrecursorCharge (), file->getPrecursorIntensity (), getAdductMass ( params.getMonoisotopicFlag () ), params.getAverageParentMonoFragments () ),
	peaks ( file, params.getMSMSPeakFilterOptions (), &parentPeak, params.getParentMassTolerance (), params.getProductMassTolerance (), params.getMonoisotopicFlag (), params.getMonoParentAverageFragments () ),
	parentMass ( parentPeak.getMass () ),
	compositionSearch ( 0 ),
	unmatchedCompositionSearch ( 0 )
{
	parentMassPlusNegTolerance = parentMass + parentTolerance;
	parentMassMinusPosTolerance = parentMass - parentTolerance;
	if ( params.getCompSearchParams ().getSearchFlag () )
		compositionSearch = new CompositionSearch ( params.getCompSearchParams (), peaks );
	else {
		if ( params.getBiemannParameters ().get_immonium () ) {
			unmatchedCompositionSearch = new UnmatchedCompositionSearch ( peaks );
		}
	}
}
MSMSSearch::~MSMSSearch ()
{
	delete compositionSearch;
	delete unmatchedCompositionSearch;
}
bool MSMSSearch::doMatch ( const string& peptide, bool nTermPeptide, bool cTermPeptide, double molWt, TagMatchVector& tagMatch, const ScoreType& minScore )
{
	tagMatch.clear ();
	tagMatch.push_back ( TagMatch ( &parentPeak, 0, 0, PeptideSequence ( peptide ) ) );
	return true;
}
bool MSTagSearch::ESI_TRAP_CID_low_res = false;
bool MSTagSearch::ESI_Q_CID = false;
bool MSTagSearch::ESI_ETD_low_res = false;
bool MSTagSearch::ESI_ETD_high_res = false;
bool MSTagSearch::etd = false;
MassType MSTagSearch::internalNTerminusWt;
MassType MSTagSearch::nTerminusWt;
MassType MSTagSearch::cTerminusWt;
bool MSTagSearch::flagsSet = false;
bool MSTagSearch::a_nh3_flag;
bool MSTagSearch::a_h2o_flag;
bool MSTagSearch::a_flag;
bool MSTagSearch::a_h3po4_flag;
bool MSTagSearch::b_h2o_flag;
bool MSTagSearch::b_plus_h2o_flag;
bool MSTagSearch::c_ladder_flag;
bool MSTagSearch::b_nh3_flag;
bool MSTagSearch::b_soch4_flag;
bool MSTagSearch::b_h3po4_flag;
bool MSTagSearch::b_flag;
bool MSTagSearch::cPlus2DaFlag;
bool MSTagSearch::cPlus1DaFlag;
bool MSTagSearch::c_flag;
bool MSTagSearch::cp2_flag;
bool MSTagSearch::cMinus1DaFlag;
bool MSTagSearch::d_flag;

bool MSTagSearch::v_flag;
bool MSTagSearch::w_flag;
bool MSTagSearch::x_flag;
bool MSTagSearch::y_h2o_flag;
bool MSTagSearch::y_nh3_flag;
bool MSTagSearch::y_soch4_flag;
bool MSTagSearch::y_h3po4_flag;
bool MSTagSearch::y_flag;
bool MSTagSearch::Y_flag;
bool MSTagSearch::z_flag;
bool MSTagSearch::zp2_flag;
bool MSTagSearch::zPlus1Dap2_flag;
bool MSTagSearch::zPlus1DaFlag;
bool MSTagSearch::zPlus2DaFlag;
bool MSTagSearch::zPlus3DaFlag;
bool MSTagSearch::n_ladder_flag;

bool MSTagSearch::bp2_flag;
bool MSTagSearch::bp3_flag;
bool MSTagSearch::bp2_nh3_flag;
bool MSTagSearch::bp2_h2o_flag;
bool MSTagSearch::bp2_soch4_flag;
bool MSTagSearch::bp2_h3po4_flag;
bool MSTagSearch::yp2_flag;
bool MSTagSearch::yp3_flag;
bool MSTagSearch::yp2_nh3_flag;
bool MSTagSearch::yp2_h2o_flag;
bool MSTagSearch::yp2_soch4_flag;
bool MSTagSearch::yp2_h3po4_flag;

ScoreType MSTagSearch::a_nh3_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::a_h2o_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::a_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::a_h3po4_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::b_h2o_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::b_plus_h2o_score =  DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::c_ladder_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::b_nh3_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::b_soch4_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::b_h3po4_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::b_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::cPlus2DaScore = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::cPlus1DaScore = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::c_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::cMinus1DaScore = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::d_score = DEFAULT_HIT_SCORE;

ScoreType MSTagSearch::v_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::w_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::x_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::y_h2o_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::y_nh3_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::y_soch4_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::y_h3po4_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::y_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::Y_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::z_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::zPlus1DaScore = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::zPlus2DaScore = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::zPlus3DaScore = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::n_ladder_score = DEFAULT_HIT_SCORE;

ScoreType MSTagSearch::bp2_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::bp3_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::bp2_nh3_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::bp2_h2o_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::bp2_soch4_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::bp2_h3po4_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::cp2_score = DEFAULT_HIT_SCORE;

ScoreType MSTagSearch::yp2_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::yp3_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::yp2_nh3_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::yp2_h2o_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::yp2_soch4_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::yp2_h3po4_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::zp2_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::zPlus1Dap2Score = DEFAULT_HIT_SCORE;

ScoreType MSTagSearch::internal_a_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::internal_b_score = DEFAULT_HIT_SCORE;
ScoreType MSTagSearch::internal_loss_score = DEFAULT_HIT_SCORE;

ScoreType MSTagSearch::unmatched_score = DEFAULT_MISS_SCORE;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_a_scoreZ1NC			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_a_z1_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_scoreZ1NC			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_z1_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_nh3_scoreZ1NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_nh3_z1_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h2o_scoreZ1NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h2o_z1_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_soch4_scoreZ1NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_soch4_z1_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h3po4_scoreZ1NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h3po4_z1_NC" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_a_scoreZ1N			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_a_z1_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_scoreZ1N			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_z1_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_nh3_scoreZ1N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_nh3_z1_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h2o_scoreZ1N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h2o_z1_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_soch4_scoreZ1N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_soch4_z1_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h3po4_scoreZ1N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h3po4_z1_N" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_a_scoreZ1C			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_a_z1_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_scoreZ1C			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_z1_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_nh3_scoreZ1C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_nh3_z1_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h2o_scoreZ1C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h2o_z1_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_soch4_scoreZ1C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_soch4_z1_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h3po4_scoreZ1C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h3po4_z1_C" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_a_scoreZ1			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_a_z1_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_scoreZ1			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_z1_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_nh3_scoreZ1		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_nh3_z1_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h2o_scoreZ1		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h2o_z1_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_soch4_scoreZ1		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_soch4_z1_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h3po4_scoreZ1		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h3po4_z1_NO" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_a_scoreZ2NC			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_a_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_scoreZ2NC			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_scoreZ2NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_nh3_scoreZ2NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_nh3_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_nh3_scoreZ2NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_nh3_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h2o_scoreZ2NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h2o_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h2o_scoreZ2NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h2o_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_soch4_scoreZ2NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_soch4_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_soch4_scoreZ2NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_soch4_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h3po4_scoreZ2NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h3po4_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ2NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h3po4_z2_NC" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_a_scoreZ2N			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_a_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_scoreZ2N			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_scoreZ2N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_nh3_scoreZ2N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_nh3_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_nh3_scoreZ2N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_nh3_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h2o_scoreZ2N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h2o_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h2o_scoreZ2N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h2o_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_soch4_scoreZ2N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_soch4_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_soch4_scoreZ2N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_soch4_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h3po4_scoreZ2N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h3po4_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ2N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h3po4_z2_N" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_a_scoreZ2C			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_a_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_scoreZ2C			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_scoreZ2C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_nh3_scoreZ2C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_nh3_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_nh3_scoreZ2C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_nh3_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h2o_scoreZ2C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h2o_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h2o_scoreZ2C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h2o_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_soch4_scoreZ2C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_soch4_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_soch4_scoreZ2C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_soch4_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h3po4_scoreZ2C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h3po4_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ2C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h3po4_z2_C" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_a_scoreZ2			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_a_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_scoreZ2			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_scoreZ2			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_nh3_scoreZ2		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_nh3_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_nh3_scoreZ2		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_nh3_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h2o_scoreZ2		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h2o_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h2o_scoreZ2		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h2o_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_soch4_scoreZ2		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_soch4_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_soch4_scoreZ2	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_soch4_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h3po4_scoreZ2		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h3po4_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ2	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h3po4_z2_NO" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_a_scoreZ3NC			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_a_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_scoreZ3NC			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_scoreZ3NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp3_scoreZ3NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp3_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_nh3_scoreZ3NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_nh3_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_nh3_scoreZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_nh3_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h2o_scoreZ3NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h2o_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h2o_scoreZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h2o_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_soch4_scoreZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_soch4_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_soch4_scoreZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_soch4_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h3po4_scoreZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h3po4_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h3po4_z3_NC" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_a_scoreZ3N			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_a_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_scoreZ3N			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_scoreZ3N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp3_scoreZ3N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp3_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_nh3_scoreZ3N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_nh3_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_nh3_scoreZ3N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_nh3_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h2o_scoreZ3N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h2o_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h2o_scoreZ3N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h2o_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_soch4_scoreZ3N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_soch4_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_soch4_scoreZ3N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_soch4_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h3po4_scoreZ3N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h3po4_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ3N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h3po4_z3_N" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_a_scoreZ3C			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_a_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_scoreZ3C			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_scoreZ3C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp3_scoreZ3C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp3_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_nh3_scoreZ3C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_nh3_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_nh3_scoreZ3C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_nh3_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h2o_scoreZ3C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h2o_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h2o_scoreZ3C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h2o_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_soch4_scoreZ3C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_soch4_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_soch4_scoreZ3C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_soch4_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h3po4_scoreZ3C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h3po4_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ3C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h3po4_z3_C" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_a_scoreZ3			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_a_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_scoreZ3			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_scoreZ3			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp3_scoreZ3			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp3_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_nh3_scoreZ3		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_nh3_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_nh3_scoreZ3		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_nh3_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h2o_scoreZ3		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h2o_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h2o_scoreZ3		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h2o_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_soch4_scoreZ3		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_soch4_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_soch4_scoreZ3	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_soch4_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h3po4_scoreZ3		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h3po4_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ3	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h3po4_z3_NO" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_a_scoreZ4NC			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_a_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_scoreZ4NC			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_scoreZ4NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp3_scoreZ4NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp3_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_nh3_scoreZ4NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_nh3_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_nh3_scoreZ4NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_nh3_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h2o_scoreZ4NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h2o_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h2o_scoreZ4NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h2o_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_soch4_scoreZ4NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_soch4_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_soch4_scoreZ4NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_soch4_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h3po4_scoreZ4NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h3po4_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ4NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h3po4_z4_NC" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_a_scoreZ4N			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_a_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_scoreZ4N			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_scoreZ4N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp3_scoreZ4N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp3_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_nh3_scoreZ4N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_nh3_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_nh3_scoreZ4N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_nh3_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h2o_scoreZ4N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h2o_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h2o_scoreZ4N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h2o_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_soch4_scoreZ4N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_soch4_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_soch4_scoreZ4N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_soch4_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h3po4_scoreZ4N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h3po4_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ4N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h3po4_z4_N" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_a_scoreZ4C			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_a_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_scoreZ4C			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_scoreZ4C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp3_scoreZ4C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp3_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_nh3_scoreZ4C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_nh3_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_nh3_scoreZ4C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_nh3_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h2o_scoreZ4C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h2o_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h2o_scoreZ4C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h2o_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_soch4_scoreZ4C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_soch4_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_soch4_scoreZ4C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_soch4_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h3po4_scoreZ4C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h3po4_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ4C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h3po4_z4_C" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_a_scoreZ4			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_a_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_scoreZ4			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_scoreZ4			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp3_scoreZ4			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp3_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_nh3_scoreZ4		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_nh3_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_nh3_scoreZ4		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_nh3_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h2o_scoreZ4		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h2o_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h2o_scoreZ4		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h2o_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_soch4_scoreZ4		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_soch4_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_soch4_scoreZ4	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_soch4_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h3po4_scoreZ4		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h3po4_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ4	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h3po4_z4_NO" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_a_scoreZ5NC			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_a_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_scoreZ5NC			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_scoreZ5NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp3_scoreZ5NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp3_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_nh3_scoreZ5NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_nh3_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_nh3_scoreZ5NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_nh3_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h2o_scoreZ5NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h2o_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h2o_scoreZ5NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h2o_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_soch4_scoreZ5NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_soch4_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_soch4_scoreZ5NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_soch4_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h3po4_scoreZ5NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h3po4_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ5NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h3po4_z5_NC" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_a_scoreZ5N			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_a_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_scoreZ5N			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_scoreZ5N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp3_scoreZ5N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp3_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_nh3_scoreZ5N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_nh3_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_nh3_scoreZ5N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_nh3_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h2o_scoreZ5N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h2o_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h2o_scoreZ5N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h2o_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_soch4_scoreZ5N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_soch4_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_soch4_scoreZ5N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_soch4_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h3po4_scoreZ5N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h3po4_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ5N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h3po4_z5_N" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_a_scoreZ5C			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_a_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_scoreZ5C			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_scoreZ5C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp3_scoreZ5C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp3_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_nh3_scoreZ5C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_nh3_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_nh3_scoreZ5C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_nh3_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h2o_scoreZ5C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h2o_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h2o_scoreZ5C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h2o_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_soch4_scoreZ5C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_soch4_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_soch4_scoreZ5C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_soch4_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h3po4_scoreZ5C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h3po4_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ5C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h3po4_z5_C" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_a_scoreZ5			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_a_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_scoreZ5			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_scoreZ5			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp3_scoreZ5			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp3_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_nh3_scoreZ5		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_nh3_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_nh3_scoreZ5		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_nh3_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h2o_scoreZ5		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h2o_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h2o_scoreZ5		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h2o_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_soch4_scoreZ5		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_soch4_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_soch4_scoreZ5	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_soch4_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_b_h3po4_scoreZ5		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_b_h3po4_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ5	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_bp2_h3po4_z5_NO" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_scoreZ1NC			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_z1_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_nh3_scoreZ1NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_nh3_z1_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h2o_scoreZ1NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h2o_z1_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_soch4_scoreZ1NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_soch4_z1_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h3po4_scoreZ1NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h3po4_z1_NC" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_scoreZ1N			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_z1_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_nh3_scoreZ1N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_nh3_z1_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h2o_scoreZ1N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h2o_z1_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_soch4_scoreZ1N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_soch4_z1_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h3po4_scoreZ1N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h3po4_z1_N" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_scoreZ1C			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_z1_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_nh3_scoreZ1C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_nh3_z1_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h2o_scoreZ1C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h2o_z1_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_soch4_scoreZ1C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_soch4_z1_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h3po4_scoreZ1C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h3po4_z1_C" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_scoreZ1			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_z1_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_nh3_scoreZ1		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_nh3_z1_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h2o_scoreZ1		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h2o_z1_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_soch4_scoreZ1		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_soch4_z1_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h3po4_scoreZ1		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h3po4_z1_NO" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_scoreZ2NC			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_scoreZ2NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_nh3_scoreZ2NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_nh3_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_nh3_scoreZ2NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_nh3_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h2o_scoreZ2NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h2o_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h2o_scoreZ2NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h2o_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_soch4_scoreZ2NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_soch4_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_soch4_scoreZ2NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_soch4_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h3po4_scoreZ2NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h3po4_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ2NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h3po4_z2_NC" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_scoreZ2N			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_scoreZ2N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_nh3_scoreZ2N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_nh3_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_nh3_scoreZ2N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_nh3_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h2o_scoreZ2N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h2o_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h2o_scoreZ2N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h2o_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_soch4_scoreZ2N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_soch4_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_soch4_scoreZ2N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_soch4_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h3po4_scoreZ2N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h3po4_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ2N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h3po4_z2_N" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_scoreZ2C			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_scoreZ2C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_nh3_scoreZ2C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_nh3_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_nh3_scoreZ2C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_nh3_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h2o_scoreZ2C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h2o_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h2o_scoreZ2C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h2o_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_soch4_scoreZ2C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_soch4_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_soch4_scoreZ2C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_soch4_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h3po4_scoreZ2C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h3po4_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ2C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h3po4_z2_C" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_scoreZ2			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_scoreZ2			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_nh3_scoreZ2		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_nh3_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_nh3_scoreZ2		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_nh3_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h2o_scoreZ2		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h2o_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h2o_scoreZ2		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h2o_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_soch4_scoreZ2		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_soch4_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_soch4_scoreZ2	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_soch4_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h3po4_scoreZ2		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h3po4_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ2	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h3po4_z2_NO" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_scoreZ3NC			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_scoreZ3NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp3_scoreZ3NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp3_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_nh3_scoreZ3NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_nh3_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_nh3_scoreZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_nh3_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h2o_scoreZ3NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h2o_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h2o_scoreZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h2o_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_soch4_scoreZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_soch4_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_soch4_scoreZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_soch4_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h3po4_scoreZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h3po4_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h3po4_z3_NC" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_scoreZ3N			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_scoreZ3N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp3_scoreZ3N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp3_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_nh3_scoreZ3N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_nh3_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_nh3_scoreZ3N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_nh3_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h2o_scoreZ3N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h2o_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h2o_scoreZ3N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h2o_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_soch4_scoreZ3N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_soch4_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_soch4_scoreZ3N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_soch4_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h3po4_scoreZ3N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h3po4_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ3N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h3po4_z3_N" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_scoreZ3C			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_scoreZ3C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp3_scoreZ3C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp3_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_nh3_scoreZ3C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_nh3_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_nh3_scoreZ3C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_nh3_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h2o_scoreZ3C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h2o_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h2o_scoreZ3C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h2o_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_soch4_scoreZ3C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_soch4_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_soch4_scoreZ3C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_soch4_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h3po4_scoreZ3C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h3po4_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ3C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h3po4_z3_C" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_scoreZ3			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_scoreZ3			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp3_scoreZ3			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp3_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_nh3_scoreZ3		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_nh3_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_nh3_scoreZ3		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_nh3_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h2o_scoreZ3		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h2o_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h2o_scoreZ3		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h2o_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_soch4_scoreZ3		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_soch4_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_soch4_scoreZ3	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_soch4_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h3po4_scoreZ3		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h3po4_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ3	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h3po4_z3_NO" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_scoreZ4NC			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_scoreZ4NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp3_scoreZ4NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp3_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_nh3_scoreZ4NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_nh3_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_nh3_scoreZ4NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_nh3_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h2o_scoreZ4NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h2o_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h2o_scoreZ4NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h2o_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_soch4_scoreZ4NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_soch4_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_soch4_scoreZ4NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_soch4_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h3po4_scoreZ4NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h3po4_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ4NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h3po4_z4_NC" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_scoreZ4N			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_scoreZ4N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp3_scoreZ4N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp3_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_nh3_scoreZ4N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_nh3_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_nh3_scoreZ4N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_nh3_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h2o_scoreZ4N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h2o_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h2o_scoreZ4N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h2o_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_soch4_scoreZ4N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_soch4_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_soch4_scoreZ4N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_soch4_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h3po4_scoreZ4N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h3po4_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ4N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h3po4_z4_N" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_scoreZ4C			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_scoreZ4C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp3_scoreZ4C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp3_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_nh3_scoreZ4C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_nh3_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_nh3_scoreZ4C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_nh3_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h2o_scoreZ4C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h2o_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h2o_scoreZ4C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h2o_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_soch4_scoreZ4C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_soch4_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_soch4_scoreZ4C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_soch4_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h3po4_scoreZ4C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h3po4_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ4C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h3po4_z4_C" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_scoreZ4			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_scoreZ4			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp3_scoreZ4			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp3_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_nh3_scoreZ4		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_nh3_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_nh3_scoreZ4		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_nh3_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h2o_scoreZ4		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h2o_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h2o_scoreZ4		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h2o_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_soch4_scoreZ4		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_soch4_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_soch4_scoreZ4	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_soch4_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h3po4_scoreZ4		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h3po4_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ4	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h3po4_z4_NO" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_scoreZ5NC			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_scoreZ5NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp3_scoreZ5NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp3_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_nh3_scoreZ5NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_nh3_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_nh3_scoreZ5NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_nh3_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h2o_scoreZ5NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h2o_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h2o_scoreZ5NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h2o_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_soch4_scoreZ5NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_soch4_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_soch4_scoreZ5NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_soch4_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h3po4_scoreZ5NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h3po4_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ5NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h3po4_z5_NC" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_scoreZ5N			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_scoreZ5N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp3_scoreZ5N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp3_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_nh3_scoreZ5N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_nh3_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_nh3_scoreZ5N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_nh3_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h2o_scoreZ5N		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h2o_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h2o_scoreZ5N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h2o_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_soch4_scoreZ5N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_soch4_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_soch4_scoreZ5N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_soch4_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h3po4_scoreZ5N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h3po4_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ5N	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h3po4_z5_N" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_scoreZ5C			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_scoreZ5C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp3_scoreZ5C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp3_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_nh3_scoreZ5C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_nh3_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_nh3_scoreZ5C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_nh3_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h2o_scoreZ5C		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h2o_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h2o_scoreZ5C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h2o_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_soch4_scoreZ5C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_soch4_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_soch4_scoreZ5C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_soch4_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h3po4_scoreZ5C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h3po4_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ5C	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h3po4_z5_C" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_scoreZ5			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_scoreZ5			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp3_scoreZ5			= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp3_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_nh3_scoreZ5		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_nh3_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_nh3_scoreZ5		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_nh3_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h2o_scoreZ5		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h2o_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h2o_scoreZ5		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h2o_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_soch4_scoreZ5		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_soch4_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_soch4_scoreZ5	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_soch4_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_y_h3po4_scoreZ5		= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_y_h3po4_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ5	= FragmentationParams::instance ().getDoubleValue ( "ESI_TRAP_CID_low_res_yp2_h3po4_z5_NO" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_a_scoreZ1NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_z1_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_nh3_scoreZ1NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_nh3_z1_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h2o_scoreZ1NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h2o_z1_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h3po4_scoreZ1NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h3po4_z1_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_scoreZ1NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_z1_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_nh3_scoreZ1NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_nh3_z1_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h2o_scoreZ1NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h2o_z1_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_soch4_scoreZ1NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_soch4_z1_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h3po4_scoreZ1NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h3po4_z1_NC" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_a_scoreZ1N			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_z1_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_nh3_scoreZ1N		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_nh3_z1_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h2o_scoreZ1N		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h2o_z1_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h3po4_scoreZ1N	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h3po4_z1_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_scoreZ1N			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_z1_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_nh3_scoreZ1N		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_nh3_z1_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h2o_scoreZ1N		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h2o_z1_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_soch4_scoreZ1N	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_soch4_z1_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h3po4_scoreZ1N	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h3po4_z1_N" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_a_scoreZ1C			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_z1_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_nh3_scoreZ1C		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_nh3_z1_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h2o_scoreZ1C		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h2o_z1_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h3po4_scoreZ1C	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h3po4_z1_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_scoreZ1C			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_z1_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_nh3_scoreZ1C		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_nh3_z1_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h2o_scoreZ1C		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h2o_z1_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_soch4_scoreZ1C	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_soch4_z1_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h3po4_scoreZ1C	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h3po4_z1_C" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_a_scoreZ1			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_z1_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_nh3_scoreZ1		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_nh3_z1_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h2o_scoreZ1		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h2o_z1_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h3po4_scoreZ1	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h3po4_z1_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_scoreZ1			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_z1_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_nh3_scoreZ1		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_nh3_z1_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h2o_scoreZ1		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h2o_z1_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_soch4_scoreZ1	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_soch4_z1_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h3po4_scoreZ1	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h3po4_z1_NO" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_a_scoreZ2NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_nh3_scoreZ2NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_nh3_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h2o_scoreZ2NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h2o_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h3po4_scoreZ2NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h3po4_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_scoreZ2NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_nh3_scoreZ2NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_nh3_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h2o_scoreZ2NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h2o_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_soch4_scoreZ2NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_soch4_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h3po4_scoreZ2NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h3po4_z2_NC" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_a_scoreZ2N			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_nh3_scoreZ2N		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_nh3_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h2o_scoreZ2N		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h2o_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h3po4_scoreZ2N	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h3po4_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_scoreZ2N			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_nh3_scoreZ2N		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_nh3_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h2o_scoreZ2N		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h2o_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_soch4_scoreZ2N	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_soch4_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h3po4_scoreZ2N	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h3po4_z2_N" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_a_scoreZ2C			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_nh3_scoreZ2C		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_nh3_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h2o_scoreZ2C		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h2o_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h3po4_scoreZ2C	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h3po4_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_scoreZ2C			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_nh3_scoreZ2C		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_nh3_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h2o_scoreZ2C		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h2o_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_soch4_scoreZ2C	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_soch4_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h3po4_scoreZ2C	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h3po4_z2_C" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_a_scoreZ2			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_nh3_scoreZ2		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_nh3_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h2o_scoreZ2		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h2o_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h3po4_scoreZ2	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h3po4_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_scoreZ2			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_nh3_scoreZ2		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_nh3_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h2o_scoreZ2		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h2o_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_soch4_scoreZ2	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_soch4_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h3po4_scoreZ2	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h3po4_z2_NO" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_a_scoreZ3NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_nh3_scoreZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_nh3_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h2o_scoreZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h2o_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h3po4_scoreZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h3po4_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_scoreZ3NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_nh3_scoreZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_nh3_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h2o_scoreZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h2o_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_soch4_scoreZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_soch4_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h3po4_scoreZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h3po4_z3_NC" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_a_scoreZ3N			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_nh3_scoreZ3N		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_nh3_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h2o_scoreZ3N		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h2o_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h3po4_scoreZ3N	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h3po4_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_scoreZ3N			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_nh3_scoreZ3N		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_nh3_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h2o_scoreZ3N		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h2o_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_soch4_scoreZ3N	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_soch4_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h3po4_scoreZ3N	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h3po4_z3_N" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_a_scoreZ3C			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_nh3_scoreZ3C		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_nh3_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h2o_scoreZ3C		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h2o_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h3po4_scoreZ3C	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h3po4_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_scoreZ3C			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_nh3_scoreZ3C		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_nh3_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h2o_scoreZ3C		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h2o_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_soch4_scoreZ3C	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_soch4_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h3po4_scoreZ3C	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h3po4_z3_C" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_a_scoreZ3			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_nh3_scoreZ3		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_nh3_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h2o_scoreZ3		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h2o_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h3po4_scoreZ3	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h3po4_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_scoreZ3			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_nh3_scoreZ3		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_nh3_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h2o_scoreZ3		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h2o_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_soch4_scoreZ3	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_soch4_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h3po4_scoreZ3	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h3po4_z3_NO" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_a_scoreZ4NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_nh3_scoreZ4NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_nh3_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h2o_scoreZ4NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h2o_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h3po4_scoreZ4NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h3po4_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_scoreZ4NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_nh3_scoreZ4NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_nh3_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h2o_scoreZ4NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h2o_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_soch4_scoreZ4NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_soch4_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h3po4_scoreZ4NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h3po4_z4_NC" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_a_scoreZ4N			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_nh3_scoreZ4N		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_nh3_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h2o_scoreZ4N		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h2o_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h3po4_scoreZ4N	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h3po4_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_scoreZ4N			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_nh3_scoreZ4N		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_nh3_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h2o_scoreZ4N		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h2o_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_soch4_scoreZ4N	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_soch4_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h3po4_scoreZ4N	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h3po4_z4_N" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_a_scoreZ4C			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_nh3_scoreZ4C		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_nh3_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h2o_scoreZ4C		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h2o_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h3po4_scoreZ4C	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h3po4_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_scoreZ4C			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_nh3_scoreZ4C		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_nh3_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h2o_scoreZ4C		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h2o_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_soch4_scoreZ4C	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_soch4_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h3po4_scoreZ4C	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h3po4_z4_C" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_a_scoreZ4			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_nh3_scoreZ4		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_nh3_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h2o_scoreZ4		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h2o_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h3po4_scoreZ4	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h3po4_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_scoreZ4			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_nh3_scoreZ4		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_nh3_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h2o_scoreZ4		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h2o_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_soch4_scoreZ4	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_soch4_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h3po4_scoreZ4	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h3po4_z4_NO" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_a_scoreZ5NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_nh3_scoreZ5NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_nh3_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h2o_scoreZ5NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h2o_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h3po4_scoreZ5NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h3po4_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_scoreZ5NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_nh3_scoreZ5NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_nh3_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h2o_scoreZ5NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h2o_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_soch4_scoreZ5NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_soch4_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h3po4_scoreZ5NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h3po4_z5_NC" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_a_scoreZ5N			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_nh3_scoreZ5N		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_nh3_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h2o_scoreZ5N		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h2o_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h3po4_scoreZ5N	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h3po4_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_scoreZ5N			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_nh3_scoreZ5N		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_nh3_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h2o_scoreZ5N		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h2o_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_soch4_scoreZ5N	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_soch4_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h3po4_scoreZ5N	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h3po4_z5_N" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_a_scoreZ5C			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_nh3_scoreZ5C		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_nh3_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h2o_scoreZ5C		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h2o_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h3po4_scoreZ5C	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h3po4_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_scoreZ5C			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_nh3_scoreZ5C		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_nh3_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h2o_scoreZ5C		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h2o_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_soch4_scoreZ5C	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_soch4_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h3po4_scoreZ5C	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h3po4_z5_C" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_a_scoreZ5			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_nh3_scoreZ5		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_nh3_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h2o_scoreZ5		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h2o_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_a_h3po4_scoreZ5	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_a_h3po4_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_scoreZ5			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_nh3_scoreZ5		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_nh3_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h2o_scoreZ5		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h2o_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_soch4_scoreZ5	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_soch4_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_b_h3po4_scoreZ5	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_b_h3po4_z5_NO" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_y_scoreZ1NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_z1_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_nh3_scoreZ1NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_nh3_z1_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h2o_scoreZ1NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h2o_z1_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_soch4_scoreZ1NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_soch4_z1_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h3po4_scoreZ1NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h3po4_z1_NC" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_y_scoreZ1N			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_z1_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_nh3_scoreZ1N		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_nh3_z1_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h2o_scoreZ1N		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h2o_z1_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_soch4_scoreZ1N	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_soch4_z1_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h3po4_scoreZ1N	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h3po4_z1_N" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_y_scoreZ1C			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_z1_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_nh3_scoreZ1C		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_nh3_z1_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h2o_scoreZ1C		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h2o_z1_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_soch4_scoreZ1C	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_soch4_z1_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h3po4_scoreZ1C	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h3po4_z1_C" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_y_scoreZ1			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_z1_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_nh3_scoreZ1		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_nh3_z1_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h2o_scoreZ1		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h2o_z1_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_soch4_scoreZ1	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_soch4_z1_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h3po4_scoreZ1	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h3po4_z1_NO" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_y_scoreZ2NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_nh3_scoreZ2NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_nh3_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h2o_scoreZ2NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h2o_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_soch4_scoreZ2NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_soch4_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h3po4_scoreZ2NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h3po4_z2_NC" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_y_scoreZ2N			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_nh3_scoreZ2N		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_nh3_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h2o_scoreZ2N		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h2o_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_soch4_scoreZ2N	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_soch4_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h3po4_scoreZ2N	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h3po4_z2_N" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_y_scoreZ2C			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_nh3_scoreZ2C		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_nh3_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h2o_scoreZ2C		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h2o_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_soch4_scoreZ2C	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_soch4_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h3po4_scoreZ2C	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h3po4_z2_C" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_y_scoreZ2			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_nh3_scoreZ2		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_nh3_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h2o_scoreZ2		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h2o_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_soch4_scoreZ2	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_soch4_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h3po4_scoreZ2	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h3po4_z2_NO" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_y_scoreZ3NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_nh3_scoreZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_nh3_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h2o_scoreZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h2o_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_soch4_scoreZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_soch4_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h3po4_scoreZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h3po4_z3_NC" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_y_scoreZ3N			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_nh3_scoreZ3N		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_nh3_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h2o_scoreZ3N		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h2o_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_soch4_scoreZ3N	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_soch4_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h3po4_scoreZ3N	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h3po4_z3_N" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_y_scoreZ3C			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_nh3_scoreZ3C		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_nh3_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h2o_scoreZ3C		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h2o_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_soch4_scoreZ3C	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_soch4_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h3po4_scoreZ3C	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h3po4_z3_C" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_y_scoreZ3			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_nh3_scoreZ3		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_nh3_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h2o_scoreZ3		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h2o_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_soch4_scoreZ3	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_soch4_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h3po4_scoreZ3	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h3po4_z3_NO" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_y_scoreZ4NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_nh3_scoreZ4NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_nh3_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h2o_scoreZ4NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h2o_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_soch4_scoreZ4NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_soch4_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h3po4_scoreZ4NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h3po4_z4_NC" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_y_scoreZ4N			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_nh3_scoreZ4N		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_nh3_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h2o_scoreZ4N		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h2o_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_soch4_scoreZ4N	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_soch4_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h3po4_scoreZ4N	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h3po4_z4_N" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_y_scoreZ4C			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_nh3_scoreZ4C		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_nh3_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h2o_scoreZ4C		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h2o_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_soch4_scoreZ4C	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_soch4_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h3po4_scoreZ4C	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h3po4_z4_C" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_y_scoreZ4			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_nh3_scoreZ4		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_nh3_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h2o_scoreZ4		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h2o_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_soch4_scoreZ4	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_soch4_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h3po4_scoreZ4	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h3po4_z4_NO" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_y_scoreZ5NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_nh3_scoreZ5NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_nh3_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h2o_scoreZ5NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h2o_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_soch4_scoreZ5NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_soch4_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h3po4_scoreZ5NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h3po4_z5_NC" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_y_scoreZ5N			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_nh3_scoreZ5N		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_nh3_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h2o_scoreZ5N		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h2o_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_soch4_scoreZ5N	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_soch4_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h3po4_scoreZ5N	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h3po4_z5_N" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_y_scoreZ5C			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_nh3_scoreZ5C		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_nh3_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h2o_scoreZ5C		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h2o_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_soch4_scoreZ5C	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_soch4_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h3po4_scoreZ5C	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h3po4_z5_C" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_Q_CID_y_scoreZ5			= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_nh3_scoreZ5		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_nh3_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h2o_scoreZ5		= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h2o_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_soch4_scoreZ5	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_soch4_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_Q_CID_y_h3po4_scoreZ5	= FragmentationParams::instance ().getDoubleValue ( "ESI_Q_CID_y_h3po4_z5_NO" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_low_res_b_scoreZ2NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_b_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_c_scoreZ2NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_c_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_cMinus1DaZ2NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_cMinus1Da_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_y_scoreZ2NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_y_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_z_scoreZ2NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_z_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zPlus1DaZ2NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zPlus1Da_z2_NC" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_low_res_b_scoreZ2N		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_b_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_c_scoreZ2N		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_c_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_cMinus1DaZ2N		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_cMinus1Da_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_y_scoreZ2N		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_y_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_z_scoreZ2N		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_z_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zPlus1DaZ2N		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zPlus1Da_z2_N" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_low_res_b_scoreZ2C		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_b_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_c_scoreZ2C		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_c_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_cMinus1DaZ2C		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_cMinus1Da_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_y_scoreZ2C		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_y_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_z_scoreZ2C		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_z_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zPlus1DaZ2C		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zPlus1Da_z2_C" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_low_res_b_scoreZ2		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_b_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_c_scoreZ2		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_c_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_cMinus1DaZ2		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_cMinus1Da_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_y_scoreZ2		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_y_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_z_scoreZ2		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_z_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zPlus1DaZ2		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zPlus1Da_z2_NO" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_low_res_c_scoreZ3NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_c_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_cp2_scoreZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_cp2_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_cMinus1DaZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_cMinus1Da_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_y_scoreZ3NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_y_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_z_scoreZ3NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_z_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zPlus1DaZ3NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zPlus1Da_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zp2_scoreZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zp2_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zPlus1Dap2Z3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zPlus1Dap2_z3_NC" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_low_res_c_scoreZ3N		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_c_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_cp2_scoreZ3N		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_cp2_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_cMinus1DaZ3N		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_cMinus1Da_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_y_scoreZ3N		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_y_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_z_scoreZ3N		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_z_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zPlus1DaZ3N		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zPlus1Da_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zp2_scoreZ3N		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zp2_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zPlus1Dap2Z3N	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zPlus1Dap2_z3_N" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_low_res_c_scoreZ3C		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_c_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_cp2_scoreZ3C		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_cp2_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_cMinus1DaZ3C		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_cMinus1Da_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_y_scoreZ3C		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_y_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_z_scoreZ3C		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_z_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zPlus1DaZ3C		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zPlus1Da_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zp2_scoreZ3C		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zp2_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zPlus1Dap2Z3C	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zPlus1Dap2_z3_C" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_low_res_c_scoreZ3		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_c_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_cp2_scoreZ3		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_cp2_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_cMinus1DaZ3		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_cMinus1Da_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_y_scoreZ3		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_y_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_z_scoreZ3		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_z_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zPlus1DaZ3		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zPlus1Da_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zp2_scoreZ3		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zp2_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zPlus1Dap2Z3		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zPlus1Dap2_z3_NO" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_low_res_c_scoreZ4NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_c_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_cp2_scoreZ4NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_cp2_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_y_scoreZ4NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_y_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_z_scoreZ4NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_z_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zPlus1DaZ4NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zPlus1Da_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zp2_scoreZ4NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zp2_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zPlus1Dap2Z4NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zPlus1Dap2_z4_NC" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_low_res_c_scoreZ4N		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_c_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_cp2_scoreZ4N		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_cp2_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_y_scoreZ4N		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_y_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_z_scoreZ4N		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_z_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zPlus1DaZ4N		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zPlus1Da_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zp2_scoreZ4N		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zp2_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zPlus1Dap2Z4N	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zPlus1Dap2_z4_N" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_low_res_c_scoreZ4C		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_c_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_cp2_scoreZ4C		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_cp2_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_y_scoreZ4C		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_y_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_z_scoreZ4C		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_z_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zPlus1DaZ4C		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zPlus1Da_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zp2_scoreZ4C		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zp2_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zPlus1Dap2Z4C	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zPlus1Dap2_z4_C" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_low_res_c_scoreZ4		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_c_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_cp2_scoreZ4		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_cp2_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_y_scoreZ4		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_y_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_z_scoreZ4		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_z_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zPlus1DaZ4		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zPlus1Da_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zp2_scoreZ4		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zp2_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zPlus1Dap2Z4		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zPlus1Dap2_z4_NO" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_low_res_c_scoreZ5NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_c_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_cp2_scoreZ5NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_cp2_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_y_scoreZ5NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_y_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_z_scoreZ5NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_z_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zPlus1DaZ5NC		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zPlus1Da_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zp2_scoreZ5NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zp2_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zPlus1Dap2Z5NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zPlus1Dap2_z5_NC" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_low_res_c_scoreZ5N		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_c_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_cp2_scoreZ5N		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_cp2_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_y_scoreZ5N		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_y_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_z_scoreZ5N		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_z_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zPlus1DaZ5N		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zPlus1Da_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zp2_scoreZ5N		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zp2_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zPlus1Dap2Z5N	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zPlus1Dap2_z5_N" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_low_res_c_scoreZ5C		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_c_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_cp2_scoreZ5C		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_cp2_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_y_scoreZ5C		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_y_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_z_scoreZ5C		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_z_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zPlus1DaZ5C		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zPlus1Da_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zp2_scoreZ5C		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zp2_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zPlus1Dap2Z5C	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zPlus1Dap2_z5_C" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_low_res_c_scoreZ5		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_c_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_cp2_scoreZ5		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_cp2_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_y_scoreZ5		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_y_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_z_scoreZ5		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_z_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zPlus1DaZ5		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zPlus1Da_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zp2_scoreZ5		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zp2_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_low_res_zPlus1Dap2Z5		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_low_res_zPlus1Dap2_z5_NO" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_hi_res_b_scoreZ2NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_b_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_c_scoreZ2NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_c_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_cMinus1DaZ2NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_cMinus1Da_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_y_scoreZ2NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_y_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_z_scoreZ2NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_z_z2_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_zPlus1DaZ2NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_zPlus1Da_z2_NC" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_hi_res_b_scoreZ2N	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_b_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_c_scoreZ2N	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_c_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_cMinus1DaZ2N	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_cMinus1Da_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_y_scoreZ2N	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_y_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_z_scoreZ2N	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_z_z2_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_zPlus1DaZ2N	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_zPlus1Da_z2_N" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_hi_res_b_scoreZ2C	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_b_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_c_scoreZ2C	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_c_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_cMinus1DaZ2C	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_cMinus1Da_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_y_scoreZ2C	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_y_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_z_scoreZ2C	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_z_z2_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_zPlus1DaZ2C	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_zPlus1Da_z2_C" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_hi_res_b_scoreZ2		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_b_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_c_scoreZ2		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_c_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_cMinus1DaZ2	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_cMinus1Da_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_y_scoreZ2		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_y_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_z_scoreZ2		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_z_z2_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_zPlus1DaZ2	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_zPlus1Da_z2_NO" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_hi_res_c_scoreZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_c_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_cMinus1DaZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_cMinus1Da_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_y_scoreZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_y_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_z_scoreZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_z_z3_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_zPlus1DaZ3NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_zPlus1Da_z3_NC" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_hi_res_c_scoreZ3N	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_c_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_cMinus1DaZ3N	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_cMinus1Da_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_y_scoreZ3N	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_y_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_z_scoreZ3N	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_z_z3_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_zPlus1DaZ3N	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_zPlus1Da_z3_N" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_hi_res_c_scoreZ3C	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_c_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_cMinus1DaZ3C	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_cMinus1Da_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_y_scoreZ3C	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_y_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_z_scoreZ3C	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_z_z3_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_zPlus1DaZ3C	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_zPlus1Da_z3_C" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_hi_res_c_scoreZ3		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_c_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_cMinus1DaZ3	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_cMinus1Da_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_y_scoreZ3		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_y_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_z_scoreZ3		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_z_z3_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_zPlus1DaZ3	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_zPlus1Da_z3_NO" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_hi_res_c_scoreZ4NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_c_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_y_scoreZ4NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_y_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_z_scoreZ4NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_z_z4_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_zPlus1DaZ4NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_zPlus1Da_z4_NC" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_hi_res_c_scoreZ4N	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_c_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_y_scoreZ4N	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_y_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_z_scoreZ4N	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_z_z4_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_zPlus1DaZ4N	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_zPlus1Da_z4_N" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_hi_res_c_scoreZ4C	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_c_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_y_scoreZ4C	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_y_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_z_scoreZ4C	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_z_z4_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_zPlus1DaZ4C	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_zPlus1Da_z4_C" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_hi_res_c_scoreZ4		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_c_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_y_scoreZ4		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_y_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_z_scoreZ4		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_z_z4_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_zPlus1DaZ4	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_zPlus1Da_z4_NO" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_hi_res_c_scoreZ5NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_c_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_y_scoreZ5NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_y_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_z_scoreZ5NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_z_z5_NC" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_zPlus1DaZ5NC	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_zPlus1Da_z5_NC" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_hi_res_c_scoreZ5N	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_c_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_y_scoreZ5N	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_y_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_z_scoreZ5N	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_z_z5_N" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_zPlus1DaZ5N	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_zPlus1Da_z5_N" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_hi_res_c_scoreZ5C	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_c_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_y_scoreZ5C	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_y_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_z_scoreZ5C	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_z_z5_C" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_zPlus1DaZ5C	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_zPlus1Da_z5_C" ) * SCORE_TYPE_MULTIPLIER;

ScoreType MSTagSearch::ESI_ETD_hi_res_c_scoreZ5		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_c_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_y_scoreZ5		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_y_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_z_scoreZ5		= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_z_z5_NO" ) * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::ESI_ETD_hi_res_zPlus1DaZ5	= FragmentationParams::instance ().getDoubleValue ( "ESI_ETD_hi_res_zPlus1Da_z5_NO" ) * SCORE_TYPE_MULTIPLIER;

int MSTagSearch::a_nh3_index;
int MSTagSearch::a_h2o_index;
int MSTagSearch::a_index;
int MSTagSearch::a_h3po4_index;
int MSTagSearch::b_h2o_index;
int MSTagSearch::b_plus_h2o_index;
int MSTagSearch::c_ladder_index;
int MSTagSearch::b_nh3_index;
int MSTagSearch::b_soch4_index;
int MSTagSearch::b_h3po4_index;
int MSTagSearch::b_index;
int MSTagSearch::cPlus2DaIndex;
int MSTagSearch::cPlus1DaIndex;
int MSTagSearch::c_index;
int MSTagSearch::cp2_index;
int MSTagSearch::cMinus1DaIndex;

int MSTagSearch::v_index;
int MSTagSearch::x_index;
int MSTagSearch::Y_index;
int MSTagSearch::y_h2o_index;
int MSTagSearch::y_nh3_index;
int MSTagSearch::y_soch4_index;
int MSTagSearch::y_h3po4_index;
int MSTagSearch::y_index;
int MSTagSearch::z_index;
int MSTagSearch::zp2_index;
int MSTagSearch::zPlus1Dap2_index;
int MSTagSearch::zPlus1DaIndex;
int MSTagSearch::zPlus2DaIndex;
int MSTagSearch::zPlus3DaIndex;

int MSTagSearch::bp2_index;
int MSTagSearch::bp3_index;
int MSTagSearch::bp2_nh3_index;
int MSTagSearch::bp2_h2o_index;
int MSTagSearch::bp2_soch4_index;
int MSTagSearch::bp2_h3po4_index;
int MSTagSearch::yp2_index;
int MSTagSearch::yp3_index;
int MSTagSearch::yp2_nh3_index;
int MSTagSearch::yp2_h2o_index;
int MSTagSearch::yp2_soch4_index;
int MSTagSearch::yp2_h3po4_index;

int MSTagSearch::numNIonTypes;
int MSTagSearch::numCIonTypes;
int MSTagSearch::numInternalIonTypes;

bool MSTagSearch::checkInternalIons;
bool MSTagSearch::checkNIons;
bool MSTagSearch::checkCIons;
bool MSTagSearch::checkSatelliteIons;
bool MSTagSearch::checkNIonAmmoniaLoss;
bool MSTagSearch::checkNIonPosChargeBearing;
bool MSTagSearch::checkNIonPhosphorylation;
bool MSTagSearch::checkNIonM;
bool MSTagSearch::checkNIonWaterLoss;
bool MSTagSearch::checkCIonPosChargeBearing;
bool MSTagSearch::checkCIonAmmoniaLoss;
bool MSTagSearch::checkCIonPhosphorylation;
bool MSTagSearch::checkCIonM;
bool MSTagSearch::checkCIonWaterLoss;

unsigned int MSTagSearch::dIonExcludeMask;
unsigned int MSTagSearch::vIonExcludeMask;
unsigned int MSTagSearch::wIonExcludeMask;
unsigned int MSTagSearch::ammoniaLossMask;
unsigned int MSTagSearch::posChargeBearingMask;
unsigned int MSTagSearch::waterLossMask;

static unsigned int waterLossFlag [AA_ARRAY_SIZE];
static unsigned int ammoniaLossFlag [AA_ARRAY_SIZE];
static unsigned int posChargeBearingFlag [AA_ARRAY_SIZE];
static unsigned int oxidizedMFlag [AA_ARRAY_SIZE];
static unsigned int phosphorylationFlag [AA_ARRAY_SIZE];
static unsigned int dIonExcludeFlag [AA_ARRAY_SIZE];
static unsigned int vIonExcludeFlag [AA_ARRAY_SIZE];
static unsigned int wIonExcludeFlag [AA_ARRAY_SIZE];

MassType MSTagSearch::maxInternalIonMass;
double MSTagSearch::xLinkMass = 0.0;
double MSTagSearch::xLinkImmoniumMass = 0.0;

ScoreType MSTagSearch::maxPScore		= 100.0 * SCORE_TYPE_MULTIPLIER;
ScoreType MSTagSearch::pCIDScore		= DEFAULT_MISS_SCORE;
ScoreType MSTagSearch::pCIDXLScore		= DEFAULT_MISS_SCORE;
ScoreType MSTagSearch::pCIDImmScore		= DEFAULT_MISS_SCORE;
ScoreType MSTagSearch::pCIDH2OScore		= DEFAULT_MISS_SCORE;
ScoreType MSTagSearch::pCIDNH3Score		= DEFAULT_MISS_SCORE;
ScoreType MSTagSearch::pCID2H2OScore	= DEFAULT_MISS_SCORE;
ScoreType MSTagSearch::pCIDXLH2OScore	= DEFAULT_MISS_SCORE;
ScoreType MSTagSearch::pCIDXLNH3Score	= DEFAULT_MISS_SCORE;
ScoreType MSTagSearch::pCIDImmH2OScore	= DEFAULT_MISS_SCORE;
ScoreType MSTagSearch::pCIDImmNH3Score	= DEFAULT_MISS_SCORE;
ScoreType MSTagSearch::pETDScore		= DEFAULT_MISS_SCORE;
ScoreType MSTagSearch::pETDXLScore		= DEFAULT_MISS_SCORE;

bool MSTagSearch::crosslinking = false;
int MSTagSearch::initialPosCharge = 0;

MSTagSearch::MSTagSearch ( MSMSDataPoint* file, MSTagParameters& params ) :
	MSMSSearch ( file, params ),
	survHist ( getHistogramLimit ( params ) ),
	spectrumRetained ( true ),
	doublyChargedIons ( parentPeak.getCharge () >= 3 ),
	triplyChargedIons ( parentPeak.getCharge () >= 4 ),
	precursorCharge ( parentPeak.getCharge () )
{
	if ( !flagsSet ) {
		initFragTagFlags ( params.getBiemannParameters () );
		crosslinking = params.isCrosslinking ();
		if ( crosslinking ) {
			initXLVariables ( params.getLinkInfo () );
		}
	}
	if ( ESI_TRAP_CID_low_res ) {
		doublyChargedIons = parentPeak.getCharge () >= 2;
		triplyChargedIons = parentPeak.getCharge () >= 3;
	}
	initFragTags ();
	previousScore = -99999;
	previousNumUnmatchedIons = -1;
	previousNTerminusWt = -99999;
	previousCTerminusWt = -99999;
	previousSequence = "";
}
MSTagSearch::~MSTagSearch ()
{
}
size_t MSTagSearch::getHistogramLimit ( const MSTagParameters& params )
{
	size_t val;
	if ( params.getSearchName () == "batchtag" ) {
		if ( params.isRandomSearch () ) val = ExpectationParameters::instance ().getMaxUsedPeptides ();
		else val = 0;
	}
	else val = ULONG_MAX;
	return val;
}
ModificationTable* MSTagSearchAllowErrors::modificationTable = 0;
double MSTagSearchAllowErrors::maxParentError = 0.0;
double MSTagSearchAllowErrors::maxNegParentError = 0.0;
double MSTagSearchAllowErrors::maxPosParentError = 0.0;
MassType MSTagSearchAllowErrors::defaultNTerminusWT;
MassType MSTagSearchAllowErrors::defaultCTerminusWT;
bool MSTagSearchAllowErrors::massMods = false;
MSTagSearchAllowErrors::MSTagSearchAllowErrors ( MSMSDataPoint* file, MSTagParameters& params ) :
	MSTagSearch ( file, params )
{
	if ( modificationTable == 0 ) {
		resetMSTagSearchAllowErrors ( params.getModificationParameters () );
	}
	precursorTolerance = parentTolerance;
	if ( ModificationTable::getMassMods () ) {
		double mmTol = params.getMassModTolerance ();
		if ( parentTolerance < mmTol ) parentTolerance = mmTol;
		massMods = true;
	}
	parentMassPlusNegTolerance = parentMass + parentTolerance + maxNegParentError;
	parentMassMinusPosTolerance = parentMass - ( parentTolerance + maxPosParentError );
}
MSTagSearchAllowErrors::~MSTagSearchAllowErrors ()
{
}
void MSTagSearchAllowErrors::resetMSTagSearchAllowErrors ( const ModificationParameters& modificationParameters )
{
	if ( modificationTable != 0 ) delete modificationTable;
	modificationTable = new ModificationTable ( modificationParameters );
	maxParentError = modificationTable->getMaxParentError ();
	double modTableNegShift = modificationTable->getMostNegMassShift ();
	double modTablePosShift = modificationTable->getMostPosMassShift ();
	if ( modTableNegShift < 0.0 ) maxNegParentError = - modTableNegShift;
	if ( modTablePosShift > 0.0 ) maxPosParentError = modTablePosShift;
	defaultNTerminusWT = static_cast <MassType> (n_terminus_wt * MASS_TYPE_MULTIPLIER);
	defaultCTerminusWT = static_cast <MassType> (c_terminus_wt * MASS_TYPE_MULTIPLIER);
}
void MSTagSearch::initXLVariables ( const LinkInfo* linkInfo )
{
	xLinkMass			= massConvert ( linkInfo->getBridgeFormula ().c_str () );
	xLinkImmoniumMass	= massConvert ( linkInfo->getPCIDImmFormula ().c_str () );

	maxPScore		= linkInfo->getMaxPScore ();

	pCIDScore		= linkInfo->getPCIDScore ();
	pCIDXLScore		= linkInfo->getPCIDXLScore ();
	pCIDImmScore	= linkInfo->getPCIDImmScore ();
	pCIDH2OScore	= linkInfo->getPCIDH2OScore ();
	pCIDNH3Score	= linkInfo->getPCIDNH3Score ();
	pCID2H2OScore	= linkInfo->getPCID2H2OScore ();
	pCIDXLH2OScore	= linkInfo->getPCIDXLH2OScore ();
	pCIDXLNH3Score	= linkInfo->getPCIDXLNH3Score ();
	pCIDImmH2OScore	= linkInfo->getPCIDImmH2OScore ();
	pCIDImmNH3Score	= linkInfo->getPCIDImmNH3Score ();

	pETDScore		= linkInfo->getPETDScore ();
	pETDXLScore		= linkInfo->getPETDXLScore ();
}
void MSTagSearch::initFragTagFlags ( const BiemannParameters& bp )
{
	if ( bp.getInstrumentIonTypes () ) {
		ESI_TRAP_CID_low_res	= instInf->getFragName () == "ESI-TRAP-CID-low-res";
		ESI_Q_CID				= instInf->getFragName () == "ESI-Q-CID";
		ESI_ETD_low_res			= instInf->getFragName () == "ESI-ETD-low-res";
		ESI_ETD_high_res		= instInf->getFragName () == "ESI-ETD-high-res";
		etd						= instInf->getETDTypeFragmentation ();
	}
	internalNTerminusWt = static_cast <MassType> (internal_n_terminus_wt * MASS_TYPE_MULTIPLIER);
	nTerminusWt = static_cast <MassType> (n_terminus_wt * MASS_TYPE_MULTIPLIER);
	cTerminusWt = static_cast <MassType> (c_terminus_wt * MASS_TYPE_MULTIPLIER);
	a_flag			= bp.get_a ();
	a_nh3_flag		= bp.get_a_nh3 ();
	a_h2o_flag		= bp.get_a_h2o ();
	a_h3po4_flag	= bp.get_a_h3po4 ();
	b_flag			= bp.get_b ();
	b_h2o_flag		= bp.get_b_h2o ();
	b_plus_h2o_flag	= bp.get_b_plus_h2o ();
	b_nh3_flag		= bp.get_b_nh3 ();
	b_h3po4_flag	= bp.get_b_h3po4 ();
	b_soch4_flag	= bp.get_b_soch4 ();
	cPlus2DaFlag	= bp.getCPlus2Da ();
	cPlus1DaFlag	= bp.getCPlus1Da ();
	c_flag			= bp.get_c ();
	cMinus1DaFlag	= bp.getCMinus1Da ();
	c_ladder_flag	= bp.get_c_term_ladder ();
	d_flag			= bp.get_d ();

	v_flag			= bp.get_v ();
	w_flag			= bp.get_w ();
	x_flag			= bp.get_x ();
	Y_flag			= bp.get_Y ();
	y_flag			= bp.get_y ();
	y_nh3_flag		= bp.get_y_nh3 ();
	y_h2o_flag		= bp.get_y_h2o ();
	y_h3po4_flag	= bp.get_y_h3po4 ();
	y_soch4_flag	= bp.get_y_soch4 ();
	z_flag			= bp.get_z ();
	zPlus1DaFlag	= bp.getZPlus1Da ();
	zPlus2DaFlag	= bp.getZPlus2Da ();
	zPlus3DaFlag	= bp.getZPlus3Da ();
	n_ladder_flag	= bp.get_n_term_ladder ();

	bp2_flag		= bp.get_bp2 ();
	bp3_flag		= bp.get_bp3 ();
	bp2_nh3_flag	= bp.get_bp2_nh3 ();
	bp2_h2o_flag	= bp.get_bp2_h2o ();
	bp2_soch4_flag	= bp.get_bp2_soch4 ();
	bp2_h3po4_flag	= bp.get_bp2_h3po4 ();
	cp2_flag		= bp.get_cp2 ();
	yp2_flag		= bp.get_yp2 ();
	yp3_flag		= bp.get_yp3 ();
	yp2_nh3_flag	= bp.get_yp2_nh3 ();
	yp2_h2o_flag	= bp.get_yp2_h2o ();
	yp2_soch4_flag	= bp.get_yp2_soch4 ();
	yp2_h3po4_flag	= bp.get_yp2_h3po4 ();
	zp2_flag		= bp.get_zp2 ();
	zPlus1Dap2_flag	= bp.get_zPlus1Dap2 ();

	checkInternalIons = bp.get_internal ();

	numNIonTypes = 0;
	if ( a_flag )			a_index			= numNIonTypes++;
	if ( a_nh3_flag )		a_nh3_index		= numNIonTypes++;
	if ( a_h2o_flag )		a_h2o_index		= numNIonTypes++;
	if ( a_h3po4_flag )		a_h3po4_index	= numNIonTypes++;
	if ( b_flag || d_flag )	b_index			= numNIonTypes++;
	if ( b_h2o_flag )		b_h2o_index		= numNIonTypes++;
	if ( b_plus_h2o_flag || c_ladder_flag ) {
							b_plus_h2o_index = numNIonTypes++;
	}
	if ( b_nh3_flag )		b_nh3_index		= numNIonTypes++;
	if ( b_soch4_flag )		b_soch4_index	= numNIonTypes++;
	if ( b_h3po4_flag )		b_h3po4_index	= numNIonTypes++;
	if ( cPlus2DaFlag )		cPlus2DaIndex	= numNIonTypes++;
	if ( cPlus1DaFlag )		cPlus1DaIndex	= numNIonTypes++;
	if ( c_flag )			c_index			= numNIonTypes++;
	if ( cp2_flag )			cp2_index		= numNIonTypes++;
	if ( cMinus1DaFlag )	cMinus1DaIndex	= numNIonTypes++;
	if ( bp2_flag )			bp2_index		= numNIonTypes++;
	if ( bp3_flag )			bp3_index		= numNIonTypes++;
	if ( bp2_nh3_flag )		bp2_nh3_index	= numNIonTypes++;
	if ( bp2_h2o_flag )		bp2_h2o_index	= numNIonTypes++;
	if ( bp2_soch4_flag )	bp2_soch4_index	= numNIonTypes++;
	if ( bp2_h3po4_flag )	bp2_h3po4_index	= numNIonTypes++;

	numCIonTypes = 0;
	if ( v_flag )					v_index			= numCIonTypes++;
	if ( x_flag )					x_index			= numCIonTypes++;
	if ( Y_flag )					Y_index			= numCIonTypes++;
	if ( y_flag || n_ladder_flag || w_flag ) {
									y_index			= numCIonTypes++;
	}
	if ( y_nh3_flag )				y_nh3_index		= numCIonTypes++;
	if ( y_h2o_flag )				y_h2o_index		= numCIonTypes++;
	if ( y_soch4_flag )				y_soch4_index	= numCIonTypes++;
	if ( y_h3po4_flag )				y_h3po4_index	= numCIonTypes++;
	if ( z_flag )					z_index			= numCIonTypes++;
	if ( zp2_flag )					zp2_index		= numCIonTypes++;
	if ( zPlus1Dap2_flag )			zPlus1Dap2_index= numCIonTypes++;
	if ( zPlus1DaFlag )				zPlus1DaIndex	= numCIonTypes++;
	if ( zPlus2DaFlag )				zPlus2DaIndex	= numCIonTypes++;
	if ( zPlus3DaFlag )				zPlus3DaIndex	= numCIonTypes++;
	if ( yp2_flag )					yp2_index		= numCIonTypes++;
	if ( yp3_flag )					yp3_index		= numCIonTypes++;
	if ( yp2_nh3_flag )				yp2_nh3_index	= numCIonTypes++;
	if ( yp2_h2o_flag )				yp2_h2o_index	= numCIonTypes++;
	if ( yp2_soch4_flag )			yp2_soch4_index	= numCIonTypes++;
	if ( yp2_h3po4_flag )			yp2_h3po4_index	= numCIonTypes++;

	numInternalIonTypes = 0;
	if ( checkInternalIons ) {
		if ( a_flag )        numInternalIonTypes++;
		if ( b_flag )        numInternalIonTypes++;
		if ( b_h2o_flag )    numInternalIonTypes++;
		if ( b_nh3_flag )    numInternalIonTypes++;
	}
	checkSatelliteIons = d_flag || w_flag;
	checkNIons = ( numNIonTypes != 0 || d_flag );
	checkCIons = ( numCIonTypes != 0 || w_flag );

	checkNIonAmmoniaLoss = a_nh3_flag || b_nh3_flag;
	checkNIonPhosphorylation = a_h3po4_flag || b_h3po4_flag;
	checkNIonM = b_soch4_flag;
	checkNIonWaterLoss = a_h2o_flag || b_h2o_flag;
	checkNIonPosChargeBearing = b_plus_h2o_flag || d_flag || bp2_flag || bp2_nh3_flag || bp2_h2o_flag || bp2_soch4_flag || bp2_h3po4_flag || cp2_flag;

	checkCIonPosChargeBearing = w_flag || yp2_flag || yp2_nh3_flag || yp2_h2o_flag || yp2_soch4_flag || yp2_h3po4_flag || zp2_flag || zPlus1Dap2_flag;
	checkCIonWaterLoss = y_h2o_flag;
	checkCIonAmmoniaLoss = y_nh3_flag;
	checkCIonPhosphorylation = y_h3po4_flag;
	checkCIonM = y_soch4_flag;

	dIonExcludeMask = instInf->getDIonExcludeMask ();
	vIonExcludeMask = instInf->getVIonExcludeMask ();
	wIonExcludeMask = instInf->getWIonExcludeMask ();
	ammoniaLossMask = instInf->getAmmoniaLossMask ();
	posChargeBearingMask = instInf->getPosChargeBearingMask ();
	waterLossMask = instInf->getWaterLossMask ();
	maxInternalIonMass = static_cast <MassType> (instInf->getMaximumInternalIonMass () * MASS_TYPE_MULTIPLIER);

	char aa [] = "ACDEFGHIKLMNPQRSTUVWXYhmstuvwxy"; // Also defined in lu_mass_pep.cpp
	for ( unsigned int i = 0 ; aa [i] != '\0' ; i++ ) {
		int index = aa[i];
		unsigned int mask = aa_composition_mask [index];
		waterLossFlag [index]		= mask & waterLossMask;
		ammoniaLossFlag [index]		= mask & ammoniaLossMask;
		posChargeBearingFlag [index]= mask & posChargeBearingMask;
		oxidizedMFlag [index]		= mask & oxidized_m_mask;
		phosphorylationFlag [index]	= mask & phosphorylation_mask;
		dIonExcludeFlag [index]		= mask & dIonExcludeMask;
		vIonExcludeFlag [index]		= mask & vIonExcludeMask;
		wIonExcludeFlag [index]		= mask & wIonExcludeMask;
	}
	a_nh3_score = instInf->getALossScore ();
	a_h2o_score = instInf->getALossScore ();
	a_score = instInf->getAScore ();
	a_h3po4_score = instInf->getAPhosLossScore ();
	b_h2o_score = instInf->getBLossScore ();
	b_plus_h2o_score = instInf->getBPlusH2OScore ();
	c_ladder_score = instInf->getCLadderScore ();
	b_nh3_score = instInf->getBLossScore ();
	b_soch4_score = instInf->getBLossScore ();
	b_h3po4_score = instInf->getBPhosLossScore ();
	b_score = instInf->getBScore ();
	cPlus2DaScore = instInf->getCPlus2DaScore ();
	cPlus1DaScore = instInf->getCPlus1DaScore ();
	c_score = instInf->getCScore ();
	cMinus1DaScore = instInf->getCMinus1DaScore ();
	d_score = instInf->getDScore ();

	v_score = instInf->getVScore ();
	w_score = instInf->getWScore ();
	x_score = instInf->getXScore ();
	y_h2o_score = instInf->getYLossScore ();
	y_nh3_score = instInf->getYLossScore ();
	y_soch4_score = instInf->getYLossScore ();
	y_h3po4_score = instInf->getYPhosLossScore ();
	y_score = instInf->getYScore ();
	Y_score = instInf->getBigYScore ();
	z_score = instInf->getZScore ();
	zPlus1DaScore = instInf->getZPlus1DaScore ();
	zPlus2DaScore = instInf->getZPlus2DaScore ();
	zPlus3DaScore = instInf->getZPlus3DaScore ();
	n_ladder_score = instInf->getNLadderScore ();

	bp2_score = instInf->getBP2Score ();
	bp2_nh3_score = instInf->getBP2LossScore ();
	bp2_h2o_score = instInf->getBP2LossScore ();
	bp2_soch4_score = instInf->getBP2LossScore ();
	bp2_h3po4_score = instInf->getBP2PhosLossScore ();

	yp2_score = instInf->getYP2Score ();
	yp2_nh3_score = instInf->getYP2LossScore ();
	yp2_h2o_score = instInf->getYP2LossScore ();
	yp2_soch4_score = instInf->getYP2LossScore ();
	yp2_h3po4_score = instInf->getYP2PhosLossScore ();

	internal_a_score = instInf->getInternalAScore ();
	internal_b_score = instInf->getInternalBScore ();
	internal_loss_score = instInf->getInternalLossScore ();

	unmatched_score = instInf->getUnmatchedScore ();

	UnmatchedCompositionSearch::setImmoniumScore ( instInf->getImmoniumScore () );
	flagsSet = true;
}
void MSTagSearch::initFragTags ()
{
	numPeaks = peaks.size ();
	fragTolerance.resize ( numPeaks );
	nIonFragTag.resize ( numPeaks );
	cIonFragTag.resize ( numPeaks );

	maxNFragTag.resize ( numPeaks );
	maxCFragTag.resize ( numPeaks );
	minNFragTag.resize ( numPeaks );
	minCFragTag.resize ( numPeaks );

	int i, j;
	for ( i = 0 ; i < numPeaks ; i++ ) {
		nIonFragTag [i].resize ( numNIonTypes );
		cIonFragTag [i].resize ( numCIonTypes );

		double fragMass = peaks [i]->getMass ();
		double fragMassCharge2 = mOverZToMPlusH ( fragMass, 2, true );
		double fragMassCharge3 = mOverZToMPlusH ( fragMass, 3, true );
		fragTolerance [i] = static_cast <MassType> (peaks [i]->getTolerance () * MASS_TYPE_MULTIPLIER);

		if ( a_flag )			nIonFragTag [i][a_index]		= static_cast <MassType> (( fragMass + a_tag_offset ) * MASS_TYPE_MULTIPLIER);
		if ( a_h2o_flag )		nIonFragTag [i][a_h2o_index]	= static_cast <MassType> (( fragMass + a_h2o_tag_offset ) * MASS_TYPE_MULTIPLIER);
		if ( a_nh3_flag )		nIonFragTag [i][a_nh3_index]	= static_cast <MassType> (( fragMass + a_nh3_tag_offset ) * MASS_TYPE_MULTIPLIER);
		if ( a_h3po4_flag )		nIonFragTag [i][a_h3po4_index]	= static_cast <MassType> (( fragMass + a_h3po4_tag_offset ) * MASS_TYPE_MULTIPLIER);
		if ( b_flag || d_flag )	nIonFragTag [i][b_index]		= static_cast <MassType> (( fragMass + b_tag_offset ) * MASS_TYPE_MULTIPLIER);
		if ( b_h2o_flag )		nIonFragTag [i][b_h2o_index]	= static_cast <MassType> (( fragMass + b_h2o_tag_offset ) * MASS_TYPE_MULTIPLIER);
		if ( b_plus_h2o_flag || c_ladder_flag )	{
								nIonFragTag [i][b_plus_h2o_index]= static_cast <MassType> (( fragMass + b_plus_h2o_tag_offset ) * MASS_TYPE_MULTIPLIER);
		}
		if ( b_nh3_flag )		nIonFragTag [i][b_nh3_index]	= static_cast <MassType> (( fragMass + b_nh3_tag_offset ) * MASS_TYPE_MULTIPLIER);
		if ( b_soch4_flag )		nIonFragTag [i][b_soch4_index]	= static_cast <MassType> (( fragMass + b_soch4_tag_offset ) * MASS_TYPE_MULTIPLIER);
		if ( b_h3po4_flag )		nIonFragTag [i][b_h3po4_index]	= static_cast <MassType> (( fragMass + b_h3po4_tag_offset ) * MASS_TYPE_MULTIPLIER);
		if ( cPlus2DaFlag )		nIonFragTag [i][cPlus2DaIndex]	= static_cast <MassType> (( fragMass + cPlus2DaTagOffset ) * MASS_TYPE_MULTIPLIER);
		if ( cPlus1DaFlag )		nIonFragTag [i][cPlus1DaIndex]	= static_cast <MassType> (( fragMass + cPlus1DaTagOffset ) * MASS_TYPE_MULTIPLIER);
		if ( c_flag )			nIonFragTag [i][c_index]		= static_cast <MassType> (( fragMass + c_tag_offset ) * MASS_TYPE_MULTIPLIER);
		if ( cp2_flag )			nIonFragTag [i][cp2_index]		= static_cast <MassType> (( fragMassCharge2 + c_tag_offset ) * MASS_TYPE_MULTIPLIER);
		if ( cMinus1DaFlag )	nIonFragTag [i][cMinus1DaIndex]	= static_cast <MassType> (( fragMass + cMinus1DaTagOffset ) * MASS_TYPE_MULTIPLIER);
		if ( bp2_flag )			nIonFragTag [i][bp2_index]		= static_cast <MassType> (( fragMassCharge2 + b_tag_offset ) * MASS_TYPE_MULTIPLIER);
		if ( bp3_flag )			nIonFragTag [i][bp3_index]		= static_cast <MassType> (( fragMassCharge3 + b_tag_offset ) * MASS_TYPE_MULTIPLIER);
		if ( bp2_nh3_flag )		nIonFragTag [i][bp2_nh3_index]	= static_cast <MassType> (( fragMassCharge2 + b_nh3_tag_offset ) * MASS_TYPE_MULTIPLIER);
		if ( bp2_h2o_flag )		nIonFragTag [i][bp2_h2o_index]	= static_cast <MassType> (( fragMassCharge2 + b_h2o_tag_offset ) * MASS_TYPE_MULTIPLIER);
		if ( bp2_soch4_flag )	nIonFragTag [i][bp2_soch4_index]= static_cast <MassType> (( fragMassCharge2 + b_soch4_tag_offset ) * MASS_TYPE_MULTIPLIER);
		if ( bp2_h3po4_flag )	nIonFragTag [i][bp2_h3po4_index]= static_cast <MassType> (( fragMassCharge2 + b_h3po4_tag_offset ) * MASS_TYPE_MULTIPLIER);

		if ( v_flag )			cIonFragTag [i][v_index]		= static_cast <MassType> (( fragMass + v_tag_offset ) * MASS_TYPE_MULTIPLIER);
		if ( x_flag )			cIonFragTag [i][x_index]		= static_cast <MassType> (( fragMass + x_tag_offset ) * MASS_TYPE_MULTIPLIER);
		if ( y_flag || w_flag || n_ladder_flag ) {
								cIonFragTag [i][y_index]		= static_cast <MassType> (( fragMass + y_tag_offset ) * MASS_TYPE_MULTIPLIER);
		}
		if ( y_nh3_flag )		cIonFragTag [i][y_nh3_index]	= static_cast <MassType> (( fragMass + y_nh3_tag_offset ) * MASS_TYPE_MULTIPLIER);
		if ( y_h2o_flag )		cIonFragTag [i][y_h2o_index]	= static_cast <MassType> (( fragMass + y_h2o_tag_offset ) * MASS_TYPE_MULTIPLIER);
		if ( y_soch4_flag )		cIonFragTag [i][y_soch4_index]	= static_cast <MassType> (( fragMass + y_soch4_tag_offset ) * MASS_TYPE_MULTIPLIER);
		if ( y_h3po4_flag )		cIonFragTag [i][y_h3po4_index]	= static_cast <MassType> (( fragMass + y_h3po4_tag_offset ) * MASS_TYPE_MULTIPLIER);
		if ( Y_flag )			cIonFragTag [i][Y_index]		= static_cast <MassType> (( fragMass + Y_tag_offset ) * MASS_TYPE_MULTIPLIER);
		if ( z_flag )			cIonFragTag [i][z_index]		= static_cast <MassType> (( fragMass + z_tag_offset ) * MASS_TYPE_MULTIPLIER);
		if ( zp2_flag )			cIonFragTag [i][zp2_index]		= static_cast <MassType> (( fragMassCharge2 + z_tag_offset ) * MASS_TYPE_MULTIPLIER);
		if ( zPlus1Dap2_flag )	cIonFragTag [i][zPlus1Dap2_index]=static_cast <MassType> (( fragMassCharge2 + zPlus1DaTagOffset ) * MASS_TYPE_MULTIPLIER);
		if ( zPlus1DaFlag )		cIonFragTag [i][zPlus1DaIndex]	= static_cast <MassType> (( fragMass + zPlus1DaTagOffset ) * MASS_TYPE_MULTIPLIER);
		if ( zPlus2DaFlag )		cIonFragTag [i][zPlus2DaIndex]	= static_cast <MassType> (( fragMass + zPlus2DaTagOffset ) * MASS_TYPE_MULTIPLIER);
		if ( zPlus3DaFlag )		cIonFragTag [i][zPlus3DaIndex]	= static_cast <MassType> (( fragMass + zPlus3DaTagOffset ) * MASS_TYPE_MULTIPLIER);
		if ( yp2_flag )			cIonFragTag [i][yp2_index]		= static_cast <MassType> (( fragMassCharge2 + y_tag_offset ) * MASS_TYPE_MULTIPLIER);
		if ( yp3_flag )			cIonFragTag [i][yp3_index]		= static_cast <MassType> (( fragMassCharge3 + y_tag_offset ) * MASS_TYPE_MULTIPLIER);
		if ( yp2_nh3_flag )		cIonFragTag [i][yp2_nh3_index]	= static_cast <MassType> (( fragMassCharge2 + y_nh3_tag_offset ) * MASS_TYPE_MULTIPLIER);
		if ( yp2_h2o_flag )		cIonFragTag [i][yp2_h2o_index]	= static_cast <MassType> (( fragMassCharge2 + y_h2o_tag_offset ) * MASS_TYPE_MULTIPLIER);
		if ( yp2_soch4_flag )	cIonFragTag [i][yp2_soch4_index]= static_cast <MassType> (( fragMassCharge2 + y_soch4_tag_offset ) * MASS_TYPE_MULTIPLIER);
		if ( yp2_h3po4_flag )	cIonFragTag [i][yp2_h3po4_index]= static_cast <MassType> (( fragMassCharge2 + y_h3po4_tag_offset ) * MASS_TYPE_MULTIPLIER);
	}
	MassType maxDSubstituent = static_cast <MassType> (getMaxDSubstituent () * MASS_TYPE_MULTIPLIER);
	MassType maxWSubstituent = static_cast <MassType> (getMaxWSubstituent () * MASS_TYPE_MULTIPLIER);

	bool p2_flag =	bp2_flag || bp2_nh3_flag || bp2_h2o_flag || bp2_soch4_flag || bp2_h3po4_flag;
	p2_flag |=		cp2_flag;
	p2_flag |=		yp2_flag || yp2_nh3_flag || yp2_h2o_flag || yp2_soch4_flag || yp2_h3po4_flag;
	p2_flag |=		zp2_flag || zPlus1Dap2_flag;
	bool p3_flag =	bp3_flag || yp3_flag;

	for ( i = 0 ; i < numPeaks ; i++ ) {
		maxNFragTag [i] = static_cast <MassType> (0.0);
		maxCFragTag [i] = static_cast <MassType> (0.0);
		minNFragTag [i] = static_cast <MassType> (1E+9);
		minCFragTag [i] = static_cast <MassType> (1E+9);
		MassType ttol = p2_flag ? fragTolerance [i] + fragTolerance [i] : fragTolerance [i];
		ttol = p3_flag ? ttol + fragTolerance [i] : ttol;
		for ( j = 0 ; j < numNIonTypes ; j++ ) {
			maxNFragTag [i] = genMax ( maxNFragTag [i], nIonFragTag [i][j] + ttol );
			minNFragTag [i] = genMin ( minNFragTag [i], nIonFragTag [i][j] - ttol );
		}
		for ( j = 0 ; j < numCIonTypes ; j++ ) {
			maxCFragTag [i] = genMax ( maxCFragTag [i], cIonFragTag [i][j] + ttol );
			minCFragTag [i] = genMin ( minCFragTag [i], cIonFragTag [i][j] - ttol );
		}
		if ( checkSatelliteIons ) {
			maxNFragTag [i] += maxDSubstituent;
			maxCFragTag [i] += maxWSubstituent;
		}
	}

	for ( i = 0 ; i < numPeaks ; i++ ) {
		maxNFragTagInData = ( i == 0 ) ? maxNFragTag [0] : genMax ( maxNFragTagInData, maxNFragTag [i] );
		maxCFragTagInData = ( i == 0 ) ? maxCFragTag [0] : genMax ( maxCFragTagInData, maxCFragTag [i] );
		minNFragTagInData = ( i == 0 ) ? minNFragTag [0] : genMin ( minNFragTagInData, minNFragTag [i] );
		minCFragTagInData = ( i == 0 ) ? minCFragTag [0] : genMin ( minCFragTagInData, minCFragTag [i] );
	}
	ionMatched.resize ( numPeaks );
}
bool MSTagSearch::doMatch ( const string& peptide, bool nTermPeptide, bool cTermPeptide, double mol_wt, TagMatchVector& tagMatch, const ScoreType& minScore )
{
	if ( compositionSearch && compositionSearch->doCompositionSearch ( peptide ) == false );
	else {
		int unmatched = 0;
		ScoreType score = 0;
		fragmentMatch ( peptide, unmatched, score );
		if ( !survHist.addValue ( score / SCORE_TYPE_MULTIPLIER ) ) {
			spectrumRetained = false;
		}
		tagMatch.clear ();
		if ( score >= minScore ) {
			tagMatch.push_back ( TagMatch ( &parentPeak, unmatched, score, PeptideSequence ( peptide ) ) );
			return true;
		}
	}
	return false;
}
void MSTagSearchAllowErrors::resetNextMods ()
{
	if ( modificationTable ) modificationTable->resetNextMods ();
}
bool MSTagSearchAllowErrors::doMatch ( const string& peptide, bool nTermPeptide, bool cTermPeptide, double molWt, TagMatchVector& tagMatch, const ScoreType& minScore )
{
	ScoreType bestScore = -9999999;
	tagMatch.clear ();
	if ( genAbsDiff ( molWt, parentMass ) >= precursorTolerance ) {	// Handles the allow errors condition, i.e. parent mass is shifted
		double diff = parentMass - molWt;
		PeptideSequenceVector mutatedSequences;
		if ( modificationTable->getMutatedSequences ( diff - parentTolerance, diff + parentTolerance, peptide, nTermPeptide, cTermPeptide, mutatedSequences, precursorCharge ) ) {
			for ( PeptideSequenceVectorSizeType i = 0 ; i < mutatedSequences.size () ; i++ ) {
				const PeptideSequence& ps = mutatedSequences [i];
				int mmod = ps.applyModifications ( nTerminusWt, cTerminusWt, diff );
				if ( massMods && !mmod ) {		// Non-mass mods in a mass mod search
					if ( genAbsDiff ( ps.getMW (), parentMass ) >= precursorTolerance ) {
						continue;
					}
				}
				if ( compositionSearch && compositionSearch->doCompositionSearch ( ps.getSequence () ) == false );
				else {
					previousSequence = "";		// Make sure that fragmentMatch is always invoked
					int unmatched = 0;
					ScoreType score = 0;
					fragmentMatch ( ps, unmatched, score, molWt, diff );
					if ( score >= minScore ) {
						tagMatch.push_back ( TagMatch ( &parentPeak, unmatched, score, ps ) );
					}
					bestScore = genMax ( bestScore, score );
				}
				nTerminusWt = defaultNTerminusWT;
				cTerminusWt = defaultCTerminusWT;
			}
		}
	}
	else { // Handles the no errors condition, i.e. parent mass matches
		if ( compositionSearch && compositionSearch->doCompositionSearch ( peptide ) == false );
		else {
			int unmatched = 0;
			ScoreType score = 0;
			fragmentMatch ( peptide, unmatched, score );
			if ( score >= minScore ) {
				tagMatch.push_back ( TagMatch ( &parentPeak, unmatched, score, peptide ) );
			}
			bestScore = genMax ( bestScore, score );
		}
	}
	if ( bestScore != -9999999 ) {
		if ( !survHist.addValue ( bestScore / SCORE_TYPE_MULTIPLIER ) ) {
			spectrumRetained = false;
		}
	}
	return !tagMatch.empty ();
}
void MSTagSearch::xLinkMatchCID ( const PeptideSequence& ps, double molWt, double diff )
{
	double p1		= ps.getP1 ();
	double p1_h2o	= p1-h2_o;
	double p1_nh3	= p1-n_h3;
	double p1_2h2o	= p1_h2o-h2_o;
	double p1xL		= p1+xLinkMass;
	double p1xL_h2o	= p1xL-h2_o;
	double p1xL_nh3	= p1xL-n_h3;
	double p1imm	= p1xL+xLinkImmoniumMass;			// Biggest mass is the immonium mass
	double p1imm_h2o= p1imm-h2_o;
	double p1imm_nh3= p1imm-n_h3;

	double q1xL		= molWt+diff-p1+h1;
	double q1xL_h2o	= q1xL-h2_o;
	double q1xL_nh3	= q1xL-n_h3;
	double q1		= q1xL-xLinkMass;
	double q1_h2o	= q1-h2_o;
	double q1_nh3	= q1-n_h3;
	double q1_2h2o	= q1_h2o-h2_o;
	double q1imm	= q1xL+xLinkImmoniumMass;
	double q1imm_h2o= q1imm-h2_o;
	double q1imm_nh3= q1imm-n_h3;

	bool pFlag = p1 > q1;

	ScoreType totalPScore = 0;
	for ( int i = numPeaks ; i-- ; ) {
		double m = peaks [i]->getMass ();
		double t = peaks [i]->getTolerance ();
		ScoreType& matched = ionMatched [i];
		if ( pFlag ) {
			if ( m < q1_2h2o - t ) break;
			if ( m > p1imm + t ) continue;
		}
		else {
			if ( m < p1_2h2o - t ) break;
			if ( m > q1imm + t ) continue;
		}
		bool hit = true;
		if (		pCIDXLScore		&& (genAbsDiff(p1xL,m)<t		|| genAbsDiff (q1xL,m)<t ))		matched = pCIDXLScore;
		else if (	pCIDImmScore	&& (genAbsDiff(p1imm,m)<t		|| genAbsDiff (q1imm,m)<t))		matched = pCIDImmScore;
		else if (	pCIDScore		&& (genAbsDiff(p1,m)<t			|| genAbsDiff (q1,m)<t))		matched = pCIDScore;
		else if (	pCIDXLH2OScore	&& (genAbsDiff(p1xL_h2o,m)<t	|| genAbsDiff (q1xL_h2o,m)<t))	matched = pCIDXLH2OScore;
		else if (	pCIDXLNH3Score	&& (genAbsDiff(p1xL_nh3,m)<t	|| genAbsDiff (q1xL_nh3,m)<t))	matched = pCIDXLNH3Score;
		else if (	pCIDImmH2OScore	&& (genAbsDiff(p1imm_h2o,m)<t	|| genAbsDiff (q1imm_h2o,m)<t))	matched = pCIDImmH2OScore;
		else if (	pCIDImmNH3Score	&& (genAbsDiff(p1imm_nh3,m)<t	|| genAbsDiff (q1imm_nh3,m)<t))	matched = pCIDImmNH3Score;
		else if (	pCIDH2OScore	&& (genAbsDiff(p1_h2o,m)<t		|| genAbsDiff (q1_h2o,m)<t))	matched = pCIDH2OScore;
		else if (	pCIDNH3Score	&& (genAbsDiff(p1_nh3,m)<t		|| genAbsDiff (q1_nh3,m)<t))	matched = pCIDNH3Score;
		else if (	pCID2H2OScore	&& (genAbsDiff(p1_2h2o,m)<t		|| genAbsDiff (q1_2h2o,m)<t))	matched = pCID2H2OScore;
		else hit = false;
		if ( hit ) {
			totalPScore += matched;
			if ( totalPScore >= maxPScore ) {
				matched -= totalPScore - maxPScore;
				break;
			}
		}
	}
}
void MSTagSearch::xLinkMatchETD ( const PeptideSequence& ps, double molWt, double diff )
{
	double p1	= ps.getP1 ();
	double p1xL	= p1+xLinkMass+h1;
	double q1xL	= molWt+diff-p1+h2-ps.getNLossMass ();
	double q1	= q1xL-xLinkMass-h1;
	bool pFlag = p1 > q1;

	ScoreType totalPScore = 0;
	for ( int i = numPeaks ; i-- ; ) {
		double m = peaks [i]->getMass ();
		double t = peaks [i]->getTolerance ();
		ScoreType& matched = ionMatched [i];
		if ( pFlag ) {
			if ( m < q1 - t ) break;
			if ( m > p1xL + t ) continue;
		}
		else {
			if ( m < p1 - t ) break;
			if ( m > q1xL + t ) continue;
		}
		bool hit = true;
		if (		pETDXLScore && ( genAbsDiff(p1xL,m)<t	|| genAbsDiff (q1xL,m)<t ) )matched = pETDXLScore;
		else if (	pETDScore &&	(genAbsDiff(p1,m)<t		|| genAbsDiff (q1,m)<t ) )	matched = pETDScore;
		else hit = false;
		if ( hit ) {
			totalPScore += matched;
			if ( totalPScore >= maxPScore ) {
				matched -= totalPScore - maxPScore;
				break;
			}
		}
	}
}
void MSTagSearch::xLinkMatchCID2 ( const PeptideSequence& ps )
{
	double p1		= ps.getP1 ();
	double p1_h2o	= p1-h2_o;
	double p1_nh3	= p1-n_h3;
	double p1_2h2o	= p1_h2o-h2_o;
	double p1xL		= p1+xLinkMass;
	double p1xL_h2o	= p1xL-h2_o;
	double p1xL_nh3	= p1xL-n_h3;
	double p1imm	= p1xL+xLinkImmoniumMass;
	double p1imm_h2o= p1imm-h2_o;
	double p1imm_nh3= p1imm-n_h3;

	ScoreType totalPScore = 0;
	for ( int i = numPeaks ; i-- ; ) {
		double m = peaks [i]->getMass ();
		double t = peaks [i]->getTolerance ();
		ScoreType& matched = ionMatched [i];
		if ( m < p1_h2o - t ) break;
		if ( m > p1imm + t ) continue;
		bool hit = true;
		if ( pCIDXLScore && genAbsDiff(p1xL,m)<t )				matched = pCIDXLScore;
		else if ( pCIDImmScore && genAbsDiff(p1imm,m)<t )		matched = pCIDImmScore;
		else if ( pCIDScore && genAbsDiff(p1,m)<t )				matched = pCIDScore;
		else if ( pCIDXLH2OScore && genAbsDiff(p1xL_h2o,m)<t )	matched = pCIDXLH2OScore;
		else if ( pCIDXLNH3Score && genAbsDiff(p1xL_nh3,m)<t )	matched = pCIDXLNH3Score;
		else if ( pCIDImmH2OScore && genAbsDiff(p1imm_h2o,m)<t )matched = pCIDImmH2OScore;
		else if ( pCIDImmNH3Score && genAbsDiff(p1imm_nh3,m)<t )matched = pCIDImmNH3Score;
		else if ( pCIDH2OScore && genAbsDiff(p1_h2o,m)<t )		matched = pCIDH2OScore;
		else if ( pCIDNH3Score && genAbsDiff(p1_nh3,m)<t )		matched = pCIDNH3Score;
		else if ( pCID2H2OScore && genAbsDiff(p1_2h2o,m)<t )	matched = pCID2H2OScore;
		else hit = false;
		if ( hit ) {
			totalPScore += matched;
			if ( totalPScore >= maxPScore ) {
				matched -= totalPScore - maxPScore;
				break;
			}
		}
	}
}
void MSTagSearch::xLinkMatchETD2 ( const PeptideSequence& ps )
{
	double p1	= ps.getP1 ();
	double p1xL	= p1+xLinkMass+h1;
	ScoreType totalPScore = 0;
	for ( int i = numPeaks ; i-- ; ) {
		double m = peaks [i]->getMass ();
		double t = peaks [i]->getTolerance ();
		ScoreType& matched = ionMatched [i];
		if ( m < p1 - t ) break;
		if ( m > p1xL + t ) continue;
		bool hit = true;
		if ( pETDScore && genAbsDiff(p1,m)<t )			matched = pETDScore;
		else if ( pETDXLScore && genAbsDiff(p1xL,m)<t )	matched = pETDXLScore;
		else hit = false;
		if ( hit ) {
			totalPScore += matched;
			if ( totalPScore >= maxPScore ) {
				matched -= totalPScore - maxPScore;
				break;
			}
		}
	}
}
void MSTagSearch::fragmentMatch ( const string& sequence, int& numUnmatchedIons, ScoreType& score, double molWt, double diff )
{
	if ( sequence == previousSequence && previousNTerminusWt == nTerminusWt && previousCTerminusWt == cTerminusWt ) {
		numUnmatchedIons = previousNumUnmatchedIons;
		score = previousScore;
		return;
	}
	/*Composition Ion filter*/
	fill ( ionMatched.begin (), ionMatched.end (), 0 );

	fMatch ( sequence, numUnmatchedIons, score );

	previousScore = score;
	previousNumUnmatchedIons = numUnmatchedIons;
	previousNTerminusWt = nTerminusWt;
	previousCTerminusWt = cTerminusWt;
	previousSequence = sequence;
}
void MSTagSearch::fragmentMatch ( const PeptideSequence& ps, int& numUnmatchedIons, ScoreType& score, double molWt, double diff )
{
	if ( ps.getSequence () == previousSequence && previousNTerminusWt == nTerminusWt && previousCTerminusWt == cTerminusWt ) {
		numUnmatchedIons = previousNumUnmatchedIons;
		score = previousScore;
		return;
	}
	/*Composition Ion filter*/
	fill ( ionMatched.begin (), ionMatched.end (), 0 );

	if ( crosslinking ) {
		if ( ps.isMassMod () ) {
			if ( etd )	xLinkMatchETD ( ps, molWt, diff );
			else		xLinkMatchCID ( ps, molWt, diff );
			initialPosCharge = 1;
		}
		else {
			initialPosCharge = 0;
		}
	}
	fMatch ( ps.getSequence (), numUnmatchedIons, score );
	previousScore = score;
	previousNumUnmatchedIons = numUnmatchedIons;
	previousNTerminusWt = nTerminusWt;
	previousCTerminusWt = cTerminusWt;
	previousSequence = ps.getSequence ();
}
void MSTagSearch::fragmentMatch2 ( const PeptideSequence& ps, int& numUnmatchedIons, ScoreType& score, bool reset )
{
	if ( reset ) {
		fill ( ionMatched.begin (), ionMatched.end (), 0 );
	}
	if ( etd )	xLinkMatchETD2 ( ps );
	else		xLinkMatchCID2 ( ps );
	initialPosCharge = 1;
	fMatch ( ps.getSequence (), numUnmatchedIons, score );
}
void MSTagSearch::fMatch ( const string& sequence, int& numUnmatchedIons, ScoreType& score )
{
	if ( ESI_TRAP_CID_low_res ) {
		doNIonsESI_TRAP_CID_low_res ( sequence );
		doCIonsESI_TRAP_CID_low_res ( sequence );
	}
	else if ( ESI_Q_CID ) {
		if ( unmatchedCompositionSearch ) unmatchedCompositionSearch->doCompositionSearch ( sequence, ionMatched );
		doNIonsESI_Q_CID ( sequence );
		doCIonsESI_Q_CID ( sequence );
		doInternalIons ( sequence );
	}
	else if ( ESI_ETD_low_res ) {
		doNIonsESI_ETD_low_res ( sequence );
		doCIonsESI_ETD_low_res ( sequence );
	}
	else if ( ESI_ETD_high_res ) {
		doNIonsESI_ETD_high_res ( sequence );
		doCIonsESI_ETD_high_res ( sequence );
	}
	else {
		if ( unmatchedCompositionSearch ) unmatchedCompositionSearch->doCompositionSearch ( sequence, ionMatched );
		if ( checkNIons ) doNIons ( sequence );
		if ( checkCIons ) doCIons ( sequence );
		if ( checkInternalIons ) doInternalIons ( sequence );
	}
	for ( int i = 0 ; i < numPeaks ; i++ ) {
		if ( ionMatched [i] == 0 ) {
			numUnmatchedIons++;
			score += unmatched_score;
		}
		else score += ionMatched [i];
	}
}
void MSTagSearch::doNIons ( const string& sequence )
{
	int nIonAmmoniaLoss = 0;
	int nIonPosChargeBearing = initialPosCharge;
	int nIonM = 0;
	int nIonPhosphorylation = 0;
	int nIonWaterLoss = 0;
	int dIonExclude = 0;
	MassType nIon = nTerminusWt;
	int start = 0;
	int len = sequence.length () - 1;
	for ( int k = 0 ; k < len ; k++ ) {
		char aa = sequence [k];
		nIon += aaArrayMT [aa];
		if ( nIon > maxNFragTagInData ) break;

		if ( checkNIonAmmoniaLoss )			nIonAmmoniaLoss |= ammoniaLossFlag [aa];
		if ( checkNIonPosChargeBearing )	nIonPosChargeBearing |= posChargeBearingFlag [aa];
		if ( checkNIonWaterLoss )			nIonWaterLoss |= waterLossFlag [aa];
		if ( checkNIonM )					nIonM |= oxidizedMFlag [aa];
		if ( checkNIonPhosphorylation )		nIonPhosphorylation |= phosphorylationFlag [aa];
		if ( d_flag )						dIonExclude = dIonExcludeFlag [aa];
		if ( nIon < minNFragTagInData || k == 0 ) continue;
		for ( int i = start ; i < numPeaks ; i++ ) {
			if ( nIon < minNFragTag [i] ) break;
			if ( nIon > maxNFragTag [i] ) {
				start++;
				continue;
			}
			MassTypeVector& fragTag = nIonFragTag [i];
			ScoreType& matched = ionMatched [i];
			MassType fragTol = fragTolerance [i];
			MassType fragTol2 = fragTol + fragTol;
			if ( b_flag && genAbsDiff ( nIon, fragTag [b_index] ) < fragTol ) {
				matched = genMax ( b_score, matched );
			}
			else if ( doublyChargedIons && bp2_flag && nIonPosChargeBearing && genAbsDiff ( nIon, fragTag [bp2_index] ) < fragTol2 ) {
				matched = genMax ( bp2_score, matched );
			}
			else if ( a_flag && genAbsDiff ( nIon, fragTag [a_index] ) < fragTol ) {
				matched = genMax ( a_score, matched );
			}
			else if ( a_nh3_flag && nIonAmmoniaLoss && genAbsDiff ( nIon, fragTag [a_nh3_index] ) < fragTol ) {
				matched = genMax ( a_nh3_score, matched );
			}
			else if ( a_h2o_flag && nIonWaterLoss && genAbsDiff ( nIon, fragTag [a_h2o_index] ) < fragTol ) {
				matched = genMax ( a_h2o_score, matched );
			}
			else if ( b_nh3_flag && nIonAmmoniaLoss && genAbsDiff ( nIon, fragTag [b_nh3_index] ) < fragTol ) {
				matched = genMax ( b_nh3_score, matched );
			}
			else if ( doublyChargedIons && bp2_nh3_flag && nIonAmmoniaLoss && nIonPosChargeBearing && genAbsDiff ( nIon, fragTag [bp2_nh3_index] ) < fragTol2 ) {
				matched = genMax ( bp2_nh3_score, matched );
			}
			else if ( b_h2o_flag && nIonWaterLoss && genAbsDiff ( nIon, fragTag [b_h2o_index] ) < fragTol ) {
				matched = genMax ( b_h2o_score, matched );
			}
			else if ( doublyChargedIons && bp2_h2o_flag && nIonWaterLoss && nIonPosChargeBearing && genAbsDiff ( nIon, fragTag [bp2_h2o_index] ) < fragTol2 ) {
				matched = genMax ( bp2_h2o_score, matched );
			}
			else if ( b_plus_h2o_flag && k >= len - 3 && nIonPosChargeBearing && genAbsDiff ( nIon, fragTag [b_plus_h2o_index] ) < fragTol ) {
				matched = genMax ( b_plus_h2o_score, matched );
			}
			else if ( c_ladder_flag && genAbsDiff ( nIon, fragTag [b_plus_h2o_index] ) < fragTol ) {
				matched = genMax ( c_ladder_score, matched );
			}
			else if ( b_soch4_flag && nIonM && genAbsDiff ( nIon, fragTag [b_soch4_index] ) < fragTol ) {
				matched = genMax ( b_soch4_score, matched );
			}
			else if ( doublyChargedIons && bp2_soch4_flag && nIonM && nIonPosChargeBearing && genAbsDiff ( nIon, fragTag [bp2_soch4_index] ) < fragTol2 ) {
				matched = genMax ( bp2_soch4_score, matched );
			}
			else if ( a_h3po4_flag && nIonPhosphorylation && genAbsDiff ( nIon, fragTag [a_h3po4_index] ) < fragTol ) {
				matched = genMax ( a_h3po4_score, matched );
			}
			else if ( b_h3po4_flag && nIonPhosphorylation && genAbsDiff ( nIon, fragTag [b_h3po4_index] ) < fragTol ) {
				matched = genMax ( b_h3po4_score, matched );
			}
			else if ( doublyChargedIons && bp2_h3po4_flag && nIonPhosphorylation && nIonPosChargeBearing && genAbsDiff ( nIon, fragTag [bp2_h3po4_index] ) < fragTol2 ) {
				matched = genMax ( bp2_h3po4_score, matched );
			}
			else if ( c_flag && genAbsDiff ( nIon, fragTag [c_index] ) < fragTol ) {
				matched = genMax ( c_score, matched );
			}
			else if ( cMinus1DaFlag && genAbsDiff ( nIon, fragTag [cMinus1DaIndex] ) < fragTol ) {
				matched = genMax ( cMinus1DaScore, matched );
			}
			else if ( cPlus1DaFlag && genAbsDiff ( nIon, fragTag [cPlus1DaIndex] ) < fragTol ) {
				matched = genMax ( cPlus1DaScore, matched );
			}
			else if ( cPlus2DaFlag && genAbsDiff ( nIon, fragTag [cPlus2DaIndex] ) < fragTol ) {
				matched = genMax ( cPlus2DaScore, matched );
			}
			else if ( d_flag && nIonPosChargeBearing && dIonExclude == 0 && substituentDaMT [aa] ) {
				if ( genAbsDiff ( nIon, (MassType)(fragTag [b_index] - substituentDaMT [aa]) ) < fragTol ) {
					matched = genMax ( d_score, matched );
				}
				else {
					if ( substituentDbMT [aa] && genAbsDiff ( nIon, (MassType)(fragTag [b_index] - substituentDbMT [aa]) ) < fragTol ) {
						matched = genMax ( d_score, matched );
					}
				}
			}
		}
	}
}
void MSTagSearch::doNIonsESI_TRAP_CID_low_res ( const string& sequence )
{
	int nIonAmmoniaLoss = 0;
	int nIonPosChargeBearing = initialPosCharge;
	int nIonM = 0;
	int nIonPhosphorylation = 0;
	int nIonWaterLoss = 0;
	MassType nIon = nTerminusWt;
	int start = 0;
	int len = sequence.length () - 1;
	bool basicN = posChargeBearingFlag [sequence[0]] != 0;
	bool basicC = posChargeBearingFlag [sequence[len]] != 0;
	switch ( precursorCharge ) {
		case 1:
			a_score			= basicN ? ( basicC ? ESI_TRAP_CID_low_res_a_scoreZ1NC			: ESI_TRAP_CID_low_res_a_scoreZ1N )			: ( basicC ? ESI_TRAP_CID_low_res_a_scoreZ1C		: ESI_TRAP_CID_low_res_a_scoreZ1 );
			b_score			= basicN ? ( basicC ? ESI_TRAP_CID_low_res_b_scoreZ1NC			: ESI_TRAP_CID_low_res_b_scoreZ1N )			: ( basicC ? ESI_TRAP_CID_low_res_b_scoreZ1C		: ESI_TRAP_CID_low_res_b_scoreZ1 );
			b_nh3_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_b_nh3_scoreZ1NC		: ESI_TRAP_CID_low_res_b_nh3_scoreZ1N )		: ( basicC ? ESI_TRAP_CID_low_res_b_nh3_scoreZ1C	: ESI_TRAP_CID_low_res_b_nh3_scoreZ1 );
			b_h2o_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_b_h2o_scoreZ1NC		: ESI_TRAP_CID_low_res_b_h2o_scoreZ1N )		: ( basicC ? ESI_TRAP_CID_low_res_b_h2o_scoreZ1C	: ESI_TRAP_CID_low_res_b_h2o_scoreZ1 );
			b_soch4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_b_soch4_scoreZ1NC	: ESI_TRAP_CID_low_res_b_soch4_scoreZ1N )	: ( basicC ? ESI_TRAP_CID_low_res_b_soch4_scoreZ1C	: ESI_TRAP_CID_low_res_b_soch4_scoreZ1 );
			b_h3po4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_b_h3po4_scoreZ1NC	: ESI_TRAP_CID_low_res_b_h3po4_scoreZ1N )	: ( basicC ? ESI_TRAP_CID_low_res_b_h3po4_scoreZ1C	: ESI_TRAP_CID_low_res_b_h3po4_scoreZ1 );
			break;
		case 2:
			a_score			= basicN ? ( basicC ? ESI_TRAP_CID_low_res_a_scoreZ2NC			: ESI_TRAP_CID_low_res_a_scoreZ2N )			: ( basicC ? ESI_TRAP_CID_low_res_a_scoreZ2C			: ESI_TRAP_CID_low_res_a_scoreZ2 );
			b_score			= basicN ? ( basicC ? ESI_TRAP_CID_low_res_b_scoreZ2NC			: ESI_TRAP_CID_low_res_b_scoreZ2N )			: ( basicC ? ESI_TRAP_CID_low_res_b_scoreZ2C			: ESI_TRAP_CID_low_res_b_scoreZ2 );
			bp2_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_bp2_scoreZ2NC		: ESI_TRAP_CID_low_res_bp2_scoreZ2N )		: ( basicC ? ESI_TRAP_CID_low_res_bp2_scoreZ2C			: ESI_TRAP_CID_low_res_bp2_scoreZ2 );
			b_nh3_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_b_nh3_scoreZ2NC		: ESI_TRAP_CID_low_res_b_nh3_scoreZ2N )		: ( basicC ? ESI_TRAP_CID_low_res_b_nh3_scoreZ2C		: ESI_TRAP_CID_low_res_b_nh3_scoreZ2 );
			bp2_nh3_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_bp2_nh3_scoreZ2NC	: ESI_TRAP_CID_low_res_bp2_nh3_scoreZ2N )	: ( basicC ? ESI_TRAP_CID_low_res_bp2_nh3_scoreZ2C		: ESI_TRAP_CID_low_res_bp2_nh3_scoreZ2 );
			b_h2o_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_b_h2o_scoreZ2NC		: ESI_TRAP_CID_low_res_b_h2o_scoreZ2N )		: ( basicC ? ESI_TRAP_CID_low_res_b_h2o_scoreZ2C		: ESI_TRAP_CID_low_res_b_h2o_scoreZ2 );
			bp2_h2o_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_bp2_h2o_scoreZ2NC	: ESI_TRAP_CID_low_res_bp2_h2o_scoreZ2N )	: ( basicC ? ESI_TRAP_CID_low_res_bp2_h2o_scoreZ2C		: ESI_TRAP_CID_low_res_bp2_h2o_scoreZ2 );
			b_soch4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_b_soch4_scoreZ2NC	: ESI_TRAP_CID_low_res_b_soch4_scoreZ2N )	: ( basicC ? ESI_TRAP_CID_low_res_b_soch4_scoreZ2C		: ESI_TRAP_CID_low_res_b_soch4_scoreZ2 );
			bp2_soch4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_bp2_soch4_scoreZ2NC	: ESI_TRAP_CID_low_res_bp2_soch4_scoreZ2N )	: ( basicC ? ESI_TRAP_CID_low_res_bp2_soch4_scoreZ2C	: ESI_TRAP_CID_low_res_bp2_soch4_scoreZ2 );
			b_h3po4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_b_h3po4_scoreZ2NC	: ESI_TRAP_CID_low_res_b_h3po4_scoreZ2N )	: ( basicC ? ESI_TRAP_CID_low_res_b_h3po4_scoreZ2C		: ESI_TRAP_CID_low_res_b_h3po4_scoreZ2 );
			bp2_h3po4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ2NC	: ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ2N )	: ( basicC ? ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ2C	: ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ2 );
			break;
		case 3:
			a_score			= basicN ? ( basicC ? ESI_TRAP_CID_low_res_a_scoreZ3NC			: ESI_TRAP_CID_low_res_a_scoreZ3N )			: ( basicC ? ESI_TRAP_CID_low_res_a_scoreZ3C			: ESI_TRAP_CID_low_res_a_scoreZ3 );
			b_score			= basicN ? ( basicC ? ESI_TRAP_CID_low_res_b_scoreZ3NC			: ESI_TRAP_CID_low_res_b_scoreZ3N )			: ( basicC ? ESI_TRAP_CID_low_res_b_scoreZ3C			: ESI_TRAP_CID_low_res_b_scoreZ3 );
			bp2_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_bp2_scoreZ3NC		: ESI_TRAP_CID_low_res_bp2_scoreZ3N )		: ( basicC ? ESI_TRAP_CID_low_res_bp2_scoreZ3C			: ESI_TRAP_CID_low_res_bp2_scoreZ3 );
			bp3_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_bp3_scoreZ3NC		: ESI_TRAP_CID_low_res_bp3_scoreZ3N )		: ( basicC ? ESI_TRAP_CID_low_res_bp3_scoreZ3C			: ESI_TRAP_CID_low_res_bp3_scoreZ3 );
			b_nh3_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_b_nh3_scoreZ3NC		: ESI_TRAP_CID_low_res_b_nh3_scoreZ3N )		: ( basicC ? ESI_TRAP_CID_low_res_b_nh3_scoreZ3C		: ESI_TRAP_CID_low_res_b_nh3_scoreZ3 );
			bp2_nh3_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_bp2_nh3_scoreZ3NC	: ESI_TRAP_CID_low_res_bp2_nh3_scoreZ3N )	: ( basicC ? ESI_TRAP_CID_low_res_bp2_nh3_scoreZ3C		: ESI_TRAP_CID_low_res_bp2_nh3_scoreZ3 );
			b_h2o_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_b_h2o_scoreZ3NC		: ESI_TRAP_CID_low_res_b_h2o_scoreZ3N )		: ( basicC ? ESI_TRAP_CID_low_res_b_h2o_scoreZ3C		: ESI_TRAP_CID_low_res_b_h2o_scoreZ3 );
			bp2_h2o_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_bp2_h2o_scoreZ3NC	: ESI_TRAP_CID_low_res_bp2_h2o_scoreZ3N )	: ( basicC ? ESI_TRAP_CID_low_res_bp2_h2o_scoreZ3C		: ESI_TRAP_CID_low_res_bp2_h2o_scoreZ3 );
			b_soch4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_b_soch4_scoreZ3NC	: ESI_TRAP_CID_low_res_b_soch4_scoreZ3N )	: ( basicC ? ESI_TRAP_CID_low_res_b_soch4_scoreZ3C		: ESI_TRAP_CID_low_res_b_soch4_scoreZ3 );
			bp2_soch4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_bp2_soch4_scoreZ3NC	: ESI_TRAP_CID_low_res_bp2_soch4_scoreZ3N )	: ( basicC ? ESI_TRAP_CID_low_res_bp2_soch4_scoreZ3C	: ESI_TRAP_CID_low_res_bp2_soch4_scoreZ3 );
			b_h3po4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_b_h3po4_scoreZ3NC	: ESI_TRAP_CID_low_res_b_h3po4_scoreZ3N )	: ( basicC ? ESI_TRAP_CID_low_res_b_h3po4_scoreZ3C		: ESI_TRAP_CID_low_res_b_h3po4_scoreZ3 );
			bp2_h3po4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ3NC	: ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ3N )	: ( basicC ? ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ3C	: ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ3 );
			break;
		case 4:
			a_score			= basicN ? ( basicC ? ESI_TRAP_CID_low_res_a_scoreZ4NC			: ESI_TRAP_CID_low_res_a_scoreZ4N )			: ( basicC ? ESI_TRAP_CID_low_res_a_scoreZ4C			: ESI_TRAP_CID_low_res_a_scoreZ4 );
			b_score			= basicN ? ( basicC ? ESI_TRAP_CID_low_res_b_scoreZ4NC			: ESI_TRAP_CID_low_res_b_scoreZ4N )			: ( basicC ? ESI_TRAP_CID_low_res_b_scoreZ4C			: ESI_TRAP_CID_low_res_b_scoreZ4 );
			bp2_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_bp2_scoreZ4NC		: ESI_TRAP_CID_low_res_bp2_scoreZ4N )		: ( basicC ? ESI_TRAP_CID_low_res_bp2_scoreZ4C			: ESI_TRAP_CID_low_res_bp2_scoreZ4 );
			bp3_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_bp3_scoreZ4NC		: ESI_TRAP_CID_low_res_bp3_scoreZ4N )		: ( basicC ? ESI_TRAP_CID_low_res_bp3_scoreZ4C			: ESI_TRAP_CID_low_res_bp3_scoreZ4 );
			b_nh3_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_b_nh3_scoreZ4NC		: ESI_TRAP_CID_low_res_b_nh3_scoreZ4N )		: ( basicC ? ESI_TRAP_CID_low_res_b_nh3_scoreZ4C		: ESI_TRAP_CID_low_res_b_nh3_scoreZ4 );
			bp2_nh3_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_bp2_nh3_scoreZ4NC	: ESI_TRAP_CID_low_res_bp2_nh3_scoreZ4N )	: ( basicC ? ESI_TRAP_CID_low_res_bp2_nh3_scoreZ4C		: ESI_TRAP_CID_low_res_bp2_nh3_scoreZ4 );
			b_h2o_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_b_h2o_scoreZ4NC		: ESI_TRAP_CID_low_res_b_h2o_scoreZ4N )		: ( basicC ? ESI_TRAP_CID_low_res_b_h2o_scoreZ4C		: ESI_TRAP_CID_low_res_b_h2o_scoreZ4 );
			bp2_h2o_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_bp2_h2o_scoreZ4NC	: ESI_TRAP_CID_low_res_bp2_h2o_scoreZ4N )	: ( basicC ? ESI_TRAP_CID_low_res_bp2_h2o_scoreZ4C		: ESI_TRAP_CID_low_res_bp2_h2o_scoreZ4 );
			b_soch4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_b_soch4_scoreZ4NC	: ESI_TRAP_CID_low_res_b_soch4_scoreZ4N )	: ( basicC ? ESI_TRAP_CID_low_res_b_soch4_scoreZ4C		: ESI_TRAP_CID_low_res_b_soch4_scoreZ4 );
			bp2_soch4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_bp2_soch4_scoreZ4NC	: ESI_TRAP_CID_low_res_bp2_soch4_scoreZ4N )	: ( basicC ? ESI_TRAP_CID_low_res_bp2_soch4_scoreZ4C	: ESI_TRAP_CID_low_res_bp2_soch4_scoreZ4 );
			b_h3po4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_b_h3po4_scoreZ4NC	: ESI_TRAP_CID_low_res_b_h3po4_scoreZ4N )	: ( basicC ? ESI_TRAP_CID_low_res_b_h3po4_scoreZ4C		: ESI_TRAP_CID_low_res_b_h3po4_scoreZ4 );
			bp2_h3po4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ4NC	: ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ4N )	: ( basicC ? ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ4C	: ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ4 );
			break;
		default:
			a_score			= basicN ? ( basicC ? ESI_TRAP_CID_low_res_a_scoreZ5NC			: ESI_TRAP_CID_low_res_a_scoreZ5N )			: ( basicC ? ESI_TRAP_CID_low_res_a_scoreZ5C			: ESI_TRAP_CID_low_res_a_scoreZ5 );
			b_score			= basicN ? ( basicC ? ESI_TRAP_CID_low_res_b_scoreZ5NC			: ESI_TRAP_CID_low_res_b_scoreZ5N )			: ( basicC ? ESI_TRAP_CID_low_res_b_scoreZ5C			: ESI_TRAP_CID_low_res_b_scoreZ5 );
			bp2_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_bp2_scoreZ5NC		: ESI_TRAP_CID_low_res_bp2_scoreZ5N )		: ( basicC ? ESI_TRAP_CID_low_res_bp2_scoreZ5C			: ESI_TRAP_CID_low_res_bp2_scoreZ5 );
			bp3_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_bp3_scoreZ5NC		: ESI_TRAP_CID_low_res_bp3_scoreZ5N )		: ( basicC ? ESI_TRAP_CID_low_res_bp3_scoreZ5C			: ESI_TRAP_CID_low_res_bp3_scoreZ5 );
			b_nh3_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_b_nh3_scoreZ5NC		: ESI_TRAP_CID_low_res_b_nh3_scoreZ5N )		: ( basicC ? ESI_TRAP_CID_low_res_b_nh3_scoreZ5C		: ESI_TRAP_CID_low_res_b_nh3_scoreZ5 );
			bp2_nh3_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_bp2_nh3_scoreZ5NC	: ESI_TRAP_CID_low_res_bp2_nh3_scoreZ5N )	: ( basicC ? ESI_TRAP_CID_low_res_bp2_nh3_scoreZ5C		: ESI_TRAP_CID_low_res_bp2_nh3_scoreZ5 );
			b_h2o_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_b_h2o_scoreZ5NC		: ESI_TRAP_CID_low_res_b_h2o_scoreZ5N )		: ( basicC ? ESI_TRAP_CID_low_res_b_h2o_scoreZ5C		: ESI_TRAP_CID_low_res_b_h2o_scoreZ5 );
			bp2_h2o_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_bp2_h2o_scoreZ5NC	: ESI_TRAP_CID_low_res_bp2_h2o_scoreZ5N )	: ( basicC ? ESI_TRAP_CID_low_res_bp2_h2o_scoreZ5C		: ESI_TRAP_CID_low_res_bp2_h2o_scoreZ5 );
			b_soch4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_b_soch4_scoreZ5NC	: ESI_TRAP_CID_low_res_b_soch4_scoreZ5N )	: ( basicC ? ESI_TRAP_CID_low_res_b_soch4_scoreZ5C		: ESI_TRAP_CID_low_res_b_soch4_scoreZ5 );
			bp2_soch4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_bp2_soch4_scoreZ5NC	: ESI_TRAP_CID_low_res_bp2_soch4_scoreZ5N )	: ( basicC ? ESI_TRAP_CID_low_res_bp2_soch4_scoreZ5C	: ESI_TRAP_CID_low_res_bp2_soch4_scoreZ5 );
			b_h3po4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_b_h3po4_scoreZ5NC	: ESI_TRAP_CID_low_res_b_h3po4_scoreZ5N )	: ( basicC ? ESI_TRAP_CID_low_res_b_h3po4_scoreZ5C		: ESI_TRAP_CID_low_res_b_h3po4_scoreZ5 );
			bp2_h3po4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ5NC	: ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ5N )	: ( basicC ? ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ5C	: ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ5 );
			break;
	}
	for ( int k = 0 ; k < len ; k++ ) {
		char aa = sequence [k];
		nIon += aaArrayMT [aa];
		if ( nIon > maxNFragTagInData ) break;

		if ( checkNIonAmmoniaLoss )			nIonAmmoniaLoss |= ammoniaLossFlag [aa];
		if ( checkNIonPosChargeBearing && posChargeBearingFlag [aa] ) nIonPosChargeBearing++;
		if ( checkNIonWaterLoss )			nIonWaterLoss |= waterLossFlag [aa];
		if ( checkNIonM )					nIonM |= oxidizedMFlag [aa];
		if ( checkNIonPhosphorylation )		nIonPhosphorylation |= phosphorylationFlag [aa];
		if ( nIon < minNFragTagInData || k == 0 ) continue;
		for ( int i = start ; i < numPeaks ; i++ ) {
			if ( nIon < minNFragTag [i] ) break;
			if ( nIon > maxNFragTag [i] ) {
				start++;
				continue;
			}
			MassTypeVector& fragTag = nIonFragTag [i];
			ScoreType& matched = ionMatched [i];
			MassType fragTol = fragTolerance [i];
			MassType fragTol2 = fragTol + fragTol;
			if ( genAbsDiff ( nIon, fragTag [b_index] ) < fragTol ) {
				matched = genMax ( b_score, matched );
			}
			else if ( doublyChargedIons && nIonPosChargeBearing && genAbsDiff ( nIon, fragTag [bp2_index] ) < fragTol2 ) {
				matched = genMax ( bp2_score, matched );
			}
			else if ( triplyChargedIons && nIonPosChargeBearing > 1 ) {
				if ( genAbsDiff ( nIon, fragTag [bp3_index] ) < fragTol2 + fragTol ) {
					matched = genMax ( bp3_score, matched );
				}
			}
			else if ( genAbsDiff ( nIon, fragTag [a_index] ) < fragTol ) {
				matched = genMax ( a_score, matched );
			}
			else if ( nIonAmmoniaLoss && genAbsDiff ( nIon, fragTag [b_nh3_index] ) < fragTol ) {
				matched = genMax ( b_nh3_score, matched );
			}
			else if ( doublyChargedIons && nIonAmmoniaLoss && nIonPosChargeBearing && genAbsDiff ( nIon, fragTag [bp2_nh3_index] ) < fragTol2 ) {
				matched = genMax ( bp2_nh3_score, matched );
			}
			else if ( nIonWaterLoss && genAbsDiff ( nIon, fragTag [b_h2o_index] ) < fragTol ) {
				matched = genMax ( b_h2o_score, matched );
			}
			else if ( doublyChargedIons && nIonWaterLoss && nIonPosChargeBearing && genAbsDiff ( nIon, fragTag [bp2_h2o_index] ) < fragTol2 ) {
				matched = genMax ( bp2_h2o_score, matched );
			}
			else if ( nIonM && genAbsDiff ( nIon, fragTag [b_soch4_index] ) < fragTol ) {
				matched = genMax ( b_soch4_score, matched );
			}
			else if ( doublyChargedIons && nIonM && nIonPosChargeBearing && genAbsDiff ( nIon, fragTag [bp2_soch4_index] ) < fragTol2 ) {
				matched = genMax ( bp2_soch4_score, matched );
			}
			else if ( nIonPhosphorylation && genAbsDiff ( nIon, fragTag [b_h3po4_index] ) < fragTol ) {
				matched = genMax ( b_h3po4_score, matched );
			}
			else if ( doublyChargedIons && nIonPhosphorylation && nIonPosChargeBearing && genAbsDiff ( nIon, fragTag [bp2_h3po4_index] ) < fragTol2 ) {
				matched = genMax ( bp2_h3po4_score, matched );
			}
		}
	}
}
void MSTagSearch::doNIonsESI_Q_CID ( const string& sequence )
{
	int nIonAmmoniaLoss = 0;
	int nIonPosChargeBearing = initialPosCharge;
	int nIonM = 0;
	int nIonPhosphorylation = 0;
	int nIonWaterLoss = 0;
	MassType nIon = nTerminusWt;
	int start = 0;
	int len = sequence.length () - 1;
	bool basicN = posChargeBearingFlag [sequence[0]] != 0;
	bool basicC = posChargeBearingFlag [sequence[len]] != 0;
	switch ( precursorCharge ) {
		case 1:
			a_score			= basicN ? ( basicC ? ESI_Q_CID_a_scoreZ1NC			: ESI_Q_CID_a_scoreZ1N )		: ( basicC ? ESI_Q_CID_a_scoreZ1C		: ESI_Q_CID_a_scoreZ1 );
			a_nh3_score		= basicN ? ( basicC ? ESI_Q_CID_a_nh3_scoreZ1NC		: ESI_Q_CID_a_nh3_scoreZ1N )	: ( basicC ? ESI_Q_CID_a_nh3_scoreZ1C	: ESI_Q_CID_a_nh3_scoreZ1 );
			a_h2o_score		= basicN ? ( basicC ? ESI_Q_CID_a_h2o_scoreZ1NC		: ESI_Q_CID_a_h2o_scoreZ1N )	: ( basicC ? ESI_Q_CID_a_h2o_scoreZ1C	: ESI_Q_CID_a_h2o_scoreZ1 );
			a_h3po4_score	= basicN ? ( basicC ? ESI_Q_CID_a_h3po4_scoreZ1NC	: ESI_Q_CID_a_h3po4_scoreZ1N )	: ( basicC ? ESI_Q_CID_a_h3po4_scoreZ1C	: ESI_Q_CID_a_h3po4_scoreZ1 );
			b_score			= basicN ? ( basicC ? ESI_Q_CID_b_scoreZ1NC			: ESI_Q_CID_b_scoreZ1N )		: ( basicC ? ESI_Q_CID_b_scoreZ1C		: ESI_Q_CID_b_scoreZ1 );
			b_nh3_score		= basicN ? ( basicC ? ESI_Q_CID_b_nh3_scoreZ1NC		: ESI_Q_CID_b_nh3_scoreZ1N )	: ( basicC ? ESI_Q_CID_b_nh3_scoreZ1C	: ESI_Q_CID_b_nh3_scoreZ1 );
			b_h2o_score		= basicN ? ( basicC ? ESI_Q_CID_b_h2o_scoreZ1NC		: ESI_Q_CID_b_h2o_scoreZ1N )	: ( basicC ? ESI_Q_CID_b_h2o_scoreZ1C	: ESI_Q_CID_b_h2o_scoreZ1 );
			b_soch4_score	= basicN ? ( basicC ? ESI_Q_CID_b_soch4_scoreZ1NC	: ESI_Q_CID_b_soch4_scoreZ1N )	: ( basicC ? ESI_Q_CID_b_soch4_scoreZ1C	: ESI_Q_CID_b_soch4_scoreZ1 );
			b_h3po4_score	= basicN ? ( basicC ? ESI_Q_CID_b_h3po4_scoreZ1NC	: ESI_Q_CID_b_h3po4_scoreZ1N )	: ( basicC ? ESI_Q_CID_b_h3po4_scoreZ1C	: ESI_Q_CID_b_h3po4_scoreZ1 );
			break;
		case 2:
			a_score			= basicN ? ( basicC ? ESI_Q_CID_a_scoreZ2NC			: ESI_Q_CID_a_scoreZ2N )		: ( basicC ? ESI_Q_CID_a_scoreZ2C		: ESI_Q_CID_a_scoreZ2 );
			a_nh3_score		= basicN ? ( basicC ? ESI_Q_CID_a_nh3_scoreZ2NC		: ESI_Q_CID_a_nh3_scoreZ2N )	: ( basicC ? ESI_Q_CID_a_nh3_scoreZ2C	: ESI_Q_CID_a_nh3_scoreZ2 );
			a_h2o_score		= basicN ? ( basicC ? ESI_Q_CID_a_h2o_scoreZ2NC		: ESI_Q_CID_a_h2o_scoreZ2N )	: ( basicC ? ESI_Q_CID_a_h2o_scoreZ2C	: ESI_Q_CID_a_h2o_scoreZ2 );
			a_h3po4_score	= basicN ? ( basicC ? ESI_Q_CID_a_h3po4_scoreZ2NC	: ESI_Q_CID_a_h3po4_scoreZ2N )	: ( basicC ? ESI_Q_CID_a_h3po4_scoreZ2C	: ESI_Q_CID_a_h3po4_scoreZ2 );
			b_score			= basicN ? ( basicC ? ESI_Q_CID_b_scoreZ2NC			: ESI_Q_CID_b_scoreZ2N )		: ( basicC ? ESI_Q_CID_b_scoreZ2C		: ESI_Q_CID_b_scoreZ2 );
			b_nh3_score		= basicN ? ( basicC ? ESI_Q_CID_b_nh3_scoreZ2NC		: ESI_Q_CID_b_nh3_scoreZ2N )	: ( basicC ? ESI_Q_CID_b_nh3_scoreZ2C	: ESI_Q_CID_b_nh3_scoreZ2 );
			b_h2o_score		= basicN ? ( basicC ? ESI_Q_CID_b_h2o_scoreZ2NC		: ESI_Q_CID_b_h2o_scoreZ2N )	: ( basicC ? ESI_Q_CID_b_h2o_scoreZ2C	: ESI_Q_CID_b_h2o_scoreZ2 );
			b_soch4_score	= basicN ? ( basicC ? ESI_Q_CID_b_soch4_scoreZ2NC	: ESI_Q_CID_b_soch4_scoreZ2N )	: ( basicC ? ESI_Q_CID_b_soch4_scoreZ2C	: ESI_Q_CID_b_soch4_scoreZ2 );
			b_h3po4_score	= basicN ? ( basicC ? ESI_Q_CID_b_h3po4_scoreZ2NC	: ESI_Q_CID_b_h3po4_scoreZ2N )	: ( basicC ? ESI_Q_CID_b_h3po4_scoreZ2C	: ESI_Q_CID_b_h3po4_scoreZ2 );
			break;
		case 3:
			a_score			= basicN ? ( basicC ? ESI_Q_CID_a_scoreZ3NC			: ESI_Q_CID_a_scoreZ3N )		: ( basicC ? ESI_Q_CID_a_scoreZ3C		: ESI_Q_CID_a_scoreZ3 );
			a_nh3_score		= basicN ? ( basicC ? ESI_Q_CID_a_nh3_scoreZ3NC		: ESI_Q_CID_a_nh3_scoreZ3N )	: ( basicC ? ESI_Q_CID_a_nh3_scoreZ3C	: ESI_Q_CID_a_nh3_scoreZ3 );
			a_h2o_score		= basicN ? ( basicC ? ESI_Q_CID_a_h2o_scoreZ3NC		: ESI_Q_CID_a_h2o_scoreZ3N )	: ( basicC ? ESI_Q_CID_a_h2o_scoreZ3C	: ESI_Q_CID_a_h2o_scoreZ3 );
			a_h3po4_score	= basicN ? ( basicC ? ESI_Q_CID_a_h3po4_scoreZ3NC	: ESI_Q_CID_a_h3po4_scoreZ3N )	: ( basicC ? ESI_Q_CID_a_h3po4_scoreZ3C	: ESI_Q_CID_a_h3po4_scoreZ3 );
			b_score			= basicN ? ( basicC ? ESI_Q_CID_b_scoreZ3NC			: ESI_Q_CID_b_scoreZ3N )		: ( basicC ? ESI_Q_CID_b_scoreZ3C		: ESI_Q_CID_b_scoreZ3 );
			b_nh3_score		= basicN ? ( basicC ? ESI_Q_CID_b_nh3_scoreZ3NC		: ESI_Q_CID_b_nh3_scoreZ3N )	: ( basicC ? ESI_Q_CID_b_nh3_scoreZ3C	: ESI_Q_CID_b_nh3_scoreZ3 );
			b_h2o_score		= basicN ? ( basicC ? ESI_Q_CID_b_h2o_scoreZ3NC		: ESI_Q_CID_b_h2o_scoreZ3N )	: ( basicC ? ESI_Q_CID_b_h2o_scoreZ3C	: ESI_Q_CID_b_h2o_scoreZ3 );
			b_soch4_score	= basicN ? ( basicC ? ESI_Q_CID_b_soch4_scoreZ3NC	: ESI_Q_CID_b_soch4_scoreZ3N )	: ( basicC ? ESI_Q_CID_b_soch4_scoreZ3C	: ESI_Q_CID_b_soch4_scoreZ3 );
			b_h3po4_score	= basicN ? ( basicC ? ESI_Q_CID_b_h3po4_scoreZ3NC	: ESI_Q_CID_b_h3po4_scoreZ3N )	: ( basicC ? ESI_Q_CID_b_h3po4_scoreZ3C	: ESI_Q_CID_b_h3po4_scoreZ3 );
			break;
		case 4:
			a_score			= basicN ? ( basicC ? ESI_Q_CID_a_scoreZ4NC			: ESI_Q_CID_a_scoreZ4N )		: ( basicC ? ESI_Q_CID_a_scoreZ4C		: ESI_Q_CID_a_scoreZ4 );
			a_nh3_score		= basicN ? ( basicC ? ESI_Q_CID_a_nh3_scoreZ4NC		: ESI_Q_CID_a_nh3_scoreZ4N )	: ( basicC ? ESI_Q_CID_a_nh3_scoreZ4C	: ESI_Q_CID_a_nh3_scoreZ4 );
			a_h2o_score		= basicN ? ( basicC ? ESI_Q_CID_a_h2o_scoreZ4NC		: ESI_Q_CID_a_h2o_scoreZ4N )	: ( basicC ? ESI_Q_CID_a_h2o_scoreZ4C	: ESI_Q_CID_a_h2o_scoreZ4 );
			a_h3po4_score	= basicN ? ( basicC ? ESI_Q_CID_a_h3po4_scoreZ4NC	: ESI_Q_CID_a_h3po4_scoreZ4N )	: ( basicC ? ESI_Q_CID_a_h3po4_scoreZ4C	: ESI_Q_CID_a_h3po4_scoreZ4 );
			b_score			= basicN ? ( basicC ? ESI_Q_CID_b_scoreZ4NC			: ESI_Q_CID_b_scoreZ4N )		: ( basicC ? ESI_Q_CID_b_scoreZ4C		: ESI_Q_CID_b_scoreZ4 );
			b_nh3_score		= basicN ? ( basicC ? ESI_Q_CID_b_nh3_scoreZ4NC		: ESI_Q_CID_b_nh3_scoreZ4N )	: ( basicC ? ESI_Q_CID_b_nh3_scoreZ4C	: ESI_Q_CID_b_nh3_scoreZ4 );
			b_h2o_score		= basicN ? ( basicC ? ESI_Q_CID_b_h2o_scoreZ4NC		: ESI_Q_CID_b_h2o_scoreZ4N )	: ( basicC ? ESI_Q_CID_b_h2o_scoreZ4C	: ESI_Q_CID_b_h2o_scoreZ4 );
			b_soch4_score	= basicN ? ( basicC ? ESI_Q_CID_b_soch4_scoreZ4NC	: ESI_Q_CID_b_soch4_scoreZ4N )	: ( basicC ? ESI_Q_CID_b_soch4_scoreZ4C	: ESI_Q_CID_b_soch4_scoreZ4 );
			b_h3po4_score	= basicN ? ( basicC ? ESI_Q_CID_b_h3po4_scoreZ4NC	: ESI_Q_CID_b_h3po4_scoreZ4N )	: ( basicC ? ESI_Q_CID_b_h3po4_scoreZ4C	: ESI_Q_CID_b_h3po4_scoreZ4 );
			break;
		case 5:
			a_score			= basicN ? ( basicC ? ESI_Q_CID_a_scoreZ5NC			: ESI_Q_CID_a_scoreZ5N )		: ( basicC ? ESI_Q_CID_a_scoreZ5C		: ESI_Q_CID_a_scoreZ5 );
			a_nh3_score		= basicN ? ( basicC ? ESI_Q_CID_a_nh3_scoreZ5NC		: ESI_Q_CID_a_nh3_scoreZ5N )	: ( basicC ? ESI_Q_CID_a_nh3_scoreZ5C	: ESI_Q_CID_a_nh3_scoreZ5 );
			a_h2o_score		= basicN ? ( basicC ? ESI_Q_CID_a_h2o_scoreZ5NC		: ESI_Q_CID_a_h2o_scoreZ5N )	: ( basicC ? ESI_Q_CID_a_h2o_scoreZ5C	: ESI_Q_CID_a_h2o_scoreZ5 );
			a_h3po4_score	= basicN ? ( basicC ? ESI_Q_CID_a_h3po4_scoreZ5NC	: ESI_Q_CID_a_h3po4_scoreZ5N )	: ( basicC ? ESI_Q_CID_a_h3po4_scoreZ5C	: ESI_Q_CID_a_h3po4_scoreZ5 );
			b_score			= basicN ? ( basicC ? ESI_Q_CID_b_scoreZ5NC			: ESI_Q_CID_b_scoreZ5N )		: ( basicC ? ESI_Q_CID_b_scoreZ5C		: ESI_Q_CID_b_scoreZ5 );
			b_nh3_score		= basicN ? ( basicC ? ESI_Q_CID_b_nh3_scoreZ5NC		: ESI_Q_CID_b_nh3_scoreZ5N )	: ( basicC ? ESI_Q_CID_b_nh3_scoreZ5C	: ESI_Q_CID_b_nh3_scoreZ5 );
			b_h2o_score		= basicN ? ( basicC ? ESI_Q_CID_b_h2o_scoreZ5NC		: ESI_Q_CID_b_h2o_scoreZ5N )	: ( basicC ? ESI_Q_CID_b_h2o_scoreZ5C	: ESI_Q_CID_b_h2o_scoreZ5 );
			b_soch4_score	= basicN ? ( basicC ? ESI_Q_CID_b_soch4_scoreZ5NC	: ESI_Q_CID_b_soch4_scoreZ5N )	: ( basicC ? ESI_Q_CID_b_soch4_scoreZ5C	: ESI_Q_CID_b_soch4_scoreZ5 );
			b_h3po4_score	= basicN ? ( basicC ? ESI_Q_CID_b_h3po4_scoreZ5NC	: ESI_Q_CID_b_h3po4_scoreZ5N )	: ( basicC ? ESI_Q_CID_b_h3po4_scoreZ5C	: ESI_Q_CID_b_h3po4_scoreZ5 );
			break;
	}
	for ( int k = 0 ; k < len ; k++ ) {
		char aa = sequence [k];
		nIon += aaArrayMT [aa];
		if ( nIon > maxNFragTagInData ) break;

		if ( checkNIonAmmoniaLoss )			nIonAmmoniaLoss |= ammoniaLossFlag [aa];
		if ( checkNIonPosChargeBearing )	nIonPosChargeBearing |= posChargeBearingFlag [aa];
		if ( checkNIonWaterLoss )			nIonWaterLoss |= waterLossFlag [aa];
		if ( checkNIonM )					nIonM |= oxidizedMFlag [aa];
		if ( checkNIonPhosphorylation )		nIonPhosphorylation |= phosphorylationFlag [aa];
		if ( nIon < minNFragTagInData || k == 0 ) continue;
		for ( int i = start ; i < numPeaks ; i++ ) {
			if ( nIon < minNFragTag [i] ) break;
			if ( nIon > maxNFragTag [i] ) {
				start++;
				continue;
			}
			MassTypeVector& fragTag = nIonFragTag [i];
			ScoreType& matched = ionMatched [i];
			MassType fragTol = fragTolerance [i];
			if ( genAbsDiff ( nIon, fragTag [b_index] ) < fragTol ) {
				matched = genMax ( b_score, matched );
			}
			else if ( genAbsDiff ( nIon, fragTag [a_index] ) < fragTol ) {
				matched = genMax ( a_score, matched );
			}
			else if ( nIonAmmoniaLoss && genAbsDiff ( nIon, fragTag [a_nh3_index] ) < fragTol ) {
				matched = genMax ( a_nh3_score, matched );
			}
			else if ( nIonWaterLoss && genAbsDiff ( nIon, fragTag [a_h2o_index] ) < fragTol ) {
				matched = genMax ( a_h2o_score, matched );
			}
			else if ( nIonAmmoniaLoss && genAbsDiff ( nIon, fragTag [b_nh3_index] ) < fragTol ) {
				matched = genMax ( b_nh3_score, matched );
			}
			else if ( nIonWaterLoss && genAbsDiff ( nIon, fragTag [b_h2o_index] ) < fragTol ) {
				matched = genMax ( b_h2o_score, matched );
			}
			else if ( k >= len - 3 && nIonPosChargeBearing && genAbsDiff ( nIon, fragTag [b_plus_h2o_index] ) < fragTol ) {
				matched = genMax ( b_plus_h2o_score, matched );
			}
			else if ( nIonM && genAbsDiff ( nIon, fragTag [b_soch4_index] ) < fragTol ) {
				matched = genMax ( b_soch4_score, matched );
			}
			else if ( nIonPhosphorylation && genAbsDiff ( nIon, fragTag [a_h3po4_index] ) < fragTol ) {
				matched = genMax ( a_h3po4_score, matched );
			}
			else if ( nIonPhosphorylation && genAbsDiff ( nIon, fragTag [b_h3po4_index] ) < fragTol ) {
				matched = genMax ( b_h3po4_score, matched );
			}
		}
	}
}
void MSTagSearch::doNIonsESI_ETD_low_res ( const string& sequence )
{
	int nIonPosChargeBearing = initialPosCharge;
	MassType nIon = nTerminusWt;
	int start = 0;
	int len = sequence.length () - 1;
	bool basicN = posChargeBearingFlag [sequence[0]] != 0;
	bool basicC = posChargeBearingFlag [sequence[len]] != 0;
	switch ( precursorCharge ) {
		case 1:
		case 2:
			cMinus1DaScore	= basicN ? ( basicC ? ESI_ETD_low_res_cMinus1DaZ2NC : ESI_ETD_low_res_cMinus1DaZ2N ) : ( basicC ? ESI_ETD_low_res_cMinus1DaZ2C : ESI_ETD_low_res_cMinus1DaZ2 );
			c_score			= basicN ? ( basicC ? ESI_ETD_low_res_c_scoreZ2NC : ESI_ETD_low_res_c_scoreZ2N ) : ( basicC ? ESI_ETD_low_res_c_scoreZ2C : ESI_ETD_low_res_c_scoreZ2 );
			b_score			= basicN ? ( basicC ? ESI_ETD_low_res_b_scoreZ2NC : ESI_ETD_low_res_b_scoreZ2N ) : ( basicC ? ESI_ETD_low_res_b_scoreZ2C : ESI_ETD_low_res_b_scoreZ2 );
			break;
		case 3:
			cMinus1DaScore	= basicN ? ( basicC ? ESI_ETD_low_res_cMinus1DaZ3NC : ESI_ETD_low_res_cMinus1DaZ3N ) : ( basicC ? ESI_ETD_low_res_cMinus1DaZ3C : ESI_ETD_low_res_cMinus1DaZ3 );
			c_score			= basicN ? ( basicC ? ESI_ETD_low_res_c_scoreZ3NC : ESI_ETD_low_res_c_scoreZ3N ) : ( basicC ? ESI_ETD_low_res_c_scoreZ3C : ESI_ETD_low_res_c_scoreZ3 );
			cp2_score		= basicN ? ( basicC ? ESI_ETD_low_res_cp2_scoreZ3NC : ESI_ETD_low_res_cp2_scoreZ3N ) : ( basicC ? ESI_ETD_low_res_cp2_scoreZ3C : ESI_ETD_low_res_cp2_scoreZ3 );
			break;
		case 4:
			c_score			= basicN ? ( basicC ? ESI_ETD_low_res_c_scoreZ4NC : ESI_ETD_low_res_c_scoreZ4N ) : ( basicC ? ESI_ETD_low_res_c_scoreZ4C : ESI_ETD_low_res_c_scoreZ4 );
			cp2_score		= basicN ? ( basicC ? ESI_ETD_low_res_cp2_scoreZ4NC : ESI_ETD_low_res_cp2_scoreZ4N ) : ( basicC ? ESI_ETD_low_res_cp2_scoreZ4C : ESI_ETD_low_res_cp2_scoreZ4 );
			break;
		default:
			c_score			= basicN ? ( basicC ? ESI_ETD_low_res_c_scoreZ5NC : ESI_ETD_low_res_c_scoreZ5N ) : ( basicC ? ESI_ETD_low_res_c_scoreZ5C : ESI_ETD_low_res_c_scoreZ5 );
			cp2_score		= basicN ? ( basicC ? ESI_ETD_low_res_cp2_scoreZ5NC : ESI_ETD_low_res_cp2_scoreZ5N ) : ( basicC ? ESI_ETD_low_res_cp2_scoreZ5C : ESI_ETD_low_res_cp2_scoreZ5 );
			break;
	}
	if ( precursorCharge <= 2 ) {
		for ( int k = 0 ; k < len ; k++ ) {
			bool pCheck = ( k == len-1 || sequence [k+1] != 'P' );
			char aa = sequence [k];
			nIon += aaArrayMT [aa];
			if ( nIon > maxNFragTagInData ) break;
			if ( nIon < minNFragTagInData || k == 0 ) continue;
			for ( int i = start ; i < numPeaks ; i++ ) {
				if ( nIon < minNFragTag [i] ) break;
				if ( nIon > maxNFragTag [i] ) {
					start++;
					continue;
				}
				MassTypeVector& fragTag = nIonFragTag [i];
				ScoreType& matched = ionMatched [i];
				MassType fragTol = fragTolerance [i];
				if ( pCheck && genAbsDiff ( nIon, fragTag [c_index] ) < fragTol ) {
					matched = genMax ( c_score, matched );
				}
				else if ( pCheck && genAbsDiff ( nIon, fragTag [cMinus1DaIndex] ) < fragTol ) {
					matched = genMax ( cMinus1DaScore, matched );
				}
				else if ( b_score && genAbsDiff ( nIon, fragTag [b_index] ) < fragTol ) {
					matched = genMax ( b_score, matched );
				}
			}
		}
	}
	else if ( precursorCharge == 3 ) {
		for ( int k = 0 ; k < len ; k++ ) {
			bool pCheck = ( k == len-1 || sequence [k+1] != 'P' );
			char aa = sequence [k];
			nIon += aaArrayMT [aa];
			if ( nIon > maxNFragTagInData ) break;

			if ( posChargeBearingFlag [aa] ) nIonPosChargeBearing++;
			if ( nIon < minNFragTagInData || k == 0 ) continue;
			for ( int i = start ; i < numPeaks ; i++ ) {
				if ( nIon < minNFragTag [i] ) break;
				if ( nIon > maxNFragTag [i] ) {
					start++;
					continue;
				}
				if ( pCheck ) {
					MassTypeVector& fragTag = nIonFragTag [i];
					ScoreType& matched = ionMatched [i];
					MassType fragTol = fragTolerance [i];
					MassType fragTol2 = fragTol + fragTol;
					if ( genAbsDiff ( nIon, fragTag [c_index] ) < fragTol ) {
						matched = genMax ( c_score, matched );
					}
					else if ( nIonPosChargeBearing && cp2_score && genAbsDiff ( nIon, fragTag [cp2_index] ) < fragTol2 ) {
						matched = genMax ( cp2_score, matched );
					}
					else if ( cMinus1DaScore && genAbsDiff ( nIon, fragTag [cMinus1DaIndex] ) < fragTol ) {
						matched = genMax ( cMinus1DaScore, matched );
					}
				}
			}
		}
	}
	else {
		for ( int k = 0 ; k < len ; k++ ) {
			bool pCheck = ( k == len-1 || sequence [k+1] != 'P' );
			char aa = sequence [k];
			nIon += aaArrayMT [aa];
			if ( nIon > maxNFragTagInData ) break;

			if ( posChargeBearingFlag [aa] ) nIonPosChargeBearing++;
			if ( nIon < minNFragTagInData || k == 0 ) continue;
			for ( int i = start ; i < numPeaks ; i++ ) {
				if ( nIon < minNFragTag [i] ) break;
				if ( nIon > maxNFragTag [i] ) {
					start++;
					continue;
				}
				if ( pCheck ) {
					MassTypeVector& fragTag = nIonFragTag [i];
					ScoreType& matched = ionMatched [i];
					MassType fragTol = fragTolerance [i];
					if ( genAbsDiff ( nIon, fragTag [c_index] ) < fragTol ) {
						matched = genMax ( c_score, matched );
					}
					else if ( nIonPosChargeBearing && cp2_score ) {
						MassType fragTol2 = fragTol + fragTol;
						if ( genAbsDiff ( nIon, fragTag [cp2_index] ) < fragTol2 ) {
							matched = genMax ( cp2_score, matched );
						}
					}
				}
			}
		}
	}
}
void MSTagSearch::doNIonsESI_ETD_high_res ( const string& sequence )
{
	MassType nIon = nTerminusWt;
	int start = 0;
	int len = sequence.length () - 1;
	bool basicN = posChargeBearingFlag [sequence[0]] != 0;
	bool basicC = posChargeBearingFlag [sequence[len]] != 0;
	switch ( precursorCharge ) {
		case 1:
		case 2:
			cMinus1DaScore	= basicN ? ( basicC ? ESI_ETD_hi_res_cMinus1DaZ2NC : ESI_ETD_hi_res_cMinus1DaZ2N ) : ( basicC ? ESI_ETD_hi_res_cMinus1DaZ2C : ESI_ETD_hi_res_cMinus1DaZ2 );
			c_score			= basicN ? ( basicC ? ESI_ETD_hi_res_c_scoreZ2NC : ESI_ETD_hi_res_c_scoreZ2N ) : ( basicC ? ESI_ETD_hi_res_c_scoreZ2C : ESI_ETD_hi_res_c_scoreZ2 );
			b_score			= basicN ? ( basicC ? ESI_ETD_hi_res_b_scoreZ2NC : ESI_ETD_hi_res_b_scoreZ2N ) : ( basicC ? ESI_ETD_hi_res_b_scoreZ2C : ESI_ETD_hi_res_b_scoreZ2 );
			break;
		case 3:
			cMinus1DaScore	= basicN ? ( basicC ? ESI_ETD_hi_res_cMinus1DaZ3NC : ESI_ETD_hi_res_cMinus1DaZ3N ) : ( basicC ? ESI_ETD_hi_res_cMinus1DaZ3C : ESI_ETD_hi_res_cMinus1DaZ3 );
			c_score			= basicN ? ( basicC ? ESI_ETD_hi_res_c_scoreZ3NC : ESI_ETD_hi_res_c_scoreZ3N ) : ( basicC ? ESI_ETD_hi_res_c_scoreZ3C : ESI_ETD_hi_res_c_scoreZ3 );
			break;
		case 4:
			c_score			= basicN ? ( basicC ? ESI_ETD_hi_res_c_scoreZ4NC : ESI_ETD_hi_res_c_scoreZ4N ) : ( basicC ? ESI_ETD_hi_res_c_scoreZ4C : ESI_ETD_hi_res_c_scoreZ4 );
			break;
		default:
			c_score			= basicN ? ( basicC ? ESI_ETD_hi_res_c_scoreZ5NC : ESI_ETD_hi_res_c_scoreZ5N ) : ( basicC ? ESI_ETD_hi_res_c_scoreZ5C : ESI_ETD_hi_res_c_scoreZ5 );
			break;
	}
	for ( int k = 0 ; k < len ; k++ ) {
		bool pCheck = ( k == len-1 || sequence [k+1] != 'P' );
		char aa = sequence [k];
		nIon += aaArrayMT [aa];
		if ( nIon > maxNFragTagInData ) break;
		if ( nIon < minNFragTagInData || k == 0 ) continue;
		for ( int i = start ; i < numPeaks ; i++ ) {
			if ( nIon < minNFragTag [i] ) break;
			if ( nIon > maxNFragTag [i] ) {
				start++;
				continue;
			}
			MassTypeVector& fragTag = nIonFragTag [i];
			ScoreType& matched = ionMatched [i];
			MassType fragTol = fragTolerance [i];
			if ( pCheck && genAbsDiff ( nIon, fragTag [c_index] ) < fragTol ) {
				matched = genMax ( c_score, matched );
			}
			else if ( precursorCharge <= 3 && pCheck && cMinus1DaScore && genAbsDiff ( nIon, fragTag [cMinus1DaIndex] ) < fragTol ) {
				matched = genMax ( cMinus1DaScore, matched );
			}
			else if ( precursorCharge <= 2 && b_score && genAbsDiff ( nIon, fragTag [b_index] ) < fragTol ) {
				matched = genMax ( b_score, matched );
			}
		}
	}
}
void MSTagSearch::doCIons ( const string& sequence )
{
	int cIonPosChargeBearing = initialPosCharge;
	int cIonWaterLoss = 0;
	int cIonAmmoniaLoss = 0;
	int cIonM = 0;
	int cIonPhosphorylation = 0;
	int vIonExclude = 0;
	int wIonExclude = 0;
	MassType cIon = cTerminusWt;

	int start = 0;
	int len = sequence.length ();
	for ( int k = len ; --k ; ) {
		char aa = sequence [k];
		cIon += aaArrayMT [aa];
		if ( cIon > maxCFragTagInData ) break;

		if ( checkCIonPosChargeBearing )	cIonPosChargeBearing |= posChargeBearingFlag [aa];
		if ( checkCIonAmmoniaLoss )			cIonAmmoniaLoss |= ammoniaLossFlag [aa];
		if ( checkCIonM )					cIonM |= oxidizedMFlag [aa];
		if ( checkCIonPhosphorylation )		cIonPhosphorylation |= phosphorylationFlag [aa];
		if ( checkCIonWaterLoss )			cIonWaterLoss |= waterLossFlag [aa];
		if ( v_flag )						vIonExclude = vIonExcludeFlag [sequence [k-1]];
		if ( w_flag )						wIonExclude = wIonExcludeFlag [aa];
		if ( cIon < minCFragTagInData ) continue;
		for ( int i = start ; i < numPeaks ; i++ ) {
			if ( cIon < minCFragTag [i] ) break;
			if ( cIon > maxCFragTag [i] ) {
				start++;
				continue;
			}
			MassTypeVector& fragTag = cIonFragTag [i];
			ScoreType& matched = ionMatched [i];
			MassType fragTol = fragTolerance [i];
			MassType fragTol2 = fragTol + fragTol;

			if ( y_flag && genAbsDiff ( cIon, fragTag [y_index] ) < fragTol ) {
				matched = genMax ( y_score, matched );
			}
			if ( n_ladder_flag && genAbsDiff ( cIon, fragTag [y_index] ) < fragTol ) {
				matched = genMax ( n_ladder_score, matched );
			}
			else if ( doublyChargedIons && yp2_flag && cIonPosChargeBearing && genAbsDiff ( cIon, fragTag [yp2_index] ) < fragTol2 ) {
				matched = genMax ( yp2_score, matched );
			}
			else if ( y_nh3_flag && cIonAmmoniaLoss && genAbsDiff ( cIon, fragTag [y_nh3_index] ) < fragTol ) {
				matched = genMax ( y_nh3_score, matched );
			}
			else if ( doublyChargedIons && yp2_nh3_flag && cIonAmmoniaLoss && cIonPosChargeBearing && genAbsDiff ( cIon, fragTag [yp2_nh3_index] ) < fragTol2 ) {
				matched = genMax ( yp2_nh3_score, matched );
			}
			else if ( y_h2o_flag && cIonWaterLoss && genAbsDiff ( cIon, fragTag [y_h2o_index] ) < fragTol ) {
				matched = genMax ( y_h2o_score, matched );
			}
			else if ( doublyChargedIons && yp2_h2o_flag && cIonWaterLoss && cIonPosChargeBearing && genAbsDiff ( cIon, fragTag [yp2_h2o_index] ) < fragTol2 ) {
				matched = genMax ( yp2_h2o_score, matched );
			}
			else if ( y_soch4_flag && cIonM && genAbsDiff ( cIon, fragTag [y_soch4_index] ) < fragTol ) {
				matched = genMax ( y_soch4_score, matched );
			} 
			else if ( doublyChargedIons && yp2_soch4_flag && cIonM && cIonPosChargeBearing && genAbsDiff ( cIon, fragTag [yp2_soch4_index] ) < fragTol2 ) {
				matched = genMax ( yp2_soch4_score, matched );
			}
			else if ( y_h3po4_flag && cIonPhosphorylation && genAbsDiff ( cIon, fragTag [y_h3po4_index] ) < fragTol ) {
				matched = genMax ( y_h3po4_score, matched );
			}
			else if ( doublyChargedIons && yp2_h3po4_flag && cIonPhosphorylation && cIonPosChargeBearing && genAbsDiff ( cIon, fragTag [yp2_h3po4_index] ) < fragTol2 ) {
				matched = genMax ( yp2_h3po4_score, matched );
			}
			else if ( x_flag && genAbsDiff ( cIon, fragTag [x_index] ) < fragTol ) {
				matched = genMax ( x_score, matched );
			}
			else if ( Y_flag && genAbsDiff ( cIon, fragTag [Y_index] ) < fragTol ) {
				matched = genMax ( Y_score, matched );
			}
			else if ( v_flag && cIonPosChargeBearing && vIonExclude == 0 && genAbsDiff ( cIon, fragTag [v_index] ) < fragTol ) {
				matched = genMax ( v_score, matched );
			}
			else if ( z_flag && genAbsDiff ( cIon, fragTag [z_index] ) < fragTol ) {
				matched = genMax ( z_score, matched );
			}
			else if ( zPlus1DaFlag && genAbsDiff ( cIon, fragTag [zPlus1DaIndex] ) < fragTol ) {
				matched = genMax ( zPlus1DaScore, matched );
			}
			else if ( zPlus2DaFlag && genAbsDiff ( cIon, fragTag [zPlus2DaIndex] ) < fragTol ) {
				matched = genMax ( zPlus2DaScore, matched );
			}
			else if ( zPlus3DaFlag && genAbsDiff ( cIon, fragTag [zPlus3DaIndex] ) < fragTol ) {
				matched = genMax ( zPlus3DaScore, matched );
			}
			else if ( w_flag && cIonPosChargeBearing && wIonExclude == false && substituentWaMT [aa] ) {
				if ( genAbsDiff ( cIon, (MassType)(fragTag [y_index] - substituentWaMT [aa]) ) < fragTol ) {
					matched = genMax ( w_score, matched );
				}
				else {
					if ( substituentWbMT [aa] && genAbsDiff ( cIon, (MassType)(fragTag [y_index] - substituentWbMT [aa]) ) < fragTol ) {
						matched = genMax ( w_score, matched );
					}
				}
			}
		}
	}
}
void MSTagSearch::doCIonsESI_TRAP_CID_low_res ( const string& sequence )
{
	int cIonPosChargeBearing = initialPosCharge;
	int cIonWaterLoss = 0;
	int cIonAmmoniaLoss = 0;
	int cIonM = 0;
	int cIonPhosphorylation = 0;
	MassType cIon = cTerminusWt;
	int start = 0;
	int len = sequence.length ();
	bool basicN = posChargeBearingFlag [sequence[0]] != 0;
	bool basicC = posChargeBearingFlag [sequence[len-1]] != 0;
	switch ( precursorCharge ) {
		case 1:
			y_score			= basicN ? ( basicC ? ESI_TRAP_CID_low_res_y_scoreZ1NC			: ESI_TRAP_CID_low_res_y_scoreZ1N )			: ( basicC ? ESI_TRAP_CID_low_res_y_scoreZ1C		: ESI_TRAP_CID_low_res_y_scoreZ1 );
			y_nh3_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_y_nh3_scoreZ1NC		: ESI_TRAP_CID_low_res_y_nh3_scoreZ1N )		: ( basicC ? ESI_TRAP_CID_low_res_y_nh3_scoreZ1C	: ESI_TRAP_CID_low_res_y_nh3_scoreZ1 );
			y_h2o_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_y_h2o_scoreZ1NC		: ESI_TRAP_CID_low_res_y_h2o_scoreZ1N )		: ( basicC ? ESI_TRAP_CID_low_res_y_h2o_scoreZ1C	: ESI_TRAP_CID_low_res_y_h2o_scoreZ1 );
			y_soch4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_y_soch4_scoreZ1NC	: ESI_TRAP_CID_low_res_y_soch4_scoreZ1N )	: ( basicC ? ESI_TRAP_CID_low_res_y_soch4_scoreZ1C	: ESI_TRAP_CID_low_res_y_soch4_scoreZ1 );
			y_h3po4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_y_h3po4_scoreZ1NC	: ESI_TRAP_CID_low_res_y_h3po4_scoreZ1N )	: ( basicC ? ESI_TRAP_CID_low_res_y_h3po4_scoreZ1C	: ESI_TRAP_CID_low_res_y_h3po4_scoreZ1 );
			break;
		case 2:
			y_score			= basicN ? ( basicC ? ESI_TRAP_CID_low_res_y_scoreZ2NC			: ESI_TRAP_CID_low_res_y_scoreZ2N )			: ( basicC ? ESI_TRAP_CID_low_res_y_scoreZ2C		: ESI_TRAP_CID_low_res_y_scoreZ2 );
			yp2_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_yp2_scoreZ2NC		: ESI_TRAP_CID_low_res_yp2_scoreZ2N )		: ( basicC ? ESI_TRAP_CID_low_res_yp2_scoreZ2C		: ESI_TRAP_CID_low_res_yp2_scoreZ2 );
			y_nh3_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_y_nh3_scoreZ2NC		: ESI_TRAP_CID_low_res_y_nh3_scoreZ2N )		: ( basicC ? ESI_TRAP_CID_low_res_y_nh3_scoreZ2C	: ESI_TRAP_CID_low_res_y_nh3_scoreZ2 );
			yp2_nh3_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_yp2_nh3_scoreZ2NC	: ESI_TRAP_CID_low_res_yp2_nh3_scoreZ2N )	: ( basicC ? ESI_TRAP_CID_low_res_yp2_nh3_scoreZ2C	: ESI_TRAP_CID_low_res_yp2_nh3_scoreZ2 );
			y_h2o_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_y_h2o_scoreZ2NC		: ESI_TRAP_CID_low_res_y_h2o_scoreZ2N )		: ( basicC ? ESI_TRAP_CID_low_res_y_h2o_scoreZ2C	: ESI_TRAP_CID_low_res_y_h2o_scoreZ2 );
			yp2_h2o_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_yp2_h2o_scoreZ2NC	: ESI_TRAP_CID_low_res_yp2_h2o_scoreZ2N )	: ( basicC ? ESI_TRAP_CID_low_res_yp2_h2o_scoreZ2C	: ESI_TRAP_CID_low_res_yp2_h2o_scoreZ2 );
			y_soch4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_y_soch4_scoreZ2NC	: ESI_TRAP_CID_low_res_y_soch4_scoreZ2N )	: ( basicC ? ESI_TRAP_CID_low_res_y_soch4_scoreZ2C	: ESI_TRAP_CID_low_res_y_soch4_scoreZ2 );
			yp2_soch4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_yp2_soch4_scoreZ2NC	: ESI_TRAP_CID_low_res_yp2_soch4_scoreZ2N )	: ( basicC ? ESI_TRAP_CID_low_res_yp2_soch4_scoreZ2C: ESI_TRAP_CID_low_res_yp2_soch4_scoreZ2 );
			y_h3po4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_y_h3po4_scoreZ2NC	: ESI_TRAP_CID_low_res_y_h3po4_scoreZ2N )	: ( basicC ? ESI_TRAP_CID_low_res_y_h3po4_scoreZ2C	: ESI_TRAP_CID_low_res_y_h3po4_scoreZ2 );
			yp2_h3po4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ2NC	: ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ2N )	: ( basicC ? ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ2C: ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ2 );
			break;
		case 3:
			y_score			= basicN ? ( basicC ? ESI_TRAP_CID_low_res_y_scoreZ3NC			: ESI_TRAP_CID_low_res_y_scoreZ3N )			: ( basicC ? ESI_TRAP_CID_low_res_y_scoreZ3C		: ESI_TRAP_CID_low_res_y_scoreZ3 );
			yp2_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_yp2_scoreZ3NC		: ESI_TRAP_CID_low_res_yp2_scoreZ3N )		: ( basicC ? ESI_TRAP_CID_low_res_yp2_scoreZ3C		: ESI_TRAP_CID_low_res_yp2_scoreZ3 );
			yp3_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_yp3_scoreZ3NC		: ESI_TRAP_CID_low_res_yp3_scoreZ3N )		: ( basicC ? ESI_TRAP_CID_low_res_yp3_scoreZ3C		: ESI_TRAP_CID_low_res_yp3_scoreZ3 );
			y_nh3_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_y_nh3_scoreZ3NC		: ESI_TRAP_CID_low_res_y_nh3_scoreZ3N )		: ( basicC ? ESI_TRAP_CID_low_res_y_nh3_scoreZ3C	: ESI_TRAP_CID_low_res_y_nh3_scoreZ3 );
			yp2_nh3_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_yp2_nh3_scoreZ3NC	: ESI_TRAP_CID_low_res_yp2_nh3_scoreZ3N )	: ( basicC ? ESI_TRAP_CID_low_res_yp2_nh3_scoreZ3C	: ESI_TRAP_CID_low_res_yp2_nh3_scoreZ3 );
			y_h2o_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_y_h2o_scoreZ3NC		: ESI_TRAP_CID_low_res_y_h2o_scoreZ3N )		: ( basicC ? ESI_TRAP_CID_low_res_y_h2o_scoreZ3C	: ESI_TRAP_CID_low_res_y_h2o_scoreZ3 );
			yp2_h2o_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_yp2_h2o_scoreZ3NC	: ESI_TRAP_CID_low_res_yp2_h2o_scoreZ3N )	: ( basicC ? ESI_TRAP_CID_low_res_yp2_h2o_scoreZ3C	: ESI_TRAP_CID_low_res_yp2_h2o_scoreZ3 );
			y_soch4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_y_soch4_scoreZ3NC	: ESI_TRAP_CID_low_res_y_soch4_scoreZ3N )	: ( basicC ? ESI_TRAP_CID_low_res_y_soch4_scoreZ3C	: ESI_TRAP_CID_low_res_y_soch4_scoreZ3 );
			yp2_soch4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_yp2_soch4_scoreZ3NC	: ESI_TRAP_CID_low_res_yp2_soch4_scoreZ3N )	: ( basicC ? ESI_TRAP_CID_low_res_yp2_soch4_scoreZ3C: ESI_TRAP_CID_low_res_yp2_soch4_scoreZ3 );
			y_h3po4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_y_h3po4_scoreZ3NC	: ESI_TRAP_CID_low_res_y_h3po4_scoreZ3N )	: ( basicC ? ESI_TRAP_CID_low_res_y_h3po4_scoreZ3C	: ESI_TRAP_CID_low_res_y_h3po4_scoreZ3 );
			yp2_h3po4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ3NC	: ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ3N )	: ( basicC ? ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ3C: ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ3 );
			break;
		case 4:
			y_score			= basicN ? ( basicC ? ESI_TRAP_CID_low_res_y_scoreZ4NC			: ESI_TRAP_CID_low_res_y_scoreZ4N )			: ( basicC ? ESI_TRAP_CID_low_res_y_scoreZ4C		: ESI_TRAP_CID_low_res_y_scoreZ4 );
			yp2_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_yp2_scoreZ4NC		: ESI_TRAP_CID_low_res_yp2_scoreZ4N )		: ( basicC ? ESI_TRAP_CID_low_res_yp2_scoreZ4C		: ESI_TRAP_CID_low_res_yp2_scoreZ4 );
			yp3_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_yp3_scoreZ4NC		: ESI_TRAP_CID_low_res_yp3_scoreZ4N )		: ( basicC ? ESI_TRAP_CID_low_res_yp3_scoreZ4C		: ESI_TRAP_CID_low_res_yp3_scoreZ4 );
			y_nh3_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_y_nh3_scoreZ4NC		: ESI_TRAP_CID_low_res_y_nh3_scoreZ4N )		: ( basicC ? ESI_TRAP_CID_low_res_y_nh3_scoreZ4C	: ESI_TRAP_CID_low_res_y_nh3_scoreZ4 );
			yp2_nh3_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_yp2_nh3_scoreZ4NC	: ESI_TRAP_CID_low_res_yp2_nh3_scoreZ4N )	: ( basicC ? ESI_TRAP_CID_low_res_yp2_nh3_scoreZ4C	: ESI_TRAP_CID_low_res_yp2_nh3_scoreZ4 );
			y_h2o_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_y_h2o_scoreZ4NC		: ESI_TRAP_CID_low_res_y_h2o_scoreZ4N )		: ( basicC ? ESI_TRAP_CID_low_res_y_h2o_scoreZ4C	: ESI_TRAP_CID_low_res_y_h2o_scoreZ4 );
			yp2_h2o_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_yp2_h2o_scoreZ4NC	: ESI_TRAP_CID_low_res_yp2_h2o_scoreZ4N )	: ( basicC ? ESI_TRAP_CID_low_res_yp2_h2o_scoreZ4C	: ESI_TRAP_CID_low_res_yp2_h2o_scoreZ4 );
			y_soch4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_y_soch4_scoreZ4NC	: ESI_TRAP_CID_low_res_y_soch4_scoreZ4N )	: ( basicC ? ESI_TRAP_CID_low_res_y_soch4_scoreZ4C	: ESI_TRAP_CID_low_res_y_soch4_scoreZ4 );
			yp2_soch4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_yp2_soch4_scoreZ4NC	: ESI_TRAP_CID_low_res_yp2_soch4_scoreZ4N )	: ( basicC ? ESI_TRAP_CID_low_res_yp2_soch4_scoreZ4C: ESI_TRAP_CID_low_res_yp2_soch4_scoreZ4 );
			y_h3po4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_y_h3po4_scoreZ4NC	: ESI_TRAP_CID_low_res_y_h3po4_scoreZ4N )	: ( basicC ? ESI_TRAP_CID_low_res_y_h3po4_scoreZ4C	: ESI_TRAP_CID_low_res_y_h3po4_scoreZ4 );
			yp2_h3po4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ4NC	: ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ4N )	: ( basicC ? ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ4C: ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ4 );
			break;
		default:
			y_score			= basicN ? ( basicC ? ESI_TRAP_CID_low_res_y_scoreZ5NC			: ESI_TRAP_CID_low_res_y_scoreZ5N )			: ( basicC ? ESI_TRAP_CID_low_res_y_scoreZ5C		: ESI_TRAP_CID_low_res_y_scoreZ5 );
			yp2_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_yp2_scoreZ5NC		: ESI_TRAP_CID_low_res_yp2_scoreZ5N )		: ( basicC ? ESI_TRAP_CID_low_res_yp2_scoreZ5C		: ESI_TRAP_CID_low_res_yp2_scoreZ5 );
			yp3_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_yp3_scoreZ5NC		: ESI_TRAP_CID_low_res_yp3_scoreZ5N )		: ( basicC ? ESI_TRAP_CID_low_res_yp3_scoreZ5C		: ESI_TRAP_CID_low_res_yp3_scoreZ5 );
			y_nh3_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_y_nh3_scoreZ5NC		: ESI_TRAP_CID_low_res_y_nh3_scoreZ5N )		: ( basicC ? ESI_TRAP_CID_low_res_y_nh3_scoreZ5C	: ESI_TRAP_CID_low_res_y_nh3_scoreZ5 );
			yp2_nh3_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_yp2_nh3_scoreZ5NC	: ESI_TRAP_CID_low_res_yp2_nh3_scoreZ5N )	: ( basicC ? ESI_TRAP_CID_low_res_yp2_nh3_scoreZ5C	: ESI_TRAP_CID_low_res_yp2_nh3_scoreZ5 );
			y_h2o_score		= basicN ? ( basicC ? ESI_TRAP_CID_low_res_y_h2o_scoreZ5NC		: ESI_TRAP_CID_low_res_y_h2o_scoreZ5N )		: ( basicC ? ESI_TRAP_CID_low_res_y_h2o_scoreZ5C	: ESI_TRAP_CID_low_res_y_h2o_scoreZ5 );
			yp2_h2o_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_yp2_h2o_scoreZ5NC	: ESI_TRAP_CID_low_res_yp2_h2o_scoreZ5N )	: ( basicC ? ESI_TRAP_CID_low_res_yp2_h2o_scoreZ5C	: ESI_TRAP_CID_low_res_yp2_h2o_scoreZ5 );
			y_soch4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_y_soch4_scoreZ5NC	: ESI_TRAP_CID_low_res_y_soch4_scoreZ5N )	: ( basicC ? ESI_TRAP_CID_low_res_y_soch4_scoreZ5C	: ESI_TRAP_CID_low_res_y_soch4_scoreZ5 );
			yp2_soch4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_yp2_soch4_scoreZ5NC	: ESI_TRAP_CID_low_res_yp2_soch4_scoreZ5N )	: ( basicC ? ESI_TRAP_CID_low_res_yp2_soch4_scoreZ5C: ESI_TRAP_CID_low_res_yp2_soch4_scoreZ5 );
			y_h3po4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_y_h3po4_scoreZ5NC	: ESI_TRAP_CID_low_res_y_h3po4_scoreZ5N )	: ( basicC ? ESI_TRAP_CID_low_res_y_h3po4_scoreZ5C	: ESI_TRAP_CID_low_res_y_h3po4_scoreZ5 );
			yp2_h3po4_score	= basicN ? ( basicC ? ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ5NC	: ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ5N )	: ( basicC ? ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ5C: ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ5 );
			break;
	}
	for ( int k = len ; --k ; ) {
		char aa = sequence [k];
		cIon += aaArrayMT [aa];
		if ( cIon > maxCFragTagInData ) break;

		if ( checkCIonPosChargeBearing && posChargeBearingFlag [aa] )	cIonPosChargeBearing++;
		if ( checkCIonAmmoniaLoss )			cIonAmmoniaLoss |= ammoniaLossFlag [aa];
		if ( checkCIonM )					cIonM |= oxidizedMFlag [aa];
		if ( checkCIonPhosphorylation )		cIonPhosphorylation |= phosphorylationFlag [aa];
		if ( checkCIonWaterLoss )			cIonWaterLoss |= waterLossFlag [aa];
		if ( cIon < minCFragTagInData ) continue;
		for ( int i = start ; i < numPeaks ; i++ ) {
			if ( cIon < minCFragTag [i] ) break;
			if ( cIon > maxCFragTag [i] ) {
				start++;
				continue;
			}
			MassTypeVector& fragTag = cIonFragTag [i];
			ScoreType& matched = ionMatched [i];
			MassType fragTol = fragTolerance [i];
			MassType fragTol2 = fragTol + fragTol;

			if ( genAbsDiff ( cIon, fragTag [y_index] ) < fragTol ) {
				matched = genMax ( y_score, matched );
			}
			else if ( doublyChargedIons && cIonPosChargeBearing && genAbsDiff ( cIon, fragTag [yp2_index] ) < fragTol2 ) {
				matched = genMax ( yp2_score, matched );
			}
			else if ( triplyChargedIons && cIonPosChargeBearing > 1 ) {
				if ( genAbsDiff ( cIon, fragTag [yp3_index] ) < fragTol2 + fragTol ) {
					matched = genMax ( yp3_score, matched );
				}
			}
			else if ( cIonAmmoniaLoss && genAbsDiff ( cIon, fragTag [y_nh3_index] ) < fragTol ) {
				matched = genMax ( y_nh3_score, matched );
			}
			else if ( doublyChargedIons && cIonAmmoniaLoss && cIonPosChargeBearing && genAbsDiff ( cIon, fragTag [yp2_nh3_index] ) < fragTol2 ) {
				matched = genMax ( yp2_nh3_score, matched );
			}
			else if ( cIonWaterLoss && genAbsDiff ( cIon, fragTag [y_h2o_index] ) < fragTol ) {
				matched = genMax ( y_h2o_score, matched );
			}
			else if ( doublyChargedIons && cIonWaterLoss && cIonPosChargeBearing && genAbsDiff ( cIon, fragTag [yp2_h2o_index] ) < fragTol2 ) {
				matched = genMax ( yp2_h2o_score, matched );
			}
			else if ( cIonM && genAbsDiff ( cIon, fragTag [y_soch4_index] ) < fragTol ) {
				matched = genMax ( y_soch4_score, matched );
			} 
			else if ( doublyChargedIons && cIonM && cIonPosChargeBearing && genAbsDiff ( cIon, fragTag [yp2_soch4_index] ) < fragTol2 ) {
				matched = genMax ( yp2_soch4_score, matched );
			}
			else if ( cIonPhosphorylation && genAbsDiff ( cIon, fragTag [y_h3po4_index] ) < fragTol ) {
				matched = genMax ( y_h3po4_score, matched );
			}
			else if ( doublyChargedIons && cIonPhosphorylation && cIonPosChargeBearing && genAbsDiff ( cIon, fragTag [yp2_h3po4_index] ) < fragTol2 ) {
				matched = genMax ( yp2_h3po4_score, matched );
			}
		}
	}
}
void MSTagSearch::doCIonsESI_Q_CID ( const string& sequence )
{
	int cIonWaterLoss = 0;
	int cIonAmmoniaLoss = 0;
	int cIonM = 0;
	int cIonPhosphorylation = 0;
	MassType cIon = cTerminusWt;
	int start = 0;
	int len = sequence.length ();
	bool basicN = posChargeBearingFlag [sequence[0]] != 0;
	bool basicC = posChargeBearingFlag [sequence[len-1]] != 0;
	switch ( precursorCharge ) {
		case 1:
			y_score			= basicN ? ( basicC ? ESI_Q_CID_y_scoreZ1NC			: ESI_Q_CID_y_scoreZ1N )		: ( basicC ? ESI_Q_CID_y_scoreZ1C		: ESI_Q_CID_y_scoreZ1 );
			y_nh3_score		= basicN ? ( basicC ? ESI_Q_CID_y_nh3_scoreZ1NC		: ESI_Q_CID_y_nh3_scoreZ1N )	: ( basicC ? ESI_Q_CID_y_nh3_scoreZ1C	: ESI_Q_CID_y_nh3_scoreZ1 );
			y_h2o_score		= basicN ? ( basicC ? ESI_Q_CID_y_h2o_scoreZ1NC		: ESI_Q_CID_y_h2o_scoreZ1N )	: ( basicC ? ESI_Q_CID_y_h2o_scoreZ1C	: ESI_Q_CID_y_h2o_scoreZ1 );
			y_soch4_score	= basicN ? ( basicC ? ESI_Q_CID_y_soch4_scoreZ1NC	: ESI_Q_CID_y_soch4_scoreZ1N )	: ( basicC ? ESI_Q_CID_y_soch4_scoreZ1C	: ESI_Q_CID_y_soch4_scoreZ1 );
			y_h3po4_score	= basicN ? ( basicC ? ESI_Q_CID_y_h3po4_scoreZ1NC	: ESI_Q_CID_y_h3po4_scoreZ1N )	: ( basicC ? ESI_Q_CID_y_h3po4_scoreZ1C	: ESI_Q_CID_y_h3po4_scoreZ1 );
			break;
		case 2:
			y_score			= basicN ? ( basicC ? ESI_Q_CID_y_scoreZ2NC			: ESI_Q_CID_y_scoreZ2N )		: ( basicC ? ESI_Q_CID_y_scoreZ2C		: ESI_Q_CID_y_scoreZ2 );
			y_nh3_score		= basicN ? ( basicC ? ESI_Q_CID_y_nh3_scoreZ2NC		: ESI_Q_CID_y_nh3_scoreZ2N )	: ( basicC ? ESI_Q_CID_y_nh3_scoreZ2C	: ESI_Q_CID_y_nh3_scoreZ2 );
			y_h2o_score		= basicN ? ( basicC ? ESI_Q_CID_y_h2o_scoreZ2NC		: ESI_Q_CID_y_h2o_scoreZ2N )	: ( basicC ? ESI_Q_CID_y_h2o_scoreZ2C	: ESI_Q_CID_y_h2o_scoreZ2 );
			y_soch4_score	= basicN ? ( basicC ? ESI_Q_CID_y_soch4_scoreZ2NC	: ESI_Q_CID_y_soch4_scoreZ2N )	: ( basicC ? ESI_Q_CID_y_soch4_scoreZ2C	: ESI_Q_CID_y_soch4_scoreZ2 );
			y_h3po4_score	= basicN ? ( basicC ? ESI_Q_CID_y_h3po4_scoreZ2NC	: ESI_Q_CID_y_h3po4_scoreZ2N )	: ( basicC ? ESI_Q_CID_y_h3po4_scoreZ2C	: ESI_Q_CID_y_h3po4_scoreZ2 );
			break;
		case 3:
			y_score			= basicN ? ( basicC ? ESI_Q_CID_y_scoreZ3NC			: ESI_Q_CID_y_scoreZ3N )		: ( basicC ? ESI_Q_CID_y_scoreZ3C		: ESI_Q_CID_y_scoreZ3 );
			y_nh3_score		= basicN ? ( basicC ? ESI_Q_CID_y_nh3_scoreZ3NC		: ESI_Q_CID_y_nh3_scoreZ3N )	: ( basicC ? ESI_Q_CID_y_nh3_scoreZ3C	: ESI_Q_CID_y_nh3_scoreZ3 );
			y_h2o_score		= basicN ? ( basicC ? ESI_Q_CID_y_h2o_scoreZ3NC		: ESI_Q_CID_y_h2o_scoreZ3N )	: ( basicC ? ESI_Q_CID_y_h2o_scoreZ3C	: ESI_Q_CID_y_h2o_scoreZ3 );
			y_soch4_score	= basicN ? ( basicC ? ESI_Q_CID_y_soch4_scoreZ3NC	: ESI_Q_CID_y_soch4_scoreZ3N )	: ( basicC ? ESI_Q_CID_y_soch4_scoreZ3C	: ESI_Q_CID_y_soch4_scoreZ3 );
			y_h3po4_score	= basicN ? ( basicC ? ESI_Q_CID_y_h3po4_scoreZ3NC	: ESI_Q_CID_y_h3po4_scoreZ3N )	: ( basicC ? ESI_Q_CID_y_h3po4_scoreZ3C	: ESI_Q_CID_y_h3po4_scoreZ3 );
			break;
		case 4:
			y_score			= basicN ? ( basicC ? ESI_Q_CID_y_scoreZ4NC			: ESI_Q_CID_y_scoreZ4N )		: ( basicC ? ESI_Q_CID_y_scoreZ4C		: ESI_Q_CID_y_scoreZ4 );
			y_nh3_score		= basicN ? ( basicC ? ESI_Q_CID_y_nh3_scoreZ4NC		: ESI_Q_CID_y_nh3_scoreZ4N )	: ( basicC ? ESI_Q_CID_y_nh3_scoreZ4C	: ESI_Q_CID_y_nh3_scoreZ4 );
			y_h2o_score		= basicN ? ( basicC ? ESI_Q_CID_y_h2o_scoreZ4NC		: ESI_Q_CID_y_h2o_scoreZ4N )	: ( basicC ? ESI_Q_CID_y_h2o_scoreZ4C	: ESI_Q_CID_y_h2o_scoreZ4 );
			y_soch4_score	= basicN ? ( basicC ? ESI_Q_CID_y_soch4_scoreZ4NC	: ESI_Q_CID_y_soch4_scoreZ4N )	: ( basicC ? ESI_Q_CID_y_soch4_scoreZ4C	: ESI_Q_CID_y_soch4_scoreZ4 );
			y_h3po4_score	= basicN ? ( basicC ? ESI_Q_CID_y_h3po4_scoreZ4NC	: ESI_Q_CID_y_h3po4_scoreZ4N )	: ( basicC ? ESI_Q_CID_y_h3po4_scoreZ4C	: ESI_Q_CID_y_h3po4_scoreZ4 );
			break;
		default:
			y_score			= basicN ? ( basicC ? ESI_Q_CID_y_scoreZ5NC			: ESI_Q_CID_y_scoreZ5N )		: ( basicC ? ESI_Q_CID_y_scoreZ5C		: ESI_Q_CID_y_scoreZ5 );
			y_nh3_score		= basicN ? ( basicC ? ESI_Q_CID_y_nh3_scoreZ5NC		: ESI_Q_CID_y_nh3_scoreZ5N )	: ( basicC ? ESI_Q_CID_y_nh3_scoreZ5C	: ESI_Q_CID_y_nh3_scoreZ5 );
			y_h2o_score		= basicN ? ( basicC ? ESI_Q_CID_y_h2o_scoreZ5NC		: ESI_Q_CID_y_h2o_scoreZ5N )	: ( basicC ? ESI_Q_CID_y_h2o_scoreZ5C	: ESI_Q_CID_y_h2o_scoreZ5 );
			y_soch4_score	= basicN ? ( basicC ? ESI_Q_CID_y_soch4_scoreZ5NC	: ESI_Q_CID_y_soch4_scoreZ5N )	: ( basicC ? ESI_Q_CID_y_soch4_scoreZ5C	: ESI_Q_CID_y_soch4_scoreZ5 );
			y_h3po4_score	= basicN ? ( basicC ? ESI_Q_CID_y_h3po4_scoreZ5NC	: ESI_Q_CID_y_h3po4_scoreZ5N )	: ( basicC ? ESI_Q_CID_y_h3po4_scoreZ5C	: ESI_Q_CID_y_h3po4_scoreZ5 );
			break;
	}
	for ( int k = len ; --k ; ) {
		char aa = sequence [k];
		cIon += aaArrayMT [aa];
		if ( cIon > maxCFragTagInData ) break;

		if ( checkCIonAmmoniaLoss )			cIonAmmoniaLoss |= ammoniaLossFlag [aa];
		if ( checkCIonM )					cIonM |= oxidizedMFlag [aa];
		if ( checkCIonPhosphorylation )		cIonPhosphorylation |= phosphorylationFlag [aa];
		if ( checkCIonWaterLoss )			cIonWaterLoss |= waterLossFlag [aa];
		if ( cIon < minCFragTagInData ) continue;
		for ( int i = start ; i < numPeaks ; i++ ) {
			if ( cIon < minCFragTag [i] ) break;
			if ( cIon > maxCFragTag [i] ) {
				start++;
				continue;
			}
			MassTypeVector& fragTag = cIonFragTag [i];
			ScoreType& matched = ionMatched [i];
			MassType fragTol = fragTolerance [i];

			if ( genAbsDiff ( cIon, fragTag [y_index] ) < fragTol ) {
				matched = genMax ( y_score, matched );
			}
			else if ( cIonAmmoniaLoss && genAbsDiff ( cIon, fragTag [y_nh3_index] ) < fragTol ) {
				matched = genMax ( y_nh3_score, matched );
			}
			else if ( cIonWaterLoss && genAbsDiff ( cIon, fragTag [y_h2o_index] ) < fragTol ) {
				matched = genMax ( y_h2o_score, matched );
			}
			else if ( cIonM && genAbsDiff ( cIon, fragTag [y_soch4_index] ) < fragTol ) {
				matched = genMax ( y_soch4_score, matched );
			} 
			else if ( cIonPhosphorylation && genAbsDiff ( cIon, fragTag [y_h3po4_index] ) < fragTol ) {
				matched = genMax ( y_h3po4_score, matched );
			}
		}
	}
}
void MSTagSearch::doCIonsESI_ETD_low_res ( const string& sequence )
{
	int cIonPosChargeBearing = initialPosCharge;
	MassType cIon = cTerminusWt;
	int start = 0;
	int len = sequence.length ();
	bool basicN = posChargeBearingFlag [sequence[0]] != 0;
	bool basicC = posChargeBearingFlag [sequence[len-1]] != 0;
	switch ( precursorCharge ) {
		case 1:
		case 2:
			y_score			= basicN ? ( basicC ? ESI_ETD_low_res_y_scoreZ2NC	: ESI_ETD_low_res_y_scoreZ2N )		: ( basicC ? ESI_ETD_low_res_y_scoreZ2C		: ESI_ETD_low_res_y_scoreZ2 );
			z_score			= basicN ? ( basicC ? ESI_ETD_low_res_z_scoreZ2NC	: ESI_ETD_low_res_z_scoreZ2N )		: ( basicC ? ESI_ETD_low_res_z_scoreZ2C		: ESI_ETD_low_res_z_scoreZ2 );
			zPlus1DaScore	= basicN ? ( basicC ? ESI_ETD_low_res_zPlus1DaZ2NC	:ESI_ETD_low_res_zPlus1DaZ2N )		: ( basicC ? ESI_ETD_low_res_zPlus1DaZ2C	: ESI_ETD_low_res_zPlus1DaZ2 );
			break;
		case 3:
			y_score			= basicN ? ( basicC ? ESI_ETD_low_res_y_scoreZ3NC	:	ESI_ETD_low_res_y_scoreZ3N )	: ( basicC ? ESI_ETD_low_res_y_scoreZ3C		: ESI_ETD_low_res_y_scoreZ3 );
			z_score			= basicN ? ( basicC ? ESI_ETD_low_res_z_scoreZ3NC	:	ESI_ETD_low_res_z_scoreZ3N )	: ( basicC ? ESI_ETD_low_res_z_scoreZ3C		: ESI_ETD_low_res_z_scoreZ3 );
			zPlus1DaScore	= basicN ? ( basicC ? ESI_ETD_low_res_zPlus1DaZ3NC	:	ESI_ETD_low_res_zPlus1DaZ3N )	: ( basicC ? ESI_ETD_low_res_zPlus1DaZ3C	: ESI_ETD_low_res_zPlus1DaZ3 );
			zp2_score		= basicN ? ( basicC ? ESI_ETD_low_res_zp2_scoreZ3NC	:	ESI_ETD_low_res_zp2_scoreZ3N )	: ( basicC ? ESI_ETD_low_res_zp2_scoreZ3C	: ESI_ETD_low_res_zp2_scoreZ3 );
			zPlus1Dap2Score	= basicN ? ( basicC ? ESI_ETD_low_res_zPlus1Dap2Z3NC:	ESI_ETD_low_res_zPlus1Dap2Z3N )	: ( basicC ? ESI_ETD_low_res_zPlus1Dap2Z3C	: ESI_ETD_low_res_zPlus1Dap2Z3 );
			break;
		case 4:
			y_score			= basicN ? ( basicC ? ESI_ETD_low_res_y_scoreZ4NC	:	ESI_ETD_low_res_y_scoreZ4N )	: ( basicC ? ESI_ETD_low_res_y_scoreZ4C		: ESI_ETD_low_res_y_scoreZ4 );
			z_score			= basicN ? ( basicC ? ESI_ETD_low_res_z_scoreZ4NC	:	ESI_ETD_low_res_z_scoreZ4N )	: ( basicC ? ESI_ETD_low_res_z_scoreZ4C		: ESI_ETD_low_res_z_scoreZ4 );
			zPlus1DaScore	= basicN ? ( basicC ? ESI_ETD_low_res_zPlus1DaZ4NC	:	ESI_ETD_low_res_zPlus1DaZ4N )	: ( basicC ? ESI_ETD_low_res_zPlus1DaZ4C	: ESI_ETD_low_res_zPlus1DaZ4 );
			zp2_score		= basicN ? ( basicC ? ESI_ETD_low_res_zp2_scoreZ4NC	:	ESI_ETD_low_res_zp2_scoreZ4N )	: ( basicC ? ESI_ETD_low_res_zp2_scoreZ4C	: ESI_ETD_low_res_zp2_scoreZ4 );
			zPlus1Dap2Score	= basicN ? ( basicC ? ESI_ETD_low_res_zPlus1Dap2Z4NC:	ESI_ETD_low_res_zPlus1Dap2Z4N )	: ( basicC ? ESI_ETD_low_res_zPlus1Dap2Z4C	: ESI_ETD_low_res_zPlus1Dap2Z4 );
			break;
		default:
			y_score			= basicN ? ( basicC ? ESI_ETD_low_res_y_scoreZ5NC	:	ESI_ETD_low_res_y_scoreZ5N )	: ( basicC ? ESI_ETD_low_res_y_scoreZ5C		: ESI_ETD_low_res_y_scoreZ5 );
			z_score			= basicN ? ( basicC ? ESI_ETD_low_res_z_scoreZ5NC	:	ESI_ETD_low_res_z_scoreZ5N )	: ( basicC ? ESI_ETD_low_res_z_scoreZ5C		: ESI_ETD_low_res_z_scoreZ5 );
			zPlus1DaScore	= basicN ? ( basicC ? ESI_ETD_low_res_zPlus1DaZ5NC	:	ESI_ETD_low_res_zPlus1DaZ5N )	: ( basicC ? ESI_ETD_low_res_zPlus1DaZ5C	: ESI_ETD_low_res_zPlus1DaZ5 );
			zp2_score		= basicN ? ( basicC ? ESI_ETD_low_res_zp2_scoreZ5NC :	ESI_ETD_low_res_zp2_scoreZ5N )	: ( basicC ? ESI_ETD_low_res_zp2_scoreZ5C	: ESI_ETD_low_res_zp2_scoreZ5 );
			zPlus1Dap2Score	= basicN ? ( basicC ? ESI_ETD_low_res_zPlus1Dap2Z5NC:	ESI_ETD_low_res_zPlus1Dap2Z5N )	: ( basicC ? ESI_ETD_low_res_zPlus1Dap2Z5C	: ESI_ETD_low_res_zPlus1Dap2Z5 );
			break;
	}
	if ( precursorCharge <= 2 ) {
		for ( int k = len ; --k ; ) {
			char aa = sequence [k];
			bool pCheck = aa != 'P';
			cIon += aaArrayMT [aa];
			if ( cIon > maxCFragTagInData ) break;

			if ( posChargeBearingFlag [aa] ) cIonPosChargeBearing++;
			if ( cIon < minCFragTagInData ) continue;
			if ( cIonPosChargeBearing ) {
				for ( int i = start ; i < numPeaks ; i++ ) {
					if ( cIon < minCFragTag [i] ) break;
					if ( cIon > maxCFragTag [i] ) {
						start++;
						continue;
					}
					MassTypeVector& fragTag = cIonFragTag [i];
					ScoreType& matched = ionMatched [i];
					MassType fragTol = fragTolerance [i];

					if ( pCheck && genAbsDiff ( cIon, fragTag [z_index] ) < fragTol ) {
						matched = genMax ( z_score, matched );
					}
					if ( pCheck && zPlus1DaScore && genAbsDiff ( cIon, fragTag [zPlus1DaIndex] ) < fragTol ) {		// Should both be checked as z+1 score can be less than or greater than z score
						matched = genMax ( zPlus1DaScore, matched );
					}
					else if ( y_score && genAbsDiff ( cIon, fragTag [y_index] ) < fragTol ) {
						matched = genMax ( y_score, matched );
					}
				}
			}
		}
	}
	else {
		for ( int k = len ; --k ; ) {
			char aa = sequence [k];
			bool pCheck = aa != 'P';
			cIon += aaArrayMT [aa];
			if ( cIon > maxCFragTagInData ) break;

			if ( posChargeBearingFlag [aa] ) cIonPosChargeBearing++;
			if ( cIon < minCFragTagInData ) continue;
			if ( cIonPosChargeBearing ) {
				for ( int i = start ; i < numPeaks ; i++ ) {
					if ( cIon < minCFragTag [i] ) break;
					if ( cIon > maxCFragTag [i] ) {
						start++;
						continue;
					}
					MassTypeVector& fragTag = cIonFragTag [i];
					ScoreType& matched = ionMatched [i];
					MassType fragTol = fragTolerance [i];

					if ( pCheck && genAbsDiff ( cIon, fragTag [z_index] ) < fragTol ) {
						matched = genMax ( z_score, matched );
					}
					if ( pCheck && genAbsDiff ( cIon, fragTag [zPlus1DaIndex] ) < fragTol ) {
						matched = genMax ( zPlus1DaScore, matched );
					}
					else if ( pCheck && cIonPosChargeBearing >= 2 ) {
						MassType fragTol2 = fragTol + fragTol;
						if ( genAbsDiff ( cIon, fragTag [zPlus1Dap2_index] ) < fragTol2 ) {
							matched = genMax ( zPlus1Dap2Score, matched );
						}
						else if ( genAbsDiff ( cIon, fragTag [zp2_index] ) < fragTol2 ) {
							matched = genMax ( zp2_score, matched );
						}
					}
					else if ( y_score && genAbsDiff ( cIon, fragTag [y_index] ) < fragTol ) {
						matched = genMax ( y_score, matched );
					}
				}
			}
		}
	}
}
void MSTagSearch::doCIonsESI_ETD_high_res ( const string& sequence )
{
	int cIonPosChargeBearing = initialPosCharge;
	MassType cIon = cTerminusWt;
	int start = 0;
	int len = sequence.length ();
	bool basicN = posChargeBearingFlag [sequence[0]] != 0;
	bool basicC = posChargeBearingFlag [sequence[len-1]] != 0;
	switch ( precursorCharge ) {
		case 1:
		case 2:
			y_score			= basicN ? ( basicC ? ESI_ETD_hi_res_y_scoreZ2NC	: ESI_ETD_hi_res_y_scoreZ2N )		: ( basicC ? ESI_ETD_hi_res_y_scoreZ2C		: ESI_ETD_hi_res_y_scoreZ2 );
			z_score			= basicN ? ( basicC ? ESI_ETD_hi_res_z_scoreZ2NC	: ESI_ETD_hi_res_z_scoreZ2N )		: ( basicC ? ESI_ETD_hi_res_z_scoreZ2C		: ESI_ETD_hi_res_z_scoreZ2 );
			zPlus1DaScore	= basicN ? ( basicC ? ESI_ETD_hi_res_zPlus1DaZ2NC	:ESI_ETD_hi_res_zPlus1DaZ2N )		: ( basicC ? ESI_ETD_hi_res_zPlus1DaZ2C	: ESI_ETD_hi_res_zPlus1DaZ2 );
			break;
		case 3:
			y_score			= basicN ? ( basicC ? ESI_ETD_hi_res_y_scoreZ3NC	:	ESI_ETD_hi_res_y_scoreZ3N )	: ( basicC ? ESI_ETD_hi_res_y_scoreZ3C		: ESI_ETD_hi_res_y_scoreZ3 );
			z_score			= basicN ? ( basicC ? ESI_ETD_hi_res_z_scoreZ3NC	:	ESI_ETD_hi_res_z_scoreZ3N )	: ( basicC ? ESI_ETD_hi_res_z_scoreZ3C		: ESI_ETD_hi_res_z_scoreZ3 );
			zPlus1DaScore	= basicN ? ( basicC ? ESI_ETD_hi_res_zPlus1DaZ3NC	:	ESI_ETD_hi_res_zPlus1DaZ3N )	: ( basicC ? ESI_ETD_hi_res_zPlus1DaZ3C	: ESI_ETD_hi_res_zPlus1DaZ3 );
			break;
		case 4:
			y_score			= basicN ? ( basicC ? ESI_ETD_hi_res_y_scoreZ4NC	:	ESI_ETD_hi_res_y_scoreZ4N )	: ( basicC ? ESI_ETD_hi_res_y_scoreZ4C		: ESI_ETD_hi_res_y_scoreZ4 );
			z_score			= basicN ? ( basicC ? ESI_ETD_hi_res_z_scoreZ4NC	:	ESI_ETD_hi_res_z_scoreZ4N )	: ( basicC ? ESI_ETD_hi_res_z_scoreZ4C		: ESI_ETD_hi_res_z_scoreZ4 );
			zPlus1DaScore	= basicN ? ( basicC ? ESI_ETD_hi_res_zPlus1DaZ4NC	:	ESI_ETD_hi_res_zPlus1DaZ4N )	: ( basicC ? ESI_ETD_hi_res_zPlus1DaZ4C	: ESI_ETD_hi_res_zPlus1DaZ4 );
			break;
		default:
			y_score			= basicN ? ( basicC ? ESI_ETD_hi_res_y_scoreZ5NC	:	ESI_ETD_hi_res_y_scoreZ5N )	: ( basicC ? ESI_ETD_hi_res_y_scoreZ5C		: ESI_ETD_hi_res_y_scoreZ5 );
			z_score			= basicN ? ( basicC ? ESI_ETD_hi_res_z_scoreZ5NC	:	ESI_ETD_hi_res_z_scoreZ5N )	: ( basicC ? ESI_ETD_hi_res_z_scoreZ5C		: ESI_ETD_hi_res_z_scoreZ5 );
			zPlus1DaScore	= basicN ? ( basicC ? ESI_ETD_hi_res_zPlus1DaZ5NC	:	ESI_ETD_hi_res_zPlus1DaZ5N )	: ( basicC ? ESI_ETD_hi_res_zPlus1DaZ5C	: ESI_ETD_hi_res_zPlus1DaZ5 );
			break;
	}
	for ( int k = len ; --k ; ) {
		char aa = sequence [k];
		bool pCheck = aa != 'P';
		cIon += aaArrayMT [aa];
		if ( cIon > maxCFragTagInData ) break;

		if ( posChargeBearingFlag [aa] ) cIonPosChargeBearing++;
		if ( cIon < minCFragTagInData ) continue;
		if ( cIonPosChargeBearing ) {
			for ( int i = start ; i < numPeaks ; i++ ) {
				if ( cIon < minCFragTag [i] ) break;
				if ( cIon > maxCFragTag [i] ) {
					start++;
					continue;
				}
				MassTypeVector& fragTag = cIonFragTag [i];
				ScoreType& matched = ionMatched [i];
				MassType fragTol = fragTolerance [i];

				if ( pCheck && genAbsDiff ( cIon, fragTag [z_index] ) < fragTol ) {
					matched = genMax ( z_score, matched );
				}
				if ( pCheck && zPlus1DaScore && genAbsDiff ( cIon, fragTag [zPlus1DaIndex] ) < fragTol ) {		// Should both be checked as z+1 score can be less than or greater than z score
					matched = genMax ( zPlus1DaScore, matched );
				}
				else if ( y_score && genAbsDiff ( cIon, fragTag [y_index] ) < fragTol ) {
					matched = genMax ( y_score, matched );
				}
			}
		}
	}
}
void MSTagSearch::doInternalIons ( const string& sequence )
{
	int len = sequence.length ();
	for ( int startInternal = 2 ; startInternal < len ; startInternal++ ) {
		MassType nIon = internalNTerminusWt;
		int nIonAmmoniaLoss = 0;
		int nIonWaterLoss = 0;
		int start = 0;
		for ( int k = startInternal ; k < len ; k++ ) {
			char aa = sequence[k-1];
			nIon += aaArrayMT [aa];
			if ( nIon > maxInternalIonMass ) break;
			if ( checkNIonAmmoniaLoss ) nIonAmmoniaLoss |= ammoniaLossFlag [aa];
			if ( checkNIonWaterLoss ) nIonWaterLoss |= waterLossFlag [aa];
			if ( k > startInternal ) {
				for ( int i = start ; i < numPeaks ; i++ ) {
					if ( nIon < minNFragTag [i] ) break;
					if ( nIon > maxNFragTag [i] ) {
						start++;
						continue;
					}
					MassTypeVector& fragTag = nIonFragTag [i];
					ScoreType& matched = ionMatched [i];
					MassType fragTol = fragTolerance [i];
					if ( b_flag && genAbsDiff ( nIon, fragTag [b_index] ) < fragTol ) {
						matched = genMax ( internal_b_score, matched );
					}
					else if ( b_nh3_flag && nIonAmmoniaLoss && genAbsDiff ( nIon, fragTag [b_nh3_index] ) < fragTol ) {
						matched = genMax ( internal_loss_score, matched );
					}
					else if ( a_flag && genAbsDiff ( nIon, fragTag [a_index] ) < fragTol ) {
						matched = genMax ( internal_a_score, matched );
					}
					else if ( b_h2o_flag && nIonWaterLoss && genAbsDiff ( nIon, fragTag [b_h2o_index] ) < fragTol ) {
						matched = genMax ( internal_loss_score, matched );
					}
				}
			}
		}
	}
}
ScoreType UnmatchedCompositionSearch::immonium_score = DEFAULT_HIT_SCORE;
UnmatchedCompositionSearch::UnmatchedCompositionSearch ( const PeakContainer& peaks )
{
	for ( PeakContainerSizeType i = 0 ; i < peaks.size () ; i++ ) {
		if ( peaks [i]->getCharge () == 1 ) {
			double mass = peaks [i]->getMass ();
			if ( ImmoniumInfo::inImmoniumRegion ( mass ) ) {
				ionMask.push_back ( ImmoniumInfo::massToMask ( peaks [i] ) );
			}
			else break;
		}
	}
	numIons = ionMask.size ();
}
int UnmatchedCompositionSearch::doCompositionSearch ( const string& peptide ) const
{
	int numUnmatched = 0;
	unsigned int compositionMask = string_to_mask ( peptide );
	for ( int i = 0 ; i < numIons ; i++ ) {
		if ( ( ionMask [i] & compositionMask ) == 0 ) numUnmatched++;
	}
	return numUnmatched;
}
void UnmatchedCompositionSearch::doCompositionSearch ( const string& peptide, ScoreTypeVector& ionMatched ) const
{
	unsigned int compositionMask = string_to_mask ( peptide );
	for ( int i = 0 ; i < numIons ; i++ ) {
		if ( ( ionMask [i] & compositionMask ) != 0 ) {
			ionMatched [i] = genMax ( immonium_score, ionMatched [i] );
		}
	}
}
