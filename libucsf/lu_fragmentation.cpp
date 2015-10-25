/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_fragmentation.cpp                                          *
*                                                                             *
*  Created    : May 14th 2007                                                 *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2007-2009) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lg_string.h>
#include <lgen_error.h>
#include <lu_mass.h>
#include <lu_getfil.h>
#include <lu_fragmentation.h>
#include <lu_param_list.h>
using std::string;
using std::getline;

MSMSFragmentation::MSMSFragmentation ( const string& fragName ) :
	fragName		( fragName ),
	ammoniaLoss		( "RKNQ" ),
	waterLoss		( "STED" ),
	posChargeBearing( "RHK" ),
	dIonExclude		( "FHPWY" ),
	vIonExclude		( "GP" ),
	wIonExclude		( "FHWY" ),

	unmatchedScore ( DEFAULT_MISS_SCORE ),

	immoniumScore ( DEFAULT_HIT_SCORE ),
	relatedIonScore ( DEFAULT_HIT_SCORE ),
	mScore ( DEFAULT_HIT_SCORE ),

	aScore ( DEFAULT_HIT_SCORE ),
	aLossScore ( DEFAULT_HIT_SCORE ),
	aPhosLossScore ( DEFAULT_HIT_SCORE ),
	bScore ( DEFAULT_HIT_SCORE ),
	bPlusH2OScore ( DEFAULT_HIT_SCORE ),
	bLossScore ( DEFAULT_HIT_SCORE ),
	bPhosLossScore ( DEFAULT_HIT_SCORE ),
	cLadderScore ( DEFAULT_HIT_SCORE ),
	cPlus2DaScore ( DEFAULT_HIT_SCORE ),
	cPlus1DaScore ( DEFAULT_HIT_SCORE ),
	cScore ( DEFAULT_HIT_SCORE ),
	cMinus1DaScore ( DEFAULT_HIT_SCORE ),
	dScore ( DEFAULT_HIT_SCORE ),

	vScore ( DEFAULT_HIT_SCORE ),
	wScore ( DEFAULT_HIT_SCORE ),
	xScore ( DEFAULT_HIT_SCORE ),
	nLadderScore ( DEFAULT_HIT_SCORE ),
	yScore ( DEFAULT_HIT_SCORE ),
	yLossScore ( DEFAULT_HIT_SCORE ),
	yPhosLossScore ( DEFAULT_HIT_SCORE ),
	bigYScore ( DEFAULT_HIT_SCORE ),
	zScore ( DEFAULT_HIT_SCORE ),
	zPlus1DaScore ( DEFAULT_HIT_SCORE ),
	zPlus2DaScore ( DEFAULT_HIT_SCORE ),
	zPlus3DaScore ( DEFAULT_HIT_SCORE ),

	bp2Score ( DEFAULT_HIT_SCORE ),
	bp2LossScore ( DEFAULT_HIT_SCORE ),
	bp2PhosLossScore ( DEFAULT_HIT_SCORE ),

	yp2Score ( DEFAULT_HIT_SCORE ),
	yp2LossScore ( DEFAULT_HIT_SCORE ),
	yp2PhosLossScore ( DEFAULT_HIT_SCORE ),

	internalAScore ( DEFAULT_HIT_SCORE ),
	internalBScore ( DEFAULT_HIT_SCORE ),
	internalLossScore ( DEFAULT_HIT_SCORE ),

	MH3PO4Score ( DEFAULT_HIT_SCORE ),
	MSOCH4Score ( DEFAULT_HIT_SCORE ),

	scoring ( false ),
	maximumInternalIonMass ( 700.0 ),
	chargeReducedFragmentation ( false )
{
	if ( !fragName.empty () ) {
		GenIFStream fromFile ( MsparamsDir::instance ().getParamPath ( "fragmentation.txt" ) );
		string line;
		bool flag = false;
		while ( getline ( fromFile, line ) ) {
			if ( line == fragName ) {
				ParameterList params ( fromFile );
				params.getValue ( "nh3_loss", ammoniaLoss );
				params.getValue ( "h2o_loss", waterLoss );
				params.getValue ( "pos_charge", posChargeBearing );
				params.getValue ( "d_ion_exclude", dIonExclude );
				params.getValue ( "v_ion_exclude", vIonExclude );
				params.getValue ( "w_ion_exclude", wIonExclude );
				params.getValue ( "it", ionTypes );
				params.getValue ( "max_internal_ion_mass", maximumInternalIonMass );
				params.getValue ( "charge_reduced_fragmentation", chargeReducedFragmentation );

				scoring = getScoreValue ( params, "unmatched_score", unmatchedScore );

				scoring = getScoreValue ( params, "immonium_score", immoniumScore );
				scoring = getScoreValue ( params, "related_ion_score", relatedIonScore );
				scoring = getScoreValue ( params, "m_score", mScore );

				scoring = getScoreValue ( params, "a_score", aScore );
				scoring = getScoreValue ( params, "a_loss_score", aLossScore );
				scoring = getScoreValue ( params, "a_phos_loss_score", aPhosLossScore );
				scoring = getScoreValue ( params, "b_score", bScore );
				scoring = getScoreValue ( params, "b_plus_h2o_score", bPlusH2OScore );
				scoring = getScoreValue ( params, "b_loss_score", bLossScore );
				scoring = getScoreValue ( params, "b_phos_loss_score", bPhosLossScore );
				scoring = getScoreValue ( params, "c_ladder_score", cLadderScore );
				scoring = getScoreValue ( params, "c_plus_2_score", cPlus2DaScore );
				scoring = getScoreValue ( params, "c_plus_1_score", cPlus1DaScore );
				scoring = getScoreValue ( params, "c_score", cScore );
				scoring = getScoreValue ( params, "c_minus_1_score", cMinus1DaScore );
				scoring = getScoreValue ( params, "d_score", dScore );

				scoring = getScoreValue ( params, "v_score", vScore );
				scoring = getScoreValue ( params, "w_score", wScore );
				scoring = getScoreValue ( params, "x_score", xScore );
				scoring = getScoreValue ( params, "n_ladder_score", nLadderScore );
				scoring = getScoreValue ( params, "y_score", yScore );
				scoring = getScoreValue ( params, "y_loss_score", yLossScore );
				scoring = getScoreValue ( params, "y_phos_loss_score", yPhosLossScore );
				scoring = getScoreValue ( params, "Y_score", bigYScore );
				scoring = getScoreValue ( params, "z_score", zScore );
				scoring = getScoreValue ( params, "z_plus_1_score", zPlus1DaScore );
				scoring = getScoreValue ( params, "z_plus_2_score", zPlus2DaScore );
				scoring = getScoreValue ( params, "z_plus_3_score", zPlus3DaScore );

				scoring = getScoreValue ( params, "bp2_score", bp2Score );
				scoring = getScoreValue ( params, "bp2_loss_score", bp2LossScore );
				scoring = getScoreValue ( params, "bp2_phos_loss_score", bp2PhosLossScore );

				scoring = getScoreValue ( params, "yp2_score", yp2Score );
				scoring = getScoreValue ( params, "yp2_loss_score", yp2LossScore );
				scoring = getScoreValue ( params, "yp2_phos_loss_score", yp2PhosLossScore );

				scoring = getScoreValue ( params, "internal_a_score", internalAScore );
				scoring = getScoreValue ( params, "internal_b_score", internalBScore );
				scoring = getScoreValue ( params, "internal_loss_score", internalLossScore );

				scoring = getScoreValue ( params, "mh3po4_score", MH3PO4Score );
				scoring = getScoreValue ( params, "msoch4_score", MSOCH4Score );
				flag = true;
				break;
			}
		}
		if ( !flag ) ErrorHandler::genError ()->error ( "Invalid or unspecified fragmentation name (file: fragmentation.txt).\nFunction: DiscriminantScore.\n" );
	}
	setMasks ();
}
void MSMSFragmentation::setMasks ()
{
	ammoniaLossMask		= string_to_mask ( ammoniaLoss );
	waterLossMask		= string_to_mask ( waterLoss );
	posChargeBearingMask= string_to_mask ( posChargeBearing );
	dIonExcludeMask		= string_to_mask ( dIonExclude );
	vIonExcludeMask		= string_to_mask ( vIonExclude );
	wIonExcludeMask		= string_to_mask ( wIonExclude );
}
bool MSMSFragmentation::getScoreValue ( const ParameterList& params, const string& name, ScoreType& value )
{
	double dvalue;
	bool flag = params.getValue ( name, dvalue );
	if ( flag ) {
		dvalue *= SCORE_TYPE_MULTIPLIER;
		value = static_cast <ScoreType> ( dvalue );
	}
	return flag;
}
DoubleVector MSMSFragmentation::getLossMasses ()
{
	DoubleVector dv;
	dv.push_back ( 0.0 );
	for ( int i = 0 ; i < ionTypes.size () ; i++ ) {
		string it = ionTypes [i];
		if ( isPrefix ( it, "M+" ) || isPrefix ( it, "M-" ) ) {
			string loss = it.substr ( 1 );
			if ( loss.length () > 1 && isdigit ( loss [1] ) ) {
				dv.push_back ( atof ( loss.c_str () ) );
			}
		}
	}
	return dv;
}
FragmentationParams::FragmentationParams () :
	GenNameValueStream ( MsparamsDir::instance ().getParamPath ( "fragmentation.txt" ) )
{
}
FragmentationParams& FragmentationParams::instance ( bool reset )
{
	static FragmentationParams* d = new FragmentationParams ();
	if ( reset ) {
		delete d;
		d = new FragmentationParams ();
	}
	return *d;
}
