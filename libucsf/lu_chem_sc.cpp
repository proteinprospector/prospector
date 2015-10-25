/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_chem_sc.cpp                                                *
*                                                                             *
*  Created    : February 16th 2002                                            *
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
#ifdef CHEM_SCORE
#include <cmath>
#include <lgen_reg_exp.h>
#include <lu_chem_sc.h>
#include <lu_getfil.h>
#include <lu_mod_frag.h>
#include <lu_links.h>
#include <lu_param_list.h>
#include <lu_const_mod.h>
using std::string;
using std::ostream;
using std::istringstream;

ChemScoreParameters::ChemScoreParameters ( const ParameterList* params ) :
	chemScoreFlag	( params->getBoolValue ( "chem_score", false ) ),
	metOxFactor		( params->getDoubleValue ( "met_ox_factor", 1.0 ) )
{
	StringVector constModNames;
	params->getValue ( "const_mod", constModNames );
	cysMod = "";
	for ( StringVectorSizeType i = 0 ; i < constModNames.size () ; i++ ) {
		ConstMod constMod ( constModNames [i] );
		if ( constMod.getAAList ().find ( 'C' ) != string::npos ) {
			cysMod = constMod.getLongName ();
			break;
		}
	}
	if ( chemScoreFlag ) {
		PotentialMSFragment::setChemScore ( cysMod, metOxFactor );
		PotentialLinkFragment::setChemScore ( cysMod, metOxFactor );
	}
}
void ChemScoreParameters::printHTML ( ostream& os ) const
{
	if ( chemScoreFlag ) {
		ParameterList::printHTML ( os, "Chem Score", "true" );
		if ( !cysMod.empty () ) ParameterList::printHTML ( os, "Cys Mod", cysMod );
		ParameterList::printHTML ( os, "MetOx Factor", metOxFactor );
	}
}
void ChemScoreParameters::copyToCGI ( ostream& os, const ParameterList* params )
{
	if ( params->copyToCGI ( os, "chem_score" ) ) {
		params->copyToCGI ( os, "met_ox_factor" );
	}
}
void ChemScoreParameters::copyToHiddenFormEntry ( ostream& os, const ParameterList* params )
{
	if ( params->copyToHiddenFormEntry ( os, "chem_score" ) ) {
		params->copyToHiddenFormEntry ( os, "met_ox_factor" );
	}
}
ScorePair::ScorePair ( const string& exp, int alignment, double score ) :
	regExp ( new RegularExpression ( exp ) ),
	alignment ( alignment ),
	score ( score )
{
}
ChemScore::ChemScore ( const string& cysMod, double metOxF ) :
	metOxF ( metOxF )
{
	GenNameValueStream nvs ( MsparamsDir::instance ().getParamPath ( "chem_score.txt" ) );
	double argScore = nvs.getDoubleValue ( "arg_score" );
	double lysScore = nvs.getDoubleValue ( "lys_score" );
	double pyridylethylScore = nvs.getDoubleValue ( "pyridylethyl_score" );
	double detrimentalCysScore = nvs.getDoubleValue ( "detrimental_cys_score" );
	StringVector missedCleavageRules;
	nvs.getValue ( "missed_cleavage_rule", missedCleavageRules );
	StringVector terminalPeptideAdjustments;
	nvs.getValue ( "terminal_peptide_adjustment", terminalPeptideAdjustments );

	initialSettingsSP.push_back ( ScorePair ( "R", 0, argScore ) );
	if ( !genStrcasecmp ( cysMod.c_str (), "Pyridylethyl" ) ) {
		initialSettingsSP.push_back ( ScorePair ( "C", 0, pyridylethylScore ) );
	}
	initialSettingsSP.push_back ( ScorePair ( "K", 0, lysScore ) );

	detrimentalCys = 1.0;
	if ( cysMod == "" || !genStrcasecmp ( cysMod.c_str (), "Propionamide" ) ) {
		detrimentalCys = detrimentalCysScore;
	}

	locator = new RegularExpression ( "[KR][^P]" );

	for ( StringVectorSizeType i = 0 ; i < missedCleavageRules.size () ; i++ ) {
		istringstream istr ( missedCleavageRules [i] );
		string re;
		int align;
		double sc;
		istr >> re;
		istr >> align;
		istr >> sc;
		missedCleavagesSP.push_back ( ScorePair ( re, align, sc ) );
	}
	for ( StringVectorSizeType j = 0 ; j < terminalPeptideAdjustments.size () ; j++ ) {
		istringstream istr ( terminalPeptideAdjustments [j] );
		string re;
		int align;
		double sc;
		istr >> re;
		istr >> align;
		istr >> sc;
		terminalAdjustmentSP.push_back ( ScorePair ( re, align, sc ) );
	}
}
double ChemScore::getChemScore ( const string& fragment, int numMetOx ) const
{
	double score = 1.0;
	const char* cfragment = fragment.c_str ();
	for ( ScorePairVectorSizeType i = 0 ; i < initialSettingsSP.size () ; i++ ) {
		if ( initialSettingsSP [i].getRegExp ()->isPresent ( cfragment ) ) {
			score = initialSettingsSP [i].getScore ();
			break;
		}
	}
	int numCys = gen_strcharcount ( fragment, 'C' );
	if ( numCys ) {
		score /= detrimentalCys;
	}
	int numMet = gen_strcharcount ( fragment, 'M' );
	if ( numMet ) {
		int numMetRed = numMet - numMetOx;
		if ( metOxF == 1.0 ) score /= 2.0;
		if ( metOxF > 1.0 ) score /= pow ( metOxF, numMetRed );
		if ( metOxF < 1.0 ) score *= pow ( metOxF, numMetOx );
	}

	while ( locator->isPresentMulti ( cfragment ) ) {	// For each missed cleavage
		const char* cleavagePoint = locator->getLoc1 ();
		double mcf = 1.0;
		for ( ScorePairVectorSizeType j = 0 ; j < missedCleavagesSP.size () ; j++ ) {
			RegularExpression* currentRegExp = missedCleavagesSP [j].getRegExp ();
			if ( currentRegExp->isPresentMultiOverlap ( cfragment ) ) {
				const char* matchPosn = currentRegExp->getLoc1 ();
				if ( matchPosn - cleavagePoint == missedCleavagesSP [j].getAlignment () ) {
					mcf *= missedCleavagesSP [j].getScore ();
				}
				else {
					currentRegExp->resetMulti ();
				}
			}
		}
		double bmcf = 100.0;
		double chpf = ( bmcf + mcf ) / mcf;
		score /= chpf;
	}
	return score;
}
#endif
