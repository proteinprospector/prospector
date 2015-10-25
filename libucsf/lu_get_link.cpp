/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_get_link.cpp                                               *
*                                                                             *
*  Created    : February 18th 2002                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2002-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lg_string.h>
#include <lgen_error.h>
#include <lu_getfil.h>
#include <lu_get_link.h>
#include <lu_mass_elem.h>
#include <lu_param_list.h>
#include <lu_srch_form.h>
#include <lu_usermod.h>
using std::pair;
using std::string;
using std::ostream;
using std::getline;
using std::istringstream;
using std::copy;
using std::back_inserter;

string LinkInfo::NO_LINK = "No Link";
string LinkInfo::USER_DEFINED_LINK = "User Defined Link";
bool LinkInfo::initialised = false;
StringVector LinkInfo::nameList;
StringVector LinkInfo::linkAA1List;
StringVector LinkInfo::linkAA2List;
StringVector LinkInfo::bridgeFormulaList;
StringVectorVector LinkInfo::usermodList;

ScoreTypeVector LinkInfo::maxPScoreList;
StringVector LinkInfo::pCIDImmFormulaList;
ScoreTypeVector LinkInfo::pCIDScoreList;
ScoreTypeVector LinkInfo::pCIDXLScoreList;
ScoreTypeVector LinkInfo::pCIDImmScoreList;
ScoreTypeVector LinkInfo::pCIDH2OScoreList;
ScoreTypeVector LinkInfo::pCIDNH3ScoreList;
ScoreTypeVector LinkInfo::pCID2H2OScoreList;
ScoreTypeVector LinkInfo::pCIDXLH2OScoreList;
ScoreTypeVector LinkInfo::pCIDXLNH3ScoreList;
ScoreTypeVector LinkInfo::pCIDImmH2OScoreList;
ScoreTypeVector LinkInfo::pCIDImmNH3ScoreList;
ScoreTypeVector LinkInfo::pETDScoreList;
ScoreTypeVector LinkInfo::pETDXLScoreList;

LinkInfo::LinkInfo () :
	name ( NO_LINK ),
	maxLinkMolecules ( 0 ),
	bIndex ( -1 )
{
}
LinkInfo::LinkInfo ( const ParameterList* params, bool bridgeSearch ) :
	name ( params->getStringValue( "link_search_type", NO_LINK ) ),
	maxLinkMolecules ( params->getIntValue ( "max_link_molecules", 5 ) ),
	bIndex ( -1 )
{
	if ( name == USER_DEFINED_LINK ) {
		initialiseUserInfo ( params, bridgeSearch );
	}
	else if ( name != NO_LINK ) {
		initialiseInfo ( name );
	}
}
void LinkInfo::printHTML ( ostream& os ) const
{
	if ( name != "No Link" ) {
		ParameterList::printHTML ( os, "Link Search Type", name );
		ParameterList::printHTML ( os, "Maximum Link Molecules", maxLinkMolecules );
	}
	if ( name == "User Defined Link" ) {
		ParameterList::printHTML ( os, "Bridge Formula", bridgeFormula );
		for ( StringVectorSizeType i = 0 ; i < userLabels.size () ; i++ ) {
			string sNum = gen_itoa ( i + 1 );
			ParameterList::printHTML ( os, "Modification " + sNum + " Label", userLabels [i] );
			ParameterList::printHTML ( os, "Modification " + sNum + " Formula", userFormulae [i] );
			ParameterList::printHTML ( os, "Modification " + sNum + " AAs", userAAs [i] );
		}
	}
}
void LinkInfo::initialiseCrossLinks ( const string& str, int n )
{
	StringVector sv = genGetSeparatedValues ( str, "," );
	string nStr = gen_itoa ( n );
	for ( StringVectorSizeType i = 0 ; i < sv.size () ; i++ ) {
		if ( sv [i] == "Protein N-term" )
			userMod.push_back ( new Usermod ( "Cross Link " + nStr, "", 'n', "n" ) );
		else if ( sv [i] == "Protein C-term" )
			userMod.push_back ( new Usermod ( "Cross Link " + nStr, "", 'c', "c" ) );
		else
			userMod.push_back ( new Usermod ( "Cross Link " + nStr, "", 'e', sv [i] ) );
	}
}
void LinkInfo::initialiseInfo ( const string& name ) 
{
	if ( !initialised ) initialiseLinks ();

	for ( StringVectorSizeType i = 0 ; i < nameList.size () ; i++ ) {
		if ( nameList [i] == name ) {
			string linkAA1 = linkAA1List [i];
			string linkAA2 = linkAA2List [i];
			bridgeFormula = bridgeFormulaList [i];
			linkAminoAcid1 = linkAA1;
			linkAminoAcid2 = linkAA2;
			for ( StringVectorSizeType j = 0 ; j < usermodList [i].size () ; j++ ) {
				userMod.push_back ( new Usermod ( usermodList [i][j] ) );
			}
			initialiseCrossLinks ( linkAA1, 1 );
			if ( linkAA2 != linkAA1 ) {
				bIndex = 1;
				initialiseCrossLinks ( linkAA2, 2 );
			}
			maxPScore		= maxPScoreList [i];
			pCIDImmFormula	= pCIDImmFormulaList [i];
			pCIDScore		= pCIDScoreList [i];
			pCIDXLScore		= pCIDXLScoreList [i];
			pCIDImmScore	= pCIDImmScoreList [i];
			pCIDH2OScore	= pCIDH2OScoreList [i];
			pCIDNH3Score	= pCIDNH3ScoreList [i];
			pCID2H2OScore	= pCID2H2OScoreList [i];
			pCIDXLH2OScore	= pCIDXLH2OScoreList [i];
			pCIDXLNH3Score	= pCIDXLNH3ScoreList [i];
			pCIDImmH2OScore	= pCIDImmH2OScoreList [i];
			pCIDImmNH3Score	= pCIDImmNH3ScoreList [i];
			pETDScore		= pETDScoreList [i];
			pETDXLScore		= pETDXLScoreList [i];
			return;
		}
	}
	ErrorHandler::genError ()->error ( "Invalid or unspecified link search type.\nLinkInfo::initialiseInfo.\n" );
}
void LinkInfo::initialiseUserInfo ( const ParameterList* params, bool bridgeSearch )
{
	string linkAA = params->getStringValue ( "link_aa" );
	StringVector sv = genGetSeparatedValues ( linkAA, "->" );
	string linkAA1;
	string linkAA2;
	if ( sv.size () >= 2 ) {
		linkAA1 = sv [0];
		linkAA2 = sv [1];
	}
	bridgeFormula = params->getStringValue( "bridge_composition", "" );
	linkAminoAcid1 = linkAA1;
	linkAminoAcid2 = linkAA2;

	initialiseCrossLinks ( linkAA1, 1 );
	if ( linkAA2 != linkAA1 ) {
		bIndex = 1;
		initialiseCrossLinks ( linkAA2, 2 );
	}
	if ( bridgeSearch ) {
		SetString uniqLabels;
		for ( int i = 1 ; i <= CrosslinkingForm::getNumUserMods () ; i++ ) {
			addUserMod ( params, i, uniqLabels );
		}
	}
	maxPScore		= 100.0 * SCORE_TYPE_MULTIPLIER;	// Set to a value so all p ions are counted by default
	pCIDImmFormula	= "";
	pCIDScore		= DEFAULT_MISS_SCORE;
	pCIDXLScore		= DEFAULT_MISS_SCORE;
	pCIDImmScore	= DEFAULT_MISS_SCORE;
	pCIDH2OScore	= DEFAULT_MISS_SCORE;
	pCIDNH3Score	= DEFAULT_MISS_SCORE;
	pCID2H2OScore	= DEFAULT_MISS_SCORE;
	pCIDXLH2OScore	= DEFAULT_MISS_SCORE;
	pCIDXLNH3Score	= DEFAULT_MISS_SCORE;
	pCIDImmH2OScore	= DEFAULT_MISS_SCORE;
	pCIDImmNH3Score	= DEFAULT_MISS_SCORE;
	pETDScore		= DEFAULT_MISS_SCORE;
	pETDXLScore		= DEFAULT_MISS_SCORE;
}
void LinkInfo::addUserMod ( const ParameterList* params, int num, SetString& uniqLabels )
{
	string sNum = gen_itoa ( num );
	string label = params->getStringValue ( "mod_" + sNum + "_label", "" );
	if ( !label.empty () ) {
		string comp = params->getStringValue ( "mod_" + sNum + "_composition", "" );
		if ( !comp.empty () ) {
			string aaMod = params->getStringValue ( "aa_modified_" + sNum, "" );
			pair <SetStringIterator, bool> flag = uniqLabels.insert ( label+aaMod );
			if ( !flag.second ) {
				ErrorHandler::genError ()->error ( "Mod Labels need to be unique.\n" );
			}
			if ( aaMod == "Protein N-term" )
				userMod.push_back ( new Usermod ( label, comp, 'n', "n" ) );
			else if ( aaMod == "Protein C-term" )
				userMod.push_back ( new Usermod ( label, comp, 'c', "c" ) );
			else
				userMod.push_back ( new Usermod ( label, comp, 'e', aaMod ) );
			userLabels.push_back ( label );
			userFormulae.push_back ( comp );
			userAAs.push_back ( aaMod );
		}
		else
			ErrorHandler::genError ()->error ( "No elemental composition specified for Mod " + sNum + " Label.\n" );
	}
}
void LinkInfo::initialiseLinks () 
{
	GenIFStream fromFile ( MsparamsDir::instance ().getParamPath ( "links.txt" ) );
	string line;
	bool menuItem = true;
	string l1AA;
	string l2AA;
	string bf;
	StringVector u;

	ScoreType maxPScore = 100.0 * SCORE_TYPE_MULTIPLIER;	// Set to a value so all p ions are counted by default;
	string pCIDImmFormula;
	ScoreType pCIDScore = DEFAULT_MISS_SCORE;
	ScoreType pCIDXLScore = DEFAULT_MISS_SCORE;
	ScoreType pCIDImmScore = DEFAULT_MISS_SCORE;
	ScoreType pCIDH2OScore = DEFAULT_MISS_SCORE;
	ScoreType pCIDNH3Score = DEFAULT_MISS_SCORE;
	ScoreType pCID2H2OScore = DEFAULT_MISS_SCORE;
	ScoreType pCIDXLH2OScore = DEFAULT_MISS_SCORE;
	ScoreType pCIDXLNH3Score = DEFAULT_MISS_SCORE;
	ScoreType pCIDImmH2OScore = DEFAULT_MISS_SCORE;
	ScoreType pCIDImmNH3Score = DEFAULT_MISS_SCORE;
	ScoreType pETDScore = DEFAULT_MISS_SCORE;
	ScoreType pETDXLScore = DEFAULT_MISS_SCORE;

	while ( getline ( fromFile, line ) ) {
		if ( line.length () != 0 && line [0] != '#' ) {
			if ( menuItem ) {
				nameList.push_back ( line );
				menuItem = false;
			}
			if ( line != ">" ) {
				istringstream ist ( line );
				string name;
				ist >> name;
				if ( !name.empty () ) {
					string value = gen_strtrim ( line.substr ( name.length () ) );
					if ( name == "link_aa_1" )		l1AA = value;
					if ( name == "link_aa_2" )		l2AA = value;
					if ( name == "bridge_formula" )	bf = value;
					if ( name == "usermod" )		u.push_back ( value );

					if ( name == "max_p_score" )			maxPScore		= atof ( value.c_str () ) * SCORE_TYPE_MULTIPLIER;
					if ( name == "p_cid_imm_formula" )		pCIDImmFormula	= value;
					if ( name == "p_cid_score" )			pCIDScore		= atof ( value.c_str () ) * SCORE_TYPE_MULTIPLIER;
					if ( name == "p_cid_xl_score" )				pCIDXLScore		= atof ( value.c_str () ) * SCORE_TYPE_MULTIPLIER;
					if ( name == "p_cid_imm_score" )		pCIDImmScore	= atof ( value.c_str () ) * SCORE_TYPE_MULTIPLIER;
					if ( name == "p_cid_h2o_score" )		pCIDH2OScore	= atof ( value.c_str () ) * SCORE_TYPE_MULTIPLIER;
					if ( name == "p_cid_nh3_score" )		pCIDNH3Score	= atof ( value.c_str () ) * SCORE_TYPE_MULTIPLIER;
					if ( name == "p_cid_2h2o_score" )		pCID2H2OScore	= atof ( value.c_str () ) * SCORE_TYPE_MULTIPLIER;
					if ( name == "p_cid_xl_h2o_score" )			pCIDXLH2OScore	= atof ( value.c_str () ) * SCORE_TYPE_MULTIPLIER;
					if ( name == "p_cid_xl_nh3_score" )			pCIDXLNH3Score	= atof ( value.c_str () ) * SCORE_TYPE_MULTIPLIER;
					if ( name == "p_cid_imm_h2o_score" )		pCIDImmH2OScore	= atof ( value.c_str () ) * SCORE_TYPE_MULTIPLIER;
					if ( name == "p_cid_imm_nh3_score" )		pCIDImmNH3Score	= atof ( value.c_str () ) * SCORE_TYPE_MULTIPLIER;

					if ( name == "p_etd_score" )			pETDScore		= atof ( value.c_str () ) * SCORE_TYPE_MULTIPLIER;
					if ( name == "p_etd_xl_score" )				pETDXLScore		= atof ( value.c_str () ) * SCORE_TYPE_MULTIPLIER;
				}
			}
			else {
				if ( !l1AA.empty () )	linkAA1List.push_back ( l1AA );
				else ErrorHandler::genError ()->error ( "Missing parameter link_aa_1 in parameter file links.txt.\n" );
				if ( !l2AA.empty () )	linkAA2List.push_back ( l2AA );
				else ErrorHandler::genError ()->error ( "Missing parameter link_aa_2 in parameter file links.txt.\n" );
				if ( !bf.empty () )		bridgeFormulaList.push_back ( bf );
				else ErrorHandler::genError ()->error ( "Missing parameter bridge_formula in parameter file links.txt..\n" );
				usermodList.push_back ( u );

				maxPScoreList.push_back			( maxPScore );
				pCIDImmFormulaList.push_back	( pCIDImmFormula );
				pCIDScoreList.push_back			( pCIDScore );
				pCIDXLScoreList.push_back		( pCIDXLScore );
				pCIDImmScoreList.push_back		( pCIDImmScore );
				pCIDH2OScoreList.push_back		( pCIDH2OScore );
				pCIDNH3ScoreList.push_back		( pCIDNH3Score );
				pCID2H2OScoreList.push_back		( pCID2H2OScore );
				pCIDXLH2OScoreList.push_back	( pCIDXLH2OScore );
				pCIDXLNH3ScoreList.push_back	( pCIDXLNH3Score );
				pCIDImmH2OScoreList.push_back	( pCIDImmH2OScore );
				pCIDImmNH3ScoreList.push_back	( pCIDImmNH3Score );
				pETDScoreList.push_back			( pETDScore );
				pETDXLScoreList.push_back		( pETDXLScore );

				l1AA = "";
				l2AA = "";
				bf = "";
				u.clear ();

				maxPScore		= 100.0 * SCORE_TYPE_MULTIPLIER;	// Set to a value so all p ions are counted by default
				pCIDImmFormula	= "";
				pCIDScore		= DEFAULT_MISS_SCORE;
				pCIDXLScore		= DEFAULT_MISS_SCORE;
				pCIDImmScore	= DEFAULT_MISS_SCORE;
				pCIDH2OScore	= DEFAULT_MISS_SCORE;
				pCIDNH3Score	= DEFAULT_MISS_SCORE;
				pCID2H2OScore	= DEFAULT_MISS_SCORE;
				pCIDXLH2OScore	= DEFAULT_MISS_SCORE;
				pCIDXLNH3Score	= DEFAULT_MISS_SCORE;
				pCIDImmH2OScore	= DEFAULT_MISS_SCORE;
				pCIDImmNH3Score	= DEFAULT_MISS_SCORE;
				pETDScore		= DEFAULT_MISS_SCORE;
				pETDXLScore		= DEFAULT_MISS_SCORE;

				menuItem = true;
			}
		}
	}
	initialised = true;
}
StringVector LinkInfo::getLinkFormulae ()
{
	if ( !initialised ) initialiseLinks ();
	StringVector sv;
	sv.push_back ( "" );
	copy ( bridgeFormulaList.begin (), bridgeFormulaList.end (), back_inserter ( sv ) );
	sv.push_back ( "" );
	return sv;
}
DoubleVector LinkInfo::getLinkMasses ()
{
	StringVector sv = LinkInfo::getLinkFormulae ();
	DoubleVector dv;
	for ( StringVectorSizeType i = 0 ; i < sv.size () ; i++ ) {
		dv.push_back ( formula_to_monoisotopic_mass ( sv [i].c_str () ) );
	}
	return dv;
}
StringVector LinkInfo::getNameList ()
{
	if ( !initialised ) initialiseLinks ();
	StringVector sv;
	sv.push_back ( NO_LINK );
	copy ( nameList.begin (), nameList.end (), back_inserter ( sv ) );
	sv.push_back ( USER_DEFINED_LINK );
	return sv;
}
StringVector LinkInfo::getLinkAAs ( int num )
{
	return getLinkAAs ( linkAA1List [num-1], linkAA2List [num-1] );
}
StringVector LinkInfo::getLinkAAs ( const string& str )
{
	StringVector sv = genGetSeparatedValues ( str, "->" );
	return getLinkAAs ( sv [0], sv [1] );
}
StringVector LinkInfo::getLinkAAs ( const string& s1, const string& s2 )
{
	SetString ss;
	StringVector sv1 = genGetSeparatedValues ( s1, "," );
	StringVector sv2 = genGetSeparatedValues ( s2, "," );

	for ( StringVectorSizeType i = 0 ; i < sv1.size () ; i++ )	ss.insert ( sv1 [i] );
	for ( StringVectorSizeType j = 0 ; j < sv2.size () ; j++ )	ss.insert ( sv2 [j] );

	StringVector sv;
	for ( SetStringConstIterator k = ss.begin () ; k != ss.end () ; k++ )	sv.push_back ( *k );
	return sv;
}
StringVector LinkInfo::getUsermods ( int num )
{
	return usermodList [num-1];
}
StringVector LinkInfo::getAllUsermods ()
{
	SetString ss;
	for ( StringVectorVectorSizeType i = 0 ; i < usermodList.size () ; i++ ) {
		for ( StringVectorSizeType j = 0 ; j < usermodList [i].size () ; j++ ) {
			ss.insert ( usermodList [i][j] );
		}
	}
	StringVector sv;
	for ( SetStringConstIterator k = ss.begin () ; k != ss.end () ; k++ )	sv.push_back ( *k );
	return sv;
}
void LinkInfo::setUserMods ( const AAInitInfo& aaInitInfo ) const
{
	FragModContainer::setUserMods ( userMod, aaInitInfo );
}
