/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_get_link.h                                                 *
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

#ifndef __lu_get_link_h
#define __lu_get_link_h

#include <string>
#include <lgen_define.h>
#include <lu_fragmentation.h>

class Usermod;
class ParameterList;
class AAInitInfo;

class LinkInfo {
	std::string name;
	int maxLinkMolecules;
	int bIndex;
	std::string linkAminoAcid1;
	std::string linkAminoAcid2;
	std::string bridgeFormula;
	StringVector userLabels;
	StringVector userFormulae;
	StringVector userAAs;
	std::vector <Usermod*> userMod;

	ScoreType maxPScore;
	std::string pCIDImmFormula;
	ScoreType pCIDScore;
	ScoreType pCIDXLScore;
	ScoreType pCIDImmScore;
	ScoreType pCIDH2OScore;
	ScoreType pCIDNH3Score;
	ScoreType pCID2H2OScore;
	ScoreType pCIDXLH2OScore;
	ScoreType pCIDXLNH3Score;
	ScoreType pCIDImmH2OScore;
	ScoreType pCIDImmNH3Score;
	ScoreType pETDScore;
	ScoreType pETDXLScore;

	void initialiseInfo ( const std::string& name );
	void initialiseUserInfo ( const ParameterList* params, bool bridgeSearch );
	void addUserMod ( const ParameterList* params, int num, SetString& uniqLabels );
	void initialiseCrossLinks ( const std::string& str, int n );
	static std::string NO_LINK;
	static std::string USER_DEFINED_LINK;
	static bool initialised;

	static StringVector nameList;
	static StringVector linkAA1List;
	static StringVector linkAA2List;
	static StringVector bridgeFormulaList;
	static StringVectorVector usermodList;

	static ScoreTypeVector maxPScoreList;
	static StringVector pCIDImmFormulaList;
	static ScoreTypeVector pCIDScoreList;
	static ScoreTypeVector pCIDXLScoreList;
	static ScoreTypeVector pCIDImmScoreList;
	static ScoreTypeVector pCIDH2OScoreList;
	static ScoreTypeVector pCIDNH3ScoreList;
	static ScoreTypeVector pCID2H2OScoreList;
	static ScoreTypeVector pCIDXLH2OScoreList;
	static ScoreTypeVector pCIDXLNH3ScoreList;
	static ScoreTypeVector pCIDImmH2OScoreList;
	static ScoreTypeVector pCIDImmNH3ScoreList;
	static ScoreTypeVector pETDScoreList;
	static ScoreTypeVector pETDXLScoreList;

	static void initialiseLinks ();
public:
	LinkInfo ();
	LinkInfo ( const ParameterList* params, bool bridgeSearch = false );
	void printHTML ( std::ostream& os ) const;
	std::string getName () const { return name; }
	std::string getLinkAminoAcid1 () const { return linkAminoAcid1; }
	std::string getLinkAminoAcid2 () const { return linkAminoAcid2; }
	std::string getBridgeFormula () const { return bridgeFormula; }

	ScoreType getMaxPScore () const			{ return maxPScore; }
	std::string getPCIDImmFormula () const	{ return pCIDImmFormula; }
	ScoreType getPCIDScore () const			{ return pCIDScore; }
	ScoreType getPCIDXLScore () const		{ return pCIDXLScore; }
	ScoreType getPCIDImmScore () const		{ return pCIDImmScore; }
	ScoreType getPCIDH2OScore () const		{ return pCIDH2OScore; }
	ScoreType getPCIDNH3Score () const		{ return pCIDNH3Score; }
	ScoreType getPCID2H2OScore () const		{ return pCID2H2OScore; }
	ScoreType getPCIDXLH2OScore () const	{ return pCIDXLH2OScore; }
	ScoreType getPCIDXLNH3Score () const	{ return pCIDXLNH3Score; }
	ScoreType getPCIDImmH2OScore () const	{ return pCIDImmH2OScore; }
	ScoreType getPCIDImmNH3Score () const	{ return pCIDImmNH3Score; }
	ScoreType getPETDScore () const			{ return pETDScore; }
	ScoreType getPETDXLScore () const		{ return pETDXLScore; }

	bool getNoLink () const { return name == NO_LINK || maxLinkMolecules == 0; }
	int getMaxLinkMolecules () const { return maxLinkMolecules; }
	int getBIndex () const { return bIndex; }
	static StringVector getLinkFormulae ();
	static DoubleVector getLinkMasses ();
	static StringVector getNameList ();
	static StringVector getLinkAAs ( int num );
	static StringVector getLinkAAs ( const std::string& str );
	static StringVector getLinkAAs ( const std::string& s1, const std::string& s2 );
	static StringVector getUsermods ( int num );
	static StringVector getAllUsermods ();
	void setUserMods ( const AAInitInfo& aaInitInfo ) const;
};

#endif /* ! __lu_get_link_h */
