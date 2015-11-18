/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_usermod.h                                                  *
*                                                                             *
*  Created    : May 30th 2001                                                 *
*                                                                             *
*  Purpose    : Gets the information on user defined modifications from the   *
*               params/usermod_*.txt files.                                   *
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

#ifndef __lu_usermod_h
#define __lu_usermod_h

#include <algorithm>
#include <string>
#include <iostream>
#include <lu_formula.h>

struct SingleUserMod;
class ParameterList;
class Peak;
class AAInitInfo;

typedef std::vector <SingleUserMod*> SingleUserModPtrVector;
typedef SingleUserModPtrVector::size_type SingleUserModPtrVectorSizeType;

class Usermod {
	static MapStringToInt indexMap;
	static SingleUserModPtrVector singleMods;
	static StringVector names;

	static bool initialised;
	static bool productFlag;
	static SetString glycoMods;

	std::string outputString;
	std::string formulaStr;
	ElementalFormula formula;
	unsigned int mask;
	char terminalSpecificity;
	std::string aaList;
	static void readUsermodFile ( const std::string& fName );
	static void initialiseUsermod ();
public:
	Usermod ( const std::string& outputString, const std::string& formulaStr, char terminalSpecificity, const std::string& aaList );
	Usermod ( const std::string& n );
	Usermod ( const std::string& n, bool dummy );
	static void initialiseUsermodAAInfo ( const StringVector& nStr );
	static void initialiseAllUsermodAAInfo ();
	static MapStringToInt getN15Balance ();

	std::string getOutputString () const { return outputString; }
	ElementalFormula getElementalFormula () const { return formula; }
	std::string getElementalFormulaString () const { return formulaStr; }
	unsigned int getMask () const { return mask; }
	char getTerminalSpecificity () const { return terminalSpecificity; }
	std::string getAAList () const { return aaList; }
	static bool isGlyco ( const std::string& name );
	static StringVector getNames ()
	{
		if ( initialised == false ) initialiseUsermod ();
		return names;
	}
	static StringVector getNames ( const std::string& type );
	static StringVector getMSProdNames ();
	static StringVector getSCompNames ();
	static bool isUsermod ( const std::string& s )
	{
		if ( initialised == false ) initialiseUsermod ();
		return std::find ( names.begin (), names.end (), s ) != names.end ();
	}
	static bool setProductFlag () { productFlag = true; }
	static void parseUsermodSpecificityLine ( const std::string& line, std::string& aaList, char& terminalSpecificity );
	static void getAllUsermodInfo ( StringVector& namesList, StringVector& specificityList, StringVector& elementalFormulaStrList, StringVector& typeList );
	static StringVector getTypeOptions ();
	static std::string getPhosphoSTYMenuString ();
};

class MultipleModification2 {
	std::vector <Usermod*> userMod;
	StringVector modAANames;
public:
	MultipleModification2 ( const ParameterList* pList );
	MultipleModification2 ( const ParameterList* pList, bool flag );
	void setUserMods ( const AAInitInfo& aaInitInfo ) const;
	void printHTML ( std::ostream& os ) const;
	void putCGI ( std::ostream& os ) const;
	void putHiddenFormEntry ( std::ostream& os ) const;
	void putHiddenFormJavascriptEntry ( std::ostream& os ) const;
	static void copyToCGI ( std::ostream& os, const ParameterList* params );
};

class MultipleModification {
	bool pyroglutamicAcidFlag;
	bool oxidationFlag;
	bool acetylationFlag;
	bool incompleteCysFlag;
	std::vector <Usermod*> userMod;
	StringVector modAANames;
	std::string user1Name;
	std::string user2Name;
	std::string user3Name;
	std::string user4Name;
public:
	MultipleModification ( const ParameterList* pList );
	bool getPyroglutamicAcidFlag ()	const { return pyroglutamicAcidFlag; }
	bool getOxidationFlag ()		const { return oxidationFlag; }
	bool getAcetylationFlag ()		const { return acetylationFlag; }
	bool getIncompleteCysFlag ()	const { return incompleteCysFlag; }
	bool getUser1Flag ()			const { return userMod.size () > 0; }
	bool getUser2Flag ()			const { return userMod.size () > 1; }
	bool getUser3Flag ()			const { return userMod.size () > 2; }
	bool getUser4Flag ()			const { return userMod.size () > 3; }
	bool getModificationsAllowed () const
	{
		return ( modAANames.empty () == false );
	}
	std::vector <Usermod*> getUserMods () const { return userMod; }
	int getNumUserMods () const { return userMod.size (); }
	unsigned int getUserMask ( int i ) const { return userMod [i]->getMask (); }
	void printHTML ( std::ostream& os ) const;
	void putCGI ( std::ostream& os ) const;
	void putHiddenFormEntry ( std::ostream& os ) const;
	static void copyToCGI ( std::ostream& os, const ParameterList* params );
	static void copyToHiddenFormEntry ( std::ostream& os, const ParameterList* params );
	static void copyToHiddenFormJavascriptEntry ( std::ostream& os, const ParameterList* params );
};

class PossFragMods {
	int pyroglutamicAcid;
	int oxidation;
	int acetylation;
	int incompleteCys;
	IntVector user;

	static bool pyroglutamicAcidFlag;
	static bool oxidationFlag;
	static bool acetylationFlag;
	static bool incompleteCysFlag;
	static unsigned int acetylationMask;
	static unsigned int pyroglutamicAcidMask;
	static unsigned int oxidationMask;
	static unsigned int cysMask;
	static UIntVector userMask;
	static int numUserMods;

	static bool initialised;
public:
	PossFragMods () {};
	PossFragMods ( const std::string& fragment, bool firstFragment );
	static void initialiseMasks ( const MultipleModification& multiMod );
	static unsigned int getAcetylationMask () { return acetylationMask; }
	static unsigned int getPyroglutamicAcidMask () { return pyroglutamicAcidMask; }
	static unsigned int getOxidationMask () { return oxidationMask; }
	static unsigned int getCysMask () { return cysMask; }
	static unsigned int getUser1Mask () { return userMask.empty () ? 0 : userMask [0]; }
	static unsigned int getUser2Mask () { return userMask.size () < 2 ? 0 : userMask [1]; }
	static unsigned int getUser3Mask () { return userMask.size () < 3 ? 0 : userMask [2]; }
	static unsigned int getUser4Mask () { return userMask.size () < 4 ? 0 : userMask [3]; }
};

class FragModContainer {
	bool firstFragment;
	bool lastFragment;
	int nTermIdx;
	int cTermIdx;
	std::string aasModified;
	std::set <MapIntToInt> modList;	// Set of maps of modification indexes to number of times they occur
	int modListSize;
	int consideredListSize;
	std::map <MapIntToInt, std::pair <int, IntVector> > numXLMap;
	MapIntToInt curMod;				// Maps the modification indicies to the number of times they occur
	bool curNTerm;
	bool curCTerm;
	bool curNLoss;
	static unsigned char protNTermMask;
	static unsigned char protCTermMask;
	static unsigned char pepNTermMask;
	static unsigned char pepCTermMask;
	static unsigned char nLossMask;
	static unsigned char anyMask;
	static IntVector xLinkIdx;
	static MapCharToUChar modAAs;
	static char enzTermSpec;
	static std::map <char, std::vector <std::pair <int, char> > > modMap;	// Map of amino acids to list of indexes to possible modifications
	static DoubleVector userMod;
	static DoubleVector userModMi;
	static DoubleVector userModAv;
	static double minMod;
	static double minModMi;
	static double minModAv;
	static double maxMod;
	static double maxModMi;
	static double maxModAv;
	static std::vector <ElementalFormula> userModElemForm;
	static StringVector userModOutputString;
	void getNextHit ( int level, int start );
	bool checkTermSpec ( char ts, int level );
	bool checkMultiTerm ( char ts );
	std::set <MapIntToInt>::const_iterator cur;
	static unsigned char termSpecToMask ( char ts );
public:
	FragModContainer ( const std::string& fragment, bool firstFragment, bool lastFragment );
	void first () { cur = modList.begin (); }
	bool isDone () const { return cur != modList.end (); }
	void next () { cur++; }
	MapIntToInt getModList () const { return *cur; }
	static void setUserMods ( const std::vector <Usermod*>& userMod, const AAInitInfo& aaInitInfo );
	static void setMassType ( bool flag );
	static double getMass ( const MapIntToInt& mii );
	static double getMonoisotopicMass ( const MapIntToInt& mii );
	static double getAverageMass ( const MapIntToInt& mii );
	static double getMinMod () { return minMod; }
	static double getMinModMi () { return minModMi; }
	static double getMinModAv () { return minModAv; }
	static double getMaxMod () { return maxMod; }
	static double getMaxModMi () { return maxModMi; }
	static double getMaxModAv () { return maxModAv; }
	std::pair <int, IntVector> countBridgeSites () const;
	static bool containsMods ( const MapIntToInt& mii, const MapStringToInt& mods );
	static bool isExactMod ( const MapIntToInt& mii, const MapStringToInt& mods );
	static ElementalFormula getElementalFormula ( const MapIntToInt& mii );
	static bool getOxidation ( const MapIntToInt& mii );
	static int getNumOxidation ( const MapIntToInt& mii );
	static bool getModified ( const MapIntToInt& mii );
	static void printHTML ( std::ostream& os, const MapIntToInt& mii );
	static void printDelimited ( std::ostream& os, const MapIntToInt& mii );
	static void printXML ( std::ostream& os, const MapIntToInt& mii );
};

class ExtraUserModInfo {
	std::string formula;
	std::string userAMass;
	std::string limit;
	std::string motifOffset;
	std::string motif;
public:
	ExtraUserModInfo ();
	ExtraUserModInfo ( const std::string& formula, const std::string& userAMass, const std::string& limit, const std::string& motifOffset, const std::string& motif );
	std::string getFormula () const { return formula; }
	std::string getUserAMass () const { return userAMass; }
	std::string getLimit () const { return limit; }
	std::string getMotifOffset () const { return motifOffset; }
	std::string getMotif () const { return motif; }
};

typedef std::map <PairStringString, ExtraUserModInfo> MapPairStringStringExtraUserModInfo;
typedef MapPairStringStringExtraUserModInfo::const_iterator MapPairStringStringExtraUserModInfoConstIterator;

class ExtraUserMods {
	MapPairStringStringExtraUserModInfo extraUserMods;
	ExtraUserMods ();
public:
	static ExtraUserMods& instance ();
	void addUserMods ( const ParameterList* params );
	void addUserMods2 ( const ParameterList* params, int maxUserMods );
	int size () const { extraUserMods.size (); }
	bool empty () const { extraUserMods.empty (); }
	MapPairStringStringExtraUserModInfo getExtraUserMods () const { return extraUserMods; }
	void copyToCGI ( std::ostream& os );
};

#endif /* ! __lu_usermod_h */
