/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_aa_info.h                                                  *
*                                                                             *
*  Created    : May 22nd 2001                                                 *
*                                                                             *
*  Purpose    : Gets amino acid information from the params/aa.txt file.      *                                                             *
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

#ifndef __lu_aa_info_h
#define __lu_aa_info_h

#include <map>
#include <string>
#include <lg_string.h>
#include <lgen_define.h>
#include <lgen_math.h>
#include <lu_formula.h>

typedef std::map <std::string, char, genStrcasecmpAscending> MapNoCaseStringToChar;
typedef MapNoCaseStringToChar::iterator MapNoCaseStringToCharIterator;
typedef MapNoCaseStringToChar::const_iterator MapNoCaseStringToCharConstIterator;

struct aaInfoNoElementalFormula {};

class AAInfo {
	CharVector aminoAcids;
	MapCharToString longNames;
	MapNoCaseStringToChar longNameToAA;
	MapCharToElementalFormula elementalFormula;
	MapCharToDouble monoisotopicMass;
	MapCharToDouble averageMass;
	MapCharToString substituentA;
	MapCharToString substituentB;
	MapCharToDouble pKCTerm;
	MapCharToDouble pKNTerm;
	MapCharToDouble pKAcidicSC;
	MapCharToDouble pKBasicSC;

	MapNoCaseStringToElementalFormula modElementalFormula;
	MapNoCaseStringToDouble modMonoisotopicMass;
	MapNoCaseStringToDouble modAverageMass;

	unsigned int aaConvMask;

	static double convert_pk ( const std::string& pk_string );

	static const int AMINO_ACID_LINES_PER_ENTRY;
	AAInfo ();
	bool getPlusSeparatedMod ( double& mod, const std::string& aa, const MapCharToDouble& mMap, const MapNoCaseStringToDouble& modMMap ) const;
public:
	static AAInfo& getInfo ();

	void addAminoAcid ( char aa, const ElementalFormula& formula );
	void addModAminoAcid ( char aa, const std::string& mod, const ElementalFormula& modFormula );
	double getMonoisotopicMass ( char aa ) const { return monoisotopicMass.find ( aa )->second; }
	double getModMonoisotopicMass ( const std::string& aa ) const;
	double getAverageMass ( char aa ) const { return averageMass.find ( aa )->second; }
	double getModAverageMass ( const std::string& aa ) const;
	ElementalFormula getElemental ( char aa ) const { return elementalFormula.find ( aa )->second; }
	std::string getElementalString ( char aa ) { return elementalFormula.find ( aa )->second.getFormula (); }
	std::string getModElementalString ( const std::string& aa )
	{
		MapNoCaseStringToElementalFormulaIterator i = modElementalFormula.find ( aa );
		if ( i != modElementalFormula.end () ) return i->second.getFormula ();
		else throw aaInfoNoElementalFormula ();
	}
	std::string getSubstituentA ( char aa ) const {	return substituentA.find ( aa )->second; }
	std::string getSubstituentB ( char aa ) const {	return substituentB.find ( aa )->second; }
	double getCTermPK ( char aa ) const { return pKCTerm.find ( aa )->second; }
	double getNTermPK ( char aa ) const { return pKNTerm.find ( aa )->second; }
	double getAcidicSCPK ( char aa ) const { return pKAcidicSC.find ( aa )->second; }
	double getBasicSCPK ( char aa ) const {	return pKBasicSC.find ( aa )->second; }
	StringVector aaListToElementalList ( const std::string& aaList );
	int getUserAAElementNumber ( const std::string& element, char code );
	bool isUsermod ( const std::string& aa ) const
	{
		MapNoCaseStringToDoubleConstIterator i = modMonoisotopicMass.find ( aa );
		return i != modMonoisotopicMass.end ();
	}
	std::string getLongName ( char aa ) const {	return longNames.find ( aa )->second; }
	char getAA ( const std::string& longName ) const
	{
		MapNoCaseStringToCharConstIterator iter = longNameToAA.find ( longName );
		if ( iter != longNameToAA.end () )
			return (*iter).second;
		else
			return 0;
	}
	CharVector getAminoAcids () const { return aminoAcids; }
	static std::string convert_pk ( double pk );
};

#endif /* ! __lu_aa_info_h */
