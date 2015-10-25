/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_const_mod.h                                                *
*                                                                             *
*  Created    : October 21st 2002                                             *
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

#ifndef __lu_const_mod_h
#define __lu_const_mod_h

#include <string>
#include <vector>
#include <lgen_define.h>
#include <lu_formula.h>

struct SingleMod;

typedef std::vector <SingleMod*> SingleModPtrVector;
typedef SingleModPtrVector::size_type SingleModPtrVectorSizeType;

class ConstMod {
	static MapStringToInt indexMap;
	static SingleModPtrVector singleMods;
	static SingleModPtrVector extraSCompSingleMods;
	static bool productFlag;
	static void addExtraUsermods ( SetString& nameSet, SetString& nameSet2 );
	static void readUsermodFile ( const std::string& fName, SetString& nameSet, SetString& nameSet2 );
	static void processSingleUsermod ( const std::string& longName, const std::string& formulaStr, const std::string& line3, SetString& nameSet, SetString& nameSet2 );
	static void initialiseMods ();
	static void initialiseTerm ( const std::string& filename, const std::string& menuName, const std::string& code );
	static bool initialised;
	std::string aaList;
	std::string longName;
	std::string formulaStr;
	ElementalFormula formula;
	mutable double mass;
public:
	ConstMod ( const std::string& n );
	std::string getAAList () const { return aaList; }
	std::string getLongName () const { return longName; }
	ElementalFormula getElementalFormula () const { return formula; }
	bool getValidFormula () const { return formulaStr [0] != '*'; }
	static StringVector getNames ();
	static StringVector getNames ( const std::string& type );
	static StringVector getNTermNames () { return getNames ( "n" ); }
	static StringVector getCTermNames () { return getNames ( "c" ); }
	static StringVector getNLossNames ();
	static StringVector getNonTerminalNames ();
	static bool setProductFlag () { productFlag = true; }
	void setMass ();
	double getMass () const { return mass; }
};

typedef std::map <char, ConstMod*> MapCharConstModPtr;
typedef MapCharConstModPtr::const_iterator MapCharConstModPtrConstIterator;

typedef std::map <std::string, ConstMod*> MapStringConstModPtr;
typedef MapStringConstModPtr::const_iterator MapStringConstModPtrConstIterator;

#endif /* ! __lu_const_mod_h */
