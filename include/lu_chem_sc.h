/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_chem_sc.h                                                  *
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
#ifndef __lu_chem_sc_h
#define __lu_chem_sc_h

#ifdef CHEM_SCORE
#include <string>
#include <ostream>
#include <vector>

class ParameterList;
class RegularExpression;

class ChemScoreParameters {
	bool chemScoreFlag;
	std::string cysMod;
	double metOxFactor;
public:
	ChemScoreParameters ( const ParameterList* params );
	void printHTML ( std::ostream& os ) const;
	static void copyToCGI ( std::ostream& os, const ParameterList* params );
	static void copyToHiddenFormEntry ( std::ostream& os, const ParameterList* params );
};

class ScorePair {
	RegularExpression* regExp;
	int alignment;		// Offset from cleavage site.
	double score;
public:
	ScorePair ( const std::string& exp, int alignment, double score );
	~ScorePair () {}
	RegularExpression* getRegExp () const { return regExp; }
	double getAlignment () const { return alignment; }
	double getScore () const { return score; }
};

typedef std::vector <ScorePair> ScorePairVector;
typedef ScorePairVector::size_type ScorePairVectorSizeType;

class ChemScore {
	ScorePairVector initialSettingsSP;
	double detrimentalCys;
	RegularExpression* locator;
	ScorePairVector missedCleavagesSP;
	ScorePairVector terminalAdjustmentSP;
	double metOxF;
public:
	ChemScore ( const std::string& cysMod, double metOxF );
	double getChemScore ( const std::string& fragment, int numMetOx ) const;
};
#endif

#endif /* ! __lu_chem_sc_h */
