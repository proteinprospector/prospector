/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_comp_srch.h                                                *
*                                                                             *
*  Created    : October 22nd 2001                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_comp_srch_h
#define __lu_comp_srch_h

#include <lu_program.h>
#include <lu_formula.h>

class Peak;
class MSCompParameters;
class CompositionSearch;
struct Combination_hit;
class AACalculator;
class MSIsotopeLink;

class CompHit {
	std::string sequence;
	AAFormula aaComp;
	double numPermutations;
public:
	double pepMW;
	int index;
	ElementalFormula elemComp;
	double doubleBondEquivalent;

	CompHit ( double mw, const std::string& sequence );
	CompHit ( double mw, ElementalFormula& elemComp );
	CompHit ( ElementalFormula& elemComp, double mw );
	double getNumPermutations () const { return numPermutations; }
	std::string getSequence () const { return sequence; }
	static void printHTMLHeader ( std::ostream& os, const std::string& combinationType, const MSCompParameters& params, const std::string& ionType );
	void printHTML ( std::ostream& os, const std::string& combinationType, const MSIsotopeLink& isotopeLink, const AAFormula& knownAAComposition, const Peak* parentPeak, double offset, const MSCompParameters& params );

	friend class sort_comp_hits;
};

class sort_comp_hits {
public:
	int operator () ( const CompHit& a, const CompHit& b ) const
	{
		return ( a.pepMW < b.pepMW );
	}
};
typedef std::vector <CompHit> HitsContainer;
typedef HitsContainer::size_type HitsContainerSizeType;

class MSCompSearch : public MSProgram {
	const MSCompParameters& compParams;
	ElementalFormula cTermFormula;
	CompositionSearch* compositionSearch;
	AACalculator* aaCalc;
	HitsContainer getHitsInformationAminoAcid ( Combination_hit* combinations, int numCombinations, const std::string& ionType );
	HitsContainer getHitsInformationPeptideElemental ( Combination_hit* combinations, int numCombinations, const std::string& ionType );
	HitsContainer getHitsInformationElemental ( char** elementalCombinations, int numCombinations, const std::string& ionType );
	int getNumUniqElemCompositions ( const HitsContainer& hits );
	void outputHTMLHeader ( std::ostream& os, Peak* parentPeak, double effectiveParentMass );
	void outputHTMLTable ( std::ostream& os, const MSCompParameters& params, HitsContainer& hits, Peak* parentPeak, int maxHitsExceeded, const std::string& ionType, const MSIsotopeLink& isotopeLink );
	double ionTypeToOffset ( const std::string& ionType );
	ElementalFormula adjustedElementalFormula ( const ElementalFormula& elemental_formula, const std::string& ionType );
public:
	MSCompSearch ( const MSCompParameters& params );
	void printBodyHTML ( std::ostream& os );
	void printBodyXML ( std::ostream& os );
	void printParamsBodyHTML ( std::ostream& os ) const;
};

#endif /* ! __lu_comp_srch_h */
