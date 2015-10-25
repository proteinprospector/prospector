/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_frag_info.h                                                *
*                                                                             *
*  Created    : June 20th 1996                                                *
*                                                                             *
*  Purpose    : Function to calculate information on enzyme fragments.        *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1996-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_frag_info_h
#define __lu_frag_info_h

#include <vector>
#include <string>
#include <lg_string.h>
#include <lu_formula.h>

class MSProductLink;
class AACalculator;
class EnzymeParameters;

class EnzymeFragment {
	std::string fragment;
	int startAA;
	int endAA;
	char previousAA;
	char nextAA;
	int missedCleavages;
	int maxCharge;
public:
	EnzymeFragment () {};
	EnzymeFragment ( const std::string& protein, int start, int end, int missedCleavages );
	std::string getFragment () const { return fragment; }
	bool getFirstFragment () const { return startAA == 1; }
	bool getLastFragment () const { return nextAA == '-'; }
	int getStartAA () const { return startAA; }
	int getEndAA () const { return endAA; }
	char getPreviousAA () const { return previousAA; }
	char getNextAA () const { return nextAA; }
	int getMissedCleavages () const { return missedCleavages; }
	int getMaxCharge () const { return maxCharge; }
	double getMass () const;
	int getNumPossMod ( char aa, bool nTermLinkPossible ) const;
	ElementalFormula getElementalFormula ( const AACalculator& aaCalc ) const;
	int getNumPossCys () const { return gen_strcharcount ( fragment, 'C' ); }
	bool containsAA ( char aa ) const { return gen_strcontains ( fragment, aa ); }
	bool containsAAMaskAnd ( unsigned int mask ) const;
	bool containsAAMaskOr ( unsigned int mask ) const;
	static void printHeaderHTML ( std::ostream& os );
	static void printHeaderDelimited ( std::ostream& os );
	void printHTML ( std::ostream& os, bool hideLinks, const MSProductLink& productLink, int maxCharge = 0, char aa = 0 ) const;
	void printDelimited ( std::ostream& os ) const;
	void printXML ( std::ostream& os ) const;
};

class EnzymeFragmentContainer {
	std::vector <EnzymeFragment> enzymeFragmentList;
	int protLen;
public:
	typedef std::vector <EnzymeFragment>::size_type size_type;
	EnzymeFragmentContainer ( const std::string& protein, const EnzymeParameters& enzymeParameters );
	int getProtLen () const { return protLen; }
	friend class EnzymeFragmentIterator;
};

typedef std::vector <EnzymeFragmentContainer> EnzymeFragmentContainerVector;
typedef EnzymeFragmentContainerVector::size_type EnzymeFragmentContainerVectorSizeType;
typedef EnzymeFragmentContainer::size_type EnzymeFragmentContainerSizeType;

class EnzymeFragmentIterator {
	const EnzymeFragmentContainer& info;
	EnzymeFragmentContainerSizeType cur;
public:
	EnzymeFragmentIterator ( const EnzymeFragmentContainer& info ) :
		info ( info ), cur ( 0 ) {}
	EnzymeFragment getEnzFrag () const { return info.enzymeFragmentList [cur]; }
	bool more () const { return cur < info.enzymeFragmentList.size (); }
	void advance () { cur++; }
};

#endif /* ! __lu_frag_info_h */
