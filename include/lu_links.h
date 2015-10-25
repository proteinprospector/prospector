/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_links.h                                                    *
*                                                                             *
*  Created    : February 2nd 2000                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2000-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_links_h
#define __lu_links_h

#include <vector>
#include <lu_usermod.h>
#include <lu_frag_info.h>
#include <lu_sing_fit.h>

#ifdef CHEM_SCORE
class ChemScore;
#endif
class Peak;
class LinkInfo;

class PotentialLinkFragment {
	static StringVector accessionNumbers;
	EnzymeFragment enzymeFragment;
	double mass;
	MapIntToInt modList;
	std::pair <int, IntVector> numLinkAA;
	int entryNumber;
#ifdef CHEM_SCORE
	static ChemScore* chemScore;
#endif
public:
	PotentialLinkFragment ( const EnzymeFragment& enzymeFragment, double mass, const MapIntToInt& modList, const std::pair <int, IntVector>& numLinkAA, int entryNumber = 0 );
	double getMass () const { return mass; }
	int getNumLinkAA () const { return numLinkAA.first; }
	int getNumLinkAA ( int idx ) const { return numLinkAA.second [idx]; }
	int getStartAA () const { return enzymeFragment.getStartAA (); }
	int getEndAA () const { return enzymeFragment.getEndAA (); }
	std::string getFragment () const { return enzymeFragment.getFragment (); }
	int getEntryNumber () const { return entryNumber; }
	ElementalFormula getElementalFormula ( const AACalculator& aaCalc ) const;
	static void printHTMLHeader ( std::ostream& os, bool modifications, bool multiMolecules );
	void printHTML ( std::ostream& os, bool modifications, bool multiMolecules, bool multipleFragments, const MSProductLink& productLink, int maxCharge, bool hideLinks ) const;
	void printXML ( std::ostream& os, bool modifications ) const;
	static void setAccessionNumbers ( const StringVector& sv ) { accessionNumbers = sv; }
	static int getEntryNumber ( const std::string& acc );
	static int getNumAccessionNumbers () { return accessionNumbers.size (); }
	friend class sort_fragments_by_mass;
	bool containsMods ( const MapStringToInt& mods ) const;
	bool isExactMod ( const MapStringToInt& mods ) const;
#ifdef CHEM_SCORE
	static void setChemScore ( const std::string& cysMod, double metOxF );
	static void deleteChemScore ();
#endif
};

typedef std::vector <PotentialLinkFragment> PotentialLinkFragmentVector;
typedef PotentialLinkFragmentVector::size_type PotentialLinkFragmentVectorSizeType;
typedef PotentialLinkFragmentVector::const_iterator PotentialLinkFragmentVectorConstIterator;
typedef std::vector <PotentialLinkFragmentVector> PotentialLinkFragmentVectorVector;
typedef PotentialLinkFragmentVectorVector::size_type PotentialLinkFragmentVectorVectorSizeType;

class sort_fragments_by_mass {
public:
	int operator () ( const PotentialLinkFragment& a, const PotentialLinkFragment& b ) const
	{
		return ( a.mass < b.mass );
	}
};

class LinkHits {
	PotentialLinkFragmentVector hitFragments;
	PeakMatch peakMatch;
	std::string moleculeType;
	ElementalFormula elementalFormula;
public:
	LinkHits ( const PotentialLinkFragmentVector& hitFragments, const Peak* peak, double mass, const std::string& moleculeType, const AACalculator& aaCalc, ElementalFormula& startElementalFormula );
	bool operator!= ( const LinkHits& rhs ) const;
	bool isHomologyMatch ( const LinkHits& rhs ) const;
	double getError ( const PeakMatchContext& peakMatchContext ) const;
	bool getMultipleFragments () const { return hitFragments.size () > 1; }
	void printHTMLHeader ( std::ostream& os, const PeakMatchContext& peakMatchContext, bool modifications, bool multiMolecules, bool peptideCombination );
	const PotentialLinkFragmentVector& getHitFragments () const { return hitFragments; }
	void printHTML ( std::ostream& os, const PeakMatchContext& peakMatchContext, bool modifications, bool multiMolecules, bool peptideCombination, const MSProductLink& productLink, const MSIsotopeLink& isotopeLink, bool hideLinks );
	void printXML ( std::ostream& os, const PeakMatchContext& peakMatchContext, bool modifications, bool peptideCombination );
	friend class sortLinkHits;
};
typedef std::vector <LinkHits> LinkHitsVector;
typedef LinkHitsVector::size_type LinkHitsVectorSizeType;

class sortLinkHits {
public:
	int operator () ( const LinkHits& a, const LinkHits& b ) const
	{
		if ( a.peakMatch.getMatchedMass () == b.peakMatch.getMatchedMass () ) {
			if ( a.hitFragments.size () == 1 && b.hitFragments.size () == 1 ) {
				return a.hitFragments [0].getFragment () < b.hitFragments [0].getFragment ();
			}
		}
		return ( a.peakMatch.getMatchedMass () < b.peakMatch.getMatchedMass () );
	}
};

class LinksSearch : public SingleFitSearch {
	PotentialLinkFragmentVector potentialLinkFrags;
	double ( *mass_convert ) (const char*);
	bool multiMolecules;
	bool modifications;
	bool peptideCombinationFlag;

	void initLinksSearch ( const std::vector <EnzymeFragmentContainer>& enzFragInfo, const PeakContainer& peaks, const LinkInfo* linkInfo, const AACalculator& aaCalc, const StringVector knownSequences );
	void calculatePotentialFragments ( const EnzymeFragmentContainer& enzFrag, int entryNumber, double maxMass );
protected:
	LinkHitsVector linkHits;
	bool areHits () const;
public:
	LinksSearch ( const std::vector <EnzymeFragmentContainer>& enzFragInfo, const PeakContainer& parentPeaks, const LinkInfo* linkInfo, const AACalculator& aaCalc, const StringVector knownSequences = StringVector () );
	LinksSearch ( const EnzymeFragmentContainer& enzFragInfo, const PeakContainer& parentPeaks, const LinkInfo* linkInfo, const AACalculator& aaCalc, const StringVector knownSequences = StringVector () );
	virtual ~LinksSearch ();
	void calculateAACovered ( int index, int protLen );
	void calculateNumUnique ();
	BoolDeque getPeakUsed () const { return peakUsed; }
	int calculateNumHomologyMatches ( const LinksSearch* sfs ) const;
	void printHTMLBody ( std::ostream& os, bool hideLinks );
	void printXMLBody ( std::ostream& os );
};

#endif /* ! __lu_links_h */
