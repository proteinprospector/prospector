/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_spep_srch.h                                                *
*                                                                             *
*  Created    : February 12th 2002                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2002-2009) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_spep_srch_h
#define __lu_spep_srch_h

#include <vector>
#include <lu_charge.h>
#include <lu_mut_mtrx.h>
#include <lu_links.h>

class ModificationHit {
	EnzymeFragment enzFrag;
	PeakMatch peakMatch;
	PeptideSequence mutatedSequence;
public:
	ModificationHit ( const EnzymeFragment& enzFrag, const PeakMatch& peakMatch, const PeptideSequence& mutatedSequence );
	bool operator!= ( const ModificationHit& rhs ) const;
	const EnzymeFragment& getEnzymeFragment () const { return enzFrag; }
	bool isHomologyMatch ( const ModificationHit& rhs ) const;
	void printHTMLHeader ( std::ostream& os, const PeakMatchContext& peakMatchContext, bool reportModifications );
	void printHTML ( std::ostream& os, const PeakMatchContext& peakMatchContext, const MSProductLink& productLink, bool reportModifications, bool hideLinks ) const;
	void printXML ( std::ostream& os, const PeakMatchContext& peakMatchContext ) const;
	bool getMutated () const { return mutatedSequence.getLength () != 0; }
	friend class sortModificationHits;
	friend class sortNSpecHits;
};

typedef std::vector <ModificationHit> ModificationHitVector;
typedef ModificationHitVector::size_type ModificationHitVectorSizeType;

class SinglePeptideSearch : public SingleFitSearch {
	mutable bool reportModifications;
	bool reportModificationsSet;
protected:
	std::vector <ModificationHit> hits;
	bool areHits () const;
public:
	SinglePeptideSearch ( const PeakContainer& peaks ) :
		SingleFitSearch ( peaks ),
		reportModificationsSet ( false ) {};
	virtual ~SinglePeptideSearch ();

	void calculateAACovered ( int index, int protLen );
	void calculateNumUnique ();
	void printHTMLBody ( std::ostream& os, bool hideLinks );
	void printXMLBody ( std::ostream& os );
	bool getReportModifications () const;
};

class sortModificationHits {
public:
	int operator () ( const ModificationHit& a, const ModificationHit& b ) const
	{
		if ( a.peakMatch.getMatchedMass () == b.peakMatch.getMatchedMass () ) {
			return a.mutatedSequence.getSequence () < b.mutatedSequence.getSequence ();
		}
		return ( a.peakMatch.getMatchedMass () < b.peakMatch.getMatchedMass () );
	}
};

class sortNSpecHits {
public:
	int operator () ( const ModificationHit& a, const ModificationHit& b ) const
	{
		if ( a.peakMatch.getMatchedMass () == b.peakMatch.getMatchedMass () ) {
			return a.enzFrag.getFragment () < b.enzFrag.getFragment ();
		}
		return ( a.peakMatch.getMatchedMass () < b.peakMatch.getMatchedMass () );
	}
};

class ModificationSearch : public SinglePeptideSearch {
public:
	ModificationSearch ( const EnzymeFragmentContainer& enzFrags, const ModificationTable& modificationTable, const PeakContainer& peaks );
	virtual ~ModificationSearch ();
	int calculateNumHomologyMatches ( const ModificationSearch* ls ) const;
};

class NonSpecificSearch : public SinglePeptideSearch {
public:
	NonSpecificSearch ( const char* frame, int frameLen, const PeakContainer& peaks, int maxNonSpecificHits );
	virtual ~NonSpecificSearch ();
};

#endif /* ! __lu_spep_srch_h */
