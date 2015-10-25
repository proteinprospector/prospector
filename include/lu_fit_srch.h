/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_fit_srch.h                                                 *
*                                                                             *
*  Created    : July 13th 2001                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_fit_srch_h
#define __lu_fit_srch_h

#include <vector>

#include <lu_hit.h>
#include <lu_fit_par.h>
#include <lu_db_srch.h>
#include <lu_bdg_par.h>
#include <lu_nspec_par.h>

class MSFitSearch;
class SingleFitSearch;

struct HitStats {
	double mowseScore;
	CharVector massMatched;
	CharVector mowseMissedCleavages;
	HitStats () : mowseScore ( 0.0 ) {}
};

class FitHit : public ProteinHit {
	SingleFitSearch* sfs;
	int numMatches;
	int numUnique;
	int numHomology;
	int hitNumber;
public:
	HitStats hitStats;
	FitHit ( FastaServer* fs, int ind, int drf, int orf, int numMatches, const HitStats& hitStats );
	SingleFitSearch* getSfs () const { return sfs; }
	void setSingleFitSearch ( const MSFitParameters& params, const PeakContainer& peaks );
	void setNumHomology ( int n ) { numHomology = n; }
	void setHitNumber ( int h ) { hitNumber = h; }
	int getNumUnique () const { return numUnique; }
	int getHitNumber () const { return hitNumber; }
	double getPercentMatches ( int numPeaks ) const { return 100.0 * numMatches / numPeaks; }
	void printHTMLHeader ( std::ostream& os, int numPeaks, bool scoreFlag, bool homFlag ) const;
	void printHTML ( std::ostream& os, int numPeaks, bool scoreFlag, bool homFlag ) const;
	void printCoverageHTML ( std::ostream& os ) const;
	void printXML ( std::ostream& os, int numPeaks, bool scoreFlag ) const;
	void printHTMLDetailSummary ( std::ostream& os, int i, int numPeaks ) const;
	void printHTMLDetail ( std::ostream& os, const MSFitParameters& params, const PeakContainer& peaks ) const;
	void printXMLDetail ( std::ostream& os, const MSFitParameters& params, const PeakContainer& peaks ) const;
	friend class sort_hits_by_pi;
	friend class sort_hits_by_protein_mw;
	friend class sort_fit_hits;
	friend class sort_mowse_hits;
};
typedef std::vector <FitHit> FitHitVector;
typedef FitHitVector::size_type FitHitVectorSizeType;

typedef std::vector <FitHit> FitHitVector;
typedef FitHitVector::size_type FitHitVectorSizeType;

class FitHitsContainer {
	typedef FitHit T;
	int searchNumber;
	std::string spot;
	std::vector <FitHit> hits;
	VectorPairIntIntVector mainHits;
	int numHits;
	const MSFitSearch* fitSearch;
	PeakContainer peaks;
	MSFitParameters& params;
	void calculateMainAndSupplementaryHits ( int (*numHomologyMassesCalculator) ( const FitHit&, const FitHit& ) );
	static int calculateNumHomologyMatchesMS ( const FitHit& hit, const FitHit& homHit );
	static int calculateNumHomologyMatchesLS ( const FitHit& hit, const FitHit& homHit );
public:
	FitHitsContainer ( const MSFitSearch* fitSearch, MSFitParameters& params, const std::string& spot, int searchNumber );
	int size () const { return hits.size (); }
	T& operator [] ( int index ) { return hits [index]; }
	const T& operator [] ( int index ) const { return hits [index]; }
	void sortHits ();
	void eraseHits ();
	void sortHitsByPIorMW ();
	void setSingleFitSearch ();
	void calculateMainAndSupplementaryHits ();
	void push_back ( const T& hit ) { hits.push_back ( hit ); }
	void printHTMLReport ( std::ostream& os );
	void printXMLReport ( std::ostream& os );
	void printHTMLSummary ( std::ostream& os );
	void printCoverageTable ( std::ostream& os, const VectorPairIntIntVector& hitList ) const;
	void printHTMLDetail ( std::ostream& os );
	void printXMLDetail ( std::ostream& os );
	int getIndex ( int num ) const
		{ return hits [num].getIndex (); }
	int getDBIndex ( int num ) const
		{ return hits [num].getDBIndex (); }
	std::string getAccessionNumber ( int num ) const
		{ return hits [num].getAccessionNumber (); }
	int getOpenReadingFrame ( int num ) const
		{ return hits [num].getOpenReadingFrame (); }
	int getDNAReadingFrame ( int num ) const
		{ return hits [num].getDNAReadingFrame (); }
};

class FitHits : public DatabaseHits {
	std::vector <FitHitsContainer*> hits;
	int numSearches;
public:
	FitHits ( const std::vector <MSFitSearch*>& msFitSearch, MSFitParameters& params );

	int getIndex ( int num, int dataSet = 0 ) const
		{ return hits [dataSet]->getIndex ( num ); }
	int getDBIndex ( int num, int dataSet = 0 ) const
		{ return hits [dataSet]->getDBIndex ( num ); }
	std::string getAccessionNumber ( int num, int dataSet = 0 ) const
		{ return hits [dataSet]->getAccessionNumber ( num ); }
	int getOpenReadingFrame ( int num, int dataSet = 0 ) const
		{ return hits [dataSet]->getOpenReadingFrame ( num ); }
	int getDNAReadingFrame ( int num, int dataSet = 0 ) const
		{ return hits [dataSet]->getDNAReadingFrame ( num ); }

	int size ( int dataSet = 0 ) const { return hits [dataSet]->size (); }
	void printHTMLHeader ( std::ostream& os, int num, int dataSet = 0 ) const {}
	void printHTMLHit ( std::ostream& os, int num, int dataSet = 0 ) const {}
	void printXMLHit ( std::ostream& os, int num, int dataSet = 0 ) const {}

	void sortAndRank ();
	void addHit ( const FitHit& hit, int dataSet = 0 ) { hits [dataSet]->push_back ( hit ); }

	void printHTMLReport ( std::ostream& os, int dataSet = 0 )
		{ hits [dataSet]->printHTMLReport ( os );	}
	void printXMLReport ( std::ostream& os, int dataSet = 0 )
		{ hits [dataSet]->printXMLReport ( os ); }

	FitHitsContainer* getHit ( int index ) const { return hits [index]; }
};

class FrameIterator;

class FitSearch : public DatabaseSearch {
protected:
	const MSFitParameters& fitParams;
	std::vector <MSFitSearch*> fitSearch;
	int numSearches;
	int maxHits;
	FitHits* fitHits;
	void printParamsBodyHTML ( std::ostream& os ) const { fitParams.printHTML ( os ); }
	void doSearch ();
	void doSearch ( FastaServer* fsPtr, int num );
	void calculateMowseScores ();
public:
	FitSearch ( MSFitParameters& params );
	~FitSearch ();
	void printHTMLHitsReport ( std::ostream& os );
	void printXMLHits ( std::ostream& os ) const;
	void printHTMLHitsJavascript ( std::ostream& os ) const;
};

#endif /* ! __lu_fit_srch_h */
