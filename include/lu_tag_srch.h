/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_tag_srch.h                                                 *
*                                                                             *
*  Created    : July 21st 1996                                                *
*                                                                             *
*  Purpose    : Interface code to the MS-Seq and MS-Tag search algorithms.    *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1996-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_tag_srch_h
#define __lu_tag_srch_h

#include <set>
#include <lu_hit.h>
#include <lu_db_srch.h>
#include <lu_mut_mtrx.h>
#include <lu_charge.h>
#include <lu_tag_frag.h>

class RegularExpression;
class FrameIterator;
class MSProductLink;

class TagHit : public ProteinHit {
	std::string sequence;
	char previousAA;
	char nextAA;
	int startAA;
	int unmatched;
	PeptideSequence hitSequence;
	PeakMatch peakMatch;
	ScoreType score;
	ScoreType scoreDiff;
	PeakMatch getPeakMatch ( Peak* parentPeak );
	static const ScoreType SCORE_DIFF_NOT_CALCULATED;
	static bool eValFlag;
public:
	int rank;
	TagHit ( const std::string& sequence, char previousAA, char nextAA, int startAA, const TagMatch& tagMatch );
	TagHit ( const FrameIterator& fi, const std::string& sequence, char previousAA, char nextAA, int startAA, const TagMatch& tagMatch );
	void printHTMLProteinHeader ( std::ostream& os ) const;
	void printHTMLHeader ( std::ostream& os, const std::string& parentToleranceUnits ) const;
	void printHTML ( std::ostream& os, const PeakMatchContext& peakMatchContext, const MSProductLink& productLink, double eValue ) const;
	void printHTMLSequence ( std::ostream& os, const MSProductLink& productLink ) const;
	void printDelimited ( std::ostream& os, const PeakMatchContext& peakMatchContext ) const;
	void printDelimitedContinue ( std::ostream& os ) const;
	std::string getSequence () const { return sequence; }
	char getPreviousAA () const { return previousAA; }
	char getNextAA () const { return nextAA; }
	int getStartAA () const { return startAA; }
	PeptideSequence getHitSequence () const { return hitSequence; }
	SetPairIntString getModIndicies () const { return hitSequence.getModIndicies (); }
	unsigned int getMassModAAMask () const { return hitSequence.getMassModAAMask (); }
	double getMassModRemainder () const { return peakMatch.getMatchedMass () - hitSequence.getMassMod (); }
	int getMmodIdx () const { return hitSequence.getMmodIdx (); }
	int getMmodProteinIdx () const { return startAA + hitSequence.getMmodIdx (); }
	void setMassMod ( double mm ) { hitSequence.setMassMod ( mm ); }
	int getNumUnmatched () const { return unmatched; }
	double getMassDiff () const { return peakMatch.getMassDiff (); }
	double getMatchedMass () const { return peakMatch.getMatchedMass (); }
	double getDataMOverZ () const { return peakMatch.getDataMOverZ (); }
	double getCharge () const { return peakMatch.getCharge (); }
	ScoreType getScore () const { return score; }
	void setScore ( ScoreType s ) { score = s; }
	void setScoreDiff ( ScoreType s ) { scoreDiff = s; }
	double getScoreOutput () const { return (double) score / (double) SCORE_TYPE_MULTIPLIER; }
	double getScoreDiffOutput () const { return (double) scoreDiff / (double) SCORE_TYPE_MULTIPLIER; }
	friend class sort_tag_hits;
	friend class sort_tag_hits_by_score;
	double getNLossMass () const { return hitSequence.getNLossMass (); }
	static bool setEValFlag ( bool flag ) { eValFlag = flag; }
	std::string getHitSequenceString () const { return hitSequence.getPeptideString (); }
	PairStringInt getAccessionNumberPair () const { return std::make_pair ( accessionNumber, startAA ); }
};

class sort_tag_hits {
public:
	int operator () ( const TagHit& a, const TagHit& b ) const
	{
		if ( a.unmatched == b.unmatched ) {
			return ( a.hitSequence.getSequence () < b.hitSequence.getSequence () );
		}
		return a.unmatched < b.unmatched;
	}
};

class sort_tag_hits_by_score {
public:
	int operator () ( const TagHit& a, const TagHit& b ) const
	{
		if ( a.score == b.score ) {
			int aMMod = a.hitSequence.getMmod ();
			int bMMod = b.hitSequence.getMmod ();
			if ( aMMod == bMMod )
				return fabs ( a.getMassDiff () ) < fabs ( b.getMassDiff () );
			else
				return aMMod < bMMod;
		}
		return a.score > b.score;
	}
};

class CombinedHitAccNumInfo {
	TagHit hit1;
	TagHit hit2;
public:
	CombinedHitAccNumInfo ( const TagHit& hit1, const TagHit& hit2 );
	bool getAccessionNumbersEqual () const { return hit1.getAccessionNumber () == hit2.getAccessionNumber (); }
	int getDecoyLevel () const
	{
		int decoyLevel = 0;
		if ( hit1.isDecoy () ) decoyLevel++;
		if ( hit2.isDecoy () ) decoyLevel++;
		return decoyLevel;
	}
	std::string getAccessionNumber () const { return hit1.getAccessionNumber () + hit2.getAccessionNumber (); }
	std::string getStartAA () const { return gen_itoa ( hit1.getStartAA () ) + gen_itoa ( hit2.getStartAA () ); }
	void printHTMLHit1 ( std::ostream& os ) const;
	void printHTMLHit2 ( std::ostream& os ) const;
	void printDelimitedHit1 ( std::ostream& os ) const;
	void printDelimitedHit2 ( std::ostream& os ) const;
	void printDelimitedHitContinue ( std::ostream& os ) const;
};

class SortCombinedTagAscending2 {	// string accNumber, int startAA
	public:
		bool operator () ( const CombinedHitAccNumInfo& lhs, const CombinedHitAccNumInfo& rhs ) const
		{
			if ( lhs.getDecoyLevel () == rhs.getDecoyLevel () ) {
				if ( lhs.getAccessionNumbersEqual () == rhs.getAccessionNumbersEqual () ) {
					if ( lhs.getAccessionNumber () == rhs.getAccessionNumber () ) {
						return lhs.getStartAA () < rhs.getStartAA ();
					}
					else
						return lhs.getAccessionNumber () < rhs.getAccessionNumber ();
				}
				else
					return lhs.getAccessionNumbersEqual () > rhs.getAccessionNumbersEqual ();
			}
			else
				return lhs.getDecoyLevel () < rhs.getDecoyLevel ();
		}
};

class CombinedTagHit {
	TagHit hit1;
	TagHit hit2;
	int numUnmatchedIons;
	ScoreType firstPeptideScore;
	ScoreType score;
	double error;
	std::set <CombinedHitAccNumInfo, SortCombinedTagAscending2> schani;
	void printHTMLSequence ( std::ostream& os, const MSProductLink& productLink, const LinkInfo* linkInfo ) const;
	const TagHit& initHit ( const TagHit& h1, const TagHit& h2 );
public:
	CombinedTagHit ( const TagHit& h1, int rank1, const TagHit& h2, int rank2, double linkMass, double error, const MSMSSearch* msMSSearch );
	void printHTMLHeader ( std::ostream& os, const std::string& parentToleranceUnits ) const;
	void printHTML ( std::ostream& os, int rank, const PeakMatchContext& pmc, const MSProductLink& productLink, const LinkInfo* linkInfo ) const;
	void printDelimited ( std::ostream& os, int rank ) const;
	ScoreType getFirstPeptideScore () const { return firstPeptideScore; }
	ScoreType getScore () const { return score; }
	ScoreType getScore1 () const { return hit1.getScore (); }
	ScoreType getScore2 () const { return hit2.getScore (); }
	bool getAccessionNumbersEqual () const { return hit1.getAccessionNumber () == hit2.getAccessionNumber (); }
	int getDecoyLevel () const
	{
		int decoyLevel = 0;
		if ( hit1.isDecoy () ) decoyLevel++;
		if ( hit2.isDecoy () ) decoyLevel++;
		return decoyLevel;
	}
	std::string getAccessionNumber () const { return hit1.getAccessionNumber () + hit2.getAccessionNumber (); }
	std::string getHitSequence () const { return hit1.getHitSequenceString () + hit2.getHitSequenceString (); }
	std::string getStartAA () const { return gen_itoa ( hit1.getStartAA () ) + gen_itoa ( hit2.getStartAA () ); }
	void add ( const TagHit& h1, const TagHit& h2 );
};

class SortCombinedTagAscending {	// This is to decide what gets stored as a separate entry
	public:
		bool operator () ( const CombinedTagHit* lhs, const CombinedTagHit* rhs ) const
		{
			if ( lhs->getScore () == rhs->getScore () ) {
				if ( lhs->getDecoyLevel () == rhs->getDecoyLevel () ) {
					if ( lhs->getAccessionNumbersEqual () == rhs->getAccessionNumbersEqual () ) {
						if ( lhs->getAccessionNumber () == rhs->getAccessionNumber () ) {
							if ( lhs->getHitSequence () == rhs->getHitSequence () ) {
								return lhs->getStartAA () < rhs->getStartAA ();
							}
							else
								return lhs->getHitSequence () < rhs->getHitSequence ();
						}
						else
							return lhs->getAccessionNumber () < rhs->getAccessionNumber ();
					}
					else
						return lhs->getAccessionNumbersEqual () > rhs->getAccessionNumbersEqual ();
				}
				else
					return lhs->getDecoyLevel () < rhs->getDecoyLevel ();
			}
			else
				return lhs->getScore () > rhs->getScore ();
		}
};
typedef std::set <CombinedTagHit*, SortCombinedTagAscending> CombinedTagHitSet;
typedef CombinedTagHitSet::const_iterator CombinedTagHitSetConstIterator;

class TagHitsContainer {
	typedef TagHit T;
	std::vector <T> tHits;
	int numHits;
	int numRetainedHits;
	int numSavedHits;
	const MSMSSearch* msMSSearch;
	const Peak* parentPeak;
	const PeakContainer* peaks;
	MSTagParameters& params;
	ScoreType minScore;
	CombinedTagHitSet combinedTagHits;
	void calculateCrossLinks ();
	void calculateModScores ();
	void printCrossLinksHTML ( std::ostream& os, const MSProductLink& productLink ) const;
	int getNumHitsToSave ( int numToSave );
public:
	TagHitsContainer ( const MSMSSearch* msMSSearch, MSTagParameters& params );
	T& operator [] ( int index ) { return tHits [index]; }
	void checkRegularExpression ( RegularExpression& regExp );
	void sortHits () { std::stable_sort ( tHits.begin (), tHits.end (), sort_tag_hits_by_score () ); }
	void rankHits ();
	void pruneHits ();
	void push_back ( const T& hit );
	int size () const { return tHits.size (); }
	int getNumPeaks () const { return peaks->size (); }
	int getNumSavedHits () const { return numSavedHits; }
	const Peak* getParentPeak ( int dataSet = 0 ) const { return parentPeak; }
	const PeakContainer* getPeaks ( int dataSet = 0 ) const { return peaks; }
	void printHTMLReport ( std::ostream& os, int dataSet );
	void printDelimited ( std::ostream& os );
	void printHTML ( std::ostream& os, const MSProductLink& productLink ) const;

	int getIndex ( int num ) const
		{ return tHits [num].getIndex (); }
	int getDBIndex ( int num ) const
		{ return tHits [num].getDBIndex (); }
	std::string getAccessionNumber ( int num ) const
		{ return tHits [num].getAccessionNumber (); }
	int getOpenReadingFrame ( int num ) const
		{ return tHits [num].getOpenReadingFrame (); }
	int getDNAReadingFrame ( int num ) const
		{ return tHits [num].getDNAReadingFrame (); }
	ScoreType getScore ( int i ) const { return tHits [i].getScore (); }
	ScoreType getMinScore () const { return minScore; }
};

class TagHits : public DatabaseHits {
	std::vector <TagHitsContainer*> hits;
	int numSearches;
public:
	TagHits ( const std::vector <MSMSSearch*>& msMSSearch, MSTagParameters& params );
	~TagHits ();

	TagHit getTagHit ( int num, int dataSet ) const
		{ return (*hits [dataSet]) [num]; }
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
	int getNumPeaks ( int dataSet = 0 ) const { return hits [dataSet]->getNumPeaks (); }
	int getNumSavedHits ( int dataSet = 0 ) const { return hits [dataSet]->getNumSavedHits (); }
	const Peak* getParentPeak ( int dataSet = 0 ) const { return hits [dataSet]->getParentPeak (); }
	const PeakContainer* getPeaks ( int dataSet = 0 ) const { return hits [dataSet]->getPeaks (); }
	void printHTMLHeader ( std::ostream& os, int num, int dataSet = 0 ) const {}
	void printHTMLHit ( std::ostream& os, int num, int dataSet = 0 ) const {}
	void printXMLHit ( std::ostream& os, int num, int dataSet = 0 ) const {}

	void checkRegularExpression ( RegularExpression& regExp );
	void sortAndRank ();
	void prune ();
	void addHit ( const TagHit& hit, int dataSet = 0 );

	ScoreType getScore ( int i, int j ) const { return hits [i]->getScore ( j ); }
	ScoreType getMinScore ( int index ) const { return hits [index]->getMinScore (); }
	void printHTMLReport ( std::ostream& os, int dataSet = 0 )
		{ hits [dataSet]->printHTMLReport ( os, dataSet ); }
	void printDelimited ( std::ostream& os, int dataSet = 0 )
		{ hits [dataSet]->printDelimited ( os ); }
};
typedef std::vector <TagHit> TagHitVector;
typedef TagHitVector::size_type TagHitVectorSizeType;

class TagSearch : public DatabaseSearch {
protected:
	MSTagParameters& tagParams;
	std::vector <MSMSSearch*> msMSSearch;
	int numSearches;
	int maxTagMatches;
	bool randomSearch;
	TagHits* tagHits;
	TagMatchVector tagMatch;
	unsigned int compMask;
	bool compMaskTypeAnd;
	bool compMaskTypeOr;
	void addTagHit ( int searchNumber, const TagMatchVector& tagMatch, const std::string& sequence, const FrameIterator& fi, int previousAA, int nextAA, int startAA );
	virtual void tag_search ( const FrameIterator& fi, const char* frame ) = 0;
	void printParamsBodyHTML ( std::ostream& os ) const;
	bool checkComposition ( const std::string& fragment );
	static int getPruneInterval ( int numSearches )
	{
		if ( numSearches < 1000 ) return 1000;
		else if ( numSearches < 5000 ) return 500;
		else if ( numSearches < 10000 ) return 100;
		else return 50;
	}
	static int getActualPruneInterval ( int index, int pruneInterval )
	{
		if ( index < 4 ) return 1;
		else if ( index < 20 ) return 4;
		else if ( index < 100 ) return 20;
		else if ( index < 1000 ) return genMin ( pruneInterval, 100 );
		else return pruneInterval;
	}
	bool continueSearch () const;
	bool ( *nonSpecificTermini ) ( const IntVector& cleavageIndex, int start, int end );
	bool checkD;
	static bool nonSpecificCTermini ( const IntVector& cleavageIndex, int start, int end );
	static bool nonSpecificNTermini ( const IntVector& cleavageIndex, int start, int end );
	static bool nonSpecific1Termini ( const IntVector& cleavageIndex, int start, int end );
	static bool nonSpecific2Termini ( const IntVector& cleavageIndex, int start, int end );
	static bool noEnzyme ( const IntVector& cleavageIndex, int start, int end ) { return true; }
	static char allowNonSpecificType;
	static int missedCleavages;
	void initNonSpecific ( const MSTagParameters& params );
	void doNormalSearch ( FastaServer* fsPtr, int num );
public:
	TagSearch ( MSTagParameters& params );
	virtual ~TagSearch ();
	void doSearch ();
	void printHTMLHitsReport ( std::ostream& os );
	void printSingleProteinResultsHTML ( std::ostream& os, int index ) const;
	void printSingleProteinResultsXML ( std::ostream& os, int index ) const;
	void printXMLHits ( std::ostream& os ) const;
	void printBodyXML ( std::ostream& os, bool showPreSearch );
};

TagSearch* getTagSearch ( MSTagParameters& tagParams );

class NoEnzymeSearch : public TagSearch {
	DoubleVector startMasses;
	DoubleVector endMasses;
protected:
	void tag_search ( const FrameIterator& fi, const char* frame );
public:
	NoEnzymeSearch ( MSTagParameters& params );
	~NoEnzymeSearch ();
};

class NoEnzymeIntSearch : public TagSearch {
	IntVector startMasses;
	IntVector endMasses;
protected:
	void tag_search ( const FrameIterator& fi, const char* frame );
public:
	NoEnzymeIntSearch ( MSTagParameters& params );
	~NoEnzymeIntSearch ();
};

class YesEnzymeSearch : public TagSearch {
	double startLimit;
	double endLimit;
	double cleavedLimit;
protected:
	void tag_search ( const FrameIterator& fi, const char* frame );
public:
	YesEnzymeSearch ( MSTagParameters& params );
	~YesEnzymeSearch ();
};

#endif /* ! __lu_tag_srch_h */
