/******************************************************************************
*                                                                             *
*  Program    : searchCompare                                                 *
*                                                                             *
*  Filename   : sc_xlink.h                                                    *
*                                                                             *
*  Created    : July 9th 2012                                                 *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2012-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __sc_xlink_h
#define __sc_xlink_h

#include <sc_search_res.h>

class MSProductLink;

class XLinkDecoyInfo {
	StringVector accNum;
	IntVector startAA;
	StringVector accNum2;
	IntVector startAA2;
	int pep1IntraCount;
	int pep1InterCount;
	int pep1DecoyCount;
	int pep2IntraCount;
	int pep2InterCount;
	int pep2DecoyCount;
	static std::map <std::string, ProteinInfo*> proteinInfoMap;
public:
	XLinkDecoyInfo () {}
	void add ( const std::string& a, int s, const std::string& a2, int s2 );
	void calc ( int num );
	std::string getPreviousAA1 ( int num, int reportPreviousAA ) const;
	std::string getPreviousAA2 ( int num, int reportPreviousAA ) const;
	std::string getNextAA1 ( int num, int endAA, int reportNextAA ) const;
	std::string getNextAA2 ( int num, int endAA, int reportNextAA ) const;
	int getStartAA1 ( int num ) const;
	int getStartAA2 ( int num ) const;
	bool isDecoy1 ( int num ) const;
	bool isDecoy2 ( int num ) const;
	std::string getFullAccessionNumber1 ( int num ) const;
	std::string getFullAccessionNumber2 ( int num ) const;
	void printProteinInfoHTML ( std::ostream& os, const SResLink& sresLink, int num ) const;
	void printHTMLANum1 ( std::ostream& os, int num ) const;
	void printHTMLANum2 ( std::ostream& os, int num ) const;
	void printDelimitedANum1 ( std::ostream& os, int num ) const;
	void printDelimitedANum2 ( std::ostream& os, int num ) const;
	void printDelimited1 ( std::ostream& os, int num ) const;
	void printDelimited2 ( std::ostream& os, int num ) const;
	std::string getRepeats1 ( int num ) const;
	std::string getRepeats2 ( int num ) const;
	void printRepeatsDelimited1 ( std::ostream& os, int num ) const;
	void printRepeatsDelimited2 ( std::ostream& os, int num ) const;
	static void printRepeatsDelimitedHeader1 ( std::ostream& os );
	static void printRepeatsDelimitedHeader2 ( std::ostream& os );
};

class SearchResultsCrosslinkPeptideHit {
	PPPeptideHitInfo peptideHitInfo;
	const SpecID* specID;
	const MSMSSpectrumInfo* mmsi;
	double error;
	const HitPeptide* hitPeptide1;
	const HitPeptide* hitPeptide2;
	XLinkDecoyInfo xldi;
	short searchIndex;
	double xScore1;
	int xRank1;
	double xScore2;
	int xRank2;
	double xFirstScore;
	mutable PeakFitData* quanRatio;

	static StringVector instrument;
	static std::vector <Tolerance*> parentTolerances;
	static std::vector <Tolerance*> fragmentTolerances;
	static BoolDeque spottingPlates;
	static BoolDeque spectrumNumber;
	static StringVectorVector fractionNames;
	static bool multipleFractionNames;
	static StringVectorVector rawTypes;
	static StringVector searchKey;
	static bool reportMPlusH;
	static bool reportMOverZ;
	static bool reportCharge;
	static bool reportIntensity;
	static bool reportMSMSInfo;
	static bool reportStartAA;
	static bool reportEndAA;
	static bool reportRepeats;
	static bool reportLinks;
	static bool reportError;
	static bool reportTime;

	static bool reportXLPeptide;
	static bool reportXLScore;
	static bool reportXLExpectation;
	static bool reportXLMValue;
	static bool reportXLRank;
	static bool reportXLLowScore;
	static bool reportXLLowEVal;
	static bool reportXLLowMValue;
	static bool reportXLAA;

	static int reportPreviousAA;
	static bool reportDBPeptide;
	static int reportNextAA;
	static bool extraInfo;

	std::string getPreviousAA1 () const { return xldi.getPreviousAA1 ( 0, reportPreviousAA ); }
	std::string getPreviousAA2 () const { return xldi.getPreviousAA2 ( 0, reportPreviousAA ); }
	std::string getNextAA1 () const { return xldi.getNextAA1 ( 0, getEndAA1 (), reportNextAA ); }
	std::string getNextAA2 () const { return xldi.getNextAA2 ( 0, getEndAA2 (), reportNextAA ); }
	int getEndAA1 () const { return getStartAA1 () + getDBPeptide1 ().length () - 1; }
	int getEndAA2 () const { return getStartAA2 () + getDBPeptide2 ().length () - 1; }
	std::string getDBPeptide1 () const
	{
		if ( hitPeptide1->getDatabasePeptide () == "" ) return hitPeptide1->getPeptide ();
		else return hitPeptide1->getDatabasePeptide ();
	}
	std::string getDBPeptide2 () const
	{
		if ( hitPeptide2->getDatabasePeptide () == "" ) return hitPeptide2->getPeptide ();
		else return hitPeptide2->getDatabasePeptide ();
	}
	void runMSProduct ( int i, double score, const LinkInfo* linkInfo ) const;
public:
	SearchResultsCrosslinkPeptideHit ( const SpecID* specID, const MSMSSpectrumInfo* mmsi, const PeptideSpectralInfo* psi, double error, const HitPeptide* hitPeptide1, const HitPeptide* hitPeptide2, const XLinkDecoyInfo& xldi, int searchIndex, double xScore1, int xRank1, double xScore2, int xRank2, double xFirstScore );
	static void printHeaderHTML ( std::ostream& os, const std::string& styleID, int searchIdx );
	void printProteinInfoHTML ( std::ostream& os, const SResLink& sresLink ) const;
	void printHTML ( std::ostream& os, const std::string& styleID, const MSProductLink* productLink, const SCMSTagLink& smtl, const LinkInfo* linkInfo ) const;
	void printHeaderDelimited ( std::ostream& os, int searchIdx, const PPProteinHitQuanInfo& ppphqi ) const;
	void printDelimited ( std::ostream& os, const PPProteinHitQuanInfo& ppphqi ) const;
	static void setReportLinks ( bool f ) { reportLinks = f; }
	static void setReportError ( bool f ) { reportError = f; }
	int getFraction () const { return specID->getFraction (); }
	const SpecID getSpecIDasID () const { return *specID; }
	double getSpotAsNumber () const { return specID->getSpotAsNumber (); }
	double getMPlusH () const { return mOverZToMPlusH ( getMOverZ (), getCharge (), true ); }
	double getMOverZ () const { return mmsi->getMOverZ (); }
	double getIntensity () const { return mmsi->getIntensity (); }
	int getCharge () const { return mmsi->getCharge (); }
	double getError () const { return error; }
	double getEValue () const { return peptideHitInfo.getExpectationValue (); }
	double getScore () const { return peptideHitInfo.getScore (); }
	int getStartAA1 () const { return xldi.getStartAA1 ( 0 ); }
	int getStartAA2 () const { return xldi.getStartAA2 ( 0 ); }
	int getMModPosn1 () const { return getStartAA1 () + hitPeptide1->getMModPosition () - 1; }
	int getMModPosn2 () const { return getStartAA2 () + hitPeptide2->getMModPosition () - 1; }
	bool isDecoy1 () const { return xldi.isDecoy1 ( 0 ); }
	bool isDecoy2 () const { return xldi.isDecoy2 ( 0 ); }
	std::string getAccessionNumber1 () const { return xldi.getFullAccessionNumber1 ( 0 ); }
	std::string getAccessionNumber2 () const { return xldi.getFullAccessionNumber2 ( 0 ); }
	std::string getAccessionNumbers () const
	{
		const std::string& a1 = getAccessionNumber1 ();
		const std::string& a2 = getAccessionNumber2 ();
		if ( a1 == a2 ) return a1;
		if ( a1 < a2 )	return a1 + "\t" + a2;
		else			return a2 + "\t" + a1;
	}
	bool getIntramolecular () const { return getAccessionNumber1 () == getAccessionNumber2 () && !getDoubleDecoy (); }
	bool getIntermolecular () const { return getAccessionNumber1 () != getAccessionNumber2 () && !getDecoy () && !getDoubleDecoy (); }
	bool getDecoy () const { return isDecoy1 () != isDecoy2 (); }			// Exclusive OR
	bool getDoubleDecoy () const { return isDecoy1 () && isDecoy2 (); }
	std::string getXLinkType () const
	{
		if ( getIntramolecular () )			return "Intramolecular";
		else if ( getIntermolecular () )	return "Intermolecular";
		else if ( getDecoy () )				return "Decoy";
		else if ( getDoubleDecoy () )		return "Double Decoy";
		else								return "";
	}
	std::string getSequence1 () const
	{
		return hitPeptide1->getSequence ();
	}
	std::string getSequence2 () const
	{
		return hitPeptide2->getSequence ();
	}
	static void init ();
	void setQuanResults ( const LinkInfo* linkInfo ) const;
	bool outputQuanResults ( std::ostream& os, bool area ) const;
	DoubleVector getRatios ( bool area ) const;
};
int sortXlinkPeptideHits ( const SearchResultsCrosslinkPeptideHit* a, const SearchResultsCrosslinkPeptideHit* b, int ineq );
class sortCrosslinkHitsByMOverZ {
public:
	int operator () ( const SearchResultsCrosslinkPeptideHit* a, const SearchResultsCrosslinkPeptideHit* b ) const
	{
		return sortXlinkPeptideHits ( a, b, a->getMOverZ () < b->getMOverZ () );
	}
};
class sortCrosslinkHitsByMPlusH {
public:
	int operator () ( const SearchResultsCrosslinkPeptideHit* a, const SearchResultsCrosslinkPeptideHit* b ) const
	{
		return sortXlinkPeptideHits ( a, b, a->getMPlusH () < b->getMPlusH () );
	}
};
class sortCrosslinkHitsByError {
public:
	int operator () ( const SearchResultsCrosslinkPeptideHit* a, const SearchResultsCrosslinkPeptideHit* b ) const
	{
		return sortXlinkPeptideHits ( a, b, a->getError () < b->getError () );
	}
};
class sortCrosslinkHitsByIntensity {
public:
	int operator () ( const SearchResultsCrosslinkPeptideHit* a, const SearchResultsCrosslinkPeptideHit* b ) const
	{
		return sortXlinkPeptideHits ( a, b, a->getIntensity () > b->getIntensity () );
	}
};
class sortCrosslinkHitsBySpot {
public:
	int operator () ( const SearchResultsCrosslinkPeptideHit* a, const SearchResultsCrosslinkPeptideHit* b ) const
	{
		return sortXlinkPeptideHits ( a, b, a->getSpecIDasID () < b->getSpecIDasID () );
	}
};
class sortCrosslinkHitsByRT {
public:
	int operator () ( const SearchResultsCrosslinkPeptideHit* a, const SearchResultsCrosslinkPeptideHit* b ) const
	{
		return sortXlinkPeptideHits ( a, b, a->getSpotAsNumber () < b->getSpotAsNumber () );
	}
};
class sortCrosslinkHitsByStartResidue {
public:
	int operator () ( const SearchResultsCrosslinkPeptideHit* a, const SearchResultsCrosslinkPeptideHit* b ) const
	{
		return sortXlinkPeptideHits ( a, b, a->getStartAA1 () != b->getStartAA1 () ? a->getStartAA1 () < b->getStartAA1 () : a->getStartAA2 () < b->getStartAA2 () );
	}
};
class sortCrosslinkHitsByEndResidue {
public:
	int operator () ( const SearchResultsCrosslinkPeptideHit* a, const SearchResultsCrosslinkPeptideHit* b ) const
	{
		return sortXlinkPeptideHits ( a, b, a->getStartAA1 () != b->getStartAA1 () ? a->getStartAA1 () < b->getStartAA1 () : a->getStartAA2 () < b->getStartAA2 () );
	}
};
class sortCrosslinkHitsByXLinkPosn {
public:
	int operator () ( const SearchResultsCrosslinkPeptideHit* a, const SearchResultsCrosslinkPeptideHit* b ) const
	{
		return sortXlinkPeptideHits ( a, b, a->getMModPosn1 () != b->getMModPosn1 () ? a->getMModPosn1 () < b->getMModPosn1 () : a->getMModPosn2 () < b->getMModPosn2 () );
	}
};
class sortCrosslinkHitsByScore {
public:
	int operator () ( const SearchResultsCrosslinkPeptideHit* a, const SearchResultsCrosslinkPeptideHit* b ) const
	{
		return sortXlinkPeptideHits ( a, b, a->getScore () > b->getScore () );
	}
};
class sortCrosslinkHitsByExpectationValue {
public:
	int operator () ( const SearchResultsCrosslinkPeptideHit* a, const SearchResultsCrosslinkPeptideHit* b ) const
	{
		return sortXlinkPeptideHits ( a, b, a->getEValue () < b->getEValue () );
	}
};

class SearchResultsCrosslinkProteinHit {
	std::vector <const SearchResultsCrosslinkPeptideHit*> cLinkPeptideHits;
	double bestExpectationValue;

	PPProteinHitQuanInfo ppphqi;

	static bool reportLinks;
	void quanPlot ( std::ostream& os, bool area ) const;
public:
	SearchResultsCrosslinkProteinHit ();
	void push_back ( const SearchResultsCrosslinkPeptideHit* hit )
	{
		if ( cLinkPeptideHits.empty () )	bestExpectationValue = hit->getEValue ();
		else								bestExpectationValue = genMin ( bestExpectationValue, hit->getEValue () );
		cLinkPeptideHits.push_back ( hit );
	}
	void printHTML ( std::ostream& os, const SResLink& sresLink, int searchNumber, const LinkInfo* linkInfo ) const;
	void printHeaderDelimited ( std::ostream& os, int searchIdx ) const;
	void printDelimited ( std::ostream& os ) const;
	static void init ();
	double getBestExpectationValue () const { return bestExpectationValue; }
	bool getIntermolecular () const { return cLinkPeptideHits [0]->getIntermolecular (); }
	bool getIntramolecular () const { return cLinkPeptideHits [0]->getIntramolecular (); }
	bool getDecoy () const { return cLinkPeptideHits [0]->getDecoy (); }
	bool getDoubleDecoy () const { return cLinkPeptideHits [0]->getDoubleDecoy (); }
	std::string getXLinkType () const { return cLinkPeptideHits [0]->getXLinkType (); }
	void addSpecIDs ( SetSpecID& ssid ) const;
	void quanProteinStats ( bool area );
	std::string getAccessionNumber1 () const { return cLinkPeptideHits [0]->getAccessionNumber1 (); }
	std::string getAccessionNumber2 () const { return cLinkPeptideHits [0]->getAccessionNumber2 (); }
};

class sortCrosslinkProteinHitsByBestExpectationValue {
public:
	int operator () ( const SearchResultsCrosslinkProteinHit* a, const SearchResultsCrosslinkProteinHit* b ) const
	{
		bool decoy1 = a->getDecoy ();
		bool decoy2 = b->getDecoy ();
		bool doubleDecoy1 = a->getDoubleDecoy ();
		bool doubleDecoy2 = b->getDoubleDecoy ();
		bool inter1 = a->getIntermolecular ();
		bool inter2 = b->getIntermolecular ();
		bool intra1 = a->getIntramolecular ();
		bool intra2 = b->getIntramolecular ();
		if ( (decoy1 && decoy2) || (doubleDecoy1 && doubleDecoy2) || (inter1 && inter2) || (intra1 && intra2) ) {
			return a->getBestExpectationValue () < b->getBestExpectationValue ();
		}
		else {
			int num1 = doubleDecoy1 ? 1 : ( decoy1 ? 2 : ( inter1 ? 4 : 3 ) );
			int num2 = doubleDecoy2 ? 1 : ( decoy2 ? 2 : ( inter2 ? 4 : 3 ) );
			return num1 > num2;
		}
	}
};

#endif /* ! __sc_xlink_h */
