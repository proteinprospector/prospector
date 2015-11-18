/******************************************************************************
*                                                                             *
*  Program    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_ambiguity.h                                                *
*                                                                             *
*  Created    : February 15th 2011                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2011-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_ambiguity_h
#define __lu_ambiguity_h

#include <lgen_define.h>

class ModificationAmbiguity {
	VectorVectorVectorPairStringInt vvvpsi;	// A element in vvpsi is a single modification state (ie the mods and their positions)
	VectorVectorPairStringInt modArray;
	VectorPairStringInt curMod;
	VectorPairStringInt unambigMods;
	IntVector slip;
	mutable int cur;
	void initUnambiguity ( const std::string& unambiguity, bool chopScores );
	void initAmbiguity ( const std::string& ambiguity );
	void parseSimple ( VectorVectorPairStringInt& vvpsi, const std::string& ms1 );
	void parseComplex ( VectorVectorPairStringInt& vvpsi, const std::string& ms1 );
	void getNextMod ( int level );

	StringVector mNames;
	std::string complexVar;
	int vIdx;
	VectorPairStringInt compPairs;
	void parseNextComplex ( int start, VectorVectorPairStringInt& vvpsi );

	static int getResidueIdx ( const std::string& resStr );
public:
	ModificationAmbiguity ( const std::string& mods, bool chopScores = true );
	bool getNextSequence ( std::string& s, std::string& nTerm, std::string& cTerm, std::string& nLoss ) const;
	void getUnambiguousIndexList ( const std::string& mod, IntVector& sites, IntVector& scores ) const;
};

class SpecID;

typedef std::map <const SpecID*, MapStringToString> MapSpecIDToScoreString;
typedef MapSpecIDToScoreString::iterator MapSpecIDToScoreStringIterator;
typedef MapSpecIDToScoreString::const_iterator MapSpecIDToScoreStringConstIterator;

typedef std::pair <int, std::string> PairIntString;
typedef std::map <PairIntString, double> MapPairIntStringToDouble;
typedef MapPairIntStringToDouble::iterator MapPairIntStringToDoubleIterator;
typedef MapPairIntStringToDouble::const_iterator MapPairIntStringToDoubleConstIterator;

typedef std::set <MapPairIntStringToDouble> SetMapPairIntStringToDouble;
typedef SetMapPairIntStringToDouble::iterator SetMapPairIntStringToDoubleIterator;
typedef SetMapPairIntStringToDouble::const_iterator SetMapPairIntStringToDoubleConstIterator;

typedef std::set <PairIntString> SetPairIntString;
typedef SetPairIntString::iterator SetPairIntStringIterator;
typedef SetPairIntString::const_iterator SetPairIntStringConstIterator;

typedef std::set <SetPairIntString> SetSetPairIntString;
typedef SetSetPairIntString::iterator SetSetPairIntStringIterator;
typedef SetSetPairIntString::const_iterator SetSetPairIntStringConstIterator;

class ZeroMod {
	SetSetPairIntString zm;
	IntVectorVector ivv;
	int siz1;
	int siz2;
public:
	ZeroMod () {}
	void insertIfEmpty ( const MapPairIntStringToDouble& tm );
	void insert ( const MapPairIntStringToDouble& tm );
	void deleteMod ( const PairIntString& pss );
	void print ( std::ostream& os );
	void print2 ( std::ostream& os );
	void writeNext ( std::ostream& os, int start, int end, int level );
};

class SCModInfo {
	//static vector <std::map <pair <string, PairStringInt >, pair <const SpecID*, double> > > bSS;
	static std::vector <MapSpecIDToScoreString> msiss;
	static double threshold;
	static StringVectorVector cMods;
	static StringVector cNTerm;
	static StringVector cCTerm;
	StringVector topDBPeptide;
	DoubleVector topDBScore;
	bool eval;
	std::vector <MapPairIntStringToDouble> topMods;
	std::vector <ZeroMod> zeroMods;
	static MapPairIntStringToDouble getMods ( const std::string& nTerm, const std::string& peptide, const std::string& cTerm, const std::string& neutralLoss );
	static std::string addStartAAToModString ( const std::string& modString, int startAA, int len );
	std::string removeCMods ( const std::string& peptide );
public:
	SCModInfo () {}
	void addHit ( const std::string& nTerm, const std::string& peptide, const std::string& cTerm, const std::string& neutralLoss, double score, double expectation, int numSpectra );
	MapStringToString getScoreString ();
	void compress ();
	void add ( const SpecID* spID );
	static std::string getConstModsString ( const std::string& dbPeptide, int searchIdx );
	static std::string getModsString ( int searchNumber, const SpecID* spID, const std::string& dbPeptide, int startAA = 0 );
	static std::string getAllModsString ( int searchNumber, const SpecID* spID, const std::string& dbPeptide, int startAA = 0 );
	static void init ( int searchNumber, const StringVector& constMods, double threshold );
	//static std::string getAmbiguityString ( int searchNumber, const SpecID* spID, const std::string& dbPeptide );
	//static std::string getAmbiguityString ( const std::string& modsString );
	static void getAmbiguityString ( const std::string& modsString, std::string& ambString, std::string& unambString, bool chopScores = true );	// Gets an ambiguity string for passing to MS-Product
};
class ProbabilityAmbiguity {
	IntVector siteList;
	int siteListLength;
	IntVector sequence;
	int numSites;
	IntVectorVector hits;
	std::string modStr;
	void getNext ( int level );
public:
	ProbabilityAmbiguity ( const std::string& probStr, const std::string& mod, double probLimit );
	std::string getModStr () const { return modStr; }
};

struct SiteInfo {
	int site;
	char aa;
	int slip;
	int index;
	SiteInfo () :
		slip ( -2 ),
		index ( -2 ) {}
	SiteInfo ( int site, char aa, int slip, int index ) :
		site ( site ),
		aa ( aa ),
		slip ( slip ),
		index ( index ) {}

	bool getIndex ( int& idx ) const;

	static void printDelimitedHeaderSLIP ( std::ostream& os, int idx );
	void printDelimitedSLIP ( std::ostream& os ) const;
	bool printDelimitedAA ( std::ostream& os ) const;
	static void printHTMLHeaderSLIP ( std::ostream& os, int idx, int rowspan );
	void printHTMLSLIP ( std::ostream& os ) const;
	bool printHTMLAA ( std::ostream& os ) const;
};

struct GlycoSiteInfo : public SiteInfo {
	std::string mod;
	GlycoSiteInfo ( int site, char aa, int slip, int index, const std::string& mod ) :
		SiteInfo ( site, aa, slip, index ),
		mod ( mod ) {}
};

typedef std::vector <SiteInfo> SiteInfoVector;
typedef std::vector <SiteInfoVector> SiteInfoVectorVector;

typedef std::vector <GlycoSiteInfo> GlycoSiteInfoVector;

class ParameterList;

class SiteScores {
protected:
	static SetString glycoMods;
	static MapStringToMapCharToInt msmci;
	static int idx;
	static bool initialised;
	VectorMapIntToPairIntInt siteScores;
public:
	SiteScores ();
	static void init ( const ParameterList* pList );
	static void init ( const StringVector& mods, bool flag );
	static void init ( const StringVector& mods );
	static void initMod ( const std::string& mod, char aa );
	void add ( const std::string& peptide, const std::string& mods, int start, int index );
	void clear ();
	void getSiteInfo ( StringVector& vModString, SiteInfoVectorVector& vvSiteInfo, GlycoSiteInfoVector& vGlycoSiteInfo ) const;
};

typedef std::vector <SiteScores> SiteScoresVector;

typedef std::map <int, SiteInfoVector> MapIntToSiteInfoVector;
typedef MapIntToSiteInfoVector::iterator MapIntToSiteInfoVectorIterator;
typedef MapIntToSiteInfoVector::const_iterator MapIntToSiteInfoVectorConstIterator;

typedef std::map <std::string, MapIntToSiteInfoVector> MapStringToMapIntToSiteInfoVector;
typedef MapStringToMapIntToSiteInfoVector::iterator MapStringToMapIntToSiteInfoVectorIterator;
typedef MapStringToMapIntToSiteInfoVector::const_iterator MapStringToMapIntToSiteInfoVectorConstIterator;

typedef std::map <std::string, SiteInfoVector> MapStringToSiteInfoVector;
typedef MapStringToSiteInfoVector::iterator MapStringToSiteInfoVectorIterator;
typedef MapStringToSiteInfoVector::const_iterator MapStringToSiteInfoVectorConstIterator;

typedef std::map <int, MapStringToSiteInfoVector> MapIntToMapStringToSiteInfoVector;
typedef MapIntToMapStringToSiteInfoVector::iterator MapIntToMapStringToSiteInfoVectorIterator;
typedef MapIntToMapStringToSiteInfoVector::const_iterator MapIntToMapStringToSiteInfoVectorConstIterator;

#endif /* ! __lu_ambiguity_h */
