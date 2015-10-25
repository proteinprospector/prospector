/******************************************************************************
*                                                                             *
*  Program    : searchCompare                                                 *
*                                                                             *
*  Filename   : sc_search_res.h                                               *
*                                                                             *
*  Created    : March 21st 2003                                               *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2003-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __sc_search_res_h
#define __sc_search_res_h

#include <lg_string.h>
#include <lu_aa_calc.h>
#include <lu_acc_num.h>
#include <lu_coverage.h>
#include <lu_disc_sc.h>
#include <lu_histogram.h>
#include <lu_mass_conv.h>
#include <lu_acc_link.h>
#include <lu_pros_form.h>
#include <lu_spec_id.h>
#include <lu_tol.h>
#include <lu_xml.h>
#include <lr_main.h>
#include <sc_params.h>

class TaxonomyMatch;
class FastaServer;
class UpdatingJavascriptMessage;
class LinkInfo;

typedef std::vector <FastaServer*> FastaServerPtrVector;

class ProteinInfo {
	static bool reportUniprotID;
	static bool reportGeneName;
	static bool reportAccession;
	static bool reportLength;
	static bool reportMW;
	static bool reportPI;
	static bool reportSpecies;
	static bool reportName;
	static bool reportLinks;
	static std::vector <TaxonomyMatch*> taxMatch;
	static bool taxCheck;
	static int preferredMatchLength;
	int databaseIndex;
	int indexNum;
	mutable std::string proteinSeq;
	StringVector accessionNumbers;
	StringVector accessionInfo;
	StringVector species;
	StringVector names;
	StringVector uniprotIDs;

	StringVector organismNames;
	StringVector geneNames;
	StringVector proteinExistences;
	StringVector sequenceVersion;

	int dbaseIndex;
	std::string aNum;
	std::string acc;
	int dnaReadingFrame;
	int openReadingFrame;
	double proteinMW;
	double proteinPI;
	int length;
	static IntVector dbIdx;
	static MapStringToInt dbNameMap;
	static FastaServerPtrVector fs;
	static std::vector <AccessionNumberMap*> am;
	static StringVector dbase;
	static std::vector <AccessionNumberLinkInfo> anli;
	static std::vector <AccessionNumberLinkInfo> upidli;
	static StringVector tempDirs;
	void setAnumInfo ( const std::string& str );
	static void getUserProteinIndex ( int& index, const std::string& db, const std::string& possDB, const std::string& diskDB );
public:
	ProteinInfo ();
	ProteinInfo ( const std::string& aNum );
	void setFields ();
	int getPreferredSpeciesIndex () const;
	std::string getAccessionInfo () const { return accessionInfo [getPreferredSpeciesIndex ()]; }
	std::string getSpecies () const { return species [getPreferredSpeciesIndex ()]; }
	std::string getUniprotID () const { return uniprotIDs [getPreferredSpeciesIndex ()]; }

	std::string getOrganismName () const	{ return organismNames [getPreferredSpeciesIndex ()]; }
	std::string getGeneName () const		{ return geneNames [getPreferredSpeciesIndex ()]; }
	std::string getProteinExistence () const{ return proteinExistences [getPreferredSpeciesIndex ()]; }
	std::string getSequenceVersion () const	{ return sequenceVersion [getPreferredSpeciesIndex ()]; }

	std::string getFirstAccessionNumber () const { return accessionNumbers [0]; }

	std::string getName () const { return names [getPreferredSpeciesIndex ()]; }
	std::string getProteinSequence () const;
	double getProteinMW () const { return proteinMW; }
	double getProteinPI () const { return proteinPI; }
	int getLength () const { return length; }
	std::string getFullAccessionNumber () const { return aNum; }
	std::string getAcc () const { return acc; }
	bool isDecoy () const;
	static bool isDecoy ( const std::string& a );
	static int getColspan ();
	static void printHTMLHeader ( std::ostream& os, int rowspan );
	static void printHTMLANumHeader ( std::ostream& os, int rowspan );
	void printHTML ( std::ostream& os, bool empty ) const;
	void printHTMLANum ( std::ostream& os, bool empty ) const;
	void printHTMLLines ( std::ostream& os ) const;
	static void printDelimitedHeader ( std::ostream& os );
	static void printDelimitedANumHeader ( std::ostream& os );
	void printDelimited ( std::ostream& os ) const;
	void printDelimitedANum ( std::ostream& os ) const;
	static IntVector initialise ( const ParameterList* params );
	static void initialiseAccessionNumberLink ( std::ostream& os );
	static void setReportUniprotID ( bool f ) { reportUniprotID = f; }
	static void setReportGeneName ( bool f ) { reportGeneName = f; }
	static void setReportAccession ( bool f ) { reportAccession = f; }
	static void setReportLength ( bool f ) { reportLength = f; }
	static void setReportMW ( bool f ) { reportMW = f; }
	static void setReportPI ( bool f ) { reportPI = f; }
	static void setReportSpecies ( bool f ) { reportSpecies = f; }
	static void setReportName ( bool f ) { reportName = f; }
	static void setReportLinks ( bool f ) { reportLinks = f; }
	static bool getReport ()
	{
		return reportGeneName || reportUniprotID || reportAccession || reportLength || reportMW || reportPI || reportSpecies || reportName;
	}
	static bool getReportAccession () { return reportAccession;	}
	static bool getReportUniprotID () { return reportUniprotID; }
	static bool getReportGeneName () { return reportGeneName; }
	static void setTaxonomyMatch ( const StringVector& preferredSpecies );
	static std::vector <TaxonomyMatch*> getTaxonomyMatch () { return taxMatch; }
	static bool getTaxonomyCheck () { return taxCheck; }
	void initProteinSequence () const;
	std::string getPreviousAA ( int startAA, int numAA ) const;
	std::string getNextAA ( int endAA, int numAA ) const;
	static time_t getDatabaseTime ();
	static std::string getDatabasePath ();
	static std::string getDatabaseName ( int i );
	static std::string getDatabasePath ( int i );
	static int fsSize ();
	static bool getDNADatabase ();
	static int getNumEntries ();
	static int getNumEntries ( int i );
	static StringVectorSizeType getNumDatabases () { return dbase.size (); }
	static PairIntString getANumPair ( const std::string& aNum );
	static StringVector getDBSearchList1 ();
	static BoolDeque getDBSearchFlagList1 ();
	static bool getNonDecoyFlag ( int idx );
	static int getDBIndex ( const std::string& db );
	int getDatabaseIndex () const { return dbaseIndex; }
	std::string getDatabaseMZIdentMLRef () const;
	static void addDatabaseMZIdentMLInfo ( VectorXMLOutputItemPtr& subItems );
	static void resetANumMap ();
	static void deleteTempDirs ();
};

class PeakFitData;

class MSMSSpectrumInfo {
	double mOverZ;
	float intensity;
	short charge;
	short numPeaks;
public:
	MSMSSpectrumInfo () {}
	MSMSSpectrumInfo ( double mOverZ, int charge, double intensity, int numPeaks );
	double getMOverZ () const { return mOverZ; }
	int getCharge () const { return charge; }
	double getIntensity () const { return intensity; }
	int getNumPeaks () const { return numPeaks; }
};
typedef std::map <std::string, MSMSSpectrumInfo*> MapSpecIDMSMSSpectrumInfo;
typedef std::map <std::string, MapSpecIDMSMSSpectrumInfo> MapIDSpecIDMSMSSpectrumInfo;

class SearchResults;
typedef std::vector <SearchResults*> SearchResultsPtrVector;
typedef SearchResultsPtrVector::size_type SearchResultsPtrVectorSizeType;

class HitPeptide {
	static MapStringToDouble msd;
	std::string databasePeptide;
	std::string nTerm;
	std::string peptide;
	std::string cTerm;
	std::string neutralLoss;
	bool mmod;
public:
	HitPeptide ( const std::string& nTerm, const std::string& peptide, const std::string& cTerm, const std::string& neutralLoss );
	std::string getDatabasePeptide () const
	{
		if ( databasePeptide.empty () )	return peptide;
		else							return databasePeptide;
	}
	void getModMassesAndIndicies ( VectorPairIntDouble& vpid ) const;
	void getModMassesIndiciesAndString ( VectorPairIntPairStringDouble& vpid ) const;
	double getModMass ( const std::string& s ) const;
	double getModNTermMass () const { return getModMass ( nTerm ); }
	std::string getNTerm () const { return nTerm; }
	std::string getPeptide () const { return peptide; }
	std::string getBlibPeptide () const;
	void getBlibMods ( VectorPairIntDouble& vpid ) const;
	double getModCTermMass () const { return getModMass ( cTerm ); }
	std::string getCTerm () const { return cTerm; }
	std::string getNeutralLoss () const { return neutralLoss; }
	std::string getModString () const;
	int getPeptideLength () const
	{
		return gen_strstriptags2 ( peptide, '(', ')' ).length ();
	}
	bool getMmod () const { return mmod; }
	double getMModValue () const;
	int getMModPosition () const;
	std::string getCommandLineNVPair ( int i ) const;
	std::string getSequence () const
	{
		return nTerm + '-' + peptide + '-' + cTerm + '+' + neutralLoss;
	}
};
class PeptidePosition {
	std::string accessionNumber;
	const HitPeptide* hitPeptide;
	const SpecID* specID;
	const MSMSSpectrumInfo* mmsi;
	double error;
	mutable PeakFitData* quanRatio;
	int startAA;
	mutable bool firstOccurence;
	short searchIndex;
	static double defaultResolution;
	static int numIsotopePeaks;
	static StringVector instrument;
	static BoolDeque chargeReducedFragmentation;
	static StringVectorVector fractionNames;
	static BoolDeque spottingPlates;
	static BoolDeque spectrumNumber;
	static StringVectorVector rawTypes;
	static VectorConstParameterListPtr params;
	static StringVector searchKey;
	static std::vector <Tolerance*> parentTolerances;
	static StringVector sysErrorStr;
	static DoubleVectorVector offsets;
	static std::vector <Tolerance*> fragmentTolerances;
	static AACalculator* aaCalc;
	static std::string compMaskType;
	static unsigned int compMask;
	static MapStringToUInt compMaskMap;
	static bool multipleErrorUnits;
	static bool multipleFractionNames;
	static bool spottingPlatesFlag;
	static bool spectrumNumberFlag;
	static bool quanMultiNormalFlag;
	static bool enzymeInit;
	static bool reportCheckboxes;
	static bool reportSearchNumber;
	static bool reportMPlusH;
	static bool reportMOverZ;
	static bool reportCharge;
	static bool reportIntensity;
	static bool reportMPlusHCalc;
	static bool reportMOverZCalc;
	static bool reportError;
	static bool reportDBPeptide;
	static std::string peptideModType;
	static bool reportProteinMods;
	static bool reportTime;
	static bool reportMSMSInfo;
	static bool reportStartAA;
	static bool reportEndAA;
	static int reportPreviousAA;
	static int reportNextAA;
	static bool reportMissedCleavages;
	static bool reportLength;
	static bool reportComposition;
	static bool reportMModValue;
	static bool reportLinks;
	static bool runMSProductFlag;
	static double rtIntervalStart;
	static double rtIntervalEnd;
	friend class SortPeptidePositionAscending2;
	friend class SortPeptidePositionAscending3;
	friend class PeptidePositionQuan;
	unsigned int peptideToMask () const;
	void setSpectrumNumber ( const SpecID* spID, int searchIndex );
public:
	PeptidePosition ( const std::string& accessionNumber, const HitPeptide* hitPeptide, double error, const SpecID* spID, const MSMSSpectrumInfo* mmsi, int startAA, int searchIndex );
	PeptidePosition ( const SpecID* spID, const MSMSSpectrumInfo* mmsi, int searchIndex );
	void setSubsequentOccurence () const { firstOccurence = false; }
	void setQuanResults () const;
	bool checkComposition () const;
	std::string getAccessionNumber () const { return accessionNumber; }
	bool isDecoyHit () const;
	int getFraction () const { return specID->getFraction (); }
	std::string getSpot () const { return specID->getID (); }
	std::string getMSMSInfo () const { return specID->getMSMSInfo (); }
	int getRun () const { return specID->getRun (); }
	double getSpotAsNumber () const { return specID->getSpotAsNumber (); }
	static StringVector getInstrument () { return instrument; }
	static std::string getInstrument ( int i ) { return instrument [i]; }
	std::string getFractionName () const { return fractionNames [searchIndex][specID->getFraction ()-1]; }
	std::string getSpecID () const { return specID->getSpecID (); }
	const SpecID getSpecIDasID () const { return *specID; }
	bool hasSpecID () const { return specID != 0; }
	bool hasHitPeptide () const { return hitPeptide != 0; }
	std::string getFullSpecID () const { return specID->getFullSpecID (); }
	std::string getDBPeptide () const
	{
		if ( hitPeptide == 0 ) return "";
		else {
			if ( hitPeptide->getDatabasePeptide () == "" ) return hitPeptide->getPeptide ();
			else return hitPeptide->getDatabasePeptide ();
		}
	}
	std::string getNextAA ( const ProteinInfo& proteinInfo ) const
	{
		if ( hitPeptide == 0 ) return "";
		else return proteinInfo.getNextAA ( getEndAA (), genMax ( reportNextAA, 1 ) );
	}
	std::string getPrevAA ( const ProteinInfo& proteinInfo ) const
	{
		if ( hitPeptide == 0 ) return "";
		else return proteinInfo.getPreviousAA ( startAA, genMax ( reportPreviousAA, 1 ) );
	}
	std::string getModString () const { return hitPeptide == 0 ? "" : hitPeptide->getModString (); }
	void getModMassesAndIndicies ( VectorPairIntDouble& vpid ) const { hitPeptide->getModMassesAndIndicies ( vpid ); }
	void getModMassesIndiciesAndString ( VectorPairIntPairStringDouble& vpid ) const { hitPeptide->getModMassesIndiciesAndString ( vpid ); }
	std::string getNTerm () const { return hitPeptide == 0 ? "" : hitPeptide->getNTerm (); }
	double getModNTermMass () const { return hitPeptide == 0 ? 0.0 : hitPeptide->getModNTermMass (); }
	std::string getPeptide () const { return hitPeptide == 0 ? "" : hitPeptide->getPeptide (); }
	std::string getBlibPeptide () const { return hitPeptide == 0 ? "" : hitPeptide->getBlibPeptide (); }
	void getBlibMods ( VectorPairIntDouble& vpid ) const
	{
		if ( hitPeptide != 0 ) hitPeptide->getBlibMods ( vpid );
	}
	int getPeptideLength () const { return hitPeptide == 0 ? 0 : hitPeptide->getPeptideLength (); }
	std::string getCTerm () const { return hitPeptide == 0 ? "" : hitPeptide->getCTerm (); }
	double getModCTermMass () const { return hitPeptide == 0 ? 0.0 : hitPeptide->getModCTermMass (); }
	std::string getNeutralLoss () const { return hitPeptide == 0 ? "" : hitPeptide->getNeutralLoss (); }
	bool getMmod () const { return hitPeptide->getMmod (); }
	double getMModValue () const { return hitPeptide == 0 ? 0.0 : hitPeptide->getMModValue (); }
	std::string getSequence () const;
	std::string getSequencePlusMods () const;
	std::string getPrintedSequence () const;
	int getNumTolTerm ( const ProteinInfo& proteinInfo ) const;
	std::string getMissedCleavages () const;
	double getMOverZ () const { return mmsi->getMOverZ (); }
	int getCharge () const { return mmsi->getCharge (); }
	int getMaxFragmentCharge ( int searchNumber ) const
	{
		return mmsi->getCharge ();
	}
	double getM () const { return mOverZToM ( getMOverZ (), getCharge (), true ); }
	double getMPlusH () const { return mOverZToMPlusH ( getMOverZ (), getCharge (), true ); }
	double getMCalc ( int searchNumber ) const { return parentTolerances [searchNumber]->getActualMass ( getM (), error, getCharge () ); }
	double getMPlusHCalc ( int searchNumber ) const { return parentTolerances [searchNumber]->getActualMass ( getMPlusH (), error, getCharge () ); }
	double getMOverZCalc ( int searchNumber ) const { return mPlusHToMOverZ ( getMPlusHCalc ( searchNumber ), getCharge (), true ); }
	double getDaError ( int searchNumber ) const { return getMOverZ () - getMOverZCalc ( searchNumber ); }
	double getIntensity () const { return mmsi->getIntensity (); }
	int getNumPeaks () const { return mmsi->getNumPeaks (); }
	double getPeptideError () const
	{
		if ( getMmod () )	return PeptidePosition::INVALID_ERROR;
		else				return error;
	}
	double getError () const { return error; }
	int getStartAA () const { return startAA; }
	int getEndAA () const { return startAA + getDBPeptide ().length () - 1; }
	bool getFirstOccurence () const { return firstOccurence; }
	short getSearchIndex () const { return searchIndex; }
	static int getColspan ( int searchNumber, bool tabDelim = false, bool showTimes = true );
	static int getColspanPeak1 ();
	static int getColspanPeak2 ( int searchNumber, bool tabDelim );
	static void printHeaderHTML ( std::ostream& os, int searchNumber, const std::string& styleID, bool showTimes = true );
	static void printHeaderHTMLPeak1 ( std::ostream& os, const std::string& styleID, int rowspan = 0 );
	static void printHeaderHTMLPeak2 ( std::ostream& os, int searchNumber, const std::string& styleID, int rowspan = 0 );
	void printHTML ( std::ostream& os, const ProteinInfo& proteinInfo, const std::string& styleID, const SCMSTagLink& smtl, const std::string& id, bool joint, bool showTimes = true ) const;
	void printHTMLPeak1 ( std::ostream& os, const std::string& styleID ) const;
	void printHTMLPeak2 ( std::ostream& os, const std::string& styleID, const SCMSTagLink& smtl, bool joint ) const;
	static void printHeaderDelimited ( std::ostream& os, int searchNumber, bool showTimes = true );
	static void printHeaderDelimitedPeak1 ( std::ostream& os );
	static void printHeaderDelimitedPeak2 ( std::ostream& os, int searchNumber );
	void printDelimited ( std::ostream& os, const ProteinInfo& proteinInfo, bool joint, bool showTimes = true ) const;
	void printDelimitedPeak1 ( std::ostream& os ) const;
	void printDelimitedPeak2 ( std::ostream& os, bool joint ) const;
	bool outputQuanResults ( std::ostream& os, const std::string& searchName, int numRepeats, bool area ) const;
	DoubleVector getIntensityRatios () const;
	DoubleVector getAreaRatios () const;
	static void printParamsHTML (std::ostream& os, int searchNumber );
	static void printHTMLEmpty ( std::ostream& os, int searchNumber, const std::string& styleID, bool showTimes = true );
	static void printHTMLEmptyPeak1 ( std::ostream& os, const std::string& styleID );
	static void printHTMLEmptyPeak2 ( std::ostream& os, int searchNumber, const std::string& styleID );
	static void printDelimitedEmpty ( std::ostream& os, int searchNumber, bool showTimes = true );
	static void printDelimitedEmptyPeak1 ( std::ostream& os );
	static void printDelimitedEmptyPeak2 ( std::ostream& os, int searchNumber );
	static void initialiseParams ( const ParameterList* p );
	static void initialiseDataSetInfo ( double timeWindowStart, double timeWindowEnd );
	static void initialise ( const SearchResultsPtrVector& searchResults );
	static void initialiseComposition ( const StringVector& compIons, const StringVector& massCompIons, const std::string& maskType );
	static void setReportCheckboxes	( bool f )	{ reportCheckboxes = f; }
	static void setReportMPlusH		( bool f )	{ reportMPlusH = f; }
	static void setReportMOverZ		( bool f )	{ reportMOverZ = f; }
	static void setReportCharge		( bool f )	{ reportCharge = f; }
	static void setReportIntensity	( bool f )	{ reportIntensity = f; }
	static void setReportSearchNumber( bool f )	{ reportSearchNumber = f; }
	static void setReportMPlusHCalc	( bool f )	{ reportMPlusHCalc = f; }
	static void setReportMOverZCalc	( bool f )	{ reportMOverZCalc = f; }
	static void setReportError		( bool f )	{ reportError = f; }
	static void setReportDBPeptide	( bool f )	{ reportDBPeptide = f; }
	static void setPeptideModType	( const std::string& s )	{ peptideModType = s; }
	static void setReportProteinMods( bool f )	{ reportProteinMods = f; }
	static void setReportTime		( bool f )	{ reportTime = f; }
	static void setReportMSMSInfo	( bool f )	{ reportMSMSInfo = f; }
	static void setReportStartAA	( bool f )	{ reportStartAA = f; }
	static void setReportEndAA		( bool f )	{ reportEndAA = f; }
	static void setReportPreviousAA	( int num )	{ reportPreviousAA = num; }
	static void setReportNextAA		( int num )	{ reportNextAA = num; }
	static void setReportMissedCleavages( bool f )	{ reportMissedCleavages = f; }
	static void setReportLength		( bool f )	{ reportLength = f; }
	static void setReportComposition( bool f )	{ reportComposition = f; }
	static void setReportMModValue	( bool f )	{ reportMModValue = f; }
	static void setReportLinks		( bool f )	{ reportLinks = f; }
	static void setRunMSProductFlag	( bool f )	{ runMSProductFlag = f; }

	static std::vector <Tolerance*> getParentTolerances () { return parentTolerances; }
	static std::vector <Tolerance*> getFragmentTolerances () { return fragmentTolerances; }
	static BoolDeque getSpottingPlates () { return spottingPlates; }
	static BoolDeque getSpectrumNumber () { return spectrumNumber; }
	static StringVectorVector getFractionNames () { return fractionNames; }
	static bool getMultipleFractionNames () { return multipleFractionNames; }
	static StringVectorVector getRawTypes () { return rawTypes; }

	static bool getReportMPlusH () { return reportMPlusH; }
	static bool getReportMOverZ () { return reportMOverZ; }
	static bool getReportCharge () { return reportCharge; }
	static bool getReportIntensity () { return reportIntensity; }
	static bool getReportMSMSInfo () { return reportMSMSInfo; }
	static bool getReportStartAA () { return reportStartAA; }
	static bool getReportEndAA () { return reportEndAA; }
	static bool getReportError () { return reportError; }
	static bool getReportTime () { return reportTime; }
	static bool getReportLinks () { return reportLinks; }
	static int getReportPreviousAA () { return reportPreviousAA; }
	static bool getReportDBPeptide () { return reportDBPeptide; }
	static int getReportNextAA () { return reportNextAA; }

	static const double INVALID_ERROR;
	static const ParameterList* getParams0 ()		{ return params [0]; }
	static const ParameterList* getParams ( int i )	{ return params [i]; }
	static const VectorConstParameterListPtr getParams () { return params; }
	static std::string getSearchKey ( int i )		{ return searchKey [i]; }
	static StringVector getSearchKey () { return searchKey; }
	static Tolerance* getParentTolerance ( int i )	{ return parentTolerances [i]; }
	static Tolerance* getFragmentTolerance ( int i )	{ return fragmentTolerances [i]; }
	static int getNumQuanRatioColumns ();
	static bool getReportComposition () { return reportComposition; }
	static bool getReportMModValue () { return reportMModValue; }
	static bool getRunMSProductFlag () { return runMSProductFlag; }
	static double getRTIntervalStart () { return rtIntervalStart; }
	static double getRTIntervalEnd () { return rtIntervalEnd; }
	static std::string getSysErrorStr ( int i ) { return sysErrorStr [i]; }
	void runMSProduct ( int i, double score ) const;
};
typedef std::vector <PeptidePosition> PeptidePositionVector;
typedef PeptidePositionVector::size_type PeptidePositionVectorSizeType;

#ifdef SEARCH_RES_MAIN
#define SEARCH_RES_EXTERN
#else
#define SEARCH_RES_EXTERN extern
#endif
SEARCH_RES_EXTERN bool sresProt;
SEARCH_RES_EXTERN bool sresTime;
SEARCH_RES_EXTERN bool sresXLinks;
SEARCH_RES_EXTERN bool sresFPR;
SEARCH_RES_EXTERN bool sresKeepReplicates;
SEARCH_RES_EXTERN bool sresKeepCharges;
SEARCH_RES_EXTERN bool sresKeepTimeReplicates;
SEARCH_RES_EXTERN bool sresSingleProject;
SEARCH_RES_EXTERN bool sresMergedFlag;
SEARCH_RES_EXTERN bool sresMainAndSupplementary;
SEARCH_RES_EXTERN UpdatingJavascriptMessage* ujm;
SEARCH_RES_EXTERN bool sresViewer;

class SearchResultsProteinInfo {
protected:
	double score;
	double bestScore;
	int numUnique;
	int peptideCount;
	static bool reportScore;
	static bool reportNumUnique;
	static bool reportPeptideCount;
	static bool reportBestScore;
	static bool reportCoverage;
	static bool reportBestDiscScore;
	static bool reportBestExpectVal;
public:
	SearchResultsProteinInfo ( double score, double bestScore, int numUnique, int peptideCount ) :
		score ( score ),
		bestScore ( bestScore ),
		numUnique ( numUnique ),
		peptideCount ( peptideCount ) {}
	virtual ~SearchResultsProteinInfo ();
	virtual int getColspan () const;
	virtual void printHeaderHTML ( std::ostream& os, int index, const std::string& styleID ) const;
	virtual void printHTML ( std::ostream& os, bool empty, const std::string& styleID ) const;
	virtual void printHeaderDelimited ( std::ostream& os, int index ) const;
	virtual void printDelimited ( std::ostream& os ) const;
	double getScore () const { return score; }
	double getBestScore () const { return bestScore; }
	int getNumUnique () const { return numUnique; }
	int getPeptideCount () const { return peptideCount; }
	virtual CoverageMap getCoverageMap () const { return CoverageMap (); }
	virtual double getBestExpectationValue () const { return -1; }
	virtual double getBestDiscriminantScore () const { return DiscriminantScore::MIN_DISC_SCORE; }
	virtual double getTotalDiscriminantScore () const { return DiscriminantScore::MIN_DISC_SCORE; }
	virtual void setAreaRatios ( const DoubleVector& m, const DoubleVector& q1, const DoubleVector& q2, const DoubleVector& mean, const DoubleVector& stDev, const IntVector& num ) {}
	virtual void setIntensityRatios ( const DoubleVector& m, const DoubleVector& q1, const DoubleVector& q2, const DoubleVector& mean, const DoubleVector& stDev, const IntVector& num ) {}
	void setScore ( double s ) { score = s; }
	void setBestScore ( double bs ) { bestScore = bs; }
	void setNumUnique ( int nu ) { numUnique = nu; }
	void setPeptideCount ( int pc ) { peptideCount = pc; }
	static void setParams ( const ParameterList* params );
};
typedef std::vector <const SearchResultsProteinInfo*> ConstSearchResultsProteinInfoPtrVector;
typedef ConstSearchResultsProteinInfoPtrVector::size_type ConstSearchResultsProteinInfoPtrVectorSizeType;

class PeptideSpectralInfo {
	int unmatched;
	int numPeaks;
	int rank;
	double score;
	double expectation;
	double scoreDiff;
	int numPrecursor;
	double a;
	double b;
	static const double NO_SCORE_DIFF;
	static bool linearTailFitExpectation;
	static bool noExpectation;
public:
	PeptideSpectralInfo ( int unmatched, int numPeaks, int rank, double score, double expectation, const std::string& scDiff, int numPrecursor, double a = 0.0, double b = 0.0 );
	int getMatched () const { return numPeaks - unmatched; }
	int getUnmatched () const { return unmatched; }
	int getNumPeaks () const { return numPeaks; }
	int getRank () const { return rank; }
	double getScore () const { return score; }
	double getExpectation () const { return expectation; }
	double getPValue () const
	{
		if ( expectation == -1.0 ) return -1.0;		// PValue invalid
		else return expectation / numPrecursor;
	}
	double getMascotScore () const
	{
		double pv = getPValue ();
		if ( pv == -1.0 ) return INVALID_MASCOT_SCORE;
		else return - 10.0 * log10 ( pv );
	}

	double getExpectationValue ( double score ) const;
	double getPValue ( double score ) const;
	double getMascotScore ( double score ) const;

	double getScoreDifference () const { return scoreDiff == NO_SCORE_DIFF ? score : scoreDiff; }
	int getNumPrecursor () const { return numPrecursor; }
	double getA () const { return a; }
	double getB () const { return b; }
	static void setExpectationFlag ( const std::string& expectationCalculationType );
	static const double INVALID_MASCOT_SCORE;
};

class PPPeptideHitInfo {
	const PeptideSpectralInfo* psi;
	double discriminantScore;
	int repeats;
	static const double NO_SCORE_DIFF;
	static bool reportUnmatched;
	static bool reportNumPks;
	static bool reportRank;
	static bool reportScore;
	static bool reportScoreDiff;
	static bool reportExpectation;
	static bool reportPValue;
	static bool reportMValue;
	static bool reportNumPrecursor;
	static bool reportGradient;
	static bool reportOffset;
	static bool reportDiscScore;
	static bool reportRepeats;
public:
	PPPeptideHitInfo ();
	PPPeptideHitInfo ( const PeptideSpectralInfo* psi );
	virtual ~PPPeptideHitInfo ();
	double getScore () const { return psi->getScore (); }
	double getScoreDifference () const { return psi->getScoreDifference (); }
	int getRepeats () const	{ return repeats; }
	int getMatched () const	{ return psi->getMatched (); }
	void setRepeats ( int r ) { repeats = r; }
	void setDiscriminantScore ( double s ) { discriminantScore = s; }
	double getDiscriminantScore () const { return discriminantScore; }
	double getExpectationValue () const { return psi->getExpectation (); }
	double getPValue () const { return psi->getPValue (); }
	double getMascotScore () const { return psi->getMascotScore (); }
	double getExpectationValue ( double score ) const { return psi->getExpectationValue ( score ); }
	double getPValue ( double score ) const { return psi->getPValue ( score ); }
	double getMascotScore ( double score ) const { return psi->getMascotScore ( score ); }
	int getNumPeaks () const { return psi->getNumPeaks (); }
	int getRank () const { return psi->getRank (); }
	int getNumPrecursor () const { return psi->getNumPrecursor (); }
	int getColspan () const;
	bool notEmpty () const { return psi->getScore () != 0.0; }
	static void printHeaderHTML ( std::ostream& os, const std::string& styleID );
	void printHTML ( std::ostream& os, const std::string& styleID, int colspan = 0, int rowspan = 0 ) const;
	static void printHeaderDelimited ( std::ostream& os );
	void printDelimited ( std::ostream& os ) const;
	static void setReportUnmatched ( bool f )	{ reportUnmatched = f; }
	static void setReportNumPks ( bool f )		{ reportNumPks = f; }
	static void setReportRank ( bool f )		{ reportRank = f; }
	static void setReportScore ( bool f )		{ reportScore = f; }
	static void setReportScoreDiff ( bool f )	{ reportScoreDiff = f; }
	static void setReportExpectation ( bool f )	{ reportExpectation = f; }
	static void setReportPValue ( bool f )		{ reportPValue = f; }
	static void setReportMValue ( bool f )		{ reportMValue = f; }
	static void setReportNumPrecursor ( bool f ){ reportNumPrecursor = f; }
	static void setReportGradient ( bool f )	{ reportGradient = f; }
	static void setReportOffset ( bool f )		{ reportOffset = f; }
	static void setReportDiscScore ( bool f )	{ reportDiscScore = f; }
	static void setReportRepeats ( bool f )		{ reportRepeats = f; }

	static bool getReportRepeats () { return reportRepeats; }
};

typedef std::vector <const PPPeptideHitInfo*> ConstPPPeptideHitInfoPtrVector;
typedef ConstPPPeptideHitInfoPtrVector::size_type ConstPPPeptideHitInfoPtrVectorSizeType;

class SearchResultsPeptideHit {
protected:
	PPPeptideHitInfo peptideHitInfo;
	PeptidePosition peptidePosition;
public:
	SearchResultsPeptideHit ( const SpecID* spID, const MSMSSpectrumInfo* mmsi, int searchIndex );
	SearchResultsPeptideHit ( const SpecID* specID, const MSMSSpectrumInfo* mmsi, const PeptideSpectralInfo* psi, double error, const HitPeptide* hitPeptide, int startAA, const std::string& accNum, int searchIndex );
	void printHTML ( std::ostream& os, const ProteinInfo& proteinInfo, int searchNumber, const std::string& styleID, const SCMSTagLink& smtl, const std::string& id, bool joint, bool showTimes = true ) const;
	void printHTMLPeak1 ( std::ostream& os, const std::string& styleID ) const;
	void printHTMLPeak2 ( std::ostream& os, int searchNumber, const std::string& styleID, const SCMSTagLink& smtl, bool joint ) const;
	void printDelimited ( std::ostream& os, const ProteinInfo& proteinInfo, int searchNumber ) const;
	static void printDelimitedEmpty ( std::ostream& os, int searchNumber );
	const PPPeptideHitInfo* getPeptideHitInfo () const { return &peptideHitInfo; }
	double getScore () const { return peptideHitInfo.getScore (); }
	double getScoreDifference () const { return peptideHitInfo.getScoreDifference (); }
	void setDiscriminantScore ( double ds ) { peptideHitInfo.setDiscriminantScore ( ds ); }
	double getExpectationValue () const { return peptideHitInfo.getExpectationValue (); }
	double getPValue () const { return peptideHitInfo.getPValue (); }
	double getMascotScore () const { return peptideHitInfo.getMascotScore (); }
	int getRank () const { return peptideHitInfo.getRank (); }
	int getNumPeaks () const { return peptideHitInfo.getNumPeaks (); }
	int getMatched () const { return peptideHitInfo.getMatched (); }
	void setRepeats ( int repeats ) { peptideHitInfo.setRepeats ( repeats ); }
	double getDiscriminantScore () const { return peptideHitInfo.getDiscriminantScore (); }
	bool notEmpty () const { return peptideHitInfo.notEmpty (); }
	int getRepeats () const { return peptideHitInfo.getRepeats (); }

	const PeptidePosition* getPeptidePosition () const { return &peptidePosition; };
	std::string getHitSequence () const { return peptidePosition.getSequence (); }
	std::string getSequencePlusMods () const { return peptidePosition.getSequencePlusMods (); }
	std::string getPrintedSequence () const { return peptidePosition.getPrintedSequence (); }
	bool checkComposition () const { return peptidePosition.checkComposition () ; }
	std::string getPeptide () const { return peptidePosition.getPeptide (); }
	std::string getSequence () const { return peptidePosition.getSequence (); }
	void getModMassesAndIndicies ( VectorPairIntDouble& vpid ) const { peptidePosition.getModMassesAndIndicies ( vpid ); }
	void getModMassesIndiciesAndString ( VectorPairIntPairStringDouble& vpid ) const { peptidePosition.getModMassesIndiciesAndString ( vpid ); }
	double getModNTermMass () const { return peptidePosition.getModNTermMass (); }
	double getModCTermMass () const { return peptidePosition.getModCTermMass (); }
	std::string getModNTerm () const { return peptidePosition.getNTerm (); }
	std::string getModCTerm () const { return peptidePosition.getCTerm (); }
	std::string getDBPeptide () const { return peptidePosition.getDBPeptide (); }
	std::string getBlibPeptide () const { return peptidePosition.getBlibPeptide (); }
	void getBlibMods ( VectorPairIntDouble& vpid ) const { peptidePosition.getBlibMods ( vpid ); }
	std::string getPrevAA ( const ProteinInfo& proteinInfo ) const { return peptidePosition.getPrevAA ( proteinInfo ); }
	std::string getNextAA ( const ProteinInfo& proteinInfo ) const { return peptidePosition.getNextAA ( proteinInfo ); }
	int getNumTolTerm ( const ProteinInfo& proteinInfo ) const { return peptidePosition.getNumTolTerm ( proteinInfo ); }
	std::string getMissedCleavages () const { return peptidePosition.getMissedCleavages (); }
	void setSubsequentOccurence () const { peptidePosition.setSubsequentOccurence (); }
	std::string getAccessionNumber () const { return peptidePosition.getAccessionNumber (); }
	std::string getFullSpecID () const { return peptidePosition.getFullSpecID (); }
	bool hasSpecID () const { return peptidePosition.hasSpecID (); }
	bool hasHitPeptide () const { return peptidePosition.hasHitPeptide (); }
	const SpecID getSpecIDasID () const { return peptidePosition.getSpecIDasID (); }
	double getM () const { return peptidePosition.getM (); }
	double getMOverZ () const { return peptidePosition.getMOverZ (); }
	double getMOverZCalc ( int searchNumber ) const { return peptidePosition.getMOverZCalc ( searchNumber ); }
	double getMPlusH () const { return peptidePosition.getMPlusH (); }
	double getMCalc ( int searchNumber ) const { return peptidePosition.getMCalc ( searchNumber ); }
	double getDaError ( int searchNumber ) const { return peptidePosition.getDaError ( searchNumber ); }
	int getCharge () const { return peptidePosition.getCharge (); }
	double getPeptideError () const { return peptidePosition.getPeptideError (); }
	double getError () const { return peptidePosition.getError (); }
	double getIntensity () const { return peptidePosition.getIntensity (); }
	double getMModValue () const { return peptidePosition.getMModValue (); }
	short getSearchIndex () const { return peptidePosition.getSearchIndex (); }
	int getStartAA () const { return peptidePosition.getStartAA (); }
	int getEndAA () const { return peptidePosition.getEndAA (); }
	double getSpotAsNumber () const { return peptidePosition.getSpotAsNumber (); }
	std::string getFractionName () const { return peptidePosition.getFractionName (); }
	std::string getSpot () const { return peptidePosition.getSpot (); }
	std::string getSpecID () const { return peptidePosition.getSpecID (); }
	int getFraction () const { return peptidePosition.getFraction (); }
	std::string getMSMSInfo () const { return peptidePosition.getMSMSInfo (); }
	bool getFirstOccurence () const { return peptidePosition.getFirstOccurence (); }
};
typedef std::vector <SearchResultsPeptideHit*> SearchResultsPeptideHitPtrVector;
typedef SearchResultsPeptideHitPtrVector::const_iterator SearchResultsPeptideHitPtrVectorConstIterator;
typedef SearchResultsPeptideHitPtrVector::size_type SearchResultsPeptideHitPtrVectorSizeType;
typedef std::pair <SearchResultsPeptideHitPtrVectorConstIterator, SearchResultsPeptideHitPtrVectorConstIterator> PairSearchResultsPeptideHitPtrVectorConstIterator;

typedef std::vector <const SearchResultsPeptideHit*> SearchResultsPeptideHitConstPtrVector;
typedef SearchResultsPeptideHitConstPtrVector::const_iterator SearchResultsPeptideHitConstPtrVectorConstIterator;
typedef SearchResultsPeptideHitConstPtrVector::size_type SearchResultsPeptideHitConstPtrVectorSizeType;

class SortPeptidePositionAscendingMerged {
	public:
		bool operator () ( const SearchResultsPeptideHit* lhs, const SearchResultsPeptideHit* rhs ) const
		{
			if ( sresKeepTimeReplicates ) {
				if ( lhs->getSpecIDasID () == rhs->getSpecIDasID () ) {
					if ( lhs->getHitSequence () == rhs->getHitSequence () )
						return lhs->getFirstOccurence () > rhs->getFirstOccurence ();
					else
						return lhs->getHitSequence () < rhs->getHitSequence ();
				}
				else
					return lhs->getSpecIDasID () < rhs->getSpecIDasID ();
			}
			else {
				if ( lhs->getSpecIDasID () == rhs->getSpecIDasID () ) {
					if ( lhs->getFirstOccurence () == rhs->getFirstOccurence () )
						return lhs->getDiscriminantScore () > rhs->getDiscriminantScore ();
					else
						return lhs->getFirstOccurence () > rhs->getFirstOccurence ();
				}
				else
					return lhs->getSpecIDasID () < rhs->getSpecIDasID ();
			}
		}
};
class SortPeptidePositionAscendingSeparated {
	public:
		bool operator () ( const SearchResultsPeptideHit* lhs, const SearchResultsPeptideHit* rhs ) const
		{
			if ( lhs->getSpecIDasID () == rhs->getSpecIDasID () ) {
				if ( lhs->getSearchIndex () == rhs->getSearchIndex () ) {
					if ( lhs->getSequencePlusMods () == rhs->getSequencePlusMods () )
						return lhs->getFirstOccurence () > rhs->getFirstOccurence ();
					else
						return lhs->getSequencePlusMods () < rhs->getSequencePlusMods ();
				}
				else 
					return lhs->getSearchIndex () < rhs->getSearchIndex ();
			}
			else
				return lhs->getSpecIDasID () < rhs->getSpecIDasID ();
		}
};
class SortPeptidePositionAscending2 {	// This is to decide what gets stored as a separate entry
	public:
		bool operator () ( const SearchResultsPeptideHit* lhs, const SearchResultsPeptideHit* rhs ) const
		{
			if ( sresKeepReplicates ) {
				if ( lhs->getAccessionNumber () == rhs->getAccessionNumber () ) {
					if ( lhs->getHitSequence () == rhs->getHitSequence () ) {
						return lhs->getSpecIDasID () < rhs->getSpecIDasID ();
					}
					else return lhs->getHitSequence () < rhs->getHitSequence ();
				}
				else return lhs->getAccessionNumber () < rhs->getAccessionNumber ();
			}
			else if ( sresKeepCharges ) {
				if ( lhs->getAccessionNumber () == rhs->getAccessionNumber () ) {
					if ( lhs->getHitSequence () == rhs->getHitSequence () ) {
						return lhs->getCharge () < rhs->getCharge ();
					}
					else return lhs->getHitSequence () < rhs->getHitSequence ();
				}
				else return lhs->getAccessionNumber () < rhs->getAccessionNumber ();
			}
			else {
				if ( lhs->getAccessionNumber () == rhs->getAccessionNumber () ) {
					return lhs->getHitSequence () < rhs->getHitSequence ();
				}
				else return lhs->getAccessionNumber () < rhs->getAccessionNumber ();
			}
		}
};
class SortPeptidePositionAscending3 {	// As basis for getting raw data
	public:
		bool operator () ( const PeptidePosition* lhs, const PeptidePosition* rhs ) const
		{
			return lhs->getSpecIDasID () < rhs->getSpecIDasID ();
		}
};

typedef std::set <SearchResultsPeptideHit*, SortPeptidePositionAscending2> SetSearchResultsPeptideHit;
typedef SetSearchResultsPeptideHit::const_iterator SetSearchResultsPeptideHitConstIterator;
typedef std::pair <SetSearchResultsPeptideHitConstIterator, SetSearchResultsPeptideHitConstIterator> PairSearchResultsPeptideHitCIter;
typedef std::map <std::string, PairSearchResultsPeptideHitCIter> MapStringToPairSearchResultsPeptideHitCIter;
typedef std::map <std::string, SearchResultsPeptideHitPtrVector> MapAccNumSearchResultsPeptideHitPtrVector;
typedef MapAccNumSearchResultsPeptideHitPtrVector::iterator MapAccNumSearchResultsPeptideHitPtrVectorIterator;
typedef MapAccNumSearchResultsPeptideHitPtrVector::const_iterator MapAccNumSearchResultsPeptideHitPtrVectorConstIterator;
typedef std::map <std::string, MapAccNumSearchResultsPeptideHitPtrVector> MapIDMapAccNumSearchResultsPeptideHitPtrVector;

class SearchResultsProteinHit {
protected:
	const SearchResultsProteinInfo* proteinHitInfo;
	std::string accessionNumber;
public:
	SearchResultsProteinHit ( const SearchResultsProteinInfo* proteinHitInfo, const std::string& accessionNumber );
	virtual ~SearchResultsProteinHit ();
	std::string getAccessionNumber () const { return accessionNumber; }
	const SearchResultsProteinInfo* getProteinHitInfo () const { return proteinHitInfo; }
	virtual int getNumMSMSHits () const = 0;
	virtual SearchResultsPeptideHit* getMSMSHit ( int i ) const = 0;
	double getBestExpectationValue () const { return proteinHitInfo->getBestExpectationValue (); }
	double getBestDiscriminantScore () const { return proteinHitInfo->getBestDiscriminantScore (); }
	double getTotalDiscriminantScore () const { return proteinHitInfo->getTotalDiscriminantScore (); }
	double getScore () const { return proteinHitInfo->getScore (); }
};
typedef std::vector <SearchResultsProteinHit*> SearchResultsProteinHitPtrVector;
typedef SearchResultsProteinHitPtrVector::size_type SearchResultsProteinHitPtrVectorSizeType;

typedef std::map <std::string, const SearchResultsProteinInfo*> MapAccessionNumberToSearchResultsProteinInfo;
typedef MapAccessionNumberToSearchResultsProteinInfo::const_iterator MapAccessionNumberToSearchResultsProteinInfoConstIterator;

typedef std::pair <const PeptidePosition*, const PPPeptideHitInfo*> PairPeptidePositionPPPeptideHitInfo;

typedef std::map <std::string, SearchResultsProteinHitPtrVector> MapIDSearchResultsProteinHitPtrVector;
typedef std::map <std::string, MapAccessionNumberToSearchResultsProteinInfo> MapIDMapAccessionNumberToSearchResultsProteinInfo;
typedef std::map <std::string, SetSearchResultsPeptideHit> MapIDSetSearchResultsPeptideHit;
typedef std::map <std::string, MapStringToPairSearchResultsPeptideHitCIter> MapIDMapStringToPairSearchResultsPeptideHitCIter;

class DatabaseResults {
	StringVector database;
	IntVector numDatabaseEntries;
	IntVector numSearchedEntries;
public:
	DatabaseResults ( const StringVector& database, std::istream& istr );
	bool isConcat () const;
	void printHTML ( std::ostream& os ) const;
};

class PPProteinHitQuanInfo {
	DoubleVector medianAreaRatio;
	DoubleVector lowQAreaRatio;
	DoubleVector highQAreaRatio;
	DoubleVector meanAreaRatio;
	DoubleVector stDevAreaRatio;
	IntVector numAreaRatio;

	DoubleVector medianIntensityRatio;
	DoubleVector lowQIntensityRatio;
	DoubleVector highQIntensityRatio;
	DoubleVector meanIntensityRatio;
	DoubleVector stDevIntensityRatio;
	IntVector numIntensityRatio;

	DoubleVectorSizeType numQuanPeaks;

	static bool reportArea;
	static bool reportIntensity;
	static bool reportMedian;
	static bool reportIQR;
	static bool reportMean;
	static bool reportStDev;
	static bool reportNum;

	static double numStdDev;
	static std::string stdevMinusHTMLHeader;
	static std::string stdevPlusHTMLHeader;
	static std::string stdevMinusDelimHeader;
	static std::string stdevPlusDelimHeader;
public:
	void initQuan ();
	std::string getRatioStr ( int i ) const;
	void setAreaRatios ( const DoubleVector& med, const DoubleVector& q1, const DoubleVector& q2, const DoubleVector& mean, const DoubleVector& stDev, const IntVector& num )
	{
		medianAreaRatio = med;
		lowQAreaRatio = q1;
		highQAreaRatio = q2;
		meanAreaRatio = mean;
		stDevAreaRatio = stDev;
		numAreaRatio = num;
	}
	void setIntensityRatios ( const DoubleVector& med, const DoubleVector& q1, const DoubleVector& q2, const DoubleVector& mean, const DoubleVector& stDev, const IntVector& num )
	{
		medianIntensityRatio = med;
		lowQIntensityRatio = q1;
		highQIntensityRatio = q2;
		meanIntensityRatio = mean;
		stDevIntensityRatio = stDev;
		numIntensityRatio = num;
	}
	int getColspan () const;
	void printHeaderHTML ( std::ostream& os, const std::string& styleID ) const;
	void printHTMLEmpty ( std::ostream& os, const std::string& styleID ) const;
	void printHTML ( std::ostream& os, const std::string& styleID ) const;
	void printHeaderDelimited ( std::ostream& os ) const;
	void printDelimitedEmpty ( std::ostream& os ) const;
	void printDelimited ( std::ostream& os ) const;
	static bool getQuan ()
	{
		return ( reportArea || reportIntensity ) && ( reportMedian || reportIQR || reportMean || reportStDev || reportNum );
	}
	static void setQuanParams ( const ParameterList* params );
};

class PPProteinHitInfo : public SearchResultsProteinInfo {
	double totalDiscriminantScore;
	double bestDiscriminantScore;
	double bestExpectationValue;
	CoverageMap coverageMap;
	std::string accessionNumber;

	PPProteinHitQuanInfo ppphqi;
	std::string getRatioStr ( int i ) const;
public:
	PPProteinHitInfo ();
	PPProteinHitInfo ( const std::string& accessionNumber );
	PPProteinHitInfo ( const StringVector& str );
	int getNumUnique () const { return numUnique; }
	int getPeptideCount () const { return peptideCount; }
	void setTotalDiscriminantScore ( double tds ) { totalDiscriminantScore = tds; }
	void setBestDiscriminantScore ( double bs ) { bestDiscriminantScore = bs; }
	void setBestExpectationValue ( double be ) { bestExpectationValue = be; }
	CoverageMap getCoverageMap () const { return coverageMap; }
	double getTotalDiscriminantScore () const { return totalDiscriminantScore; }
	double getBestDiscriminantScore () const { return bestDiscriminantScore; }
	double getBestExpectationValue () const { return bestExpectationValue; }
	std::string getAccessionNumber () const { return accessionNumber; }
	void setCoverageMap ( const CoverageMap& c ) { coverageMap = c; }
	void setAreaRatios ( const DoubleVector& med, const DoubleVector& q1, const DoubleVector& q2, const DoubleVector& mean, const DoubleVector& stDev, const IntVector& num )
	{
		ppphqi.setAreaRatios ( med, q1, q2, mean, stDev, num );
	}
	void setIntensityRatios ( const DoubleVector& med, const DoubleVector& q1, const DoubleVector& q2, const DoubleVector& mean, const DoubleVector& stDev, const IntVector& num )
	{
		ppphqi.setIntensityRatios ( med, q1, q2, mean, stDev, num );
	}
	int getColspan () const;
	void printHeaderHTML ( std::ostream& os, int index, const std::string& styleID ) const;
	void printHTML ( std::ostream& os, bool empty, const std::string& styleID ) const;
	void printHeaderDelimited ( std::ostream& os, int index ) const;
	void printDelimited ( std::ostream& os ) const;
};

typedef std::map <std::string, std::pair<std::pair<double,double>, std::string> > MapSpecIDBestDiscriminantScore;
typedef MapSpecIDBestDiscriminantScore::const_iterator MapSpecIDBestDiscriminantScoreConstIterator;

typedef std::map <std::pair<std::string, std::string>, std::pair<double, int> > MapSpecIDAndPeptideDiscriminantScore;
typedef MapSpecIDAndPeptideDiscriminantScore::iterator MapSpecIDAndPeptideDiscriminantScoreIterator;
typedef MapSpecIDAndPeptideDiscriminantScore::const_iterator MapSpecIDAndPeptideDiscriminantScoreConstIterator;

class PPXMLProteinHit : public SearchResultsProteinHit {
	std::vector <SearchResultsPeptideHit*>& tagHits;
	PPProteinHitInfo* getPPProteinHitInfo ()
	{
		SearchResultsProteinInfo* srpi = const_cast <SearchResultsProteinInfo*> (proteinHitInfo);
		return static_cast <PPProteinHitInfo*> (srpi);
	};
public:
	PPXMLProteinHit ( const std::string& aNum, std::vector <SearchResultsPeptideHit*>& th );
	SearchResultsPeptideHit* getMSMSHit ( int i ) const { return tagHits [i]; }
	int getNumMSMSHits () const { return tagHits.size (); }
	void calculateDiscriminantScores ( MapSpecIDAndPeptideDiscriminantScore& dsMap, MapSpecIDBestDiscriminantScore& bestScores, bool setRepeats, const DiscriminantScore& discScore, double maxPeptideEValue );
	void filterPeptides ( bool bestDiscrimOnly, bool noExpectation, double minPeptideScore, double minBestDiscScore, double maxPeptideEValue, MapSpecIDBestDiscriminantScore& bestScores );
	void filterPeptidesIntersection ( const std::string& id, int numMatches );
	void removePeptides ( const VectorPairStringString& removePeps );
	void calculateStats ();
	void setProteinStats ();
};
typedef std::vector <PPXMLProteinHit*> PPXMLProteinHitPtrVector;
typedef PPXMLProteinHitPtrVector::size_type PPXMLProteinHitPtrVectorSizeType;

class ExpectationValueInfo {
	double gradient;
	double offset;
	int numSpectra;
public:
	ExpectationValueInfo () {}
	ExpectationValueInfo ( double gradient, double offset, int numSpectra ) :
		gradient ( gradient ), offset ( offset ), numSpectra ( numSpectra ) {}
	double getGradient () const { return gradient; }
	double getOffset () const { return offset; }
};
typedef std::map <std::string, ExpectationValueInfo> MapIDExpectationValueInfo;
typedef std::map <std::string, SearchResultsProteinHit*> MapAccNumSearchResultsProteinHitPtr;
typedef MapAccNumSearchResultsProteinHitPtr::const_iterator MapAccNumSearchResultsProteinHitPtrConstIterator;
typedef std::map <std::string, MapAccNumSearchResultsProteinHitPtr> MapIDAccNumSearchResultsProteinHitPtrConstIterator;

class SearchResults;
class SearchResultsCrosslinkProteinHit;
class SearchResultsCrosslinkPeptideHit;

typedef std::vector <SearchResults*> VectorSearchResultsPtr;
typedef std::vector <SearchResultsCrosslinkPeptideHit*>::iterator SearchResultsCrosslinkPeptideHitPtrVectorIterator;

class SResLink;

class SearchResults {
	std::string projectName;
	std::string resultsName;
	std::string fname;
	bool discScoreGraph;
	double modificationScoreThreshold;
	std::string discFilename;
	bool discFilenameExists;
	std::string searchEndTime;
	std::string searchTime;
	bool noExpectation;
	static const std::string defaultID;
	MapIDSearchResultsProteinHitPtrVector proteinHits;
	MapIDMapAccessionNumberToSearchResultsProteinInfo proteinMap;
	SearchResultsProteinInfo* emptyProteinHit;

	MapStringToInt peptideHits;
	MapIDSetSearchResultsPeptideHit peptideMap;

	MapIDSpecIDMSMSSpectrumInfo midmsi;
	bool spottingPlate;
	StringVector fractionNames;
	StringVector rawTypes;
	IntVector numMSSpectra;
	ToleranceInfo* parentToleranceInfo;
	double sysError;
	std::string sysErrorStr;
	DoubleVector offsets;
	ToleranceInfo* fragmentToleranceInfo;

	DatabaseResults* databaseResults;
	SetString idSet;

	std::vector <SearchResultsCrosslinkProteinHit*> cLinkProteinHits;
	std::vector <SearchResultsCrosslinkPeptideHit*> cLinkPeptideHits;
	LinkInfo* linkInfo;

	static VectorSearchResultsPtr vsr;
	static int numSearches;
	static bool noMergeExpectation;
	static MapIDAccNumSearchResultsProteinHitPtrConstIterator compProtHits;
	static MapIDMapAccNumSearchResultsPeptideHitPtrVector compPepHits;
	friend class CheckTagHitsIntersection;
	void setProteinMap ( const std::string& id = defaultID );
	void setPeptideHitsAndMap ( const std::string& id = defaultID );
	typedef std::map <std::string, std::vector <SearchResultsPeptideHit*> > MapAccNoAndVectorSearchResultsPeptideHit;
	typedef MapAccNoAndVectorSearchResultsPeptideHit::const_iterator MapAccNoAndVectorSearchResultsPeptideHitConstIterator;
	typedef std::map <std::string, MapAccNoAndVectorSearchResultsPeptideHit > MapIDMapAccNoAndVectorSearchResultsPeptideHit;
	typedef MapIDMapAccNoAndVectorSearchResultsPeptideHit::const_iterator MapIDMapAccNoAndVectorSearchResultsPeptideHitConstIterator;
	Histogram histogram;
	double getEval ( double score, double a, double b, int numSpectra, bool linearTailFitExpectation ) const;
	bool getEValues ( const std::string& fname, MapIDExpectationValueInfo& mievi );
	std::string getMapTagHits ( MapIDMapAccNoAndVectorSearchResultsPeptideHit& mapTagHits, const std::string& fname, bool multisample, const SetInt& idFilterSet, double minPeptideScore, double maxPeptideEValue, double xlMinLowScore, double xlMinScoreDiff, double xlMaxLowExpectation );
	void readBlock ( std::string& spec, MapIDMapAccNoAndVectorSearchResultsPeptideHit& mapTagHits, bool multisample, const SetInt& idFilterSet, bool idFilter, bool eValueFlag, MapIDExpectationValueInfo& mievi, const std::string& expName, bool linearTailFitExpectation, double minPeptideScore, double maxPeptideEValue, double xlMinLowScore, double xlMinScoreDiff, double xlMaxLowExpectation, const IntVector& iv );
	void readXLinkBlock ( std::string& spec, std::string& s, std::string::size_type& start, std::string::size_type& end, double maxScore, double a, double b, int numSpectra, bool linearTailFitExpectation, double minPeptideScore, double maxPeptideEValue, double xlMinLowScore, double xlMinScoreDiff, double xlMaxLowExpectation, bool idFilter, SpecID* spID, const std::string& specID, int numPeaks, MSMSSpectrumInfo* mmsi );
	static void processHits ( SearchResultsProteinHitPtrVector& pHits, bool noExpectation, double minBestDiscScore, MapSpecIDBestDiscriminantScore& bestScores, MapSpecIDAndPeptideDiscriminantScore& dsMap, DiscriminantScore& discScore, const SearchCompareParams& params, int fileIndex, const std::string& id );
	void createHistogram ( const MapSpecIDAndPeptideDiscriminantScore& dsMap );
	void getPeptideFromResults ( std::string& spec, std::string::size_type& start, std::string& nTerm, std::string& peptide, std::string& cTerm, std::string& neutralLoss );
	static void sortCLinkPeptideLines ( const std::string& sortType, const SearchResultsCrosslinkPeptideHitPtrVectorIterator& begin, const SearchResultsCrosslinkPeptideHitPtrVectorIterator& end );
	static std::string getFullAccNum ( const IntVector& iv, const std::string& a );
public:
	SearchResults ( const std::string& projectName, const std::string& resultName, const std::string& fname, const SearchCompareParams& params, int fileIndex, const std::string& searchEndTime, const std::string& searchTime );
	void drawHistogram ( std::ostream& os ) const;
	StringVector getAccessionNumbers ( const std::string& id = defaultID );
	static PairSearchResultsPeptideHitPtrVectorConstIterator getSearchPeptidePositionIters ( const std::string& accessionNumber, const std::string& id );
	int getNumProteins ( const std::string& id = defaultID ) { return proteinHits [id].size (); }
	int getNumPeptides ( const std::string& id = defaultID ) { return peptideHits [id]; }
	const SearchResultsProteinInfo* getSearchResultsProteinInfo ( const std::string& aNum = "", const std::string& id = "" );
	static const SearchResultsProteinInfo* getSearchResultsProteinInfo2 ( const std::string& aNum, const std::string& id );
	SearchResultsPeptideHit* getPPPeptideHitInfoPair ( SearchResultsPeptideHit* srph = 0, const std::string& id = "" );
	void printDatabaseInfoHTML ( std::ostream& os ) const { databaseResults->printHTML ( os ); }
	StringVector getIDList ();
	static const std::string getDefaultID () { return defaultID; }
	std::string getProjectName () const { return projectName; }
	std::string getResultsName () const { return resultsName; }
	std::string getFName () const { return fname; }
	MSMSSpectrumInfo* getMSMSSpectrumInfo ( const std::string& specID, const std::string& id = defaultID )
		{
			return midmsi [id][specID];
		}
	MapSpecIDMSMSSpectrumInfo::const_iterator getMSMSSpectrumInfoBegin ( const std::string& id = defaultID )
		{ return midmsi [id].begin (); }
	MapSpecIDMSMSSpectrumInfo::const_iterator getMSMSSpectrumInfoEnd ( const std::string& id = defaultID )
		{ return midmsi [id].end (); }
	int getMSMSSpectrumInfoSize ( const std::string& id = defaultID )
		{ return midmsi [id].size (); }
	std::string getSearchEndTime () const { return searchEndTime; }
	std::string getSearchTime () const { return searchTime; }
	StringVector getFractionNames () const { return fractionNames; }
	bool getSpottingPlate () const { return spottingPlate; }
	StringVector getRawTypes () const { return rawTypes; }
	int getNumFractions () const { return fractionNames.size (); }
	std::string getFractionName ( int i ) const { return fractionNames [i]; }
	bool getMSData () const { return !numMSSpectra.empty (); }
	Tolerance* getParentTolerance () const { return parentToleranceInfo->getTolerance (); }
	double getSysError () const { return sysError; }
	std::string getSysErrorStr () const { return sysErrorStr; }
	DoubleVector getOffsets () const { return offsets; }
	Tolerance* getFragmentTolerance () const { return fragmentToleranceInfo->getTolerance (); }
	bool isConcat () { return databaseResults->isConcat (); }
	static void mergeResults ( const SearchCompareParams& params );
	static StringVector getMergedAccNumbers ( const std::string& id = defaultID );
	void printCLinkProteinHitsLinks ( std::ostream& os, int idx ) const;
	void printCLinkProteinHitsHTML ( std::ostream& os, const SResLink& sresLink, int idx ) const;
	void printCLinkProteinHitsDelimited ( std::ostream& os, int idx ) const;
	const std::vector <SearchResultsCrosslinkProteinHit*>& getCLinksHits () const { return cLinkProteinHits; }
	void sortCLinkPeptideLines ( const std::string& sortType );
	void setQuanCLink ( int idx );
};
#endif /* ! __sc_search_res_h */
