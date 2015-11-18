/******************************************************************************
*                                                                             *
*  Program    : searchCompare                                                 *
*                                                                             *
*  Filename   : sc_sres_rep.h                                                 *
*                                                                             *
*  Created    : March 27th 2003                                               *
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

#ifndef __sc_sres_rep_h
#define __sc_sres_rep_h

#include <limits>
#include <vector>
#include <lgen_define.h>
#include <lu_ambiguity.h>
#include <lu_species.h>
#include <sc_search_res.h>
#include <sc_sres_link.h>

class SearchResultsProteinLine;

class SearchResultsProteinLine {
	StringVector idStr;
	std::string accessionNumber;
	ProteinInfo proteinInfo;
	int numUnique;
	ConstSearchResultsProteinInfoPtrVector searchResultsInfo;
	int numHomology;
	static bool reportNumber;
	static bool reportLinks;
	friend class SortProteinReportByTotalDiscriminantScore;
	static std::string styleID1;
	static std::string styleID2;
public:
	SearchResultsProteinLine ( const std::string& accessionNumber, const ProteinInfo& proteinInfo, int numUnique, const ConstSearchResultsProteinInfoPtrVector& searchResultsInfo );
	std::vector <SearchResultsProteinLine*> supSRL;
	std::string getAccessionNumber () const { return accessionNumber; }
	std::string getIDStr () const { return idStr.empty () ? "" : idStr [0]; }
	StringVector getIDStrVec () const { return idStr; }
	std::string getIDStrVecOutput () const
	{
		std::string idOut;
		for ( StringVectorSizeType i = 0 ; i < idStr.size () ; i++ ) {
			idOut += idStr [i];
			if ( i != idStr.size () - 1 ) idOut += ", ";
		}
		return idOut;
	}
	int getNumHomology () const { return numHomology; }
	int getNumUnique () const { return numUnique; }
	std::string getProteinSequence () const { return proteinInfo.getProteinSequence (); }
	std::string getSpecies () const { return proteinInfo.getSpecies (); }
	std::string getUniprotID () { return proteinInfo.getUniprotID (); }
	std::string getGeneName () { return proteinInfo.getGeneName (); }
	std::string getName () const { return proteinInfo.getName (); }
	int getLength () const { return proteinInfo.getLength (); }
	std::string getAcc () const { return proteinInfo.getAcc (); }
	std::string getDatabaseMZIdentMLRef () const { return proteinInfo.getDatabaseMZIdentMLRef (); }
	void setNumHomology ( int n ) { numHomology = n; }
	void setIDStr ( const std::string& id ) { idStr.push_back ( id ); }
	void printHTMLHeader ( std::ostream& os, const StringVector& searchNames, bool reportUniqPeps ) const;
	void printHTMLHeader2 ( std::ostream& os, const StringVector& searchNames, bool reportUniqPeps ) const;
	void printHTMLHeader3 ( std::ostream& os ) const;
	void printHTML ( std::ostream& os, const SResLink& sresLink, const std::string& id, bool reportUniqPeps ) const;
	void printDelimitedHeader ( std::ostream& os, bool ID, const StringVector& searchNames, bool reportUniqPeps ) const;
	void printDelimited ( std::ostream& os, const std::string& id, bool reportUniqPeps ) const;
	static void setReportNumber ( bool f ) { reportNumber = f; }
	static void setReportLinks ( bool f ) { reportLinks = f; }
	static bool getReportLinks () { return reportLinks; }
};
typedef std::vector <SearchResultsProteinLine*> SearchResultsProteinLinePtrVector;
typedef SearchResultsProteinLinePtrVector::size_type SearchResultsProteinLinePtrVectorSizeType;

class SearchResultsPeptideLine {
	ProteinInfo proteinInfo;
	SearchResultsPeptideHitConstPtrVector srph;
	ConstSearchResultsProteinInfoPtrVector srpi;
	static bool reportNumber;
	static bool reportLinks;
	friend inline int sortBySpecID ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b );
	friend inline int sortByDiscriminantScore ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b );
	friend inline int sortByExpectationValue ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b );
	friend inline int sortByMOverZ ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b );
	friend inline int sortByMPlusH ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b );
	friend inline int sortByError ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b );
	friend inline int sortByCharge ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b );
	friend inline int sortByIntensity ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b );
	friend inline int sortByStartAA ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b );
	friend inline int sortByEndAA ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b );
	friend inline int sortByMassMod ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b );
	friend inline int sortByPeptideScore ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b );
	friend inline int sortBySpot ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b );
	static std::string styleID1;
	static std::string styleID2;
public:
	SearchResultsPeptideLine ( const ProteinInfo& proteinInfo, const SearchResultsPeptideHitConstPtrVector& srph, const ConstSearchResultsProteinInfoPtrVector& srpi );
	SearchResultsPeptideLine ( const SearchResultsPeptideHitConstPtrVector& srph, const ConstSearchResultsProteinInfoPtrVector& srpi );
	void printProteinSequenceHTML ( std::ostream& os, bool coverage ) const;

	std::string getFullAccessionNumber () const { return proteinInfo.getFullAccessionNumber (); }
	std::string getAcc () const { return proteinInfo.getAcc (); }
	std::string getDatabaseMZIdentMLRef () const { return proteinInfo.getDatabaseMZIdentMLRef (); }
	std::string getAccessionInfo () const { return proteinInfo.getAccessionInfo (); }
	std::string getUniprotID () const { return proteinInfo.getUniprotID (); }
	std::string getGeneName () const { return proteinInfo.getGeneName (); }
	std::string getName () const { return proteinInfo.getName (); }
	double getProteinMW () const { return proteinInfo.getProteinMW (); }
	double getProteinPI () const { return proteinInfo.getProteinPI (); }

	double getModNTermMass () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( hasHitPeptide ( i ) ) return srph [i]->getModNTermMass ();
		}
		return 0.0;
	}
	double getModCTermMass () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( hasHitPeptide ( i ) ) return srph [i]->getModCTermMass ();
		}
		return 0.0;
	}
	std::string getModNTerm () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( hasHitPeptide ( i ) ) return srph [i]->getModNTerm ();
		}
		return "";
	}
	std::string getModCTerm () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( hasHitPeptide ( i ) ) return srph [i]->getModCTerm ();
		}
		return "";
	}
	void getModMassesAndIndicies ( VectorPairIntDouble& vpid ) const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( hasHitPeptide ( i ) ) srph [i]->getModMassesAndIndicies ( vpid );
		}
	}
	void getModMassesIndiciesAndString ( VectorPairIntPairStringDouble& vpid ) const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( hasHitPeptide ( i ) ) srph [i]->getModMassesIndiciesAndString ( vpid );
		}
	}
	std::string getDBPeptide () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( hasHitPeptide ( i ) ) return srph [i]->getDBPeptide ();
		}
		return "";
	}
	std::string getBlibPeptide () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( hasHitPeptide ( i ) ) return srph [i]->getBlibPeptide ();
		}
		return "";
	}
	void getBlibMods ( VectorPairIntDouble& vpid ) const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( hasHitPeptide ( i ) ) {
				srph [i]->getBlibMods ( vpid );
				return;
			}
		}
	}
	std::string getPeptide () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( hasHitPeptide ( i ) ) return srph [i]->getPeptide ();
		}
		return "";
	}
	std::string getPrintedSequence () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( hasHitPeptide ( i ) ) return srph [i]->getPrintedSequence ();
		}
		return "";
	}
	std::string getThePrevAA () const
	{
		std::string prevAA = getPrevAA ();
		return prevAA.substr ( prevAA.length () - 1 );
	}
	std::string getTheNextAA () const
	{
		std::string nextAA = getNextAA ();
		return nextAA.substr ( 0, 1 );
	}
	std::string getPrevAA () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( hasHitPeptide ( i ) ) return srph [i]->getPrevAA ( proteinInfo );
		}
		return "-";
	}
	std::string getNextAA () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( hasHitPeptide ( i ) ) return srph [i]->getNextAA ( proteinInfo );
		}
		return "-";
	}
	int getStartAA () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( hasHitPeptide ( i ) ) return srph [i]->getStartAA ();
		}
		return 0;
	}
	int getEndAA () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( hasHitPeptide ( i ) ) return srph [i]->getEndAA ();
		}
		return 0;
	}
	std::string getMissedCleavages () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( hasHitPeptide ( i ) ) return srph [i]->getMissedCleavages ();
		}
		return "-";
	}
	int getNumTolTerm () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( hasHitPeptide ( i ) ) return srph [i]->getNumTolTerm ( proteinInfo );
		}
		return -1;
	}
	std::string getSpecID () const		// Should only be use in Peptide Time reports where the times are equal
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( hasSpecID ( i ) ) return srph [i]->getSpecID ();
		}
		return "";
	}
	short getSearchIndex ( int i ) const { return srph [i]->getSearchIndex (); }
	short getSearchIndex0 () const { return srph [0]->getSearchIndex (); }
	std::string getFullSpecID () const		// Should only be use in Peptide Time reports where the times are equal
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( hasSpecID ( i ) ) return srph [i]->getFullSpecID ();
		}
		return "";
	}
	std::string getMSMSInfo () const		// Should only be use in Peptide Time reports where the times are equal
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( hasSpecID ( i ) ) return srph [i]->getMSMSInfo ();
		}
		return "";
	}
	std::string getFullSpecID ( int i ) const		// Should only be use in Peptide Time reports where the times are equal
	{
		if ( hasSpecID ( i ) ) return srph [i]->getFullSpecID ();
		return "";
	}
	int getFraction () const	// This function is only used when there is a single compared results file
	{
		return srph [0]->getFraction ();
	}
	int getCharge () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( hasSpecID ( i ) ) return srph [i]->getCharge ();
		}
		return 0;
	}
	double getM () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( hasSpecID ( i ) ) return srph [i]->getM ();
		}
		return 0.0;
	}
	double getMCalc () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( hasSpecID ( i ) ) return srph [i]->getMCalc ( i );
		}
		return 0.0;
	}
	double getMOverZ () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( hasSpecID ( i ) ) return srph [i]->getMOverZ ();
		}
		return 0.0;
	}
	double getMOverZCalc () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( hasSpecID ( i ) ) return srph [i]->getMOverZCalc ( i );
		}
		return 0.0;
	}
	double getDaError () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( hasSpecID ( i ) ) return srph [i]->getDaError ( i );
		}
		return 0.0;
	}
	double getSpotAsNumber () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( hasSpecID ( i ) ) return srph [i]->getSpotAsNumber ();
		}
		return 0;
	}
	std::string getFractionName () const	// This function is only used when there is a single compared results file
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( hasSpecID ( i ) ) return srph [i]->getFractionName ();
		}
		return "";
	}
	bool hasSpecID ( int i ) const { return srph [i]->hasSpecID (); }
	bool hasHitPeptide ( int i ) const { return srph [i]->hasHitPeptide (); }
	bool isPeptideHitInfo ( int i ) const { return srph [i]->notEmpty (); }
	void setAreaRatios ( int i, const DoubleVector& m, const DoubleVector& q1, const DoubleVector& q2, const DoubleVector& mean, const DoubleVector& stDev, const IntVector& num )
	{
		if ( i == -1 )	const_cast <SearchResultsProteinInfo*> (srpi [0])->setAreaRatios ( m, q1, q2, mean, stDev, num );
		else			const_cast <SearchResultsProteinInfo*> (srpi [i])->setAreaRatios ( m, q1, q2, mean, stDev, num );
	}
	void setIntensityRatios ( int i, const DoubleVector& m, const DoubleVector& q1, const DoubleVector& q2, const DoubleVector& mean, const DoubleVector& stDev, const IntVector& num )
	{
		if ( i == -1 )	const_cast <SearchResultsProteinInfo*> (srpi [0])->setIntensityRatios ( m, q1, q2, mean, stDev, num );
		else			const_cast <SearchResultsProteinInfo*> (srpi [i])->setIntensityRatios ( m, q1, q2, mean, stDev, num );
	}
	const PeptidePosition* getPeptidePosition ( int index ) const
	{
		const PeptidePosition* pp = srph [index]->getPeptidePosition ();
		if ( pp->hasSpecID () )	return pp;
		else					return 0;
	}
	bool isDecoyHit () const
	{
		bool flag = true;
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( isDecoyHit ( i ) == false ) flag = false;
		}
		return flag;
	}
	bool isDecoyHit ( int i ) const
	{
		return hasSpecID ( i ) ? srph [i]->getPeptidePosition ()->isDecoyHit () : false; 
	}
	double getMModValue ( int i ) const { return srph [i]->getMModValue (); }
	int getRepeats () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( srph [i]->notEmpty () ) return srph [i]->getRepeats ();
		}
	}
	int getNumPeaks () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( srph [i]->notEmpty () ) return srph [i]->getNumPeaks ();
		}
		return 0;
	}
	int getRank () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( srph [i]->notEmpty () ) return srph [i]->getRank ();
		}
		return 0;
	}
	double getExpectation () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( srph [i]->notEmpty () ) return srph [i]->getExpectationValue ();
		}
		return -1.0;
	}
	double getPValue () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( srph [i]->notEmpty () ) return srph [i]->getPValue ();
		}
		return -1.0;
	}
	double getMascotScore () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( srph [i]->notEmpty () ) return srph [i]->getMascotScore ();
		}
		return PeptideSpectralInfo::INVALID_MASCOT_SCORE;
	}
	double getDiscriminantScore () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( srph [i]->notEmpty () ) return srph [i]->getDiscriminantScore ();
		}
		return -10.0;
	}
	double getScore () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( srph [i]->notEmpty () ) return srph [i]->getScore ();
		}
		return 0.0;
	}
	double getScoreDifference () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( srph [i]->notEmpty () ) return srph [i]->getScoreDifference ();
		}
		return 0.0;
	}
	int getMatched () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( srph [i]->notEmpty () ) return srph [i]->getMatched ();
		}
	}
	int getRepeats ( int index ) const { return srph [index]->getPeptideHitInfo ()->getRepeats (); }
	double getExpectationValue ( int index ) const { return srph [index]->getPeptideHitInfo ()->getExpectationValue (); }
	double getMascotScore ( int index ) const { return srph [index]->getPeptideHitInfo ()->getMascotScore (); }
	double getPeptideScore ( int i ) const { return srph [i]->getPeptideHitInfo ()->getScore (); }
	double getDiscriminantScore ( int i ) const { return srph [i]->getDiscriminantScore (); }
	double getPeptideError () const		// This function is only used when there is a single compared results file
	{
		return srph [0]->getPeptideError ();
	}
	std::string getHitSequence () const { return srph [0]->getHitSequence (); }	// Currently only used for merged hits
	bool getCompositionFlag () const
	{
		for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {				// Find the first time in the line
			if ( hasSpecID ( i ) ) return srph [i]->checkComposition ();
		}
	}
	void printHTMLTimeHeader ( std::ostream& os, const StringVector& searchNames ) const;
	void printHTMLHeader ( std::ostream& os, const StringVector& searchNames ) const;
	void printHTMLHeader2 ( std::ostream& os, const StringVector& searchNames ) const;
	void printHTMLHeader3 ( std::ostream& os ) const;
	void printHTMLProteinHeader ( std::ostream& os, const StringVector& searchNames ) const;
	void printHTML ( std::ostream& os, const SCMSTagLink& smtl, const std::string& id ) const;
	void printHTML2 ( std::ostream& os, const SCMSTagLink& smtl, const std::string& id ) const;
	void printHTML2 ( std::ostream& os, const SCMSTagLink& smtl, const std::string& id, int i ) const;
	static void printHTMLEmpty ( std::ostream& os, int i );
	void printHTML ( std::ostream& os, int pline, bool empty, bool aNumEmpty, const SCMSTagLink& smtl ) const;
	void printHTMLProtein ( std::ostream& os, const StringVector& searchNames ) const;
	void printDelimitedHeader ( std::ostream& os, const StringVector& searchNames, bool ID, bool reportUniqPeps ) const;
	void printDelimitedHeader2 ( std::ostream& os, bool ID, bool reportUniqPeps ) const;
	void printDelimitedHeader3 ( std::ostream& os, const StringVector& searchNames ) const;
	void printDelimitedHeader4 ( std::ostream& os ) const;
	void printDelimitedHeader5 ( std::ostream& os, bool ID, bool reportUniqPeps ) const;
	void printDelimitedHeader6 ( std::ostream& os ) const;
	void printDelimitedHeader7 ( std::ostream& os ) const;
	void printDelimited ( std::ostream& os, const std::string& idStr, int numHomology, const std::string& idStr2, bool reportUniqPeps ) const;
	void printDelimited2 ( std::ostream& os, const std::string& idStr, int numHomology, const std::string& idStr2, bool reportUniqPeps ) const;
	void printDelimited3 ( std::ostream& os, const std::string& idStr, int numHomology, const std::string& idStr2, bool reportUniqPeps ) const;
	void printDelimited4 ( std::ostream& os ) const;
	void printDelimited4 ( std::ostream& os, int i ) const;
	void printDelimited5 ( std::ostream& os ) const;
	void printDelimitedEmpty4 ( std::ostream& os, int i );
	bool outputQuanResults ( std::ostream& os, const StringVector& searchNames, bool area ) const;
	DoubleVectorVector getRatios ( bool area ) const;
	static void setReportNumber ( bool f ) { reportNumber = f; }
	static void setReportLinks ( bool f ) { reportLinks = f; }
	void addSiteScores ( SiteScoresVector& siteScores, int line ) const;
};
typedef std::vector <SearchResultsPeptideLine*> SearchResultsPeptideLinePtrVector;
typedef SearchResultsPeptideLinePtrVector::size_type SearchResultsPeptideLinePtrVectorSizeType;
typedef SearchResultsPeptideLinePtrVector::iterator SearchResultsPeptideLinePtrVectorIterator;


//***************************
//***************************
//***************************
//***************************
inline int sortBySearchIndex ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b )
{
	int aIdx = a->getSearchIndex ( 0 );
	int bIdx = b->getSearchIndex ( 0 );
	if ( aIdx == bIdx ) return 0;
	return aIdx < bIdx ? -1 : 1;
}
inline int sortBySpecID ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b )
{
	SpecID minA;
	bool minASet = false;
	for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < a->srph.size () ; i++ ) {
		if ( a->hasSpecID ( i ) ) {
			if ( minASet ) minA = genMin ( minA, a->srph [i]->getSpecIDasID () ); 
			else {
				minA = a->srph [i]->getSpecIDasID ();
				minASet = true;
			}
		}
	}
	SpecID minB;
	bool minBSet = false;
	for ( SearchResultsPeptideHitConstPtrVectorSizeType j = 0 ; j < b->srph.size () ; j++ ) {
		if ( b->hasSpecID ( j ) ) {
			if ( minBSet ) minB = genMin ( minB, b->srph [j]->getSpecIDasID () );
			else {
				minB = b->srph [j]->getSpecIDasID ();
				minBSet = true;
			}
		}
	}
	if ( minA == minB ) return 0;
	return minA < minB ? -1 : 1;
}
inline int sortBySpot ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b )
{
	double minA = std::numeric_limits<double>::max();
	for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < a->srph.size () ; i++ ) {
		if ( a->hasSpecID ( i ) ) minA = genMin ( minA, a->srph [i]->getSpotAsNumber () );
	}
	double minB = std::numeric_limits<double>::max();
	for ( SearchResultsPeptideHitConstPtrVectorSizeType j = 0 ; j < b->srph.size () ; j++ ) {
		if ( b->hasSpecID ( j ) ) minB = genMin ( minB, b->srph [j]->getSpotAsNumber () );
	}
	if ( minA == minB ) return 0;
	return minA < minB ? -1 : 1;
}
inline int sortByPeptideScore ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b )
{
	double maxA = -std::numeric_limits<double>::max();
	for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < a->srph.size () ; i++ ) {
		if ( a->isPeptideHitInfo ( i ) ) maxA = genMax ( maxA, a->getPeptideScore ( i ) );
	}
	double maxB = -std::numeric_limits<double>::max();
	for ( SearchResultsPeptideHitConstPtrVectorSizeType j = 0 ; j < b->srph.size () ; j++ ) {
		if ( b->isPeptideHitInfo ( j ) ) maxB = genMax ( maxB, b->getPeptideScore ( j ) );
	}
	if ( maxA == maxB ) return 0;
	return maxA < maxB ? -1 : 1;
}
inline int sortByDiscriminantScore ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b )
{
	double maxA = -std::numeric_limits<double>::max();
	for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < a->srph.size () ; i++ ) {
		if ( a->isPeptideHitInfo ( i ) ) maxA = genMax ( maxA, a->getDiscriminantScore ( i ) );
	}
	double maxB = -std::numeric_limits<double>::max();
	for ( SearchResultsPeptideHitConstPtrVectorSizeType j = 0 ; j < b->srph.size () ; j++ ) {
		if ( b->isPeptideHitInfo ( j ) ) maxB = genMax ( maxB, b->getDiscriminantScore ( j ) );
	}
	if ( maxA == maxB ) return 0;
	return maxA < maxB ? -1 : 1;
}
inline int sortByExpectationValue ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b )
{
	double minA = std::numeric_limits<double>::max();
	for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < a->srph.size () ; i++ ) {
		if ( a->isPeptideHitInfo ( i ) ) minA = genMin ( minA, a->getExpectationValue ( i ) );
	}
	double minB = std::numeric_limits<double>::max();
	for ( SearchResultsPeptideHitConstPtrVectorSizeType j = 0 ; j < b->srph.size () ; j++ ) {
		if ( b->isPeptideHitInfo ( j ) ) minB = genMin ( minB, b->getExpectationValue ( j ) );
	}
	if ( minA == minB ) return 0;
	return minA < minB ? -1 : 1;
}
inline int sortByMOverZ ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b )
{
	double minA = std::numeric_limits<double>::max();
	for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < a->srph.size () ; i++ ) {
		if ( a->hasSpecID ( i ) ) minA = genMin ( minA, a->srph [i]->getMOverZ () );
	}
	double minB = std::numeric_limits<double>::max();
	for ( SearchResultsPeptideHitConstPtrVectorSizeType j = 0 ; j < b->srph.size () ; j++ ) {
		if ( b->hasSpecID ( j ) ) minB = genMin ( minB, b->srph [j]->getMOverZ () );
	}
	if ( minA == minB ) return 0;
	return minA < minB ? -1 : 1;
}
inline int sortByMPlusH ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b )
{
	double minA = std::numeric_limits<double>::max();
	for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < a->srph.size () ; i++ ) {
		if ( a->hasSpecID ( i ) ) minA = genMin ( minA, a->srph [i]->getMPlusH () );
	}
	double minB = std::numeric_limits<double>::max();
	for ( SearchResultsPeptideHitConstPtrVectorSizeType j = 0 ; j < b->srph.size () ; j++ ) {
		if ( b->hasSpecID ( j ) ) minB = genMin ( minB, b->srph [j]->getMPlusH () );
	}
	if ( minA == minB ) return 0;
	return minA < minB ? -1 : 1;
}
inline int sortByCharge ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b )
{
	int minA = std::numeric_limits<int>::max();
	for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < a->srph.size () ; i++ ) {
		if ( a->hasSpecID ( i ) ) minA = genMin ( minA, a->srph [i]->getCharge () );
	}
	int minB = std::numeric_limits<int>::max();
	for ( SearchResultsPeptideHitConstPtrVectorSizeType j = 0 ; j < b->srph.size () ; j++ ) {
		if ( b->hasSpecID ( j ) ) minB = genMin ( minB, b->srph [j]->getCharge () );
	}
	if ( minA == minB ) return 0;
	return minA < minB ? -1 : 1;
}
inline int sortByError ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b )
{
	double minA = std::numeric_limits<double>::max();
	for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < a->srph.size () ; i++ ) {
		if ( a->hasSpecID ( i ) ) minA = genMin ( minA, a->srph [i]->getError () );
	}
	double minB = std::numeric_limits<double>::max();
	for ( SearchResultsPeptideHitConstPtrVectorSizeType j = 0 ; j < b->srph.size () ; j++ ) {
		if ( b->hasSpecID ( j ) ) minB = genMin ( minB, b->srph [j]->getError () );
	}
	if ( minA == minB ) return 0;
	return minA < minB ? -1 : 1;
}
inline int sortByIntensity ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b )
{
	double maxA = -std::numeric_limits<double>::max();
	for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < a->srph.size () ; i++ ) {
		if ( a->hasSpecID ( i ) ) maxA = genMax ( maxA, a->srph [i]->getIntensity () );
	}
	double maxB = -std::numeric_limits<double>::max();
	for ( SearchResultsPeptideHitConstPtrVectorSizeType j = 0 ; j < b->srph.size () ; j++ ) {
		if ( b->hasSpecID ( j ) ) maxB = genMax ( maxB, b->srph [j]->getIntensity () );
	}
	if ( maxA == maxB ) return 0;
	return maxA < maxB ? -1 : 1;
}
inline int sortByStartAA ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b )
{
	int minA = std::numeric_limits<int>::max();
	for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < a->srph.size () ; i++ ) {
		if ( a->hasSpecID ( i ) ) minA = genMin ( minA, a->srph [i]->getStartAA () );
	}
	int minB = std::numeric_limits<int>::max();
	for ( SearchResultsPeptideHitConstPtrVectorSizeType j = 0 ; j < b->srph.size () ; j++ ) {
		if ( b->hasSpecID ( j ) ) minB = genMin ( minB, b->srph [j]->getStartAA () );
	}
	if ( minA == minB ) return 0;
	return minA < minB ? -1 : 1;
}
inline int sortByEndAA ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b )
{
	int minA = std::numeric_limits<int>::max();
	for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < a->srph.size () ; i++ ) {
		if ( a->hasSpecID ( i ) ) minA = genMin ( minA, a->srph [i]->getEndAA () );
	}
	int minB = std::numeric_limits<int>::max();
	for ( SearchResultsPeptideHitConstPtrVectorSizeType j = 0 ; j < b->srph.size () ; j++ ) {
		if ( b->hasSpecID ( j ) ) minB = genMin ( minB, b->srph [j]->getEndAA () );
	}
	if ( minA == minB ) return 0;
	return minA < minB ? -1 : 1;
}
inline int sortByMassMod ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b )
{
	double minA = std::numeric_limits<double>::max();
	for ( SearchResultsPeptideHitConstPtrVectorSizeType i = 0 ; i < a->srph.size () ; i++ ) {
		if ( a->hasSpecID ( i ) ) {
			double mm = a->getMModValue ( i );
			if ( mm ) minA = genMin ( minA, mm );
		}
	}
	double minB = std::numeric_limits<double>::max();
	for ( SearchResultsPeptideHitConstPtrVectorSizeType j = 0 ; j < b->srph.size () ; j++ ) {
		if ( b->hasSpecID ( j ) ) {
			double mm = b->getMModValue ( j );
			if ( mm ) minB = genMin ( minB, mm );
		}
	}
	if ( minA == minB ) return 0;
	return minA < minB ? -1 : 1;
}
class sortPeptideReportByMOverZ {
public:
	int operator () ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b )
	{
		return sortByMOverZ ( a, b ) == -1;
	}
};

class sortPeptideReportByMPlusH {
public:
	int operator () ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b )
	{
		return sortByMPlusH ( a, b ) == -1;
	}
};

class sortPeptideReportByError {
public:
	int operator () ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b )
	{
		return sortByError ( a, b ) == -1;
	}
};

class sortPeptideReportByIntensity {
public:
	int operator () ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b ) const
	{
		return sortByIntensity ( a, b ) == 1;
	}
};

class sortPeptideReportBySpot {
public:
	int operator () ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b ) const
	{
		return sortBySpecID ( a, b ) == -1;
	}
};
class sortPeptideReportByStartResidue {
public:
	int operator () ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b ) const
	{
		return sortByStartAA ( a, b ) == -1;
	}
};
class sortPeptideReportByEndResidue {
public:
	int operator () ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b ) const
	{
		return sortByEndAA ( a, b ) == -1;
	}
};
class sortPeptideReportByPeptideScore {
public:
	int operator () ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b ) const
	{
		return sortByPeptideScore ( a, b ) == 1;
	}
};
class sortPeptideReportByDiscriminantScore {
public:
	int operator () ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b ) const
	{
		return sortByDiscriminantScore ( a, b ) == 1;
	}
};
class sortPeptideReportByExpectationValue {
public:
	int operator () ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b ) const
	{
		return sortByExpectationValue ( a, b ) == -1;
	}
};
class sortPeptideReportByMassMod {
public:
	int operator () ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b ) const
	{
		return sortByMassMod ( a, b ) == -1;
	}
};
class sortPeptideReportByTime {		// Only used by Time report
public:
	int operator () ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b ) const
	{
		if ( !sresSingleProject ) {
			int f0 = sortBySearchIndex ( a, b );
			if ( f0 ) {
				return f0 == -1;
			}
		}
		int f1 = sortBySpecID ( a, b );
		if ( f1 == 0 && sresMainAndSupplementary ) {
			int f2 = sortByDiscriminantScore ( a, b );
			return f2 == 1;
		}
		else return f1 == -1;
	}
};
class sortPeptideReportByRT {
public:
	int operator () ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b ) const
	{
		int f1 = sortBySpot ( a, b );
		if ( f1 == 0 ) {
			int f2 = sortByDiscriminantScore ( a, b );
			return f2 == 1;
		}
		else return f1 == -1;
	}
};
class sortPeptideTimeReportByMOverZ {
public:
	int operator () ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b ) const
	{
		int f1 = sortByMOverZ ( a, b );
		if ( f1 == 0 ) {
			int f2 = sortByDiscriminantScore ( a, b );
			return f2 == 1;
		}
		else return f1 == -1;
	}
};
class sortPeptideTimeReportByError {
public:
	int operator () ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b ) const
	{
		int f1 = sortByError ( a, b );
		if ( f1 == 0 ) {
			int f2 = sortByDiscriminantScore ( a, b );
			return f2 == 1;
		}
		else return f1 == -1;
	}
};
class sortPeptideTimeReportByMPlusH {
public:
	int operator () ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b ) const
	{
		int f1 = sortByMPlusH ( a, b );
		if ( f1 == 0 ) {
			int f2 = sortByDiscriminantScore ( a, b );
			return f2 == 1;
		}
		else return f1 == -1;
	}
};
class sortPeptideTimeReportByChargeAndMPlusH {
public:
	int operator () ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b ) const
	{
		int f1 = sortByCharge ( a, b );
		if ( f1 == 0 ) {
			int f2 = sortByMPlusH ( a, b );
			if ( f2 == 0 ) {
				int f3 = sortByDiscriminantScore ( a, b );
				return f3 == 1;
			}
			else return f2 == -1;
		}
		else return f1 == -1;
	}
};
class sortPeptideTimeReportByIntensity {
public:
	int operator () ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b ) const
	{
		int f1 = sortByIntensity ( a, b );
		if ( f1 == 0 ) {
			int f2 = sortByDiscriminantScore ( a, b );
			return f2 == 1;
		}
		else return f1 == 1;
	}
};

class SearchResultsReport {
protected:
	std::vector <SearchResults*> searchResults;
	StringVector fullSearchNames;
	StringVector searchNames;
	std::string reportHitsType;
	std::string id;
public:
	SearchResultsReport ( const std::vector <SearchResults*>& sr, const std::string& reportHitsType, const std::string& id );
	virtual ~SearchResultsReport ();
	void printIDHTML ( std::ostream& os, const std::string& id ) const;
	void printDelimitedReport ( std::ostream& os, int i, bool last, const std::string& outputDirectory, const std::string& outputFilename, bool& delimHeaderPrinted ) const;
	void printHistogramDelimited ( std::ostream& os ) const;
	void printHistogramHTML ( std::ostream& os ) const;
	bool isRawForwarding () const;
#ifdef MYSQL_DATABASE
	virtual void printReportHeader ( std::ostream& os ) const {}
	virtual void printReportFooter ( std::ostream& os ) const {}
#endif
	virtual void printHTML ( std::ostream& os ) const = 0;
	virtual void printDelimited ( std::ostream& os ) const = 0;
	virtual bool printDelimitedHeader ( std::ostream& os ) const = 0;
	virtual void printMGF ( std::ostream& os, const std::string& outputDirectory ) const = 0;
	virtual void printPepXML ( std::ostream& os, const std::string& outputDirectory, const std::string& outputFilename ) const = 0;
	virtual void printMZIdentML ( std::ostream& os, const std::string& outputDirectory, const std::string& outputFilename ) const = 0;
	virtual void writeBiblioSpec ( std::ostream& os, const std::string& outputDirectory, const std::string& outputFilename, const std::string& id, bool norm ) const = 0;
	virtual void printCrosslinksDelimited ( std::ostream& os ) const = 0;
};

class SearchResultsProteinReport : public SearchResultsReport {
	std::string reportHomologousProteins;
	void calculateMainAndSupplementaryHits ();
	void makeProteinLines ( const std::vector <SearchResultsProteinLine*>& mainSRL );
	bool checkHomologyMatch ( std::vector <SearchResultsProteinLine*>& mainSRL, SearchResultsProteinLine* srprotl );
	int getNumHomologyMatches ( const std::string& aNum, const std::string& mainANum );
	PairStringVectorStringVector getAccessionNumberListPair () const;
	static void getReportedProteinHits ( const std::string& reportedHitsType, const StringVectorVector& aNums, StringVector& reportedHits );
	static void printBatchTagAccNumberList ( std::ostream& os, const StringVector& aNum );
protected:
	static bool reportUniqPeps;
	std::vector <SearchResultsProteinLine*> srprotl;
	std::pair <StringVector, StringVector> fullAccList;
public:
	SearchResultsProteinReport ( const std::vector <SearchResults*>& sr, bool remove, const MapStringToStringVector& aNumList, const std::string& reportHitsType, const std::string& reportHomologousProteins, const std::string& id );
	void printHTML ( std::ostream& os ) const;
	void printDelimited ( std::ostream& os ) const;
	bool printDelimitedHeader ( std::ostream& os ) const;
	void printMGF ( std::ostream& os, const std::string& outputDirectory ) const {}
	void printPepXML ( std::ostream& os, const std::string& outputDirectory, const std::string& outputFilename ) const {}
	void printMZIdentML ( std::ostream& os, const std::string& outputDirectory, const std::string& outputFilename ) const {}
	void writeBiblioSpec ( std::ostream& os, const std::string& outputDirectory, const std::string& outputFilename, const std::string& id, bool norm ) const {}
	void printCrosslinksDelimited ( std::ostream& os ) const {}
};

class ErrorHistogram;
class SearchResultsPeptideReport : public SearchResultsProteinReport {
	std::vector <SearchResultsPeptideLine*> srpepl;
	ErrorHistogram* errorHistogram;
	XYData mModData;
	bool errorFlag;
	DoubleVector mean;
	DoubleVector sdev;
	IntVector size;
	double aveError;
	double sdevError;
	void doErrorHistogram ( std::vector <SearchResultsPeptideLine*>& pepLine );
	void doMModHistogram ( std::vector <SearchResultsPeptideLine*>& pepLine );
	void addUnmatchedSpectra ();
	void formResultsLines ( const std::string& sortType, const std::string& sortType2 );
	void formResultsLines2 ( bool eraseNonUnique );
	void setQuanResultsCLink ();
	void setQuanResults ();
	void addSearchResultPeptideLines ( const std::string& prevSpecID, std::vector <const SearchResultsPeptideHit*>& visrph );
	void printLink ( std::ostream& os, const std::string& fraction, const std::string& outputPath, const std::string suffix ) const;
	void printArchiveLink ( std::ostream& os, const std::string& projectName, const std::string& outputPath ) const;
	void printMZIdentML_SequenceCollection ( std::ostream& ost ) const;
	void printMZIdentML_AnalysisCollection ( std::ostream& ost ) const;
	void printMZIdentML_AnalysisProtocolCollection ( std::ostream& ost ) const;
	void printMZIdentML_DataCollection ( std::ostream& ost ) const;
	bool createMods ( MapStringToMapIntToSiteInfoVector& msmivpii, MapIntToMapStringToSiteInfoVector& mimssiv, int num ) const;
	void printHTMLModsHeader ( std::ostream& os ) const;
	void printHTMLModsRow ( std::ostream& os, const std::string& mod, int site, const SiteInfoVector& siv, const SCMSTagLink& smtl ) const;
	void printHTMLMods ( std::ostream& os, int num, const SCMSTagLink& smtl ) const;
	void printDelimitedModsRow ( std::ostream& os, const std::string& mod, int site, const SiteInfoVector& siv, const std::string& idStr, const std::string& idStr2, int numHomology ) const;
	void printDelimitedMods ( std::ostream& os, int num, const std::string& idStr, int numHomology, const std::string& id, bool reportUniqPeps ) const;
	int getBestSLIPIndex ( const SiteInfoVector& siv ) const;
public:
	SearchResultsPeptideReport ( const std::vector <SearchResults*>& sr, bool remove, const MapStringToStringVector& aNumList, const std::string& sortType, const std::string& sortType2, const std::string& reportHitsType, const std::string& reportHomologousProteins, const std::string& id );
	SearchResultsPeptideReport ( const std::vector <SearchResults*>& sr, bool remove, const MapStringToStringVector& aNumList, const std::string& sortType, const std::string& sortType2, bool unmatchedSpectra, const std::string& reportHitsType, const std::string& reportHomologousProteins, const std::string& id, bool eraseNonUnique );
	SearchResultsPeptideReport ( const std::vector <SearchResults*>& sr, const std::string& id );
	StringVector getAccessionNumbers2 () const;
#ifdef MYSQL_DATABASE
	void printReportHeader ( std::ostream& os ) const;
	void printReportFooter ( std::ostream& os ) const;
#endif
	void printHTML ( std::ostream& os ) const;
	void printHTMLHeader ( std::ostream& os ) const;
	void printHTMLTimeTable ( std::ostream& os, const SCMSTagLink& smtl ) const;
	void printHTMLPeptideTables ( std::ostream& os, const SCMSTagLink& smtl, const SResLink& sresLink ) const;
	void printHTMLCalibration ( std::ostream& os ) const;
	void expectationDensityPlot ( std::ostream& os ) const;
	void quanPlot ( std::ostream& os, int n, bool area ) const;
	void quanProteinStats ( int n, bool area ) const;
	void printMGF ( std::ostream& os, const std::string& outputDirectory ) const;
	void printPepXML ( std::ostream& os, const std::string& outputDirectory, const std::string& outputFilename ) const;
	void printMZIdentML ( std::ostream& os, const std::string& outputDirectory, const std::string& outputFilename ) const;
	void writeBiblioSpec ( std::ostream& os, const std::string& outputDirectory, const std::string& outputFilename, const std::string& id, bool norm ) const;
	void writeBiblioSpecDB ( const std::string& actualPath, const std::string& outputName, bool norm ) const;
	void getRTRanges ( MapStringToPairIntDouble& cdFirstRT, MapStringToPairIntDouble& cdLastRT ) const;
	void printDelimited ( std::ostream& os ) const;
	bool printDelimitedHeader ( std::ostream& os ) const;
	void sortPeptideLines ( const std::string& sortType, const SearchResultsPeptideLinePtrVectorIterator& begin, const SearchResultsPeptideLinePtrVectorIterator& end );
	void sortPeptideTimesLines ( const std::string& sortType, const SearchResultsPeptideLinePtrVectorIterator& begin, const SearchResultsPeptideLinePtrVectorIterator& end );
	void printCrosslinksHTML ( std::ostream& os, const SResLink& sresLink ) const;
	void printCrosslinksDelimited ( std::ostream& os ) const;
};


#endif /* ! __sc_sres_rep_h */
