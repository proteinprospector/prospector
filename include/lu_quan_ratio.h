/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_quan_ratio.h                                               *
*                                                                             *
*  Created    : January 23rd 2003                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2003-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_quan_ratio_h
#define __lu_quan_ratio_h

#include <stdexcept>
#include <nr.h>
#include <lgen_define.h>
#include <lu_formula.h>
#include <lu_parent.h>

class QuantitationInfo {
	StringVector name;
	MapStringToStringVector singQuanInfo;
	void initialise ();
	QuantitationInfo ();
public:
	static QuantitationInfo& instance ();
	StringVector getQuanInfo ( const std::string& n )
	{
		MapStringToStringVectorConstIterator cur = singQuanInfo.find ( n );
		if ( cur != singQuanInfo.end () ) return (*cur).second;
		else {
			std::string err ( "The entry " );
			err += n;
			err += " does not occur in the parameter file quan.txt.";
			throw std::runtime_error ( err );
		}
	}
	StringVector getNames () const { return name; }
};

class LinkInfo;

class QuanPeptide {
	std::string peptide1;
	std::string nterm1;
	std::string cterm1;
	std::string nloss;
	std::string peptide2;
	std::string nterm2;
	std::string cterm2;
	std::string linkName;
	static LinkInfo* linkInfo;
public:
	QuanPeptide ( const std::string& peptide1, const std::string& nterm1, const std::string& cterm1, const std::string& nloss, const std::string& peptide2 = "", const std::string& nterm2 = "", const std::string& cterm2 = "", const std::string& linkName = "" );
	std::string getPeptide1 () const { return peptide1; }
	std::string getNTerm1 () const { return nterm1; }
	std::string getCTerm1 () const { return cterm1; }
	std::string getNLoss () const { return nloss; }
	std::string getPeptide2 () const { return peptide2; }
	std::string getNTerm2 () const { return nterm2; }
	std::string getCTerm2 () const { return cterm2; }
	std::string getLinkName () const { return linkName; }
	bool isXLink () const { return !peptide2.empty (); }
	std::string getBridgeFormula () const;
	static void setLinkInfo ( LinkInfo* li );
};

class QuanModInfo {
	std::string name;
	std::string name1;
	std::string name2;
	ElementalFormula diffFormula;
	double monoMass;
	static std::string convertPSIlabelToFormula ( const std::string& modStr );
public:
	QuanModInfo ( const std::string& aa, const std::string& name );
	int getNumMods ( const QuanPeptide& qPeptide );
	ElementalFormula getDiffFormula ( int numMods ) const;
	double getMonoMass () const { return monoMass; }
	double getMonoMass ( int numMods ) const { return numMods * monoMass; }
};

typedef std::vector <QuanModInfo> VectorQuanModInfo;
typedef VectorQuanModInfo::size_type VectorQuanModInfoSizeType;
typedef std::map <std::string, VectorQuanModInfo> QuanModInfoMap;
typedef QuanModInfoMap::iterator QuanModInfoMapIterator;
typedef QuanModInfoMap::const_iterator QuanModInfoMapConstIterator;

class QuantitationRatio : public QuantitationData {
	int charge;
	int numPeaks;

	bool quanPep;

	IntVectorVector ok;				// The previous BoolDeque here caused a memory leak
	DoubleVectorVector mOverZ;
	DoubleVectorVector intensity;
	DoubleVectorVector snr;
	DoubleVectorVector fwhm;
	DoubleVectorVector resolution;
	DoubleVectorVector area;

	StringVector formulaStrings;

	static QuanModInfoMap quanMods;
	static StringVector mod;
	static bool o18Flag;
	static bool n15Flag;
	static int numQuanStates;
	static int numRatios;
	void init2 ( const DoubleVectorVectorVector& coeff, ElementalFormulaVector& formulae );
	void init ( const XYData& xyData, const ElementalFormula* ef, bool efLink, double resolution, double mOZ, const QuanPeptide& qPeptide );
	void printHTMLMOverZ ( std::ostream& os, DoubleVectorVectorSizeType peakNumber, const std::string& styleID ) const;
	static int getNumPotentialMods ( const QuanPeptide& qPeptide, const std::string& aas );
	static void getModInfo ( const std::string& modLine, std::string& mName, StringVector& aas );
	static void printHTMLHeader1 ( std::ostream& os, const std::string& str, const std::string& styleID );
	static void printHTMLHeader2 ( std::ostream& os, const std::string& str, const std::string& styleID );
	static void printDelimitedHeader1 ( std::ostream& os, const std::string& str );
	static void printDelimitedHeader2 ( std::ostream& os, const std::string& str );
	void printHTMLLine1 ( std::ostream& os, const DoubleVectorVector& dvv, DoubleVectorVectorSizeType peakNumber, const std::string& styleID ) const;
	void printHTMLLine2 ( std::ostream& os, const DoubleVectorVector& dvv, DoubleVectorVectorSizeType peakNumber, const std::string& styleID ) const;
	void printHTMLLine3 ( std::ostream& os, const DoubleVectorVector& dvv, DoubleVectorVectorSizeType peakNumber, const std::string& styleID ) const;
	void printHTMLLine4 ( std::ostream& os, const DoubleVectorVector& ratio, const CharVectorVector& ratioType, DoubleVectorVectorSizeType peakNumber, const std::string& styleID ) const;
	void printDelimitedLine1 ( std::ostream& os, const DoubleVectorVector& dvv, DoubleVectorVectorSizeType peakNumber ) const;
	void printDelimitedLine2 ( std::ostream& os, const DoubleVectorVector& dvv, DoubleVectorVectorSizeType peakNumber ) const;
	void printDelimitedLine3 ( std::ostream& os, const DoubleVectorVector& dvv, DoubleVectorVectorSizeType peakNumber ) const;
	void printDelimitedLine4 ( std::ostream& os, const DoubleVectorVector& ratio, const CharVectorVector& ratioType, DoubleVectorVectorSizeType peakNumber ) const;
public:
	QuantitationRatio ( const XYData& xyData, double mOZ, int charge, const QuanPeptide& qPeptide, const ElementalFormula* ef, bool efLink, double resolution, int numPeaks );
	void printHTML ( std::ostream& os ) const;
	static void setQuanResidue ( const std::string& quanType );

	static std::string getUnmodifiedPeptide ( const std::string& peptide );
	static bool getQuanMasses ( const QuanPeptide& qPeptide, std::vector <ElementalFormula>& efDelta );
	static bool getQuanMasses ( const ElementalFormula* ef, std::vector <ElementalFormula>& efDelta );
	static bool getQuanMassesN15 ( const QuanPeptide& qPeptide, const ElementalFormula* ef, std::vector <ElementalFormula>& efDelta );
	static bool getDataRange ( const QuanPeptide& qPeptide, double mOverZ, int charge, double& startMass, double& endMass );

	static int getColspan ( DoubleVectorVectorSizeType peakNumber );
	static void printHTMLHeader ( std::ostream& os, DoubleVectorVectorSizeType peakNumber, const std::string& styleID );
	static void printDelimitedHeader ( std::ostream& os, DoubleVectorVectorSizeType peakNumber );
	void printHTMLLine ( std::ostream& os, DoubleVectorVectorSizeType peakNumber, const std::string& styleID ) const;
	void printDelimitedLine ( std::ostream& os, DoubleVectorVectorSizeType peakNumber ) const;
	static void printHTMLBlankLine ( std::ostream& os, DoubleVectorVectorSizeType peakNumber, const std::string& styleID );
	static void printDelimitedBlankLine ( std::ostream& os, DoubleVectorVectorSizeType peakNumber );
	bool outputQuanResults ( std::ostream& os, const std::string& searchName, int numRepeats, bool area ) const;
	DoubleVector getAreaRatios () const;
	DoubleVector getIntensityRatios () const;
	static int getNumRatios () { return numRatios; }
	static int getNumRatios ( const std::string& quanType );
	static std::string getRatioString ( int i );
	static int getNumPeaks ( int numPeaks );
};

#endif /* ! __lu_quan_ratio_h */
