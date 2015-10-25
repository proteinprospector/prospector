/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_mass_frag.h                                                *
*                                                                             *
*  Created    : May 30th 2001                                                 *
*                                                                             *
*  Purpose    : Calculates the masses of fragment ions for a peptide          *
*               sequence.                                                     *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_mass_frag_h
#define __lu_mass_frag_h

#include <vector>
#include <nr.h>

class ParameterList;
class PeakMatchContext;
class PeakContainer;
class LinkInfo;

class LossInfo {
	double offset;
	std::string label;
	int maxNumber;
public:
	LossInfo ( double offset, const std::string& label, int maxNumber ) :
		offset ( offset ), label ( label ), maxNumber ( maxNumber ) {}
	double getOffset () const { return offset; }
	std::string getLabel () const { return label; }
	int getMaxNumber () const { return maxNumber; }
};

class IonType;
struct BasicLossInfo;
typedef std::vector <BasicLossInfo> BasicLossInfoVector;

class BiemannParameters {
	StringVector it;
	bool instrumentIonTypes;
	std::string maxLosses;
	bool multiZInternal;
	bool countPosCharges;
	IonType* immIon;
	IonType* mIon;
	std::vector <IonType*> nIonList;
	std::vector <IonType*> cIonList;
	std::vector <IonType*> internalIonList;
	std::vector <IonType*> intactLossIonList;

	bool nIonLossType;
	std::string label;
	double startOffset;
	std::vector <LossInfo> lossInfo;
	std::vector <unsigned int> masks;
	int numMasks;
	IntVector numHits;
	StringVector lab;
	bool allowIncorrectCharge;
	void getLossPermutations ( const std::string& ionLabel, bool nIonFlag, double offset, const std::string& sequence, const BasicLossInfoVector& bliv );
	void getNextLossPermutation ( int level );
	
	bool inList ( const StringVector& it, const std::string& s1 ) const;
	bool inList ( const StringVector& it, const std::string& s1, const std::string& s2 ) const;
	BiemannParameters& operator= ( BiemannParameters& rhs );
	BiemannParameters ( const BiemannParameters& rhs );
public:
	BiemannParameters ( const ParameterList* params );
	~BiemannParameters ();
	bool getInstrumentIonTypes () const { return instrumentIonTypes; }
	bool getCountPosCharges () const { return countPosCharges; }
	StringVector initIonType ( const ParameterList* params );
	std::string getIonSeries () const { return it [0]; }
	int getNIonListSize () const { return nIonList.size (); }
	int getCIonListSize () const { return cIonList.size (); }
	int getInternalIonListSize () const { return internalIonList.size (); }
	int getIntactLossIonListSize () const { return intactLossIonList.size (); }
	IonType* getImmIon () const { return immIon; }
	IonType* getMIon () const { return mIon; }
	IonType* getNIonList ( int index ) const { return nIonList [index]; }
	IonType* getCIonList ( int index ) const { return cIonList [index]; }
	IonType* getInternalIonList ( int index ) const { return internalIonList [index]; }
	IonType* getIntactLossIonList ( int index ) const { return intactLossIonList [index]; }
	bool get_immonium () const;
	bool get_m () const;
	bool get_a_h3po4 () const;
	bool get_a_soch4 () const;
	bool get_a_h2o () const;
	bool get_a_nh3 () const;
	bool get_a () const;
	bool get_b_h3po4 () const;
	bool get_b_soch4 () const;
	bool get_b_h2o () const;
	bool get_b_nh3 () const;
	bool get_b () const;
	bool get_b_plus_h2o () const;
	bool get_c_term_ladder () const;
	bool getCPlus2Da () const;
	bool getCPlus1Da () const;
	bool get_c () const;
	bool getCMinus1Da () const;
	bool get_d () const;

	bool get_v () const;
	bool get_w () const;
	bool get_x () const;
	bool get_n_term_ladder () const;
	bool get_y () const;
	bool get_y_h2o () const;
	bool get_y_nh3 () const;
	bool get_y_soch4 () const;
	bool get_y_h3po4 () const;
	bool get_Y () const;
	bool get_z () const;
	bool getZPlus1Da () const;
	bool getZPlus2Da () const;
	bool getZPlus3Da () const;

	bool get_bp2 () const;
	bool get_bp3 () const;
	bool get_bp2_h3po4 () const;
	bool get_bp2_soch4 () const;
	bool get_bp2_h2o () const;
	bool get_bp2_nh3 () const;
	bool get_cp2 () const;
	bool get_yp2 () const;
	bool get_yp3 () const;
	bool get_yp2_h3po4 () const;
	bool get_yp2_soch4 () const;
	bool get_yp2_h2o () const;
	bool get_yp2_nh3 () const;
	bool get_zp2 () const;
	bool get_zPlus1Dap2 () const;

	StringVector get_M_Losses () const;
	bool get_M () const;

	bool get_internal () const;

	bool getAllowIncorrectCharge () const { return allowIncorrectCharge; }

	void printHTML ( std::ostream& os ) const;
};

class ImmoniumIonSeries;
class RelatedIonSeries;
class IonSeries;
class InternalIonSeries;
class PIonSeries;
class TheoreticalSpectrum;
class SpectrumMatch;
typedef std::vector <IonSeries*> IonSeriesList;
typedef IonSeriesList::size_type IonSeriesListSizeType;
typedef std::vector <InternalIonSeries*> InternalIonSeriesList;
typedef InternalIonSeriesList::size_type InternalIonSeriesListSizeType;
typedef std::vector <PIonSeries*> PIonSeriesList;
typedef PIonSeriesList::size_type PIonSeriesListSizeType;

class PepInfo {
	StringVector pep;
	int pepLen;
	double nTermMass;
	std::string nTermStr;
	double cTermMass;
	std::string cTermStr;
	bool show;
	friend class SInfo;
	friend class BiemannFragments;
public:
	PepInfo ( const StringVector& pep, double nTermMass, const std::string& nTermStr, double cTermMass, const std::string& cTermStr, bool show );
	double peptideFormulaToPIon () const;
	double peptideFormulaToMolecularWeight () const;
	void printNTerm ( std::ostream& os ) const;
	void printCTerm ( std::ostream& os ) const;
	void printSequence ( std::ostream& os, int idx ) const;
	void printSequence ( std::ostream& os ) const;
	int getLength () const { return pepLen; }
};

class SInfo {
	const PepInfo* pInfo;
	int precision;
	RelatedIonSeries* relatedIonSeries;
	IonSeries* mIonSeries;
	IonSeriesList nIonSeriesList;
	IonSeriesList cIonSeriesList;
	IonSeriesList intactLossIonSeriesList;
	InternalIonSeriesList internalIonSeriesList;
	PIonSeriesList pIonList;
	TheoreticalSpectrum* theoreticalSpectrum;
	const LinkInfo* linkInfo;
	SetInt mainNTermIonIndicies;
	SetInt mainCTermIonIndicies;
	friend class BiemannFragments;
	void printInternalIons ( std::ostream& os, const SpectrumMatch* spectrumMatch, int sNum ) const;
	void printVerticalSequence ( std::ostream& os, const SpectrumMatch* spectrumMatch, int sNum ) const;
	void printHorizontalSequence ( std::ostream& os, const SpectrumMatch* spectrumMatch, int sNum ) const;
	void addXLinkIons ( int mxZ );
public:
	SInfo ( const PepInfo* pInfo, const BiemannParameters* bp, int maxReportedCharge, int precision, bool countPosCharges, const LinkInfo* linkInfo );
	~SInfo ();
	void printHTML ( std::ostream& os, const SpectrumMatch* spectrumMatch, int sNum ) const;
	void printXML ( std::ostream& os ) const;
};

typedef std::vector <SInfo*> SInfoPtrVector;
typedef SInfoPtrVector::size_type SInfoPtrVectorSizeType;

class BiemannFragments {
	SInfoPtrVector sInfo;
	const std::vector <PeakContainer*> dataPeaks;
	std::vector <SpectrumMatch*> spectrumMatch;
public:
	BiemannFragments ( const std::vector <PepInfo*> pInfo, const BiemannParameters* bp, const std::vector <PeakContainer*> dataPeaks, const IntVector& maxReportedCharge, bool calibrate = false, double calTol = 0.0, bool alternative = false, const LinkInfo* linkInfo = 0 );
	~BiemannFragments ();
	static double peptideFormulaToMolecularWeight ( const StringVector& pep, double nTermWt, double cTermWt );
	void printHTML ( std::ostream& os ) const;
	void printXML ( std::ostream& os, bool discriminating ) const;
	void printApplet ( std::ostream& os, const std::vector <XYData>& vXYData, bool plotCentroids, bool discriminating ) const;
	void printErrorChart ( std::ostream& os ) const;
	void printSpectrumMatchHTML ( std::ostream& os, bool discriminating ) const;
	void printSpectrumMatchTabDelimitedText ( std::ostream& os, bool discriminating ) const;
	bool getSpectrumMatch () const { return !spectrumMatch.empty (); }
	StringVector getPeptide ( SInfoPtrVectorSizeType idx ) const { return idx < sInfo.size () ? sInfo [idx]->pInfo->pep : StringVector (); }
};

void setMSProductPrecursorCharge ( int ch );

#endif /* ! __lu_mass_frag_h */
