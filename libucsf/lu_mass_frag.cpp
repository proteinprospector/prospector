/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_mass_frag.cpp                                              *
*                                                                             *
*  Created    : July 20th 1996                                                *
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
*  Copyright (1996-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef VIS_C
#pragma warning( disable : 4503 )	// Disable warnings about decorated name length being too long
#endif
#include <iterator>
#include <lg_string.h>
#include <lgen_math.h>
#include <lu_getfil.h>
#include <lu_delim.h>
#include <lu_inst.h>
#include <lu_get_link.h>
#include <lu_immonium.h>
#include <lu_mass_elem.h>
#include <lu_mass_frag.h>
#include <lu_html.h>
#include <lu_mass.h>
#include <lu_r_plot.h>
#include <lu_table.h>
#include <lu_mass_seq.h>
#include <lu_param_list.h>
#include <lu_app_gr.h>
#include <lu_charge.h>
using std::pair;
using std::make_pair;
using std::set;
using std::vector;
using std::string;
using std::stringstream;
using std::ostream;
using std::fill;
using std::find;
using std::sort;
using std::endl;
using std::copy;
using std::ostream_iterator;
using std::ostringstream;
using std::showpos;
using std::noshowpos;
using std::getline;
using std::stable_sort;

double calculateIonMOverZ ( double mass, int charge = 1 )
{
	int j = charge - 1;
	return ( mass + ( j * cation_wt ) - ELECTRON_REST_MASS ) / charge;
}

class B1Info {
	set <string> b1NTermSet;
	B1Info ();
public:
	static B1Info& instance ();
	bool isB1 ( const string& s );
};

typedef unsigned int CompositionMask;

class InternalIon {
	StringVector sequence;
	double mOverZ;
public:
	InternalIon ( const StringVector& sequence, double mOverZ );
	StringVector getSequence () const { return sequence; }
	double getMOverZ () const { return mOverZ; }
	friend class SortInternalIons;
};

class SortInternalIons {
public:
	int operator () ( const InternalIon& a, const InternalIon& b ) const
	{
		return ( a.mOverZ < b.mOverZ );
	}
};

class PeptideInfo {
	StringVector peptideFormula;
	int pepLen;
	double mPlusCation;
	DoubleVector aaMass;
	DoubleVector seriesMass;
	UIntVector seriesCompMask;
	IntVector seriesMaxZ;
	int maxZ;
	vector <InternalIon> internalIonList;
	void calculateAAMasses ();
	void calculateSeriesMasses ( double terminusWt );
	void calculateSeriesCompMasks ();
	static int getMaxCharge ( int maxReportedCharge, bool countPosCharges, int posChargeBearing, bool ntFlag );
	void calculateSeriesMaxZ ( int maxReportedCharge, bool countPosCharges, bool ntFlag );
	void calculateInternalFragments ( double terminusWt );
	static bool chargeReducedFragmentation;
	static const int PRECURSOR_CHARGE_NOT_SET;
	static int precursorCharge;
	static int maximumLossZ;
public:
	PeptideInfo ( const StringVector& peptideFormula, bool ntFlag, int maxCharge, double mPlusCation, double terminusWt, bool calcInternal, bool countPosCharges );
	int getMaxZ () const { return maxZ; }
	int getIntactLossMaxZ ( double offset ) const;
	int getPepLen () const { return pepLen; }
	double getMPlusCation () const { return mPlusCation; }
	double getM () const { return mPlusCation - cation_wt; }
	double getAAMass ( int index ) const { return aaMass [index]; }
	double getSeriesMass ( int index ) const { return seriesMass [index]; }
	CompositionMask getSeriesCompMask ( int index ) const { return seriesCompMask [index]; }
	char getSeriesAA ( int index ) const { return peptideFormula [index] [0]; }
	int getSeriesMaxZ ( int index ) const { return seriesMaxZ [index]; }
	int getInternalMaxZ ( int index ) const;
	double getInternalSeriesMass ( int index  ) const { return internalIonList [index].getMOverZ (); }
	StringVector getInternalSeriesSequence ( int index ) const { return internalIonList [index].getSequence (); }
	int getNumInternalIons () const { return internalIonList.size (); }
	static void setChargeReducedFragmentation ( bool f ) { chargeReducedFragmentation = f; }
	static bool getChargeReducedFragmentation () { return chargeReducedFragmentation; }
	static void setPrecursorCharge ( int ch ) { precursorCharge = ch; }
};
const int PeptideInfo::PRECURSOR_CHARGE_NOT_SET = -99999;
int PeptideInfo::precursorCharge = PRECURSOR_CHARGE_NOT_SET;
int PeptideInfo::maximumLossZ = PRECURSOR_CHARGE_NOT_SET;

class IonType {
protected:
	string label;
	string subLabel;
	CompositionMask mask;
	double offset;
	int startGap;
	int endGap;
	bool singleChargeOnly;
public:
	IonType ( const string& label, const string& subLabel, CompositionMask mask, double offset, int startGap, int endGap, bool singleChargeOnly );
	virtual ~IonType ();
	double getOffset () const { return offset; }
	string getLabel () const { return label; }
	string getSubLabel () const { return subLabel; }
	int getStartGap () const { return startGap; }
	int getEndGap () const { return endGap; }
	void modifyStartGap ( int g ) { startGap += g; }
	void modifyEndGap ( int g ) { endGap += g; }
	CompositionMask getMask () const { return mask; }
	bool getSingleChargeOnly () const { return singleChargeOnly; }
	virtual double calculateIonMass ( const PeptideInfo& peptideInfo, int index, int charge ) const;
	virtual BoolDeque getSeriesFlag ( const PeptideInfo& peptideInfo ) const;
};

typedef vector <IonType*> IonTypePtrVector;
typedef IonTypePtrVector::size_type IonTypePtrVectorSizeType;

class StandardIonType : public IonType {
public:
	StandardIonType ( const string& label, double offset, int startGap, int endGap, bool singleChargeOnly );
	virtual ~StandardIonType ();
};

class LossIonType : public IonType {
	UIntVector masks;
	int numMasks;
	IntVector numHits;
public:
	LossIonType ( const string& label, const string& subLabel, const UIntVector& masks, const IntVector& numHits, double offset, int startGap, int endGap, bool singleChargeOnly );
	virtual ~LossIonType ();
	BoolDeque getSeriesFlag ( const PeptideInfo& peptideInfo ) const;
};

class BPlusWaterIonType : public IonType {
public:
	BPlusWaterIonType ();
	virtual ~BPlusWaterIonType ();
	BoolDeque getSeriesFlag ( const PeptideInfo& peptideInfo ) const;
};

class CPlus2DaIonType : public IonType {
public:
	CPlus2DaIonType ();
	virtual ~CPlus2DaIonType ();
	BoolDeque getSeriesFlag ( const PeptideInfo& peptideInfo ) const;
};

class CPlus1DaIonType : public IonType {
public:
	CPlus1DaIonType ();
	virtual ~CPlus1DaIonType ();
	BoolDeque getSeriesFlag ( const PeptideInfo& peptideInfo ) const;
};

class CIonType : public IonType {
public:
	CIonType ();
	virtual ~CIonType ();
	BoolDeque getSeriesFlag ( const PeptideInfo& peptideInfo ) const;
};

class CMinus1DaIonType : public IonType {
public:
	CMinus1DaIonType ();
	virtual ~CMinus1DaIonType ();
	BoolDeque getSeriesFlag ( const PeptideInfo& peptideInfo ) const;
};

class ZIonType : public IonType {
public:
	ZIonType ();
	virtual ~ZIonType ();
	BoolDeque getSeriesFlag ( const PeptideInfo& peptideInfo ) const;
};

class ZPlus1DaIonType : public IonType {
public:
	ZPlus1DaIonType ();
	virtual ~ZPlus1DaIonType ();
	BoolDeque getSeriesFlag ( const PeptideInfo& peptideInfo ) const;
};

class ZPlus2DaIonType : public IonType {
public:
	ZPlus2DaIonType ();
	virtual ~ZPlus2DaIonType ();
	BoolDeque getSeriesFlag ( const PeptideInfo& peptideInfo ) const;
};

class ZPlus3DaIonType : public IonType {
public:
	ZPlus3DaIonType ();
	virtual ~ZPlus3DaIonType ();
	BoolDeque getSeriesFlag ( const PeptideInfo& peptideInfo ) const;
};

class SatelliteIonType : public IonType {
public:
	SatelliteIonType ( const string& label, CompositionMask mask, double offset );
	virtual ~SatelliteIonType ();
	virtual BoolDeque getSeriesFlag ( const PeptideInfo& peptideInfo ) const;
	virtual double calculateIonMass ( const PeptideInfo& peptideInfo, int index, int charge ) const = 0;
};

class VIonType : public SatelliteIonType {
public:
	VIonType ();
	virtual ~VIonType ();
	double calculateIonMass ( const PeptideInfo& peptideInfo, int index, int charge ) const;
};

class DWIonType : public SatelliteIonType {
	const double* substituentList;
public:
	DWIonType ( const string& label, CompositionMask mask, double offset, const double* substituentList );
	virtual ~DWIonType ();
	double calculateIonMass ( const PeptideInfo& peptideInfo, int index, int charge ) const;
};

class MIonType : public IonType {
public:
	MIonType ();
	virtual ~MIonType ();
	double calculateIonMass ( const PeptideInfo& peptideInfo, int index, int charge ) const;
};

class ImmoniumIonType : public IonType {
public:
	ImmoniumIonType ();
	virtual ~ImmoniumIonType ();
	double calculateIonMass ( const PeptideInfo& peptideInfo, int index, int charge ) const;
};

class IntactLossIonType : public IonType {
	static const int PRECURSOR_CHARGE_NOT_SET;
	static int precursorCharge;
public:
	IntactLossIonType ( const string& subLabel, CompositionMask mask, double offset );
	virtual ~IntactLossIonType ();
	double calculateIonMass ( const PeptideInfo& peptideInfo, int index, int charge ) const;
	static void setPrecursorCharge ( int ch ) { precursorCharge = ch; }
};
const int IntactLossIonType::PRECURSOR_CHARGE_NOT_SET = -99999;
int IntactLossIonType::precursorCharge = PRECURSOR_CHARGE_NOT_SET;

void setMSProductPrecursorCharge ( int ch )
{
	IntactLossIonType::setPrecursorCharge ( ch );
	PeptideInfo::setPrecursorCharge ( ch );
}

class InternalIonType : public IonType {
public:
	InternalIonType ( const string& subLabel, double offset, bool singleChargeOnly );
	virtual ~InternalIonType ();
	BoolDeque getSeriesFlag ( const PeptideInfo& peptideInfo ) const;
	double calculateIonMass ( const PeptideInfo& peptideInfo, int index, int charge ) const;
};

class InternalLossIonType : public InternalIonType {
	UIntVector masks;
	int numMasks;
	IntVector numHits;
public:
	InternalLossIonType ( const string& subLabel, const UIntVector& masks, const IntVector& numHits, double offset, bool singleChargeOnly );
	virtual ~InternalLossIonType ();
	BoolDeque getSeriesFlag ( const PeptideInfo& peptideInfo ) const;
};

class FragmentPeakLabel {
	string label;
	string subLabel;
	int number;
	int charge;
	static string getLabel ( const StringVector& lab );
public:
	FragmentPeakLabel ( const string& label, const string& subLabel, int number, int charge = 1 );
	FragmentPeakLabel ( const StringVector& label, const string& subLabel, int number, int charge = 1 );
	void printHTML ( ostream& os ) const;
	void printXML ( ostream& os ) const;
	bool isPrecursorIon () const { return label == "MH"; }
};

class FragmentPeak : public Peak {
	FragmentPeakLabel fragmentPeakLabel;
	int category;
public:
	FragmentPeak ( const FragmentPeakLabel& fragmentPeakLabel, int category, double mOverZ, int charge, double intensity );
	bool isMatch ( double mass, int charge ) const;
	void printMassHTML ( ostream& os ) const;
	void printLabelHTML ( ostream& os ) const;
	void printLabelXML ( ostream& os ) const;
	string getLabelXML () const;
	int getCategoryHTML () const;
	string getLabelHTML () const;
	friend class SortFragmentPeaks;
	bool isPrecursorIon () const { return fragmentPeakLabel.isPrecursorIon (); }
};

class SortFragmentPeaks {
public:
	int operator () ( const FragmentPeak& a, const FragmentPeak& b ) const
	{
		return ( a.mOverZ < b.mOverZ );
	}
};
typedef vector <FragmentPeak> FragmentPeakVector;
typedef FragmentPeakVector::size_type FragmentPeakVectorSizeType;

typedef vector <FragmentPeakVector> FragmentPeakVectorVector;
typedef FragmentPeakVectorVector::size_type FragmentPeakVectorVectorSizeType;

class FragmentPeakMatch {
	const Peak* dataPeak;
	FragmentPeakVectorVector hitInfo;
	int siz;
	bool multiplePeptides;
	pair <MapStringToInt, bool> getLabelMap () const;
	static string getColor ( int col );
	void printLabel ( ostream& os, const PeakMatchContext& pmc, const FragmentPeak& hi ) const;
	void printColoredLabel ( ostream& os, const PeakMatchContext& pmc, const FragmentPeak& hi, int col ) const;
	void printDelimitedLine ( ostream& os, const PeakMatchContext& pmc, bool printPk = true, const string& idx = "", const string& label = "", double err = 0.0 ) const;
	void printXMLMatchInfo ( ostream& os, const FragmentPeak& hi, double err, int errSF, double inten, int intenSF ) const;
public:
	FragmentPeakMatch ( const Peak* dataPeak, const FragmentPeakVectorVector& hitInfo, int siz );
	double getMOverZ () const { return dataPeak->getMOverZ (); }
	double getIntensity () const { return dataPeak->getIntensity (); }
	bool isPrecursorIon () const
	{
		for ( FragmentPeakVectorVectorSizeType i = 0 ; i < hitInfo.size () ; i++ ) {
			for ( FragmentPeakVectorSizeType j = 0 ; j < hitInfo [i].size () ; j++ ) {
				if ( hitInfo [i][j].isPrecursorIon () ) return true;
			}
		}
		return false;
	}
	bool match () const;
	bool match ( int sNum ) const;
	bool isMatch ( double mass, int charge, int sNum ) const;
	void calibrate ( XYData& xydata, const PeakMatchContext& pmc ) const;
	void printApplet ( LabelledCatagorizedGraphData& graphData, const PeakMatchContext& pmc, bool discriminating, double normalizationFactor = 0.0 ) const;
	void printErrorChart ( ostream& os, const PeakMatchContext& pmc, int sNum ) const;
	void printErrorChart ( ostream& os, const PeakMatchContext& pmc ) const;
	void printMOverZHTML ( ostream& os, const PeakMatchContext& pmc ) const;
	string getMOverZHTML ( const PeakMatchContext& pmc ) const;
	void printMassHTML ( ostream& os, const PeakMatchContext& pmc ) const;
	string getMassHTML ( const PeakMatchContext& pmc ) const;
	string getMatchHTML ( const PeakMatchContext& pmc, bool discriminating ) const;
	void printTabDelimitedText ( ostream& os, const PeakMatchContext& pmc, bool discriminating ) const;
	void printXML ( ostream& os, const PeakMatchContext& pmc, bool discriminating ) const;
};

typedef vector <FragmentPeakMatch> FragmentPeakMatchVector;
typedef FragmentPeakMatchVector::size_type FragmentPeakMatchVectorSizeType;

class IonSeries {
protected:
	DoubleVectorVector ions;
	const IonType* ionType;
	const PeptideInfo& peptideInfo;
	BoolDeque seriesFlag;
	BoolDeque printCol;
	bool printSeries;
	int getMaxCharge ( int index );
	void calculateMasses ( int start, int end );
	void calculateIntactLossMasses ( int n );
	void calculatePrintCol ();
	static double minMOverZ;
	static double maxMOverZ;
public:
	IonSeries ( const IonType* ionType, const PeptideInfo& peptideInfo );
	virtual void printHTMLEmptyCells ( ostream& os, bool noMaxCharge = true ) const;
	virtual void printHTMLHeader ( ostream& os, bool noMaxCharge = true ) const;
	virtual void printHTML ( ostream& os, int precision, const SpectrumMatch* sm, int sNum, int idx, bool noMaxCharge = true ) const;
	virtual void printHTML ( ostream& os, bool ascending, int precision, const SpectrumMatch* sm, int sNum ) const;
	virtual void printXML ( ostream& os, bool ascending, int precision ) const;
	static void printIonMassHTML ( ostream& os, double mOverZ, int charge, int precision, bool hit );
	static void printIonMassXML ( ostream& os, double mOverZ, int charge, int number, int precision );
	static void setMOverZRange ( double min, double max )
	{
		minMOverZ = min;
		maxMOverZ = max;
	}
	virtual void getPeaks ( FragmentPeakVector& fpc ) const;
};

class ImmoniumIonSeries : public IonSeries {
public:
	ImmoniumIonSeries ( const IonType* ionType, const PeptideInfo& peptideInfo ) :
		IonSeries ( ionType, peptideInfo ) {}
	void printHTML ( ostream& os, bool ascending, int precision, const SpectrumMatch* sm, int sNum ) const;
	void getPeaks ( FragmentPeakVector& fpc ) const;
};

class IntactLossIonSeries : public IonSeries {
public:
	IntactLossIonSeries ( const IonType* ionType, const PeptideInfo& peptideInfo ) :
		IonSeries ( ionType, peptideInfo ) {}
	void printHTML ( ostream& os, bool ascending, int precision, const SpectrumMatch* sm, int sNum ) const;
	void getPeaks ( FragmentPeakVector& fpc ) const;
};

class InternalIonSeries {
	DoubleVectorVector ions;
	const IonType* ionType;
	const PeptideInfo peptideInfo;
	BoolDeque seriesFlag;
	void calculateMasses ( int start, int end );
public:
	InternalIonSeries ( const IonType* ionType, const PeptideInfo& peptideInfo );
	StringVector getSequence ( int index ) const { return peptideInfo.getInternalSeriesSequence ( index ); }
	double getMass ( int index, int z ) const { return ions [z-1][index]; }
	int getSize () const { return peptideInfo.getNumInternalIons (); }
	void printLabel ( ostream& os ) const;
	void getPeaks ( FragmentPeakVector& fpc ) const;
	int getMaxCharge ( int index );
};

class PIonSeries {
	DoubleVector ions;
	string label;
	string subLabel;
	int category;
	int maxCharge;
public:
	PIonSeries ( const string& label, const string& subLabel, int category, double mass, int maxCharge );
	void printHTML ( ostream& os, int precision, const SpectrumMatch* sm, int sNum ) const;
	void getPeaks ( FragmentPeakVector& fpc ) const;
};

class RelatedIonSeries {
	const PeptideInfo peptideInfo;
	DoubleVectorVector ions;
public:
	RelatedIonSeries ( const PeptideInfo& peptideInfo );
	void printHTMLEmptyCells ( ostream& os ) const;
	void printHTMLHeader ( ostream& os ) const;
	void printHTML ( ostream& os, int precision, const SpectrumMatch* sm, int sNum, int idx ) const;
	void printHTML ( ostream& os, int precision, const SpectrumMatch* sm, int sNum ) const;
	void printXML ( ostream& os, int precision ) const;
	static void printIonMassXML ( ostream& os, double mOverZ, int number, int precision );
	void getPeaks ( FragmentPeakVector& fpc ) const;
};

typedef vector <IntactLossIonSeries*> IntactLossIonSeriesList;
typedef IntactLossIonSeriesList::size_type IntactLossIonSeriesListSizeType;

class SpectrumMatch;

class TheoreticalSpectrum {
	static bool xLink;
	int precision;
	FragmentPeakVector peaks;
public:
	TheoreticalSpectrum ( int precision );
	void add ( const PIonSeriesList& pIonSeriesList );
	void add ( const RelatedIonSeries* ionSeries );
	void add ( const IonSeries* ionSeries );
	void add ( const IonSeriesList& ionSeriesList );
	void add ( const InternalIonSeriesList& ionSeriesList );
	void sortPeaks ();
	void printPeakTable ( ostream& os, int numColumns, const SpectrumMatch* sm, int sNum ) const;
	void printApplet ( ostream& os, const PeakContainer& dataPeaks ) const;
	void reportIonMatches ( ostream& os, const PeakContainer& dataPeaks ) const;
	friend class SpectrumMatch;
};

class SpectrumMatch {
	BoolDeque set;
	vector <FragmentPeakMatch> matches;
	double maxIntensity;
	bool sequenceFlag;
	int numMatched;
	double percentMatchedIntensity;
	double percentMatchedIntensity2;
	IntVector numMatchedArray;
	DoubleVector percentMatchedIntensityArray;
	DoubleVector percentMatchedIntensityArray2;
	int numPeptides;
	void init ( const vector <TheoreticalSpectrum*> tsv, const PeakContainer* dataPeaks, bool allowIncorrectCharge, const PeakMatchContext& pmc );
	void calculateMatchStats ();
	void calculateMatchStatsAlternate ();
	void printStats ( ostream& os, int nm, double pmi, double pmi2 ) const;
	void initSetList ( const vector <TheoreticalSpectrum*> tsv );
public:
	SpectrumMatch ( const vector <TheoreticalSpectrum*> tsv, const PeakContainer* dataPeaks, bool allowIncorrectCharge, bool calFlag, double calTol, bool alternative );
	void calibrate ( const PeakMatchContext& pmc, double& gradient, double& offset ) const;
	bool isMatch ( double mass, int charge, int sNum ) const;
	void printApplet ( ostream& os, const PeakMatchContext& pmc, const vector <XYData>& vXYData, bool plotCentroids, bool discriminating ) const;
	void printStats ( ostream& os ) const;
	void printErrorChart ( ostream& os, const PeakMatchContext& pmc ) const;
	void printHTML ( ostream& os, const PeakMatchContext& pmc, bool discriminating ) const;
	void printTabDelimitedText ( ostream& os, const PeakMatchContext& pmc, bool discriminating ) const;
	void printXML ( ostream& os, const PeakMatchContext& pmc, bool discriminating ) const;
};

static const int MASS_CATEGORY = 0;
static const int SERIES_CATEGORY = 1;
static const int LOSS_CATEGORY = 2;
static const int IMMONIUM_CATEGORY = 3;
static const int INTERNAL_CATEGORY = 4;
static const int XLINK_CATEGORY = 5;

IonType::IonType ( const string& label, const string& subLabel, CompositionMask mask, double offset, int startGap, int endGap, bool singleChargeOnly ) :
	label ( label ), subLabel ( subLabel ),
	mask ( mask ), offset ( offset ),
	startGap ( startGap ), endGap ( endGap ),
	singleChargeOnly ( singleChargeOnly )
{
}
IonType::~IonType () {}
double IonType::calculateIonMass ( const PeptideInfo& peptideInfo, int index, int charge ) const
{
	return calculateIonMOverZ ( peptideInfo.getSeriesMass ( index ) - offset, charge );
}
StandardIonType::StandardIonType ( const string& label, double offset, int startGap, int endGap, bool singleChargeOnly ) :
	IonType ( label, "", 0, offset, startGap, endGap, singleChargeOnly )
{
}
StandardIonType::~StandardIonType () {}
BoolDeque IonType::getSeriesFlag ( const PeptideInfo& peptideInfo ) const
{
	return BoolDeque ( peptideInfo.getPepLen (), true );
}
LossIonType::LossIonType ( const string& label, const string& subLabel, const UIntVector& masks, const IntVector& numHits, double offset, int startGap, int endGap, bool singleChargeOnly ) :
	IonType ( label, subLabel, 0, offset, startGap, endGap, singleChargeOnly ),
	masks ( masks ),
	numMasks ( masks.size () ),
	numHits ( numHits )
{
}
LossIonType::~LossIonType () {}
BoolDeque LossIonType::getSeriesFlag ( const PeptideInfo& peptideInfo ) const
{
	int pepLen = peptideInfo.getPepLen ();
	BoolDeque seriesFlag ( pepLen );
	IntVector maskCount (numMasks);
	fill ( maskCount.begin (), maskCount.end (), 0 );
	for ( int i = 0 ; i < pepLen ; i++ ) {
		int num = 0;
		for ( int j = 0 ; j < numMasks ; j++ ) {
			if ( peptideInfo.getSeriesCompMask ( i ) & masks [j] ) maskCount [j] += 1;
			if ( ( maskCount [j] >= numHits [j] ) ) num++;
		}
		seriesFlag [i] = ( num == numMasks );
	}
	return seriesFlag;
}

BPlusWaterIonType::BPlusWaterIonType () :
	IonType ( "b", "+H<sub>2</sub>O", 0, b_plus_h2o_tag_offset, 1, 1, false )
{
}
BPlusWaterIonType::~BPlusWaterIonType () {}
BoolDeque BPlusWaterIonType::getSeriesFlag ( const PeptideInfo& peptideInfo ) const
{
	int pepLen = peptideInfo.getPepLen ();
	BoolDeque seriesFlag ( pepLen );
	bool posChargeBearing = false;
	for ( int i = 0 ; i < pepLen ; i++ ) {
		if ( peptideInfo.getSeriesCompMask ( i ) & instInf->getPosChargeBearingMask () ) posChargeBearing = true;
		seriesFlag [i] = ( ( i == pepLen - 1 || i == pepLen - 2 || i == pepLen - 3 ) && posChargeBearing );
	}
	return seriesFlag;
}
CPlus2DaIonType::CPlus2DaIonType () :
	IonType ( "c+2", "", 0, cPlus2DaTagOffset, 0, 1, false )
{
}
CPlus2DaIonType::~CPlus2DaIonType () {}
BoolDeque CPlus2DaIonType::getSeriesFlag ( const PeptideInfo& peptideInfo ) const
{
	BoolDeque seriesFlag;
	int pepLen = peptideInfo.getPepLen ();
	if ( pepLen ) {
		seriesFlag.resize ( pepLen );
		for ( int i = 0 ; i < pepLen - 1 ; i++ ) {
			seriesFlag [i] = ( peptideInfo.getSeriesAA ( i+1 ) != 'P' );
		}
		seriesFlag [pepLen-1] = true;
	}
	return seriesFlag;
}
CPlus1DaIonType::CPlus1DaIonType () :
	IonType ( "c+1", "", 0, cPlus1DaTagOffset, 0, 1, false )
{
}
CPlus1DaIonType::~CPlus1DaIonType () {}
BoolDeque CPlus1DaIonType::getSeriesFlag ( const PeptideInfo& peptideInfo ) const
{
	BoolDeque seriesFlag;
	int pepLen = peptideInfo.getPepLen ();
	if ( pepLen ) {
		seriesFlag.resize ( pepLen );
		for ( int i = 0 ; i < pepLen - 1 ; i++ ) {
			seriesFlag [i] = ( peptideInfo.getSeriesAA ( i+1 ) != 'P' );
		}
		seriesFlag [pepLen-1] = true;
	}
	return seriesFlag;
}
CIonType::CIonType () :
	IonType ( "c", "", 0, c_tag_offset, 0, 1, false )
{
}
CIonType::~CIonType () {}
BoolDeque CIonType::getSeriesFlag ( const PeptideInfo& peptideInfo ) const
{
	BoolDeque seriesFlag;
	int pepLen = peptideInfo.getPepLen ();
	if ( pepLen ) {
		seriesFlag.resize ( pepLen );
		for ( int i = 0 ; i < pepLen - 1 ; i++ ) {
			seriesFlag [i] = ( peptideInfo.getSeriesAA ( i+1 ) != 'P' );
		}
		seriesFlag [pepLen-1] = true;
	}
	return seriesFlag;
}
CMinus1DaIonType::CMinus1DaIonType () :
	IonType ( "c-1", "", 0, cMinus1DaTagOffset, 0, 1, false )
{
}
CMinus1DaIonType::~CMinus1DaIonType () {}
BoolDeque CMinus1DaIonType::getSeriesFlag ( const PeptideInfo& peptideInfo ) const
{
	BoolDeque seriesFlag;
	int pepLen = peptideInfo.getPepLen ();
	if ( pepLen ) {
		seriesFlag.resize ( pepLen );
		for ( int i = 0 ; i < pepLen - 1 ; i++ ) {
			seriesFlag [i] = ( peptideInfo.getSeriesAA ( i+1 ) != 'P' );
		}
		seriesFlag [pepLen-1] = true;
	}
	return seriesFlag;
}
ZIonType::ZIonType () :
	IonType ( "z", "", 0, z_tag_offset, 0, 1, false )
{
}
ZIonType::~ZIonType () {}
BoolDeque ZIonType::getSeriesFlag ( const PeptideInfo& peptideInfo ) const
{
	int pepLen = peptideInfo.getPepLen ();
	BoolDeque seriesFlag ( pepLen );
	for ( int i = 0 ; i < pepLen ; i++ ) {
		seriesFlag [i] = ( peptideInfo.getSeriesAA ( i ) != 'P' );
	}
	return seriesFlag;
}
ZPlus1DaIonType::ZPlus1DaIonType () :
	IonType ( "z+1", "", 0, zPlus1DaTagOffset, 0, 1, false )
{
}
ZPlus1DaIonType::~ZPlus1DaIonType () {}
BoolDeque ZPlus1DaIonType::getSeriesFlag ( const PeptideInfo& peptideInfo ) const
{
	int pepLen = peptideInfo.getPepLen ();
	BoolDeque seriesFlag ( pepLen );
	for ( int i = 0 ; i < pepLen ; i++ ) {
		seriesFlag [i] = ( peptideInfo.getSeriesAA ( i ) != 'P' );
	}
	return seriesFlag;
}
ZPlus2DaIonType::ZPlus2DaIonType () :
	IonType ( "z+2", "", 0, zPlus2DaTagOffset, 0, 1, false )
{
}
ZPlus2DaIonType::~ZPlus2DaIonType () {}
BoolDeque ZPlus2DaIonType::getSeriesFlag ( const PeptideInfo& peptideInfo ) const
{
	int pepLen = peptideInfo.getPepLen ();
	BoolDeque seriesFlag ( pepLen );
	for ( int i = 0 ; i < pepLen ; i++ ) {
		seriesFlag [i] = ( peptideInfo.getSeriesAA ( i ) != 'P' );
	}
	return seriesFlag;
}
ZPlus3DaIonType::ZPlus3DaIonType () :
	IonType ( "z+3", "", 0, zPlus3DaTagOffset, 0, 1, false )
{
}
ZPlus3DaIonType::~ZPlus3DaIonType () {}
BoolDeque ZPlus3DaIonType::getSeriesFlag ( const PeptideInfo& peptideInfo ) const
{
	int pepLen = peptideInfo.getPepLen ();
	BoolDeque seriesFlag ( pepLen );
	for ( int i = 0 ; i < pepLen ; i++ ) {
		seriesFlag [i] = ( peptideInfo.getSeriesAA ( i ) != 'P' );
	}
	return seriesFlag;
}
SatelliteIonType::SatelliteIonType ( const string& label, CompositionMask mask, double offset ) :
	IonType ( label, "", mask, offset, 0, 0, true )
{
}
SatelliteIonType::~SatelliteIonType () {}
BoolDeque SatelliteIonType::getSeriesFlag ( const PeptideInfo& peptideInfo ) const
{
	int pepLen = peptideInfo.getPepLen ();
	BoolDeque seriesFlag ( pepLen );
	bool posChargeBearing = false;
	for ( int i = 0 ; i < pepLen ; i++ ) {
		CompositionMask compMask = peptideInfo.getSeriesCompMask ( i );
		seriesFlag [i] = ( compMask & mask ) == false  && posChargeBearing;
		if ( compMask & instInf->getPosChargeBearingMask () ) posChargeBearing = true;
	}
	return seriesFlag;
}
DWIonType::DWIonType ( const string& label, CompositionMask mask, double offset, const double* substituentList ) :
	SatelliteIonType ( label, mask, offset ),
	substituentList ( substituentList )
{
}
DWIonType::~DWIonType () {}
double DWIonType::calculateIonMass ( const PeptideInfo& peptideInfo, int index, int charge ) const
{
	double substituent = substituentList [peptideInfo.getSeriesAA ( index )];
	if ( substituent )
		return calculateIonMOverZ ( peptideInfo.getSeriesMass ( index ) + substituent - offset );
	else
		return 0.0;
}
VIonType::VIonType () :
	SatelliteIonType ( "v", instInf->getVIonExcludeMask (), v_tag_offset )
{
}
VIonType::~VIonType () {}
double VIonType::calculateIonMass ( const PeptideInfo& peptideInfo, int index, int charge ) const
{
	return calculateIonMOverZ ( peptideInfo.getSeriesMass ( index - 1 ) - offset );
}
MIonType::MIonType () :
	IonType ( "m", "", 0, c2_h2_n_o, 0, 0, true )
{
}
MIonType::~MIonType () {}
double MIonType::calculateIonMass ( const PeptideInfo& peptideInfo, int index, int charge ) const
{
	return peptideInfo.getMPlusCation () - peptideInfo.getAAMass ( index ) + offset;
}
ImmoniumIonType::ImmoniumIonType () :
	IonType ( "Immonium", "", 0, - c_o + h1, 0, 0, true )
{
}
ImmoniumIonType::~ImmoniumIonType () {}
double ImmoniumIonType::calculateIonMass ( const PeptideInfo& peptideInfo, int index, int charge ) const
{
	return calculateIonMOverZ ( peptideInfo.getAAMass ( index ) + offset, charge );
}
IntactLossIonType::IntactLossIonType ( const string& subLabel, CompositionMask mask, double offset ) :
	IonType ( "MH", subLabel, mask, offset, -1, -1, false )
{
}
IntactLossIonType::~IntactLossIonType () {}
double IntactLossIonType::calculateIonMass ( const PeptideInfo& peptideInfo, int index, int charge ) const
{
	int j = charge - 1;
	if ( precursorCharge != PRECURSOR_CHARGE_NOT_SET && PeptideInfo::getChargeReducedFragmentation () ) {
		j = precursorCharge - 1;
	}
	return ( peptideInfo.getMPlusCation () - ( index * offset ) + ( j * cation_wt ) ) / charge;
}
InternalIonType::InternalIonType ( const string& subLabel, double offset, bool singleChargeOnly ) :
	IonType ( "Internal", subLabel, 0, offset, 0, 0, singleChargeOnly )
{
}
InternalIonType::~InternalIonType () {}
BoolDeque InternalIonType::getSeriesFlag ( const PeptideInfo& peptideInfo ) const
{
	return BoolDeque ( peptideInfo.getNumInternalIons (), true );
}
double InternalIonType::calculateIonMass ( const PeptideInfo& peptideInfo, int index, int charge ) const
{
	return calculateIonMOverZ ( peptideInfo.getInternalSeriesMass ( index ) - offset, charge );
}
InternalLossIonType::InternalLossIonType ( const string& subLabel, const UIntVector& masks, const IntVector& numHits, double offset, bool singleChargeOnly ) :
	InternalIonType ( subLabel, offset, singleChargeOnly ),
	masks ( masks ),
	numMasks ( masks.size () ),
	numHits ( numHits )
{
}
InternalLossIonType::~InternalLossIonType () {}
BoolDeque InternalLossIonType::getSeriesFlag ( const PeptideInfo& peptideInfo ) const
{
	int size = peptideInfo.getNumInternalIons ();
	BoolDeque seriesFlag ( size );
	IntVector maskCount (numMasks);
	for ( int i = 0 ; i < size ; i++ ) {
		const StringVector& sequence = peptideInfo.getInternalSeriesSequence ( i );
		int num = 0;
		fill ( maskCount.begin (), maskCount.end (), 0 );
		for ( int j = 0 ; j < numMasks ; j++ ) {
			for ( StringVectorSizeType k = 0 ; k < sequence.size () ; k++ ) {
				if ( aa_composition_mask [sequence [k][0]] & masks [j] ) maskCount [j] += 1;
				if ( maskCount [j] >= numHits [j] ) {
					num++;
					break;
				}
			}
		}
		seriesFlag [i] = ( num == numMasks );
	}
	return seriesFlag;
}
InternalIon::InternalIon ( const StringVector& sequence, double mOverZ ) :
	sequence ( sequence ),
	mOverZ ( mOverZ )
{
}
bool PeptideInfo::chargeReducedFragmentation = false;
PeptideInfo::PeptideInfo ( const StringVector& peptideFormula, bool ntFlag, int maxReportedCharge, double mPlusCation, double terminusWt, bool calcInternal, bool countPosCharges ) :
	peptideFormula ( peptideFormula ),
	pepLen ( peptideFormula.size () ),
	mPlusCation ( mPlusCation )
{
	calculateAAMasses ();
	calculateSeriesMasses ( terminusWt );
	calculateSeriesCompMasks ();
	int posChargeBearing = 0;
	if ( countPosCharges ) {
		posChargeBearing = sumBasicResidues ( peptideFormula );
	}
	maxZ = getMaxCharge ( maxReportedCharge, countPosCharges, posChargeBearing, ntFlag );
	calculateSeriesMaxZ ( maxReportedCharge, countPosCharges, ntFlag );
	if ( calcInternal ) calculateInternalFragments ( terminusWt );
}
void PeptideInfo::calculateAAMasses ()
{
	for ( int i = 0 ; i < pepLen ; i++ ) {
		if ( peptideFormula [i].length () > 1 )
			aaMass.push_back ( getAminoAcidWt ( peptideFormula [i] ) );
		else
			aaMass.push_back ( amino_acid_wt [peptideFormula [i][0]] );
	}
}
void PeptideInfo::calculateSeriesMasses ( double terminusWt )
{
	double mass = terminusWt;
	for ( int i = 0 ; i < pepLen ; i++ ) {
		if ( peptideFormula [i].length () > 1 )
			mass += getAminoAcidWt ( peptideFormula [i] );
		else
			mass += amino_acid_wt [peptideFormula [i][0]];
		seriesMass.push_back ( mass );
	}
}
void PeptideInfo::calculateSeriesCompMasks ()
{
	for ( int i = 0 ; i < pepLen ; i++ ) {
		string aa = peptideFormula [i];
		if ( !genStrcasecmp ( aa, "M(Oxidation)" ) )	seriesCompMask.push_back ( aa_composition_mask ['m'] );
		else if ( !genStrcasecmp ( aa, "S(Phospho)" ) )	seriesCompMask.push_back ( aa_composition_mask ['s'] );
		else if ( !genStrcasecmp ( aa, "T(Phospho)" ) )	seriesCompMask.push_back ( aa_composition_mask ['t'] );
		else if ( !genStrcasecmp ( aa, "Y(Phospho)" ) ) seriesCompMask.push_back ( aa_composition_mask ['y'] );
		else seriesCompMask.push_back ( aa_composition_mask [aa[0]] );
	}
}
int PeptideInfo::getIntactLossMaxZ ( double offset ) const
{
	int retZ = maxZ;
	if ( maximumLossZ != PRECURSOR_CHARGE_NOT_SET && offset <= 2 ) {
		retZ = maximumLossZ;
	}
	return retZ;
}
int PeptideInfo::getMaxCharge ( int maxReportedCharge, bool countPosCharges, int posChargeBearing, bool ntFlag )
{
	int maximumZ;
	if ( countPosCharges ) {
		if ( ntFlag ) maximumZ = posChargeBearing + 1;
		else maximumZ = chargeReducedFragmentation ? posChargeBearing : posChargeBearing + 1;
		if ( maxReportedCharge ) {
			maximumZ = genMin ( maximumZ, maxReportedCharge );
		}
	}
	else {
		if ( maxReportedCharge )
			maximumZ = maxReportedCharge;
		else
			ErrorHandler::genError ()->error ( "If you set the maximum charge to No Limit you must count the basic residues.\n" );
	}
	if ( precursorCharge != PRECURSOR_CHARGE_NOT_SET && chargeReducedFragmentation ) {
		maximumLossZ = maxReportedCharge;
		maximumZ = genMin ( maximumZ, precursorCharge - 1 );
	}
	if ( maximumZ == 0 ) maximumZ = 1;		// 1 is the minimum charge
	return maximumZ;
}
void PeptideInfo::calculateSeriesMaxZ ( int maxReportedCharge, bool countPosCharges, bool ntFlag )
{
	int posChargeBearing = 0;
	for ( int i = 0 ; i < pepLen ; i++ ) {
		if ( countPosCharges ) {
			if ( seriesCompMask [i] & instInf->getPosChargeBearingMask () ) posChargeBearing++;
			else {
				const string& f = peptideFormula [i];
				if ( f.size () > 8 && isPrefix ( f.substr ( 2 ), "Cation" ) ) posChargeBearing++;
			}
		}
		seriesMaxZ.push_back ( getMaxCharge ( maxReportedCharge, countPosCharges, posChargeBearing, ntFlag ) );
	}
}
void PeptideInfo::calculateInternalFragments ( double terminusWt )
{
	set <StringVector> uniqueIntIons;
	int maxInternalLength = pepLen - 2;
	for ( int i = 1 ; i < maxInternalLength ; i++ ) {
		double offset = internal_n_terminus_wt - seriesMass [i-1];
		for ( int j = 1 ; j <= maxInternalLength - i ; j++ ) {
			StringVector sv ( peptideFormula.begin () + i, peptideFormula.begin () + i + j + 1 );
			pair <std::set <StringVector>::iterator, bool> flag = uniqueIntIons.insert ( sv );
			if ( flag.second ) {	// Make sure same internal ion doesn't get inserted twice
				internalIonList.push_back ( InternalIon ( sv, seriesMass [i+j] + offset ) );
			}
		}
	}
	sort ( internalIonList.begin (), internalIonList.end (), SortInternalIons () );
}
int PeptideInfo::getInternalMaxZ ( int index ) const
{
	int posChargeBearing = 0;
	const StringVector& pep = internalIonList [index].getSequence ();
	for ( StringVectorSizeType i = 0 ; i < pep.size () ; i++ ) {
		if ( aa_composition_mask [pep [i][0]] & instInf->getPosChargeBearingMask () ) posChargeBearing++;
	}
	return genMin ( posChargeBearing + 1, maxZ );
}
struct BasicLossInfo {
	string html;
	double offset;
	unsigned int mask;
	BasicLossInfo ( const string& formula, double offset, unsigned int mask );
};
BasicLossInfo::BasicLossInfo ( const string& formula, double offset, unsigned int mask ) :
	html ( getElementalFormulaHTML ( formula ) ),
	offset ( offset ),
	mask ( mask )
{
}
typedef vector <BasicLossInfo> BasicLossInfoVector;

BiemannParameters::BiemannParameters ( const ParameterList* params ) :
	immIon ( 0 ),
	mIon ( 0 ),
	instrumentIonTypes ( params->getBoolValue ( "use_instrument_ion_types", false ) ),
	it ( initIonType ( params ) ),
	maxLosses ( params->getStringValue ( "max_losses", "1" ) ),
	multiZInternal ( params->getBoolValue ( "multi_z_internal", false ) ),
	countPosCharges ( params->getStringValue ( "count_pos_z", "Count Basic AA" ) == "Count Basic AA" ),
	allowIncorrectCharge ( instInf->getAllowIncorrectCharge () )
{
	StringVector s = initSequence ( params->getStrippedStringValue ( "sequence", "" ) );
	int seqLength = s.size ();
	string sequence;
	if ( maxLosses != "1" ) {
		for ( StringVectorSizeType i = 0 ; i < s.size () ; i++ ) {
			const string& str = s [i];
			if ( str.length () == 1 )							sequence += str [0];
			else if ( !genStrcasecmp ( str, "S(Phospho)" ) )	sequence += 's';
			else if ( !genStrcasecmp ( str, "T(Phospho)" ) )	sequence += 't';
			else if ( !genStrcasecmp ( str, "Y(Phospho)" ) )	sequence += 'y';
			else if ( !genStrcasecmp ( str, "M(Oxidation)" ) )	sequence += 'm';
			else sequence += str [0];
		}
	}
	string da_html			= "d<sub>a</sub>";
	string db_html			= "d<sub>b</sub>";
	string wa_html			= "w<sub>a</sub>";
	string wb_html			= "w<sub>b</sub>";

	if ( get_immonium () ) {
		immIon = new ImmoniumIonType ();
	}
	if ( get_m () )				mIon = new MIonType ();

	BasicLossInfoVector blivA;
	if ( get_a_h3po4 () )	blivA.push_back ( BasicLossInfo ( "H3PO4", h3_p_o4, phosphorylation_mask ) );
	if ( get_a_soch4 () )	blivA.push_back ( BasicLossInfo ( "SOCH4", s_o_c_h4, oxidized_m_mask ) );
	if ( get_a_h2o () )		blivA.push_back ( BasicLossInfo ( "H2O", h2_o, instInf->getWaterLossMask () ) );
	if ( get_a_nh3 () )		blivA.push_back ( BasicLossInfo ( "NH3", n_h3, instInf->getAmmoniaLossMask () ) );
	getLossPermutations ( "a", true, a_tag_offset, sequence, blivA );
	if ( get_a () )				nIonList.push_back ( new StandardIonType ( "a", a_tag_offset, 1, 1, false ) );

	BasicLossInfoVector blivB;
	if ( get_b_h3po4 () )	blivB.push_back ( BasicLossInfo ( "H3PO4", h3_p_o4, phosphorylation_mask ) );
	if ( get_b_soch4 () )	blivB.push_back ( BasicLossInfo ( "SOCH4", s_o_c_h4, oxidized_m_mask ) );
	if ( get_b_h2o () )		blivB.push_back ( BasicLossInfo ( "H2O", h2_o, instInf->getWaterLossMask () ) );
	if ( get_b_nh3 () )		blivB.push_back ( BasicLossInfo ( "NH3", n_h3, instInf->getAmmoniaLossMask () ) );
	getLossPermutations ( "b", true, b_tag_offset, sequence, blivB );
	if ( get_b () )				nIonList.push_back ( new StandardIonType ( "b", 0.0, 1, 1, false ) );

	if ( get_b_plus_h2o () )	nIonList.push_back ( new BPlusWaterIonType () );
	if ( get_c_term_ladder () )	nIonList.push_back ( new StandardIonType ( "c<sub>ladder</sub>", b_plus_h2o_tag_offset, 1, 1, false ) ); 
	if ( getCPlus2Da () )		nIonList.push_back ( new CPlus2DaIonType () );
	if ( getCPlus1Da () )		nIonList.push_back ( new CPlus1DaIonType () );
	if ( get_c () )				nIonList.push_back ( new CIonType () );
	if ( getCMinus1Da () )		nIonList.push_back ( new CMinus1DaIonType () );
	if ( get_d () )				nIonList.push_back ( new DWIonType ( da_html, instInf->getDIonExcludeMask (), d_tag_offset, substituent_da ) ); 
	if ( get_d () )				nIonList.push_back ( new DWIonType ( db_html, instInf->getDIonExcludeMask (), d_tag_offset, substituent_db ) ); 

	if ( get_v () )				cIonList.push_back ( new VIonType () );
	if ( get_w () )				cIonList.push_back ( new DWIonType ( wa_html, instInf->getWIonExcludeMask (), w_tag_offset, substituent_wa ) ); 
	if ( get_w () )				cIonList.push_back ( new DWIonType ( wb_html, instInf->getWIonExcludeMask (), w_tag_offset, substituent_wb ) ); 

	if ( get_x () )				cIonList.push_back ( new StandardIonType ( "x", x_tag_offset, 0, 1, true ) );
	if ( get_n_term_ladder () )	cIonList.push_back ( new StandardIonType ( "n<sub>ladder</sub>", y_tag_offset, 0, 1, false ) );
	if ( get_y () )				cIonList.push_back ( new StandardIonType ( "y", y_tag_offset, 0, 1, false ) );
	BasicLossInfoVector blivY;
	if ( get_y_h3po4 () )	blivY.push_back ( BasicLossInfo ( "H3PO4", h3_p_o4, phosphorylation_mask ) );
	if ( get_y_soch4 () )	blivY.push_back ( BasicLossInfo ( "SOCH4", s_o_c_h4, oxidized_m_mask ) );
	if ( get_y_h2o () )		blivY.push_back ( BasicLossInfo ( "H2O", h2_o, instInf->getWaterLossMask () ) );
	if ( get_y_nh3 () )		blivY.push_back ( BasicLossInfo ( "NH3", n_h3, instInf->getAmmoniaLossMask () ) );
	getLossPermutations ( "y", false, y_tag_offset, sequence, blivY );

	if ( get_Y () )				cIonList.push_back ( new StandardIonType ( "Y", Y_tag_offset, 0, 1, true ) );
	if ( get_z () )				cIonList.push_back ( new ZIonType () );
	if ( getZPlus1Da () )		cIonList.push_back ( new ZPlus1DaIonType () );
	if ( getZPlus2Da () )		cIonList.push_back ( new ZPlus2DaIonType () );
	if ( getZPlus3Da () )		cIonList.push_back ( new ZPlus3DaIonType () );

	if ( get_internal () && seqLength >= 4 && seqLength <= params->getIntValue ( "max_internal_len", 200 ) ) {
		internalIonList.push_back ( new InternalIonType ( "", b_tag_offset, !multiZInternal ) );
		internalIonList.push_back ( new InternalIonType ( "-28", a_tag_offset, !multiZInternal ) );
		string intSequence = sequence.length () < 3 ? "" : string ( sequence, 1, sequence.length () - 2 );
		BasicLossInfoVector blivI;
		blivI.push_back ( BasicLossInfo ( "H2O", h2_o, instInf->getWaterLossMask () ) );
		blivI.push_back ( BasicLossInfo ( "NH3", n_h3, instInf->getAmmoniaLossMask () ) );
		getLossPermutations ( "", true, b_tag_offset, intSequence, blivI );
	}
	StringVector sv = get_M_Losses ();
	for ( int j = 0 ; j < sv.size () ; j++ ) {
		string loss = sv [j];
		if ( loss == "H3PO4" )		intactLossIonList.push_back ( new IntactLossIonType ( getElementalFormulaHTML ( "H3PO4" ), phosphorylation_mask, h3_p_o4 ) );
		else if ( loss == "SOCH4" )	intactLossIonList.push_back ( new IntactLossIonType ( getElementalFormulaHTML ( "SOCH4" ), oxidized_m_mask, s_o_c_h4 ) );
		else if ( loss == "H2O" )	intactLossIonList.push_back ( new IntactLossIonType ( getElementalFormulaHTML ( "H2O" ), instInf->getWaterLossMask (), h2_o ) );
		else if ( loss == "NH3" )	intactLossIonList.push_back ( new IntactLossIonType ( getElementalFormulaHTML ( "NH3" ), instInf->getAmmoniaLossMask (), n_h3 ) );
		else if ( loss == "" )		intactLossIonList.push_back ( new IntactLossIonType ( "", 0xffffffff, 0.0 ) );
		else {
			double offset = isPrefix ( loss, "+" ) ? -atof ( loss.c_str () ) : atof ( loss.c_str () );
			intactLossIonList.push_back ( new IntactLossIonType ( loss, 0xffffffff, offset ) );
		}
	}
}
BiemannParameters::~BiemannParameters ()
{
	delete immIon;
	delete mIon;

	IonTypePtrVectorSizeType i;
	for ( i = 0 ; i < nIonList.size () ; i++ ) {
		delete nIonList [i];
	}
	for ( i = 0 ; i < cIonList.size () ; i++ ) {
		delete cIonList [i];
	}
	for ( i = 0 ; i < internalIonList.size () ; i++ ) {
		delete internalIonList [i];
	}
	for ( i = 0 ; i < intactLossIonList.size () ; i++ ) {
		delete intactLossIonList [i];
	}
}
void BiemannParameters::getLossPermutations ( const string& ionLabel, bool nIonFlag, double offset, const string& sequence, const BasicLossInfoVector& bliv )
{
	int maxNumber = 1;
	nIonLossType = nIonFlag;
	label = ionLabel;
	startOffset = offset;
	masks.clear ();
	lossInfo.clear ();
	for ( StringVectorSizeType i = 0 ; i < bliv.size () ; i++ ) {
		masks.push_back ( bliv [i].mask );
		if ( sequence != "" ) maxNumber = sumResidues ( sequence, bliv [i].mask );
		lossInfo.push_back ( LossInfo ( bliv [i].offset, bliv [i].html, maxNumber ) );
	}
	numMasks = masks.size ();
	numHits.resize (numMasks);
	lab.resize (numMasks);
	if ( !lossInfo.empty () ) {
		getNextLossPermutation ( 0 );
	}
}
void BiemannParameters::getNextLossPermutation ( int level )
{
	for ( int i = 0 ; i <= lossInfo [level].getMaxNumber () ; i++ ) {	// Iterate through masks
		numHits [level] = i;
		lab [level] = ( i == 0 ) ? "" : ( i == 1 ) ? "-" + lossInfo [level].getLabel () : "-" + gen_itoa ( i ) + lossInfo [level].getLabel ();
		if ( level + 1 < numMasks )
			getNextLossPermutation ( level + 1 );
		else {
			string subLabel;
			for ( int j = 0 ; j < numMasks ; j++ ) {
				subLabel += lab [j];
			}
			if ( subLabel != "" ) {
				double offset = startOffset;
				int totalHits = 0;
				for ( int k = 0 ; k < numMasks ; k++ ) {
					offset += numHits [k] * lossInfo [k].getOffset ();
					totalHits += numHits [k];
				}
				if ( maxLosses == "All" || totalHits <= atoi ( maxLosses.c_str () ) ) {
					if ( label == "" ) {
						internalIonList.push_back ( new InternalLossIonType ( subLabel, masks, numHits, offset, !multiZInternal ) );
					}
					else {
						if ( nIonLossType ) {
							nIonList.push_back ( new LossIonType ( label, subLabel, masks, numHits, offset, 1, 1, false ) );
						}
						else {
							cIonList.push_back ( new LossIonType ( label, subLabel, masks, numHits, offset, 0, 1, false ) );
						}
					}
				}
			}
		}
	}
}
StringVector BiemannParameters::initIonType ( const ParameterList* params )
{
	StringVector ionType;
	params->getValue ( "it", ionType );
	params->getValue ( "ion_type", ionType );
	bool useInstrumentIonTypes = params->getBoolValue ( "use_instrument_ion_types", false );
	if ( useInstrumentIonTypes ) {
		ionType = instInf->getIonTypes ();
	}
	return ionType;
}
bool BiemannParameters::get_immonium () const { return inList ( it, "i" ); }
bool BiemannParameters::get_m () const { return inList ( it, "m" ); }
bool BiemannParameters::get_a_h3po4 () const { return inList ( it, "a-H3PO4" ) || inList ( it, "a", "P" ); }
bool BiemannParameters::get_a_soch4 () const { return inList ( it, "a-SOCH4" ) || inList ( it, "a", "S" ); }
bool BiemannParameters::get_a_h2o () const { return inList ( it, "a-H2O" ); }
bool BiemannParameters::get_a_nh3 () const { return inList ( it, "a-NH3" ) || inList ( it, "a", "n" ); }
bool BiemannParameters::get_a () const { return inList ( it, "a" ); }
bool BiemannParameters::get_b_h3po4 () const { return inList ( it, "b-H3PO4" ) || inList ( it, "b", "P" ); }
bool BiemannParameters::get_b_soch4 () const { return inList ( it, "b-SOCH4" ) || inList ( it, "b", "S" ); }
bool BiemannParameters::get_b_h2o () const { return inList ( it, "b-H2O" ) || inList ( it, "b", "h" ); }
bool BiemannParameters::get_b_nh3 () const { return inList ( it, "b-NH3" ) || inList ( it, "b", "n" ); }
bool BiemannParameters::get_b () const { return inList ( it, "b" ); }
bool BiemannParameters::get_b_plus_h2o () const { return inList ( it, "b+H2O" ) || inList ( it, "b", "B" ); }
bool BiemannParameters::get_c_term_ladder () const { return inList ( it, "C" ); }
bool BiemannParameters::getCPlus2Da () const { return inList ( it, "c+2" ); }
bool BiemannParameters::getCPlus1Da () const { return inList ( it, "c+1" ); }
bool BiemannParameters::get_c () const { return inList ( it, "c" ); }
bool BiemannParameters::getCMinus1Da () const { return inList ( it, "c-1" ); }
bool BiemannParameters::get_d () const { return inList ( it, "d" ); }

bool BiemannParameters::get_v () const { return inList ( it, "v" ); }
bool BiemannParameters::get_w () const { return inList ( it, "w" ); }
bool BiemannParameters::get_x () const { return inList ( it, "x" ); }
bool BiemannParameters::get_n_term_ladder () const { return inList ( it, "N" ); }
bool BiemannParameters::get_y () const { return inList ( it, "y" ); }
bool BiemannParameters::get_y_h2o () const { return inList ( it, "y-H2O" ) || inList ( it, "y", "h" ); }
bool BiemannParameters::get_y_nh3 () const { return inList ( it, "y-NH3" ) || inList ( it, "y", "n" ); }
bool BiemannParameters::get_y_soch4 () const { return inList ( it, "y-SOCH4" ) || inList ( it, "y", "S" ); }
bool BiemannParameters::get_y_h3po4 () const { return inList ( it, "y-H3PO4" ) || inList ( it, "y", "P" ); }
bool BiemannParameters::get_Y () const { return inList ( it, "Y" ); }
bool BiemannParameters::get_z () const { return inList ( it, "z" ); }
bool BiemannParameters::getZPlus1Da () const { return inList ( it, "z+1" ); }
bool BiemannParameters::getZPlus2Da () const { return inList ( it, "z+2" ); }
bool BiemannParameters::getZPlus3Da () const { return inList ( it, "z+3" ); }

bool BiemannParameters::get_bp2 ()		const { return inList ( it, "bp2" ); }
bool BiemannParameters::get_bp3 ()		const { return inList ( it, "bp3" ); }
bool BiemannParameters::get_bp2_h3po4 ()const { return inList ( it, "bp2-H3PO4" ); }
bool BiemannParameters::get_bp2_soch4 ()const { return inList ( it, "bp2-SOCH4" ); }
bool BiemannParameters::get_bp2_h2o ()	const { return inList ( it, "bp2-H2O" ); }
bool BiemannParameters::get_bp2_nh3 ()	const { return inList ( it, "bp2-NH3" ); }
bool BiemannParameters::get_cp2 ()		const { return inList ( it, "cp2" ); }
bool BiemannParameters::get_yp2 ()		const { return inList ( it, "yp2" ); }
bool BiemannParameters::get_yp3 ()		const { return inList ( it, "yp3" ); }
bool BiemannParameters::get_yp2_h3po4 ()const { return inList ( it, "yp2-H3PO4" ); }
bool BiemannParameters::get_yp2_soch4 ()const { return inList ( it, "yp2-SOCH4" ); }
bool BiemannParameters::get_yp2_h2o ()	const { return inList ( it, "yp2-H2O" ); }
bool BiemannParameters::get_yp2_nh3 ()	const { return inList ( it, "yp2-NH3" ); }
bool BiemannParameters::get_zp2 ()		const { return inList ( it, "zp2" ); }
bool BiemannParameters::get_zPlus1Dap2 () const { return inList ( it, "z+1p2" ); }

StringVector BiemannParameters::get_M_Losses () const
{
	StringVector sv;
	bool Pset = false;
	bool Sset = false;
	bool hset = false;
	bool nset = false;
	sv.push_back ( "" );
	for ( StringVectorSizeType i = 0 ; i < it.size () ; i++ ) {
		string s = it [i];
		if ( s == "P" ) s = "M-H3PO4";
		if ( s == "S" ) s = "M-SOCH4";
		if ( s == "h" ) s = "M-H2O";
		if ( s == "n" ) s = "M-NH3";
		if ( isPrefix ( s, "M-" ) ) {
			string str = s.substr ( 2 );
			if ( str == "H3PO4" ) {
				if ( Pset ) continue;
				else Pset = true;
			}
			if ( str == "SOCH4" ) {
				if ( Sset ) continue;
				else Sset = true;
			}
			if ( str == "H2O" ) {
				if ( hset ) continue;
				else hset = true;
			}
			if ( str == "NH3" ) {
				if ( nset ) continue;
				else nset = true;
			}
			sv.push_back ( str );
		}
		if ( isPrefix ( s, "M+" ) ) sv.push_back ( s.substr ( 1 ) );
	}
	return sv;
}
bool BiemannParameters::get_M () const { return true; }

bool BiemannParameters::get_internal () const { return inList ( it, "I" ); }

bool BiemannParameters::inList ( const StringVector& it, const string& s1 ) const
{
	return find ( it.begin (), it.end (), s1 ) != it.end ();
}
bool BiemannParameters::inList ( const StringVector& it, const string& s1, const string& s2 ) const
{
	return find ( it.begin (), it.end (), s1 ) != it.end () && find ( it.begin (), it.end (), s2 ) != it.end ();
}
void BiemannParameters::printHTML ( ostream& os ) const
{
	if ( maxLosses != "1" ) ParameterList::printHTML ( os, "Max Losses", maxLosses );
	ParameterList::printHTMLContainer ( os, "Ion Types Considered", it );
}
FragmentPeakLabel::FragmentPeakLabel ( const string& label, const string& subLabel, int number, int charge ) :
	label ( label ),
	subLabel ( subLabel ),
	number ( number ),
	charge ( charge )
{
}
FragmentPeakLabel::FragmentPeakLabel ( const StringVector& label, const string& subLabel, int number, int charge ) :
	label ( getLabel ( label ) ),
	subLabel ( subLabel ),
	number ( number ),
	charge ( charge )
{
}
string FragmentPeakLabel::getLabel ( const StringVector& lab )
{
	string label;
	for ( StringSizeType i = 0 ; i < lab.size () ; i++ )
		label += lab [i];
	return label;
}
void FragmentPeakLabel::printHTML ( ostream& os ) const
{
	os << label;
	if ( number != -1 ) {
		os << "<sub>";
		os << number;
		os << "</sub>";
	}
	os << subLabel;
	print_charge_superscript ( os, charge );
}
void FragmentPeakLabel::printXML ( ostream& os ) const
{
	os << gen_strstriptags ( label );
	os << gen_strstriptags ( subLabel );
	os << ',';
	if ( number != -1 ) {
		os << number;
	}
	os << ',';
	if ( charge >= 1 ) os << "+";
	if ( charge < 0 ) os << "-";
	os << charge;
}
FragmentPeak::FragmentPeak ( const FragmentPeakLabel& fragmentPeakLabel, int category, double mOverZ, int charge, double intensity ) :
	Peak ( mOverZ, 0.0, charge, intensity, cation_wt ),
	category ( category ),
	fragmentPeakLabel ( fragmentPeakLabel )
{
}
bool FragmentPeak::isMatch ( double mass, int charge ) const
{
	if ( mass == getMOverZ () && charge == getCharge () ) return true;
	else return false;
}
void FragmentPeak::printLabelHTML ( ostream& os ) const
{
	fragmentPeakLabel.printHTML ( os );
}
void FragmentPeak::printLabelXML ( ostream& os ) const
{
	fragmentPeakLabel.printXML ( os );
}
string FragmentPeak::getLabelXML () const
{
	ostringstream ostr;
	fragmentPeakLabel.printXML ( ostr );
	return ostr.str ();
}
string FragmentPeak::getLabelHTML () const
{
	stringstream ss;
	fragmentPeakLabel.printHTML ( ss );
	return ss.str ();
}
int FragmentPeak::getCategoryHTML () const
{
	return category;
}
double IonSeries::minMOverZ = 0.0;
double IonSeries::maxMOverZ = 0.0;
IonSeries::IonSeries ( const IonType* ionType, const PeptideInfo& peptideInfo ) :
	ionType ( ionType ),
	peptideInfo ( peptideInfo ),
	printCol ( peptideInfo.getMaxZ (), false ),
	printSeries ( false )
{
	int maxCharge =  peptideInfo.getMaxZ ();
	int length =  peptideInfo.getPepLen ();
	ions.resize ( maxCharge );
	for ( DoubleVectorVectorSizeType i = 0 ; i < ions.size () ; i++ ) {
		ions [i].resize ( length );
		fill ( ions [i].begin (), ions [i].end (), 0.0 );
	}
	seriesFlag = ionType->getSeriesFlag ( peptideInfo );
	if ( ionType->getStartGap () == -1 ) {
		int n = 0;
		for ( int j = 0 ; j < length ; j++ ) {
			if ( peptideInfo.getSeriesCompMask ( j ) & ionType->getMask () ) n++;
		}
		if ( n > 1 ) n = 1;
		calculateIntactLossMasses ( n );
	}
	else
		calculateMasses ( ionType->getStartGap (), length - ionType->getEndGap () );
	calculatePrintCol ();
}
void IonSeries::calculatePrintCol ()
{
	for ( DoubleVectorVectorSizeType i = 0 ; i < ions.size () ; i++ ) {
		DoubleVectorSizeType size = ions [i].size ();
		for ( DoubleVectorSizeType j = 0 ; j < size ; j++ ) {
			if ( ions [i][j] != 0.0 ) {
				printCol [i] = true;
				printSeries = true;
			}
		}
	}
}
int IonSeries::getMaxCharge ( int index )
{
	if ( ionType->getSingleChargeOnly () ) return 1;
	else return peptideInfo.getSeriesMaxZ ( index );
}
void IonSeries::calculateMasses ( int start, int end )
{
	for ( int i = start ; i < end ; i++ ) {
		if ( seriesFlag [i] ) {
			for ( int charge = 1 ; charge <= getMaxCharge ( i ) ; charge++ ) {
				ions [charge-1][i] = ionType->calculateIonMass ( peptideInfo, i, charge );
			}
		}
	}
}
void IonSeries::calculateIntactLossMasses ( int n )
{
	double offset = ionType->getOffset ();
	int maxCharge = peptideInfo.getIntactLossMaxZ ( offset );
	if ( ions.size () < maxCharge ) {
		ions.resize ( maxCharge );
		for ( int z = 0 ; z < maxCharge ; z++ ) {
			ions [z].resize ( n );
		}
	}
	for ( int i = 0 ; i < n ; i++ ) {
		for ( int charge = 1 ; charge <= maxCharge ; charge++ ) {
			ions [charge-1][i] = ionType->calculateIonMass ( peptideInfo, i+1, charge );
		}
	}
}
void IonSeries::printHTMLHeader ( ostream& os, bool noMaxCharge ) const
{
	for ( DoubleVectorVectorSizeType i = 0 ; i < ions.size () ; i++ ) {
		if ( !noMaxCharge && i > 1 ) break;
		if ( printCol [i] ) {
			tableHeaderStart ( os, "", "center" );
				os << ionType->getLabel () << ionType->getSubLabel ();
				if ( i+1 != 1 ) {
					os << "<sup>" << showpos << (int)i+1 << noshowpos << "</sup>";
				}
				os << endl;
			tableHeaderEnd ( os );
		}
	}
}
void IonSeries::printHTMLEmptyCells ( ostream& os, bool noMaxCharge ) const
{
	for ( DoubleVectorVectorSizeType i = 0 ; i < ions.size () ; i++ ) {
		if ( !noMaxCharge && i > 1 ) break;
		if ( printCol [i] ) tableEmptyCell ( os );
	}
}
void IonSeries::printHTML ( ostream& os, int precision, const SpectrumMatch* sm, int sNum, int idx, bool noMaxCharge ) const		// Vertical sequence
{
	for ( DoubleVectorVectorSizeType i = 0 ; i < ions.size () ; i++ ) {
		if ( !noMaxCharge && i > 1 ) break;
		if ( printCol [i] ) {
			bool hit = sm ? sm->isMatch ( ions [i][idx], i+1, sNum ) : false;
			printIonMassHTML ( os, ions [i][idx], i+1, precision, hit );
		}
	}
}
void IonSeries::printHTML ( ostream& os, bool ascending, int precision, const SpectrumMatch* sm, int sNum ) const	// Horizontal sequence
{
	DoubleVectorVectorSizeType i;
	for ( i = 0 ; i < ions.size () ; i++ ) {
		if ( printCol [i] ) {
			tableRowStart ( os );
				tableDataStart ( os, "", "left" );
					os << ionType->getLabel () << ionType->getSubLabel ();
					if ( i+1 != 1 ) {
						os << "<sup>" << showpos << (int)i+1 << noshowpos << "</sup>";
					}
					os << endl;
				tableDataEnd ( os );

				DoubleVectorSizeType size = ions [i].size ();
				DoubleVectorSizeType j;
				if ( ascending ) {
					for ( j = 0 ; j < size ; j++ ) {
						bool hit = sm ? sm->isMatch ( ions [i][j], i+1, sNum ) : false;
						printIonMassHTML ( os, ions [i][j], i+1, precision, hit );
					}
				}
				else {
					for ( j = size ; j-- ; ) {
						bool hit = sm ? sm->isMatch ( ions [i][j], i+1, sNum ) : false;
						printIonMassHTML ( os, ions [i][j], i+1, precision, hit );
					}
				}
			tableRowEnd ( os );
		}
	}
}
void IonSeries::printXML ( ostream& os, bool ascending, int precision ) const
{
	if ( printSeries ) {
		os << "<ion_series>" << endl;
		ParameterList::printXML ( os, "name", gen_strstriptags ( ionType->getLabel () + ionType->getSubLabel () ) );
		for ( DoubleVectorVectorSizeType i = 0 ; i < ions.size () ; i++ ) {
			if ( printCol [i] ) {
				DoubleVectorSizeType size = ions [i].size ();
				if ( ascending ) {
					for ( DoubleVectorSizeType j = 0 ; j < size ; j++ )
						printIonMassXML ( os, ions [i][j], i+1, j+1, precision );
				}
				else {
					for ( DoubleVectorSizeType j = size ; j-- ; )
						printIonMassXML ( os, ions [i][j], i+1, j+1, precision );
				}
			}
		}
		os << "</ion_series>" << endl;
	}
}
void IonSeries::printIonMassHTML ( ostream& os, double mOverZ, int charge, int precision, bool hit )
{
	static const char* styleID [] = { "", "", "msp_z_2", "msp_z_3", "msp_z_n" };
	if ( mOverZ ) {
		int idx;
		if ( charge < 2 ) idx = 1;
		else if ( charge == 2 ) idx = 2;
		else if ( charge == 3 ) idx = 3;
		else idx = 4;
		tableDataStart ( os, styleID [idx], "right" );
		bool outOfRange = minMOverZ > 0.0 && ( mOverZ < minMOverZ || mOverZ > maxMOverZ );
		if ( outOfRange ) os << "<i><font size=\"-2\">";
		if ( hit ) os << "<font color=\"#FF0000\">";
		genPrint ( os, mOverZ, precision );
		//print_charge_superscript ( os, charge );
		if ( hit ) os << "</font>";
		if ( outOfRange ) os << "</font></i>";
		os << endl;
		tableDataEnd ( os );
	}
	else {
		tableDataStart ( os, "", "center" );
			os << "---" << endl;
		tableDataEnd ( os );
	}
}
void IonSeries::printIonMassXML ( ostream& os, double mOverZ, int charge, int number, int precision )
{
	if ( mOverZ ) {
		os << "<i" << number << ">" << endl;
			os << "<mz>";
				genPrint ( os, mOverZ, precision );
				os << ",";
				os << charge;
			os << "</mz>";
			os << endl;
		os << "</i" << number << ">" << endl;
	}
	else {
		os << "<i" << number << " />" << endl;
	}
}
void IntactLossIonSeries::printHTML ( ostream& os, bool ascending, int precision, const SpectrumMatch* sm, int sNum ) const
{
	bool printRow = false;
	for ( DoubleVectorVectorSizeType i = 0 ; i < ions.size () ; i++ ) {
		for ( DoubleVectorSizeType j = ions [i].size () ; j-- ; ) {
			if ( ions [i][j] != 0.0 ) printRow = true;
		}
	}
	if ( printRow ) {
		tableRowStart ( os );
		os << "<th align=\"center\">";
		os << ionType->getLabel ();
		if ( !ionType->getSubLabel ().empty () ) {
			if ( !isPrefix ( ionType->getSubLabel (), "+" ) ) os << "-";
			os << ionType->getSubLabel ();
		}
		os << "</th>" << endl;
		for ( DoubleVectorVectorSizeType i = 0 ; i < ions.size () ; i++ ) {
			DoubleVectorSizeType size = ions [i].size ();
			for ( DoubleVectorSizeType j = size ; j-- ; ) {
				if ( ions [i][j] != 0.0 ) {
					bool hit = sm ? sm->isMatch ( ions [i][j], i+1, sNum ) : false;
					printIonMassHTML ( os, ions [i][j], i+1, precision, hit );
				}
			}
		}
		tableRowEnd ( os );
	}
}
void IonSeries::getPeaks ( FragmentPeakVector& fpc ) const
{
	for ( DoubleVectorVectorSizeType i = 0 ; i < ions.size () ; i++ ) {
		for ( DoubleVectorSizeType j = 0 ; j < ions [i].size () ; j++ ) {
			double mOverZ = ions [i][j];
			if ( mOverZ ) {
				FragmentPeakLabel fpl ( ionType->getLabel (), ionType->getSubLabel (), j+1, i+1 );
				int category = ( ionType->getSubLabel () == "" ) ? SERIES_CATEGORY : LOSS_CATEGORY;
				fpc.push_back ( FragmentPeak ( fpl, category, mOverZ, i+1, 100.0 ) );
			}
		}
	}
}
void IntactLossIonSeries::getPeaks ( FragmentPeakVector& fpc ) const
{
	for ( DoubleVectorVectorSizeType i = 0 ; i < ions.size () ; i++ ) {
		int loss = 1;
		for ( DoubleVectorSizeType j = 0 ; j < ions [i].size () ; j++ ) {
			double mOverZ = ions [i][j];
			if ( mOverZ ) {
				string subLabel = ionType->getSubLabel ();
				if ( !subLabel.empty () ) {
					string minus = !isPrefix ( subLabel, "+" ) ? "-" : "";
					if ( loss > 1 ) subLabel = minus + gen_itoa ( loss ) + subLabel;
					else subLabel = minus + subLabel;
				}
				FragmentPeakLabel fpl ( ionType->getLabel (), subLabel, -1, i+1 );
				int category = subLabel.empty () ? SERIES_CATEGORY : LOSS_CATEGORY;
				fpc.push_back ( FragmentPeak ( fpl, category, mOverZ, i+1, 100.0 ) );
				loss++;
			}
		}
	}
}
void ImmoniumIonSeries::printHTML ( ostream& os, bool ascending, int precision, const SpectrumMatch* sm, int sNum ) const
{
	tableRowStart ( os );
		tableCell ( os, "Immonium Ions" );
		for ( DoubleVectorSizeType j = 0 ; j < ions [0].size () ; j++ ) {
			tableDataStart ( os, "", "right" );
				genPrint ( os, ions [0][j], precision );
				os << endl;
			tableDataEnd ( os );
		}
	tableRowEnd ( os );
}
void ImmoniumIonSeries::getPeaks ( FragmentPeakVector& fpc ) const
{
	for ( DoubleVectorSizeType j = 0 ; j < ions [0].size () ; j++ ) {
		double mOverZ = ions [0][j];
		if ( mOverZ ) {
			FragmentPeakLabel fpl ( string ( 1, peptideInfo.getSeriesAA ( j ) ), "", -1, 1 );
			fpc.push_back ( FragmentPeak ( fpl, IMMONIUM_CATEGORY, mOverZ, 1, 100.0 ) );
		}
	}
}
PIonSeries::PIonSeries ( const string& label, const string& subLabel, int category, double mass, int maxCharge ) :
	label ( label ),
	subLabel ( subLabel ),
	category ( category ),
	maxCharge ( maxCharge )
{
	for ( int ch = 1 ; ch <= maxCharge ; ch++ ) {
		ions.push_back ( calculateIonMOverZ ( mass,	ch ) );
	}
}
void PIonSeries::printHTML ( ostream& os, int precision, const SpectrumMatch* sm, int sNum ) const
{
	tableRowStart ( os );
		os << "<th align=\"center\">" << label << subLabel << "</th>" << endl;
		for ( int ch = 1 ; ch <= maxCharge ; ch++ ) {
			int i = ch - 1;
			bool hit = sm ? sm->isMatch ( ions [i], ch, sNum ) : false;
			IonSeries::printIonMassHTML ( os, ions [i], ch, precision, hit );
		}
	tableRowEnd ( os );
}
void PIonSeries::getPeaks ( FragmentPeakVector& fpc ) const
{
	for ( int ch = 1 ; ch <= maxCharge ; ch++ ) {
		int i = ch - 1;
		fpc.push_back ( FragmentPeak ( FragmentPeakLabel ( label, subLabel, -1, ch ), category, ions [i], ch, 100.0 ) );
	}
}

RelatedIonSeries::RelatedIonSeries ( const PeptideInfo& peptideInfo ) :
	peptideInfo ( peptideInfo )
{
	int length = peptideInfo.getPepLen ();
	ions.resize ( length );
	for ( int i = 0 ; i < length ; i++ ) {
		ions [i] = ImmoniumInfo::getRelatedIons ( peptideInfo.getSeriesAA ( i ) );
	}
}
void RelatedIonSeries::printHTMLEmptyCells ( ostream& os ) const
{
	tableEmptyCell ( os ); 
}
void RelatedIonSeries::printHTMLHeader ( ostream& os ) const
{
	tableHeader ( os, "Low Mass" ); 
}
void RelatedIonSeries::printHTML ( ostream& os, int precision, const SpectrumMatch* sm, int sNum, int idx ) const
{
	DoubleVectorSizeType len = ions [idx].size ();
	if ( len ) {
		tableDataStart ( os, "", "center" );
			for ( DoubleVectorSizeType i = 0 ; i < len ; i++ ) {
				bool hit = sm ? sm->isMatch ( ions [idx][i], 1, sNum ) : false;
				if ( hit ) os << "<font color=\"#FF0000\">";
				genPrint ( os, ions [idx][i], precision );
				if ( hit ) os << "</font>";
				if ( i != len - 1 ) os << "<br />" << endl;
			}
			os << endl;
		tableDataEnd ( os );
	}
	else {
		tableDataStart ( os, "", "center" );
			os << "---" << endl;
		tableDataEnd ( os );
	}
}
void RelatedIonSeries::printHTML ( ostream& os, int precision, const SpectrumMatch* sm, int sNum ) const
{
	tableRowStart ( os );
	printHTMLHeader ( os );
	for ( DoubleVectorVectorSizeType i = 0 ; i < ions.size () ; i++ ) {
		DoubleVectorSizeType len = ions [i].size ();
		if ( len ) {
			tableDataStart ( os, "", "right" );
				for ( DoubleVectorSizeType j = 0 ; j < len ; j++ ) {
					bool hit = sm ? sm->isMatch ( ions [i][j], 1, sNum ) : false;
					if ( hit ) os << "<font color=\"#FF0000\">";
					genPrint ( os, ions [i][j], precision );
					if ( hit ) os << "</font>";
					if ( j != len - 1 ) os << "<br />" << endl;
				}
				os << endl;
			tableDataEnd ( os );
		}
		else {
			tableDataStart ( os, "", "center" );
				os << "---" << endl;
			tableDataEnd ( os );
		}
	}
	tableRowEnd ( os );
}
void RelatedIonSeries::printXML ( ostream& os, int precision ) const
{
	os << "<ion_series>" << endl;
		ParameterList::printXML ( os, "name", "Immonium and Related" );
		for ( DoubleVectorVectorSizeType i = 0 ; i < ions.size () ; i++ ) {
			DoubleVectorSizeType len = ions [i].size ();
			if ( len ) {
				os << "<i" << i+1 << ">" << endl;
					for ( DoubleVectorSizeType j = 0 ; j < len ; j++ ) {
						printIonMassXML ( os, ions [i][j], j+1, precision );
					}
				os << "</i" << i+1 << ">" << endl;
			}
			else {
				os << "<i" << i+1 << " />" << endl;
			}
		}
	os << "</ion_series>" << endl;
}
void RelatedIonSeries::printIonMassXML ( ostream& os, double mOverZ, int number, int precision )
{
	if ( mOverZ ) {
		os << "<mz>";
		genPrint ( os, mOverZ, precision );
		os << ",";
		os << 1;
		os << "</mz>";
		os << endl;
	}
}
void RelatedIonSeries::getPeaks ( FragmentPeakVector& fpc ) const
{
	set <char> uniqueAA;
	for ( DoubleVectorVectorSizeType i = 0 ; i < ions.size () ; i++ ) {			// For each amino acid
		char aa = peptideInfo.getSeriesAA ( i );
		pair <SetCharIterator, bool> flag = uniqueAA.insert ( aa );
		if ( flag.second ) {
			for ( DoubleVectorSizeType j = 0 ; j < ions [i].size () ; j++ ) {
				FragmentPeakLabel fpl ( string ( 1, aa ), "", -1, 1 );
				fpc.push_back ( FragmentPeak ( fpl, IMMONIUM_CATEGORY, ions [i][j], 1, 100.0 ) );
			}
		}
	}
}
InternalIonSeries::InternalIonSeries ( const IonType* ionType, const PeptideInfo& peptideInfo ) :
	ionType ( ionType ),
	peptideInfo ( peptideInfo )
{
	int maxCharge = ionType->getSingleChargeOnly () ? 1 : peptideInfo.getMaxZ ();
	int length =  peptideInfo.getNumInternalIons ();
	ions.resize ( maxCharge );
	for ( DoubleVectorVectorSizeType i = 0 ; i < ions.size () ; i++ ) {
		ions [i].resize ( length );
		fill ( ions [i].begin (), ions [i].end (), 0.0 );
	}
	seriesFlag = ionType->getSeriesFlag ( peptideInfo );
	calculateMasses ( 0, length );
}
void InternalIonSeries::calculateMasses ( int start, int end )
{
	for ( int i = start ; i < end ; i++ ) {
		if ( seriesFlag [i] ) {
			for ( int charge = 1 ; charge <= getMaxCharge ( i ) ; charge++ ) {
				ions [charge-1][i] = ionType->calculateIonMass ( peptideInfo, i, charge );
			}
		}
	}
}
int InternalIonSeries::getMaxCharge ( int index )
{
	if ( ionType->getSingleChargeOnly () ) return 1;
	else return peptideInfo.getInternalMaxZ ( index );
}
void InternalIonSeries::printLabel ( ostream& os ) const
{
	os << ionType->getLabel () << ionType->getSubLabel () << " ions";
}
void InternalIonSeries::getPeaks ( FragmentPeakVector& fpc ) const
{
	for ( DoubleVectorVectorSizeType i = 0 ; i < ions.size () ; i++ ) {
		for ( DoubleVectorSizeType j = 0 ; j < ions [i].size () ; j++ ) {
			double mOverZ = ions [i][j];
			if ( mOverZ ) {
				FragmentPeakLabel fpl ( peptideInfo.getInternalSeriesSequence ( j ), ionType->getSubLabel (), -1, i+1 );
				fpc.push_back ( FragmentPeak ( fpl, INTERNAL_CATEGORY, mOverZ, i+1, 100.0 ) );
			}
		}
	}
}
PepInfo::PepInfo ( const StringVector& pep, double nTermMass, const string& nTermStr, double cTermMass, const string& cTermStr, bool show ) :
	pep ( pep ),
	pepLen ( pep.size () ),
	nTermMass ( nTermMass ),
	nTermStr ( nTermStr ),
	cTermMass ( cTermMass ),
	cTermStr ( cTermStr ),
	show ( show )
{
}
double PepInfo::peptideFormulaToPIon () const
{
	double molWT = cation_wt;
	if ( !genIsNumberStart ( nTermStr [0] ) )	molWT += nTermMass;
	else										molWT += h1;

	if ( !genIsNumberStart ( cTermStr [0] ) )	molWT += cTermMass;
	else										molWT += o_h;

	for ( StringVectorSizeType i = 0 ; i < pep.size () ; i++ ) {
		const string& p = pep [i];
		if ( p.length () > 1 ) {
			if ( genIsNumberStart ( p [2] ) )
				molWT += amino_acid_wt [p[0]];
			else
				molWT += getAminoAcidWt ( p );
		}
		else
			molWT += amino_acid_wt [p[0]];
	}
	return molWT;
}
double PepInfo::peptideFormulaToMolecularWeight () const
{
	double molWT = nTermMass + cTermMass + cation_wt;

	for ( StringVectorSizeType i = 0 ; i < pep.size () ; i++ ) {
		if ( pep [i].length () > 1 )
			molWT += getAminoAcidWt ( pep [i] );
		else
			molWT += amino_acid_wt [pep [i][0]];
	}
	return molWT;
}
void PepInfo::printSequence ( ostream& os, int idx ) const
{
	tableCell ( os, pep [idx], false, true );
}
void PepInfo::printNTerm ( ostream& os ) const
{
	tableCell ( os, nTermStr, false, true );
}
void PepInfo::printCTerm ( ostream& os ) const
{
	tableCell ( os, cTermStr, false, true );
}
void PepInfo::printSequence ( ostream& os ) const
{
	int i;

	tableEmptyRow ( os );

	tableRowStart ( os );
		tableEmptyCell ( os );
		for ( i = 0 ; i < pepLen ; i++ ) {
			tableHeaderStart ( os, "", "right" );
				os << i+1 << endl;
			tableHeaderEnd ( os );
		}
	tableRowEnd ( os );

	tableRowStart ( os );
		tableHeaderStart ( os, "", "right" );
			os << nTermStr;
			os << " - ";
			os << endl;
		tableHeaderEnd ( os );
		for ( i = 0 ; i < pepLen ; i++ ) {
			tableCell ( os, pep [i], false, true );
		}
		tableHeaderStart ( os );
			os << " - ";
			os << cTermStr;
			os << endl;
		tableHeaderEnd ( os );
	tableRowEnd ( os );

	tableRowStart ( os );
		tableEmptyCell ( os );
		for ( i = 0 ; i < pepLen ; i++ ) {
			tableHeaderStart ( os, "", "left" );
				os << pepLen - i << endl;
			tableHeaderEnd ( os );
		}
	tableRowEnd ( os );

	tableEmptyRow ( os );
}
SInfo::SInfo ( const PepInfo* pInfo, const BiemannParameters* bp, int maxReportedCharge, int precision, bool countPosCharges, const LinkInfo* linkInfo ) :
	pInfo ( pInfo ),
	precision ( precision ),
	relatedIonSeries ( 0 ),
	mIonSeries ( 0 ),
	theoreticalSpectrum ( 0 ),
	linkInfo ( linkInfo )
{
	if ( pInfo->show ) {
		int i;
		double mPlusCation = pInfo->peptideFormulaToMolecularWeight ();
		bool intIonsFlag = ( bp->getInternalIonListSize () > 0 );
		PeptideInfo nPeptideInfo ( pInfo->pep, true, maxReportedCharge, mPlusCation, pInfo->nTermMass, intIonsFlag, countPosCharges );
		if ( bp->getImmIon () ) {
			relatedIonSeries = new RelatedIonSeries ( nPeptideInfo );
		}
		if ( bp->getMIon () )	mIonSeries = new IonSeries ( bp->getMIon (), nPeptideInfo );
		bool b1Flag = B1Info::instance ().isB1 ( pInfo->nTermStr );
		for ( i = 0 ; i < bp->getNIonListSize () ; i++ ) {
			IonType* it = bp->getNIonList (i);
			bool a_or_b = it->getLabel () == "a" || it->getLabel () == "b";
			if ( b1Flag && a_or_b ) it->modifyStartGap ( -1 );
			if ( ( it->getLabel () == "b" || it->getLabel () == "c" ) && ( it->getSubLabel ().empty () || it->getSubLabel () == "-H<sub>3</sub>PO<sub>4</sub>" ) ) {
				mainNTermIonIndicies.insert ( i );
			}
			nIonSeriesList.push_back ( new IonSeries ( it, nPeptideInfo ) );
			if ( b1Flag && a_or_b ) it->modifyStartGap ( 1 );
		}
		for ( i = 0 ; i < bp->getInternalIonListSize () ; i++ ) { 
			internalIonSeriesList.push_back ( new InternalIonSeries ( bp->getInternalIonList (i), nPeptideInfo ) );
		}
		for ( i = 0 ; i < bp->getIntactLossIonListSize () ; i++ ) { 
			intactLossIonSeriesList.push_back ( new IntactLossIonSeries ( bp->getIntactLossIonList (i), nPeptideInfo ) );
		}
		PeptideInfo cPeptideInfo ( reverseVector ( pInfo->pep ), false, maxReportedCharge, mPlusCation, pInfo->cTermMass, false, countPosCharges );
		for ( i = 0 ; i < bp->getCIonListSize () ; i++ ) {
			IonType* it = bp->getCIonList (i);
			bool y = it->getLabel () == "y";
			if ( b1Flag && y ) it->modifyEndGap ( -1 );
			if ( ( y || it->getLabel () == "z" ) && ( it->getSubLabel ().empty () || it->getSubLabel () == "-H<sub>3</sub>PO<sub>4</sub>" ) ) {
				mainCTermIonIndicies.insert ( i );
			}
			cIonSeriesList.push_back ( new IonSeries ( bp->getCIonList (i), cPeptideInfo ) );
			if ( b1Flag && y ) it->modifyEndGap ( 1 );
		}
		theoreticalSpectrum = new TheoreticalSpectrum ( precision );
		if ( relatedIonSeries )	theoreticalSpectrum->add ( relatedIonSeries );
		if ( mIonSeries )		theoreticalSpectrum->add ( mIonSeries );
		theoreticalSpectrum->add ( nIonSeriesList );
		theoreticalSpectrum->add ( cIonSeriesList );
		theoreticalSpectrum->add ( internalIonSeriesList );
		theoreticalSpectrum->add ( intactLossIonSeriesList );
		if ( linkInfo && linkInfo->getName () != "No Link" && linkInfo->getName () != "User Defined Link" ) {
			addXLinkIons ( maxReportedCharge );
		}
		theoreticalSpectrum->sortPeaks ();
	}
}
SInfo::~SInfo ()
{
	if ( pInfo->show ) {
		delete relatedIonSeries;
		delete mIonSeries;
		IonSeriesListSizeType i;
		for ( i = 0 ; i < nIonSeriesList.size () ; i++ )
			delete nIonSeriesList [i];
		for ( i = 0 ; i < cIonSeriesList.size () ; i++ )
			delete cIonSeriesList [i];
		InternalIonSeriesListSizeType j;
		for ( j = 0 ; j < internalIonSeriesList.size () ; j++ )
			delete internalIonSeriesList [j];
		IntactLossIonSeriesListSizeType k;
		for ( k = 0 ; k < intactLossIonSeriesList.size () ; k++ )
			delete intactLossIonSeriesList [k];
		delete theoreticalSpectrum;
	}
}
void SInfo::addXLinkIons ( int mxZ )
{
	static int num = 0;
	num++;
	string lab;
	if ( num == 1 ) lab = "P";
	else			lab = "P";		// was Q
	string labXL = lab + "+L";
	double bridgeMass = massConvert ( linkInfo->getBridgeFormula ().c_str () );
	double p = pInfo->peptideFormulaToPIon ();
	if ( instInf->getETDTypeFragmentation () ) {
		pIonList.push_back ( new PIonSeries ( lab,	"", XLINK_CATEGORY, p,					mxZ ) );
		pIonList.push_back ( new PIonSeries ( labXL,"", XLINK_CATEGORY, p+bridgeMass+h1,	mxZ ) );
	}
	else {
		string water = "-H<sub>2</sub>O";
		string water2 = "-2H<sub>2</sub>O";
		string ammonia = "-NH<sub>3</sub>";
		pIonList.push_back ( new PIonSeries ( lab,	"",		XLINK_CATEGORY,p,				mxZ ) );
		pIonList.push_back ( new PIonSeries ( lab,	water,	XLINK_CATEGORY,p-h2_o,			mxZ ) );
		pIonList.push_back ( new PIonSeries ( lab,	water2,	XLINK_CATEGORY,p-h2_o-h2_o,		mxZ ) );
		pIonList.push_back ( new PIonSeries ( lab,	ammonia,XLINK_CATEGORY,p-n_h3,			mxZ ) );
		pIonList.push_back ( new PIonSeries ( labXL,"",		XLINK_CATEGORY,p+bridgeMass,	mxZ ) );
		pIonList.push_back ( new PIonSeries ( labXL,water,	XLINK_CATEGORY,p+bridgeMass-h2_o,mxZ ) );
		pIonList.push_back ( new PIonSeries ( labXL,ammonia,XLINK_CATEGORY,p+bridgeMass-n_h3,mxZ ) );

		double immoniumMass = massConvert ( linkInfo->getPCIDImmFormula ().c_str () );
		pIonList.push_back ( new PIonSeries ( labXL,"K",				XLINK_CATEGORY,	p+bridgeMass+immoniumMass,	mxZ ) );
		pIonList.push_back ( new PIonSeries ( labXL,"K-H<sub>2</sub>O",XLINK_CATEGORY,	p+bridgeMass+immoniumMass-h2_o,	mxZ ) );
		pIonList.push_back ( new PIonSeries ( labXL,"K-NH<sub>3</sub>",XLINK_CATEGORY,	p+bridgeMass+immoniumMass-n_h3,	mxZ ) );
	}
	theoreticalSpectrum->add ( pIonList );
}
void SInfo::printHTML ( ostream& os, const SpectrumMatch* spectrumMatch, int sNum ) const
{
	if ( pInfo->pepLen && pInfo->show ) {
		ExpandableJavascriptBlock ejb ( "Main Sequence Ions", spectrumMatch == 0 );
		ejb.printHeader ( os );
		os << "<br />" << endl;
		os << "<br />" << endl;
		printVerticalSequence ( os, spectrumMatch, sNum );
		ejb.printFooter ( os );
		os << "<br />" << endl;

		ExpandableJavascriptBlock ejb2 ( "All Sequence Ions", spectrumMatch == 0 );
		ejb2.printHeader ( os );
		os << "<br />" << endl;
		os << "<br />" << endl;
		os << "<table cellspacing=\"4\" class=\"msprod_ion_table\">" << endl;
			for ( IonSeriesListSizeType i = 0 ; i < intactLossIonSeriesList.size () ; i++ ) {
				intactLossIonSeriesList [i]->printHTML ( os, true, precision, spectrumMatch, sNum );
			}
			for ( PIonSeriesListSizeType j = 0 ; j < pIonList.size () ; j++ ) {
				pIonList [j]->printHTML ( os, precision, spectrumMatch, sNum );
			}
		os << "</table>" << endl;
		printHorizontalSequence ( os, spectrumMatch, sNum );

		ejb2.printFooter ( os );
		os << "<br />" << endl;
		printInternalIons ( os, spectrumMatch, sNum );
		theoreticalSpectrum->printPeakTable ( os, 5, spectrumMatch, sNum );
	}
}
void SInfo::printVerticalSequence ( ostream& os, const SpectrumMatch* spectrumMatch, int sNum ) const
{
	StringVectorVector svv;
	os << "<table cellspacing=\"4\" class=\"msprod_ion_table\">" << endl;
		tableRowStart ( os );
			for ( SetIntConstIterator i1 = mainNTermIonIndicies.begin () ; i1 != mainNTermIonIndicies.end () ; i1++ ) {
				nIonSeriesList [*i1]->printHTMLHeader ( os, false );
			}
			tableEmptyNCells ( os, 3 );
			for ( SetIntConstIterator i2 = mainCTermIonIndicies.begin () ; i2 != mainCTermIonIndicies.end () ; i2++ ) {
				cIonSeriesList [*i2]->printHTMLHeader ( os, false );
			}
		tableRowEnd ( os );
		tableRowStart ( os );
			for ( SetIntConstIterator i3 = mainNTermIonIndicies.begin () ; i3 != mainNTermIonIndicies.end () ; i3++ ) {
				nIonSeriesList [*i3]->printHTMLEmptyCells ( os, false );
			}
			tableEmptyCell ( os );
			pInfo->printNTerm ( os );
			tableEmptyCell ( os );
			for ( SetIntConstIterator i4 = mainCTermIonIndicies.begin () ; i4 != mainCTermIonIndicies.end () ; i4++ ) {
				cIonSeriesList [*i4]->printHTMLEmptyCells ( os, false );
			}
		tableRowEnd ( os );
		int nrow = pInfo->getLength ();
		for ( int x = 0 ; x < nrow ; x++ ) {
			tableRowStart ( os );
				for ( SetIntConstIterator i5 = mainNTermIonIndicies.begin () ; i5 != mainNTermIonIndicies.end () ; i5++ ) {
					nIonSeriesList [*i5]->printHTML ( os, precision, spectrumMatch, sNum, x, false );
				}
				tableDataStart ( os );
					os << "<font size=\"-3\">" << x+1 << "</font>" << endl;
				tableDataEnd ( os );
				pInfo->printSequence ( os, x );
				tableDataStart ( os );
					os << "<font size=\"-3\">" << nrow-x << "</font>" << endl;
				tableDataEnd ( os );
				for ( SetIntConstIterator i6 = mainCTermIonIndicies.begin () ; i6 != mainCTermIonIndicies.end () ; i6++ ) {
					cIonSeriesList [*i6]->printHTML ( os, precision, spectrumMatch, sNum, nrow - x - 1, false );
				}
			tableRowEnd ( os );
		}
		tableRowStart ( os );
			for ( SetIntConstIterator i7 = mainNTermIonIndicies.begin () ; i7 != mainNTermIonIndicies.end () ; i7++ ) {
				nIonSeriesList [*i7]->printHTMLEmptyCells ( os, false );
			}
			tableEmptyCell ( os );
			pInfo->printCTerm ( os );
			tableEmptyCell ( os );
			for ( SetIntConstIterator i8 = mainCTermIonIndicies.begin () ; i8 != mainCTermIonIndicies.end () ; i8++ ) {
				cIonSeriesList [*i8]->printHTMLEmptyCells ( os, false );
			}
		tableRowEnd ( os );
	os << "</table>" << endl;
}
void SInfo::printHorizontalSequence ( ostream& os, const SpectrumMatch* spectrumMatch, int sNum ) const
{
	os << "<table cellspacing=\"4\" class=\"msprod_ion_table\">" << endl;
		if ( relatedIonSeries ) relatedIonSeries->printHTML ( os, precision, spectrumMatch, sNum );
		if ( mIonSeries ) mIonSeries->printHTML ( os, true, precision, spectrumMatch, sNum );
		tableEmptyRow ( os );

		tableRowStart ( os );
		tableHeaderStart ( os );
		os << "N-terminal" << endl;
		tableHeaderEnd ( os );
		tableRowEnd ( os );

		for ( IonSeriesListSizeType j = 0 ; j < nIonSeriesList.size () ; j++ ) {
			nIonSeriesList [j]->printHTML ( os, true, precision, spectrumMatch, sNum );
		}
		pInfo->printSequence ( os );
		tableRowStart ( os );
		tableHeaderStart ( os );
		os << "C-terminal" << endl;
		tableHeaderEnd ( os );
		tableRowEnd ( os );

		for ( IonSeriesListSizeType k = 0 ; k < cIonSeriesList.size () ; k++ ) {
			cIonSeriesList [k]->printHTML ( os, false, precision, spectrumMatch, sNum );
		}
	os << "</table>" << endl;
}
void SInfo::printInternalIons ( ostream& os, const SpectrumMatch* spectrumMatch, int sNum ) const
{
	if ( internalIonSeriesList.size () ) {
		ExpandableJavascriptBlock ejb2 ( "Internal Ions" );
		ejb2.printHeader ( os );
		os << "<table border=\"border\" cellspacing=\"3\" class=\"msprod_ion_table\">" << endl;
		tableRowStart ( os );
			tableHeader ( os, "Internal Sequence" );
			for ( InternalIonSeriesListSizeType j = 0 ; j < internalIonSeriesList.size () ; j++ ) {
				tableHeaderStart ( os );
					internalIonSeriesList [j]->printLabel ( os );
					os << endl;
				tableHeaderEnd ( os );
			}
		tableRowEnd ( os );
		for ( int i = 0 ; i < internalIonSeriesList [0]->getSize () ; i++ ) {
			int mc = internalIonSeriesList [0]->getMaxCharge ( i );
			for ( int charge = 1 ; charge <= mc ; charge++ ) {
				tableRowStart ( os );
					tableHeaderStart ( os, "", "left" );
						StringVector sv = internalIonSeriesList [0]->getSequence ( i );
						copy ( sv.begin (), sv.end (), ostream_iterator <string> ( os, "" ) );
						print_charge_superscript ( os, charge );
						os << endl;
					tableHeaderEnd ( os );
					for ( InternalIonSeriesListSizeType j = 0 ; j < internalIonSeriesList.size () ; j++ ) {
						bool hit = spectrumMatch ? spectrumMatch->isMatch ( internalIonSeriesList [j]->getMass ( i, charge ), charge, sNum ) : false;
						IonSeries::printIonMassHTML ( os, internalIonSeriesList [j]->getMass ( i, charge ), charge, precision, hit );
					}
				tableRowEnd ( os );
			}
		}
		os << "</table>" << endl;
		ejb2.printFooter ( os );
		os << "<br />" << endl;
	}
}
void SInfo::printXML ( ostream& os ) const
{
	if ( pInfo->pepLen && pInfo->show ) {
		os << "<theoretical_ions>" << endl;
			if ( relatedIonSeries ) relatedIonSeries->printXML ( os, precision );
			if ( mIonSeries ) mIonSeries->printXML ( os, true, precision );
			for ( IonSeriesListSizeType i = 0 ; i < nIonSeriesList.size () ; i++ ) {
				nIonSeriesList [i]->printXML ( os, true, precision );
			}
			for ( IonSeriesListSizeType j = 0 ; j < cIonSeriesList.size () ; j++ ) {
				cIonSeriesList [j]->printXML ( os, false, precision );
			}
		os << "</theoretical_ions>" << endl;
	}
}
BiemannFragments::BiemannFragments ( const vector <PepInfo*> pInfo, const BiemannParameters* bp, const vector <PeakContainer*> dataPeaks, const IntVector& maxReportedCharge, bool calibrate, double calTol, bool alternative, const LinkInfo* linkInfo ) :
	dataPeaks ( dataPeaks ),
	spectrumMatch ( 0 )
{
	PeptideInfo::setChargeReducedFragmentation ( instInf->getChargeReducedFragmentation () );
	int precision = instInf->getFragmentPeakPrecision ().getMassDecimalPlaces ();
	for ( std::vector <PepInfo*>::size_type i = 0 ; i < pInfo.size () ; i++ ) {
		sInfo.push_back ( new SInfo ( pInfo [i], bp, maxReportedCharge [i], precision, bp->getCountPosCharges (), linkInfo ) );
	}
	if ( !dataPeaks.empty () ) {
		double globalMinMass = 0.0;
		double globalMaxMass = 0.0;
		for ( int ii = 0 ; ii < dataPeaks.size () ; ii++ ) {
			double minMass = dataPeaks [ii]->getMinMassMinusTol ();
			double maxMass = dataPeaks [ii]->getMaxMassPlusTol ();
			globalMinMass = ( ii == 0 ) ? minMass : genMin ( minMass, globalMinMass );
			globalMaxMass = ( ii == 0 ) ? maxMass : genMax ( maxMass, globalMaxMass );
		}
		IonSeries::setMOverZRange ( globalMinMass, globalMaxMass );
		vector <TheoreticalSpectrum*> tsv;
		for ( std::vector <PepInfo*>::size_type j = 0 ; j < pInfo.size () ; j++ ) {
			tsv.push_back ( sInfo [j]->theoreticalSpectrum );
		}
		for ( int k = 0 ; k < dataPeaks.size () ; k++ ) {
			spectrumMatch.push_back ( new SpectrumMatch ( tsv, dataPeaks [k], bp == 0 ? false : bp->getAllowIncorrectCharge (), calibrate, calTol, alternative ) );
		}
	}
}
BiemannFragments::~BiemannFragments ()
{
	for ( int i = 0 ; i < spectrumMatch.size () ; i++ ) {
		delete spectrumMatch [i];
	}
	for ( SInfoPtrVectorSizeType j = 0 ; j < sInfo.size () ; j++ ) {
		delete sInfo [j];
	}
}
void BiemannFragments::printHTML ( ostream& os ) const
{
	if ( spectrumMatch.empty () ) {
		for ( SInfoPtrVectorSizeType j = 0 ; j < sInfo.size () ; j++ ) {
			sInfo [j]->printHTML ( os, 0, j );
		}
	}
	else {
		for ( int i = 0 ; i < spectrumMatch.size () ; i++ ) {
			for ( SInfoPtrVectorSizeType j = 0 ; j < sInfo.size () ; j++ ) {
				sInfo [j]->printHTML ( os, spectrumMatch [i], j );
			}
		}
	}
}
void BiemannFragments::printXML ( ostream& os, bool discriminating ) const
{
	if ( spectrumMatch.empty () ) {
		for ( SInfoPtrVectorSizeType i = 0 ; i < sInfo.size () ; i++ ) {
			sInfo [i]->printXML ( os );
		}
	}
	else {
		for ( int i = 0 ; i < spectrumMatch.size () ; i++ ) {
			for ( SInfoPtrVectorSizeType j = 0 ; j < sInfo.size () ; j++ ) {
				PeakMatchContext pmc ( dataPeaks [i]->getTolerance (), instInf->getFragmentPeakPrecision (), dataPeaks [i]->getNonUnitChargeData () );
				spectrumMatch [i]->printXML ( os, pmc, discriminating );
				sInfo [j]->printXML ( os );
			}
		}
	}
}
void BiemannFragments::printApplet ( ostream& os, const vector <XYData>& vXYData, bool plotCentroids, bool discriminating ) const
{
	for ( int i = 0 ; i < spectrumMatch.size () ; i++ ) {
		PeakMatchContext pmc ( dataPeaks [i]->getTolerance (), instInf->getFragmentPeakPrecision (), dataPeaks [i]->getNonUnitChargeData () );
		spectrumMatch [i]->printApplet ( os, pmc, vXYData, plotCentroids, discriminating );
		spectrumMatch [i]->printStats ( os );
	}
}
void BiemannFragments::printErrorChart ( ostream& os ) const
{
	for ( int i = 0 ; i < spectrumMatch.size () ; i++ ) {
		PeakMatchContext pmc ( dataPeaks [0]->getTolerance (), instInf->getFragmentPeakPrecision (), dataPeaks [0]->getNonUnitChargeData () );
		spectrumMatch [i]->printErrorChart ( os, pmc );
	}
}
void BiemannFragments::printSpectrumMatchHTML ( ostream& os, bool discriminating ) const
{
	for ( int i = 0 ; i < spectrumMatch.size () ; i++ ) {
		PeakMatchContext pmc ( dataPeaks [i]->getTolerance (), instInf->getFragmentPeakPrecision (), dataPeaks [i]->getNonUnitChargeData () );
		spectrumMatch [i]->printHTML ( os, pmc, discriminating );
	}
}
void BiemannFragments::printSpectrumMatchTabDelimitedText ( ostream& os, bool discriminating ) const
{
	for ( int i = 0 ; i < spectrumMatch.size () ; i++ ) {
		PeakMatchContext pmc ( dataPeaks [i]->getTolerance (), instInf->getFragmentPeakPrecision (), dataPeaks [i]->getNonUnitChargeData () );
		spectrumMatch [i]->printTabDelimitedText ( os, pmc, discriminating );
	}
}

bool TheoreticalSpectrum::xLink = false;
TheoreticalSpectrum::TheoreticalSpectrum ( int precision ) :
	precision ( precision )
{
}
void TheoreticalSpectrum::add ( const PIonSeriesList& pIonSeriesList )
{
	xLink = true;
	for ( int i = 0 ; i < pIonSeriesList.size () ; i++ ) {
		pIonSeriesList [i]->getPeaks ( peaks );
	}
}
void TheoreticalSpectrum::add ( const RelatedIonSeries* ionSeries )
{
	ionSeries->getPeaks ( peaks );
}
void TheoreticalSpectrum::add ( const IonSeries* ionSeries )
{
	ionSeries->getPeaks ( peaks );
}
void TheoreticalSpectrum::add ( const IonSeriesList& ionSeriesList )
{
	for ( IonSeriesListSizeType i = 0 ; i < ionSeriesList.size () ; i++ ) {
		ionSeriesList [i]->getPeaks ( peaks );
	}
}
void TheoreticalSpectrum::add ( const InternalIonSeriesList& ionSeriesList )
{
	for ( InternalIonSeriesListSizeType i = 0 ; i < ionSeriesList.size () ; i++ ) {
		ionSeriesList [i]->getPeaks ( peaks );
	}
}
void TheoreticalSpectrum::sortPeaks ()
{
	sort ( peaks.begin (), peaks.end (), SortFragmentPeaks () );
}
void TheoreticalSpectrum::printPeakTable ( ostream& os, int numColumns, const SpectrumMatch* sm, int sNum ) const
{
	int numPeaks = peaks.size ();
	int numRows = gen_round_up_int_divide ( numPeaks, numColumns );

	ExpandableJavascriptBlock ejb ( "Theoretical Peak Table" );
	ejb.printHeader ( os );
	os << "<table border=\"border\" cellspacing=\"3\" class=\"msprod_ion_table\">" << endl;
		for ( int i = 0 ; i < numRows ; i++ ) {
			tableRowStart ( os );
			for ( int j = 0 ; j < numColumns ; j++ ) {
				int k = i + ( j * numRows );
				if ( k < numPeaks ) {
					bool hit = sm ? sm->isMatch ( peaks [k].getMOverZ (), peaks [k].getCharge (), sNum ) : false;
					tableDataStart ( os );
						if ( hit ) os << "<font color=\"#FF0000\">";
							peaks [k].printMOverZHTML ( os, precision );
						if ( hit ) os << "</font>";
						os << endl;
					tableDataEnd ( os );
					tableHeaderStart ( os, "", "left" );
						if ( hit ) os << "<font color=\"#FF0000\">";
							peaks [k].printLabelHTML ( os );
						if ( hit ) os << "</font>";
						os << endl;
					tableHeaderEnd ( os );
				}
			}
			tableRowEnd ( os );
		}
	os << "</table>" << endl;
	ejb.printFooter ( os );
	os << "<br />" << endl;
}
FragmentPeakMatch::FragmentPeakMatch ( const Peak* dataPeak, const FragmentPeakVectorVector& hitInfo, int siz ) :
	dataPeak ( dataPeak ),
	hitInfo ( hitInfo ),
	siz ( siz ),
	multiplePeptides ( siz > 1 )
{
}
bool FragmentPeakMatch::match () const
{
	for ( FragmentPeakVectorVectorSizeType i = 0 ; i < hitInfo.size () ; i++ ) {
		if ( !hitInfo [i].empty () ) return true;
	}
	return false;
}
bool FragmentPeakMatch::match ( int sNum ) const
{
	return !hitInfo [sNum].empty ();
}
bool FragmentPeakMatch::isMatch ( double mass, int charge, int sNum ) const
{
	for ( FragmentPeakVectorSizeType j = 0 ; j < hitInfo [sNum].size ()  ; j++ ) {
		if ( hitInfo [sNum][j].isMatch ( mass, charge ) ) return true;
	}
	return false;
}
void FragmentPeakMatch::calibrate ( XYData& xydata, const PeakMatchContext& pmc ) const
{
	for ( FragmentPeakVectorSizeType i = 0 ; i < hitInfo.size () ; i++ ) {
		for ( FragmentPeakVectorSizeType j = 0 ; j < hitInfo [i].size () ; j++ ) {
			if ( hitInfo [i][j].getCategoryHTML () != 4 ) {
				double err = pmc.getError ( dataPeak->getMOverZ (), hitInfo [i][j].getMOverZ (), dataPeak->getCharge () );
				xydata.add ( dataPeak->getMOverZ (), err );
			}
		}
	}
}
pair <MapStringToInt, bool> FragmentPeakMatch::getLabelMap () const
{
	MapStringToInt msi;
	bool flag = false;
	if ( multiplePeptides ) {
		for ( FragmentPeakVectorSizeType i = 0 ; i < hitInfo.size () ; i++ ) {		// For each alterative peptide
			for ( FragmentPeakVectorSizeType j = 0 ; j < hitInfo [i].size () ; j++ ) {
				string lab = hitInfo [i][j].getLabelHTML ();
				MapStringToIntConstIterator cur = msi.find ( lab );
				if ( cur == msi.end () )	msi [lab] = 1;
				else						msi [lab] += 1;
			}
		}
		for ( MapStringToIntConstIterator k = msi.begin () ; k != msi.end () ; k++ ) {
			if ( (*k).second != siz ) flag = true;
		}
	}
	return make_pair ( msi, flag );
}
void FragmentPeakMatch::printApplet ( LabelledCatagorizedGraphData& graphData, const PeakMatchContext& pmc, bool discriminating, double normalizationFactor ) const
{
	PeakPrecision pp = pmc.getPeakPrecision ();
	double measuredMass = dataPeak->getMass ();
	stringstream label;
	stringstream category;
	int color = 0;
	pair <MapStringToInt, bool> pmsib = getLabelMap ();
	MapStringToInt& msi = pmsib.first;
	for ( FragmentPeakVectorSizeType i = hitInfo.size () ; i-- ; ) {		// For each alterative peptide
		if ( hitInfo [i].size () ) color = 8;
		for ( FragmentPeakVectorSizeType j = 0 ; j < hitInfo [i].size () ; j++ ) {	// For each matching ion
			string lab = hitInfo [i][j].getLabelHTML ();
			int cat = hitInfo [i][j].getCategoryHTML ();
			if ( multiplePeptides ) {
				MapStringToIntConstIterator cur = msi.find ( lab );
				if ( cur != msi.end () ) {
					if ( (*cur).second == siz ) {
						if ( !discriminating || pmsib.second ) {
							label << lab << '|';
							category << cat << '|';
						}
						msi.erase ( lab );
					}
					else {
						label << lab << '|';
						category << cat << '.' << i+1 << '|';
					}
				}
			}
			else { 
				label << lab << '|';
				category << cat << '|';
			}
		}
	}
	if ( color == 0 ) {
		genPrint ( label, dataPeak->getMOverZ (), pp.getMassDecimalPlaces () );
		label << '|';
		category << MASS_CATEGORY << "|";
	}
	double intensity = dataPeak->getIntensity ();
	if ( normalizationFactor > 0.0 ) intensity *= normalizationFactor;
	graphData.add ( dataPeak->getMOverZ (), intensity, label.str (), color, category.str () );
}
void FragmentPeakMatch::printErrorChart ( ostream& os, const PeakMatchContext& pmc, int sNum ) const
{
	for ( FragmentPeakVectorSizeType j = 0 ; j < hitInfo [sNum].size () ; j++ ) {
		double err = pmc.getError ( dataPeak->getMOverZ (), hitInfo [sNum][j].getMOverZ (), dataPeak->getCharge () );
		genPrint ( os, dataPeak->getMOverZ (), 4 );
		os << " ";
		genPrintSigFig ( os, err, 2 );
		os << " ";
		os << hitInfo [sNum][j].getCategoryHTML ();
		os << endl;
	}
}
void FragmentPeakMatch::printErrorChart ( ostream& os, const PeakMatchContext& pmc ) const
{
	for ( FragmentPeakVectorSizeType i = 0 ; i < hitInfo.size () ; i++ ) {
		for ( FragmentPeakVectorSizeType j = 0 ; j < hitInfo [i].size () ; j++ ) {
			double err = pmc.getError ( dataPeak->getMOverZ (), hitInfo [i][j].getMOverZ (), dataPeak->getCharge () );
			genPrint ( os, dataPeak->getMOverZ (), 4 );
			os << " ";
			genPrintSigFig ( os, err, 2 );
			os << " ";
			os << hitInfo [i][j].getCategoryHTML ();
			os << endl;
		}
	}
}
void FragmentPeakMatch::printMOverZHTML ( ostream& os, const PeakMatchContext& pmc ) const
{
	tableHeaderStart ( os );
		dataPeak->printMOverZHTML ( os, pmc.getPrecision () );
		os << endl;
	tableHeaderEnd ( os );
}
string FragmentPeakMatch::getMOverZHTML ( const PeakMatchContext& pmc ) const
{
	ostringstream ostr;
	printMOverZHTML ( ostr, pmc );
	return ostr.str ();
}
void FragmentPeakMatch::printMassHTML ( ostream& os, const PeakMatchContext& pmc ) const
{
	tableDataStart ( os );
		genPrint ( os, dataPeak->getMass (), pmc.getPrecision () );
		os << endl;
	tableDataEnd ( os );
}
string FragmentPeakMatch::getMassHTML ( const PeakMatchContext& pmc ) const
{
	ostringstream ostr;
	printMassHTML ( ostr, pmc );
	return ostr.str ();
}
string FragmentPeakMatch::getMatchHTML ( const PeakMatchContext& pmc, bool discriminating ) const
{
	ostringstream ostr;
	pair <MapStringToInt, bool> pmsib = getLabelMap ();
	MapStringToInt& msi = pmsib.first;
	for ( FragmentPeakVectorVectorSizeType i = 0 ; i < hitInfo.size () ; i++ ) {
		for ( FragmentPeakVectorSizeType j = 0 ; j < hitInfo [i].size () ; j++ ) {
			const FragmentPeak& hi = hitInfo [i][j];
			if ( multiplePeptides ) {
				MapStringToIntIterator cur = msi.find ( hi.getLabelHTML () );
				if ( cur != msi.end () ) {
					if ( (*cur).second == siz ) {
						if ( !discriminating || pmsib.second ) printColoredLabel ( ostr, pmc, hi, -1 );
						msi.erase ( cur );
					}
					else printColoredLabel ( ostr, pmc, hi, i );
				}
			}
			else printColoredLabel ( ostr, pmc, hi, 4 );
		}
	}
	return ostr.str ();
}
string FragmentPeakMatch::getColor ( int col )
{
	if ( col == 0 ) return "#007F00";	//green
	if ( col == 1 ) return "#0000FF";	//blue
	if ( col == 2 ) return "#FF00FF";	//magenta
	if ( col == 3 ) return "#00FFFF";	//cyan
	if ( col == 4 ) return "#FF0000";	//red
	if ( col == 5 ) return "#7F7F7F";	//grey
	return "#000000";					//black
}
void FragmentPeakMatch::printColoredLabel ( ostream& os, const PeakMatchContext& pmc, const FragmentPeak& hi, int col ) const
{
	os << "<font color=\"" << getColor ( col ) << "\">";
	printLabel ( os, pmc, hi );
	os << "</font>";
}
void FragmentPeakMatch::printLabel ( ostream& os, const PeakMatchContext& pmc, const FragmentPeak& hi ) const
{
	double err = pmc.getError ( dataPeak->getMOverZ (), hi.getMOverZ (), dataPeak->getCharge () );
	hi.printLabelHTML ( os );
	os << "(";
	genPrintSigFig ( os, err, 2 );
	os << ")";
	os << "<br />";
}
void FragmentPeakMatch::printTabDelimitedText ( ostream& os, const PeakMatchContext& pmc, bool discriminating ) const
{
	bool printed = false;
	pair <MapStringToInt, bool> pmsib = getLabelMap ();
	MapStringToInt& msi = pmsib.first;
	for ( FragmentPeakVectorVectorSizeType i = 0 ; i < hitInfo.size () ; i++ ) {
		for ( FragmentPeakVectorSizeType j = 0 ; j < hitInfo [i].size () ; j++ ) {
			const FragmentPeak& hi = hitInfo [i][j];
			double err = pmc.getError ( dataPeak->getMOverZ (), hi.getMOverZ (), dataPeak->getCharge () );
			if ( multiplePeptides ) {
				MapStringToIntIterator cur = msi.find ( hi.getLabelHTML () );
				if ( cur != msi.end () ) {
					if ( (*cur).second == siz ) {
						if ( !discriminating || pmsib.second ) {
							printDelimitedLine ( os, pmc, !printed, "All", hi.getLabelXML (), err );
							printed = true;
						}
						msi.erase ( cur );
					}
					else {
						printDelimitedLine ( os, pmc, !printed, gen_itoa ( (int)i+1 ), hi.getLabelXML (), err );
						printed = true;
					}
				}
			}
			else {
				printDelimitedLine ( os, pmc, !printed, "", hi.getLabelXML (), err );
				printed = true;
			}
		}
	}
	if ( !printed && !discriminating ) printDelimitedLine ( os, pmc );
}
void FragmentPeakMatch::printDelimitedLine ( ostream& os, const PeakMatchContext& pmc, bool printPk, const string& idx, const string& label, double err ) const
{
	delimitedRowStart ( os );

	if ( printPk )	dataPeak->printDelimited ( os, pmc.getPeakPrecision () );
	else			delimitedEmptyNCells ( os, 3 );
	if ( multiplePeptides ) {
		if ( !idx.empty () )	delimitedCell ( os, idx );
		else					delimitedEmptyCell ( os );
	}
	if ( !label.empty () ) {
		delimitedCell ( os, label );
		delimitedCellSigFig ( os, err, pmc.getErrorSigFig () );
	}
	else delimitedEmptyNCells ( os, 2 );

	delimitedRowEnd ( os );
}
void FragmentPeakMatch::printXML ( ostream& os, const PeakMatchContext& pmc, bool discriminating ) const
{
	bool print = true;
	if ( multiplePeptides && discriminating ) print = false;
	pair <MapStringToInt, bool> pmsib = getLabelMap ();
	MapStringToInt& msi = pmsib.first;
	ostringstream ostr;
	ostr << "<peak>" << endl;

	const PeakPrecision& pp = pmc.getPeakPrecision ();
	dataPeak->printXML ( ostr, pp );

	for ( FragmentPeakVectorVectorSizeType i = 0 ; i < hitInfo.size () ; i++ ) {
		for ( FragmentPeakVectorSizeType j = 0 ; j < hitInfo [i].size () ; j++ ) {
			const FragmentPeak& hi = hitInfo [i][j];
			double err = pmc.getError ( dataPeak->getMOverZ (), hi.getMOverZ (), dataPeak->getCharge () );

			if ( multiplePeptides ) {
				MapStringToIntIterator cur = msi.find ( hi.getLabelHTML () );
				if ( cur != msi.end () ) {
					if ( (*cur).second == siz ) {
						if ( !discriminating || pmsib.second ) {
							ostr << "<match" << j+1 << "_All>";
							printXMLMatchInfo ( ostr, hi, err, pmc.getErrorSigFig (), dataPeak->getIntensity (), pp.getIntensitySigFig () );
							ostr << "</match" << j+1 << "_All>" << endl;
							print = true;
						}
						msi.erase ( cur );
					}
					else {
						ostr << "<match" << j+1 << "_" << i+1 << ">";
						printXMLMatchInfo ( ostr, hi, err, pmc.getErrorSigFig (), dataPeak->getIntensity (), pp.getIntensitySigFig () );
						ostr << "</match" << j+1 << "_" << i+1 << ">" << endl;
						print = true;
					}
				}
			}
			else {
				ostr << "<match" << j+1 << ">";
				printXMLMatchInfo ( ostr, hi, err, pmc.getErrorSigFig (), dataPeak->getIntensity (), pp.getIntensitySigFig () );
				ostr << "</match" << j+1 << ">" << endl;
			}
		}
	}
	ostr << "</peak>" << endl;
	if ( print ) os << ostr.str ();
}
void FragmentPeakMatch::printXMLMatchInfo ( ostream& os, const FragmentPeak& hi, double err, int errSF, double inten, int intenSF ) const
{
	hi.printLabelXML ( os );
	os << ",";
	genPrintSigFig ( os, err, errSF );
	os << ",";
	genPrintSigFig ( os, inten, intenSF );
}
SpectrumMatch::SpectrumMatch ( const vector <TheoreticalSpectrum*> tsv, const PeakContainer* dataPeaks, bool allowIncorrectCharge, bool calFlag, double calTol, bool alternative ) :
	sequenceFlag ( false ),
	numPeptides ( tsv.size () ),
	maxIntensity ( 0 )
{
	initSetList ( tsv );
	PeakMatchContext pmc ( dataPeaks->getTolerance (), instInf->getFragmentPeakPrecision (), dataPeaks->getNonUnitChargeData () );
	init ( tsv, dataPeaks, allowIncorrectCharge, pmc );
	if ( alternative )	calculateMatchStatsAlternate ();
	else				calculateMatchStats ();
	if ( calFlag ) {
		double gradient;
		double offset;
		calibrate ( pmc, gradient, offset );
		const_cast <PeakContainer*> (dataPeaks)->calibrate ( gradient, offset );
		const_cast <PeakContainer*> (dataPeaks)->setToleranceValue ( calTol );
		matches.clear ();
		init ( tsv, dataPeaks, allowIncorrectCharge, pmc );
		if ( alternative )	calculateMatchStatsAlternate ();
		else				calculateMatchStats ();
	}
}
void SpectrumMatch::initSetList ( const vector <TheoreticalSpectrum*> tsv )
{
	for ( int i = 0 ; i < tsv.size () ; i++ ) {
		set.push_back ( tsv [i] != 0 );
	}
}
class SortHitInfoByAbsoluteError {
	double dataMZ;
	int dataZ;
	const PeakMatchContext& pmc;
	public:
		SortHitInfoByAbsoluteError ( double dataMZ, int dataZ, const PeakMatchContext& pmc ) :
			dataMZ ( dataMZ ),
			dataZ ( dataZ ),
			pmc ( pmc )
		{
		}
		int operator () ( const FragmentPeak& a, const FragmentPeak& b ) const
		{
			double errA = pmc.getError ( dataMZ, a.getMOverZ (), dataZ );
			double errB = pmc.getError ( dataMZ, b.getMOverZ (), dataZ );

			return fabs ( errA ) < fabs ( errB );
		}
};
void SpectrumMatch::init ( const vector <TheoreticalSpectrum*> tsv, const PeakContainer* dataPeaks, bool allowIncorrectCharge, const PeakMatchContext& pmc )
{
	for ( int m = 0 ; m < tsv.size () ; m++ ) {
		if ( tsv [m] ) sequenceFlag |= !tsv [m]->peaks.empty ();
	}
	IntVector j ( tsv.size (), 0 );
	for ( int i = 0 ; i < dataPeaks->size () ; i++ ) {
		const Peak* dp = (*dataPeaks) [i];
		vector <vector <FragmentPeak> > hitInfo ( tsv.size () );
		int siz = 0;
		for ( int m = 0 ; m < tsv.size () ; m++ ) {
			if ( tsv [m] ) {
				const FragmentPeakVector& peaks = tsv [m]->peaks;
				while ( j [m] < peaks.size () && !dp->isLowerMatch ( peaks [j[m]].getMass () ) ) {
					j [m] += 1;
				}
				for ( FragmentPeakVectorSizeType k = j [m] ; k < peaks.size () ; k++ ) {
					const FragmentPeak& fp = peaks [k];
					if ( allowIncorrectCharge && dp->getCharge () == 1 ) {
						if ( dp->isMatch ( fp.getMOverZ () ) ) {	// This code will fit peaks with incorrect charge. no good for highly charged (eg ECD) spectra
							hitInfo [m].push_back ( fp );
						}
					}
					else {
						if ( dp->isMatch ( fp.getMass () ) && dp->getCharge () == fp.getCharge () ) {
							hitInfo [m].push_back ( fp );
						}
					}
				}
				siz++;
			}
		}
		for ( int ii = 0 ; ii < hitInfo.size () ; ii++ ) {
			stable_sort ( hitInfo [ii].begin (), hitInfo [ii].end (), SortHitInfoByAbsoluteError ( dp->getMOverZ (), dp->getCharge (), pmc ) );
		}
		matches.push_back ( FragmentPeakMatch ( dp, hitInfo, siz ) );
		maxIntensity = genMax ( maxIntensity, matches.back ().getIntensity () );
	}
}
void SpectrumMatch::calculateMatchStatsAlternate ()
{
	numMatchedArray.resize ( numPeptides );
	percentMatchedIntensityArray.resize ( numPeptides );
	percentMatchedIntensityArray2.resize ( numPeptides );
	for ( int i = 0 ; i < numPeptides ; i++ ) {
		numMatchedArray [i] = 0;
		double totalIntensity = 0.0;
		double totalMatchedIntensity = 0.0;
		for ( int j = 0 ; j < matches.size () ; j++ ) {
			FragmentPeakMatch& fpm = matches [j];
			double intensity = fpm.getIntensity ();
			totalIntensity += intensity;
			if ( fpm.match ( i ) ) {
				numMatchedArray [i] += 1;
				totalMatchedIntensity += intensity;
			}
		}
		percentMatchedIntensityArray [i] = 100.0 * totalMatchedIntensity / totalIntensity;

		double totalIntensity2 = 0.0;
		double totalMatchedIntensity2 = 0.0;
		for ( int k = 0 ; k < matches.size () ; k++ ) {
			FragmentPeakMatch& fpm = matches [k];
			if ( fpm.getMOverZ () > 170 && !fpm.isPrecursorIon () ) {
				double intensity = fpm.getIntensity ();
				totalIntensity2 += intensity;
				if ( fpm.match ( i ) ) {
					totalMatchedIntensity2 += intensity;
				}
			}
		}
		percentMatchedIntensityArray2 [i] = 100.0 * totalMatchedIntensity2 / totalIntensity2;
	}
}
void SpectrumMatch::calculateMatchStats ()
{
	numMatched = 0;
	double totalIntensity = 0.0;
	double totalMatchedIntensity = 0.0;
	for ( int i = 0 ; i < matches.size () ; i++ ) {
		FragmentPeakMatch& fpm = matches [i];
		double intensity = fpm.getIntensity ();
		totalIntensity += intensity;
		if ( fpm.match () ) {
			numMatched++;
			totalMatchedIntensity += intensity;
		}
	}
	percentMatchedIntensity = 100.0 * totalMatchedIntensity / totalIntensity;

	double totalIntensity2 = 0.0;
	double totalMatchedIntensity2 = 0.0;
	for ( int j = 0 ; j < matches.size () ; j++ ) {
		FragmentPeakMatch& fpm = matches [j];
		if ( fpm.getMOverZ () > 170 && !fpm.isPrecursorIon () ) {
			double intensity = fpm.getIntensity ();
			totalIntensity2 += intensity;
			if ( fpm.match () ) {
				totalMatchedIntensity2 += intensity;
			}
		}
	}
	percentMatchedIntensity2 = 100.0 * totalMatchedIntensity2 / totalIntensity2;
}
void SpectrumMatch::printApplet ( ostream& os, const PeakMatchContext& pmc, const vector <XYData>& vXYData, bool plotCentroids, bool discriminating ) const
{
	bool xLink = TheoreticalSpectrum::xLink;
	double heightMultiplier = 1.0;
	StringVector categories;
	categories.push_back ( "Mass" );
	categories.push_back ( "Ser" );
	categories.push_back ( "Loss" );
	categories.push_back ( "Imm" );
	categories.push_back ( "Int" );
	if ( xLink ) categories.push_back ( "XL" );
	BoolDeque showCategory ( categories.size (), false );
	showCategory [1] = true;
	if ( xLink ) showCategory [5] = true;
	if ( !vXYData.empty () ) {									// Raw data also plotted
		double xMin = vXYData [0].minX ();
		double xMax = vXYData [0].maxX ();
		for ( int i = 1 ; i < vXYData.size () ; i++ ) {
			xMin = genMin ( xMin, vXYData [i].minX () );
			xMax = genMax ( xMax, vXYData [i].maxX () );
		}
		for ( int j = 0 ; j < vXYData.size () ; j++ ) {
			LabelledCatagorizedFileGraphData graphData ( categories, showCategory, vXYData [j], true );
			if ( plotCentroids ) {
				double maxRawIntensity = vXYData [j].maxY ();
				double normalizationFactor = maxRawIntensity / maxIntensity;
				for ( int k = 0 ; k < matches.size () ; k++ ) {
					matches [k].printApplet ( graphData, pmc, discriminating, normalizationFactor );
				}
			}
			SpectrumGraph s ( "pr_graph.par.txt" );
			s.drawGraph ( os, graphData, false, xMin, xMax, heightMultiplier );
			os << "<br />" << endl;
			ParameterList::printDoubleHTMLSigFig ( os, "Max Intensity", graphData.maxY (), 3 );
		}
	}
	else {
		LabelledCatagorizedGraphData graphData ( categories, showCategory );

		for ( int i = 0 ; i < matches.size () ; i++ ) {
			matches [i].printApplet ( graphData, pmc, discriminating );
		}
		if ( SpectrumGraph::getDrawGraph () ) {
			SpectrumGraph s ( "pr_graph.par.txt" );
			s.drawGraph ( os, graphData, false, 0.0, 0.0, heightMultiplier );
		}
		os << "<br />" << endl;
		ParameterList::printDoubleHTMLSigFig ( os, "Max Intensity", graphData.maxY (), 3 );
	}
}
void SpectrumMatch::printStats ( ostream& os ) const
{
	if ( sequenceFlag ) {
		if ( !numMatchedArray.empty () ) {
			for ( IntVectorSizeType i = 0 ; i < numMatchedArray.size () ; i++ ) {
				if ( set [i] ) printStats ( os, numMatchedArray [i], percentMatchedIntensityArray [i], percentMatchedIntensityArray2 [i] );
			}
		}
		else
			printStats ( os, numMatched, percentMatchedIntensity, percentMatchedIntensity2 );
	}
	else {
		ParameterList::printHTML ( os, "Num Peaks", matches.size () );
		os << "<br />" << endl;
	}
}
void SpectrumMatch::printStats ( ostream& os, int nm, double pmi, double pmi2 ) const
{
	os << "Num Matched: <b>" << nm << "/" << matches.size ();
	os << " (";
	genPrint ( os, ( matches.size () - nm ) * 100.0 / matches.size (), 1 );
	os << "% unmatched)</b>";
	os << " ";
	os << "Matched Intensity: <b>";
	genPrint ( os, pmi, 1 );
	os << "%</b>";
	os << " ";
	os << "Matched Series Intensity: <b>";
	genPrint ( os, pmi2, 1 );
	os << "%</b>";
	os << "<br />" << endl;

//os << "Peptide Score: <b>";
//genPrint ( os, - 10.0 * log10 ( gen_binomial_probability ( matches.size (), 6.0 / 100.0, nm ) ), 1 );
//os << "</b>";
//os << "<br />" << endl;
}
void SpectrumMatch::printErrorChart ( ostream& os, const PeakMatchContext& pmc ) const
{
	if ( sequenceFlag ) {
		if ( RPlot::getRFlag () ) {
			if ( !numMatchedArray.empty () ) {
				for ( IntVectorSizeType i = 0 ; i < numMatchedArray.size () ; i++ ) {
					RPlot rplot ( "msmsErrorScatter.R" );
					GenOFStream ofs ( rplot.getDataFileFullPath () );
					for ( FragmentPeakMatchVectorSizeType j = 0 ; j < matches.size () ; j++ ) {
						matches [j].printErrorChart ( ofs, pmc, i );
					}
					ofs.close ();
					rplot.printImageAndLink ( os );
				}
			}
			else {
				RPlot rplot ( "msmsErrorScatter.R" );
				GenOFStream ofs ( rplot.getDataFileFullPath () );
				for ( FragmentPeakMatchVectorSizeType i = 0 ; i < matches.size () ; i++ ) {
					matches [i].printErrorChart ( ofs, pmc );
				}
				ofs.close ();
				rplot.printImageAndLink ( os );
			}
		}
	}
}
void SpectrumMatch::calibrate ( const PeakMatchContext& pmc, double& gradient, double& offset ) const
{
	XYData xydata;
	for ( FragmentPeakMatchVectorSizeType i = 0 ; i < matches.size () ; i++ ) {
		matches [i].calibrate ( xydata, pmc );
	}
	xydata.linearRegression ( &offset, &gradient );
}
bool SpectrumMatch::isMatch ( double mass, int charge, int sNum ) const
{
	for ( FragmentPeakMatchVectorSizeType i = 0 ; i < matches.size () ; i++ ) {
		if ( matches [i].isMatch ( mass, charge, sNum ) ) return true;
	}
	return false;
}
void SpectrumMatch::printHTML ( ostream& os, const PeakMatchContext& pmc, bool discriminating ) const
{
	if ( sequenceFlag ) {
		int numPeaks = matches.size ();
		StringVector mzInfo;
		StringVector massInfo;
		StringVector matchInfo;
		for ( int ii = 0 ; ii < numPeaks ; ii++ ) {
			const string& mat = matches [ii].getMatchHTML ( pmc, discriminating );
			if ( !discriminating || !mat.empty () ) {
				mzInfo.push_back ( matches [ii].getMOverZHTML ( pmc ) );
				if ( pmc.getNonUnitChargeData () ) massInfo.push_back ( matches [ii].getMassHTML ( pmc ) );
				matchInfo.push_back ( mat );
			}
		}
		int numColumns = 10;
		int numCells = mzInfo.size ();
		int numRows = gen_round_up_int_divide ( numCells, numColumns );
		int j;

		ExpandableJavascriptBlock ejb ( "Peak Matches" );
		ejb.printHeader ( os );
		os << "<table border=\"border\" cellspacing=\"0\" class=\"msprod_ion_table\">" << endl;
			for ( int i = 0 ; i < numRows ; i++ ) {
				int min_j = i * numColumns;
				int max_j = genMin ( min_j + numColumns, numCells );
				tableRowStart ( os );
					for ( j = min_j ; j < max_j ; j++ )		os << mzInfo [j];
					os << endl;
				tableRowEnd ( os );
				if ( !massInfo.empty () ) {
					tableRowStart ( os );
						for ( j = min_j ; j < max_j ; j++ )	os << massInfo [j];
						os << endl;
					tableRowEnd ( os );
				}
				tableRowStart ( os );
					for ( j = min_j ; j < max_j ; j++ )		tableCell ( os, matchInfo [j] );
					os << endl;
				tableRowEnd ( os );
			}
		os << "</table>" << endl;
		os << "<br />" << endl;
		ejb.printFooter ( os );
		os << "<br />" << endl;
	}
}
void SpectrumMatch::printTabDelimitedText ( ostream& os, const PeakMatchContext& pmc, bool discriminating ) const
{
	delimitedRowStart ( os );
		delimitedHeader ( os, "m/z" );
		delimitedHeader ( os, "z" );
		delimitedHeader ( os, "Intensity" );
		if ( numPeptides > 1 ) delimitedHeader ( os, "Peptide #" );
		delimitedHeader ( os, "Ion Type" );
		delimitedHeader ( os, "Error" );
	delimitedRowEnd ( os );
	for ( FragmentPeakMatchVectorSizeType i = 0 ; i < matches.size () ; i++ ) {
		matches [i].printTabDelimitedText ( os, pmc, discriminating );
	}
}
void SpectrumMatch::printXML ( ostream& os, const PeakMatchContext& pmc, bool discriminating ) const
{
	os << "<ion_matches>" << endl;
	for ( FragmentPeakMatchVectorSizeType i = 0 ; i < matches.size () ; i++ ) {
		matches [i].printXML ( os, pmc, discriminating );
	}
	os << "</ion_matches>" << endl;
}
B1Info::B1Info ()
{
	GenIFStream fromFile ( MsparamsDir::instance ().getParamPath ( "b1.txt" ) );
	string line;
	while ( getline ( fromFile, line ) ) {
		if ( line.length () != 0 && line [0] != '#' ) {
			b1NTermSet.insert ( gen_strtrim ( line ) );
		}
	}
}
B1Info& B1Info::instance ()
{
	static B1Info b;
	return b;
}
bool B1Info::isB1 ( const string& s )
{
	return b1NTermSet.find ( s ) != b1NTermSet.end ();
}
