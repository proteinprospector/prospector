/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_quan_ratio.cpp                                             *
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
#ifdef RAW_DATA
#ifdef VIS_C
#pragma warning( disable : 4503 )	// Disable warnings about decorated name length being too long
#endif
#include <lg_string.h>
#include <lgen_math.h>
#include <ln_dot.h>
#include <lu_quan_ratio.h>
#include <lu_aa_calc.h>
#include <lu_html.h>
#include <lu_iso.h>
#include <lu_getfil.h>
#include <lu_get_link.h>
#include <lu_table.h>
#include <lu_delim.h>
#include <lu_mass_conv.h>
#include <lu_iso_par.h>
#include <lu_app_gr.h>
#include <lu_usermod.h>
using std::map;
using std::string;
using std::ostream;
using std::getline;
using std::vector;
using std::stable_sort;
using std::endl;
using std::runtime_error;
using std::ostringstream;
using nr::sqrtTwo;
using nr::sqrtPi;

LinkInfo* QuanPeptide::linkInfo = 0;

QuanPeptide::QuanPeptide ( const string& peptide1, const string& nterm1, const string& cterm1, const string& nloss, const string& peptide2, const string& nterm2, const string& cterm2, const string& linkName ) :
	peptide1 ( peptide1 ),
	nterm1 ( nterm1 ),
	cterm1 ( cterm1 ),
	nloss ( nloss ),
	peptide2 ( peptide2 ),
	nterm2 ( nterm2 ),
	cterm2 ( cterm2 ),
	linkName ( linkName )
{
}
string QuanPeptide::getBridgeFormula () const
{
	if ( linkInfo )	return linkInfo->getBridgeFormula ();
	else			return "";
}
void QuanPeptide::setLinkInfo ( LinkInfo* li )
{
	linkInfo = li;
}
CosSimilarity::CosSimilarity ( const string& formulaString, int charge, const DoubleVector& intensity, const DoubleVector& area ) :
	cosSimilarityIntensity ( -100.0 ),
	cosSimilarityArea ( -100.0 )
{
	IsotopePeakStats id ( formulaString, charge );
	int np = id.getNumPeaks ();
	int ind1 = id.getProbabilityIndexForIdealMonoisotopicMZ ();
	int distSize = genMax ( intensity.size (), area.size () );
	DoubleVector dv;
	for ( int i = ind1 ; i < np ; i++ ) {
		dv.push_back ( id.getProbability ( i ) );
		if ( dv.size () == distSize ) break;
	}
	if ( dv.size () > 1 ) {
		if ( intensity.size () > 1 )	cosSimilarityIntensity = cosSimilarity ( intensity, dv );
		if ( area.size () > 1 )			cosSimilarityArea = cosSimilarity ( area, dv );
	}
}
void CosSimilarity::printHTML ( ostream& os )
{
	os << "<b>Intensity Cosine Similarity:</b> ";
	genPrintSigFig ( os, cosSimilarityIntensity, 3 );
	os << "<br />" << endl;

	os << "<b>Area Cosine Similarity:</b> ";
	genPrintSigFig ( os, cosSimilarityArea, 3 );
	os << "<br />" << endl;
}
void CosSimilarity::printHTML ( ostream& os, int num )
{
	os << "<b>Intensity Cosine Similarity Pk " << num << ": </b>";
	genPrintSigFig ( os, cosSimilarityIntensity, 3 );
	os << "<br />" << endl;

	os << "<b>Area Cosine Similarity Pk " << num << ": </b>";
	genPrintSigFig ( os, cosSimilarityArea, 3 );
	os << "<br />" << endl;
}
CosSimilarityList::CosSimilarityList ( const StringVector& formulaString, int charge, const DoubleVectorVector& intensity, const DoubleVectorVector& area, int numQuanStates ) :
	csi ( numQuanStates, -100.0 ),
	csa ( numQuanStates, -100.0 )
{
	for ( int i = 0 ; i < numQuanStates ; i++ ) {
		CosSimilarity cs ( formulaString [i], charge, i < intensity.size () ? intensity [i] : DoubleVector (), i < area.size () ? area [i] : DoubleVector () );
		csi [i] = cs.getIntensity ();
		csa [i] = cs.getArea ();
	}
}
void CosSimilarityList::printHTMLIntensity ( ostream& os, const string& styleID ) const
{
	printHTML ( os, csi, styleID );
}
void CosSimilarityList::printHTMLArea ( ostream& os, const string& styleID ) const
{
	printHTML ( os, csa, styleID );
}
void CosSimilarityList::printHTML ( ostream& os, const DoubleVector& c, const string& styleID ) const
{
	for ( int i = 0 ; i < c.size () ; i++ ) {
		if ( c [i] == -100.0 )	tableEmptyCell ( os, styleID );
		else					tableCellSigFig ( os, c [i], 3, false, styleID );
	}
}
void CosSimilarityList::printDelimitedIntensity ( ostream& os ) const
{
	printDelimited ( os, csi );
}
void CosSimilarityList::printDelimitedArea ( ostream& os ) const
{
	printDelimited ( os, csa );
}
void CosSimilarityList::printDelimited ( ostream& os, const DoubleVector& c ) const
{
	for ( int i = 0 ; i < c.size () ; i++ ) {
		if ( c [i] == -100.0 )	delimitedEmptyCell ( os );
		else					delimitedCellSigFig ( os, c [i], 3 );
	}
}
class MSPurityCorrection {
	int numQuanPeaks;
	int numIsotopePeaks;
	int numCoefficients;
	DoubleVectorVector matrix;
	double** a;
	double** b;
public:
	MSPurityCorrection ( const vector <IsotopePeakStats*>& ips, int offset );
	~MSPurityCorrection ();
	void loadA () const;
	void correction ( DoubleVector& dv ) const;
};
MSPurityCorrection::MSPurityCorrection ( const vector <IsotopePeakStats*>& ips, int offset ) :
	numQuanPeaks ( ips.size () ),
	matrix ( numQuanPeaks )
{
	IntVector pIndex;
	for ( int i = 0 ; i < numQuanPeaks ; i++ ) {
		const IsotopePeakStats* iPkSt1 = ips [i];
		int ind1 = iPkSt1->getProbabilityIndexForIdealMonoisotopicMZ () + offset;
		for ( int j = 0 ; j < numQuanPeaks ; j++ ) {
			double prob = 0.0;
			if ( i == j )
				prob = iPkSt1->getProbability ( ind1 );
			else {
				const IsotopePeakStats* iPkSt2 = ips [j];
				int indimm = iPkSt2->getProbabilityIndexForIdealMonoisotopicMZ ();
				int integerDiff = floor ( ( iPkSt2->getIdealMonoisotopicMass () - iPkSt1->getIdealMonoisotopicMass () ) + 0.5 );
				int ind2 = indimm - integerDiff + offset;
				if ( ind2 >= 0 && ind2 < iPkSt2->size () )
					prob = iPkSt2->getProbability ( ind2 );
				else
					prob = 0.0;
			}
			matrix [i].push_back ( prob );
		}
	}
	a = nrmatrix ( 1, numQuanPeaks, 1, numQuanPeaks );
	loadA ();
	b = nrmatrix ( 1, numQuanPeaks, 1, 1 );
}
MSPurityCorrection::~MSPurityCorrection ()
{
	if ( numQuanPeaks ) {
		free_nrmatrix ( a, 1, numQuanPeaks, 1, numQuanPeaks );
		free_nrmatrix ( b, 1, numQuanPeaks, 1, 1 );
	}
}
void MSPurityCorrection::loadA () const
{
	for ( int i = 1 ; i <= numQuanPeaks ; i++ ) {
		for ( int j = 1 ; j <= numQuanPeaks ; j++ ) {
			a [i][j] = matrix [i-1][j-1];
		}
	}
}
void MSPurityCorrection::correction ( DoubleVector& dv ) const
{
	loadA ();			// a must be reloaded every time as it is overwritten by gaussj
	for ( int i = 0 ; i < numQuanPeaks ; i++ ) {
		b [i+1][1] = dv [i];
	}
	gaussj ( a, numQuanPeaks, b, 1 );
	for ( int j = 0 ; j < numQuanPeaks ; j++ ) {
		dv [j] = b [j+1][1];
	}
}
QuanModInfoMap QuantitationRatio::quanMods;
StringVector QuantitationRatio::mod;
bool QuantitationRatio::o18Flag = false;
bool QuantitationRatio::n15Flag = false;
int QuantitationRatio::numQuanStates = 0;
int QuantitationRatio::numRatios = 0;

QuantitationRatio::QuantitationRatio ( const XYData& xydata, double mOZ, int charge, const QuanPeptide& qPeptide, const ElementalFormula* ef, bool efLink, double resolution, int numPeaks ) :
	QuantitationData ( getNumPeaks ( numPeaks ) ),
	charge ( charge ),
	numPeaks ( getNumPeaks ( numPeaks ) )
{
	init ( xydata, ef, efLink, resolution, mOZ, qPeptide );
}
void QuantitationRatio::init ( const XYData& xyData, const ElementalFormula* ef, bool efLink, double resolution, double mOZ, const QuanPeptide& qPeptide )
{
	if ( ef == 0 ) {
		ErrorHandler::genError ()->error ( "Quantitation not possible as the elemental formula of the peptide can't be calculated.\n" );
	}
	vector <ElementalFormula> efDelta (numQuanStates);
	if ( n15Flag )	quanPep = getQuanMassesN15 ( qPeptide, ef, efDelta );
	else			quanPep = getQuanMasses ( qPeptide, efDelta );
	if ( !quanPep )	for ( int i = 0 ; i < numQuanStates ; i++ ) efDelta [i] = "";
	ColouredGraphData* graphData = 0;
	bool graphFlag = PeakFit::getGraphs ();
	if ( graphFlag ) graphData = new ColouredGraphData ( xyData );
	PeakFit::getNoise ( xyData, noiseMean, noiseStDev );
	DoubleVectorVectorVector coeff;
	double mPlusH = mOverZToMPlusH ( mOZ, charge, true );
	for ( int i = 0 ; i < numQuanStates ; i++ ) {
		int np = ( i == 0 && o18Flag ) ? numPeaks + 2 : numPeaks;
		double monoMass = mPlusHToMOverZ ( mPlusH + formula_to_monoisotopic_mass ( efDelta [i] ), charge, true );
		coeff.push_back ( PeakFit::getCoefficients ( xyData, graphData, monoMass, np, charge, resolution, noiseMean, noiseStDev, snrThreshold ) );
	}
	if ( graphFlag ) {
		PeakFit::drawGraph ( *graphData, false );
		delete graphData;
	}
	ElementalFormulaVector formulae;
	for ( int j = 0 ; j < numQuanStates ; j++ ) {
		formulae.push_back ( *ef );
		formulae.back () += efDelta [j];
		if ( ( reportFormulaString && efLink ) || ( reportPeakCSIntensity || reportPeakCSArea ) ) {
			formulaStrings.push_back ( formulae.back ().getFormula () );
		}
	}
	init2 ( coeff, formulae );
}
void QuantitationRatio::init2 ( const DoubleVectorVectorVector& coeff, ElementalFormulaVector& formulae )
{
	ok.resize ( numQuanStates ); 
	mOverZ.resize ( numQuanStates ); 
	lightHeavyIntRatioType.resize ( numRatios );
	lightHeavyAreaRatioType.resize ( numRatios );
	lightHeavyIntRatio.resize ( numRatios );
	lightHeavyAreaRatio.resize ( numRatios );
	DoubleVectorVector width (numQuanStates);
	DoubleVectorVector intensity1 (numQuanStates);
	DoubleVectorVector snr1 (numQuanStates);
	DoubleVectorVector fwhm1 (numQuanStates);
	DoubleVectorVector resolution1 (numQuanStates);
	DoubleVectorVector area1 (numQuanStates);
	for ( int i = 0 ; i < numQuanStates ; i++ ) {			
		const DoubleVectorVector& coeff1 = coeff [i];
		for ( int j = 0 ; j < numPeaks ; j++ ) {
			bool okFlag = ( o18Flag && i == 0 ) ? ( j + 2 <  coeff1.size () ) : ( j < coeff1.size () );
			ok [i].push_back ( ( okFlag ? 1 : 0 ) );
			intensity1 [i].push_back	( ( okFlag ? coeff1 [j][0] : 0.0 ) );
			mOverZ [i].push_back	( ( okFlag ? coeff1 [j][1] : 0.0 ) );
			width [i].push_back		( ( okFlag ? coeff1 [j][2] : 0.0 ) );
			snr1 [i].push_back		( ( okFlag ? (intensity1 [i][j] / noiseStDev) : 0.0 ) );
			fwhm1 [i].push_back		( ( okFlag ? (2.35 * width [i][j] / sqrtTwo) : 0.0 ) );
			resolution1 [i].push_back( ( okFlag ? (mOverZ [i][j] / fwhm1 [i][j]) : 0.0 ) );
			area1 [i].push_back		( ( okFlag ? (intensity1 [i][j] * width [i][j] * sqrtPi) : 0.0 ) );
		}
	}
	for ( int ii = 1 ; ii < numQuanStates ; ii++ ) {
		lightHeavyIntRatioType [ii-1].resize ( numPeaks, RATIO_NOT_CALCULATED );
		lightHeavyAreaRatioType [ii-1].resize ( numPeaks, RATIO_NOT_CALCULATED );
		lightHeavyIntRatio [ii-1].resize ( numPeaks, 0.0 );
		lightHeavyAreaRatio [ii-1].resize ( numPeaks, 0.0 );
		for ( int jj = 0 ; jj < numPeaks ; jj++ ) {
			if ( ok [0][jj] && ok [ii][jj] ) {
				double snrLow = snr1 [0][jj];
				double snrHigh = snr1 [ii][jj];
				if ( snrLow < snrThreshold && snrHigh < snrThreshold ) {
					lightHeavyIntRatioType [ii-1][jj] = RATIO_NOT_CALCULATED;
					lightHeavyAreaRatioType [ii-1][jj] = RATIO_NOT_CALCULATED;
				}
				else {
					if ( snrLow < snrThreshold ) {
						lightHeavyIntRatioType [ii-1][jj] = RATIO_LESS_THAN;
						lightHeavyAreaRatioType [ii-1][jj] = RATIO_LESS_THAN;
						intensity1 [0][jj] = noiseStDev * snrThreshold;
						area1 [0][jj] = intensity1 [0][jj] * width [ii][jj] * sqrtPi;	// Use width for other peak
					}
					else if ( snrHigh < snrThreshold ) {
						lightHeavyIntRatioType [ii-1][jj] = RATIO_GREATER_THAN;
						lightHeavyAreaRatioType [ii-1][jj] = RATIO_GREATER_THAN;
						intensity1 [ii][jj] = noiseStDev * snrThreshold;
						area1 [ii][jj] = intensity1 [ii][jj] * width [0][jj] * sqrtPi;	// Use width for other peak
					}
					else {
						lightHeavyIntRatioType [ii-1][jj] = RATIO_OK;
						lightHeavyAreaRatioType [ii-1][jj] = RATIO_OK;
					}
				}
			}
		}
	}
	DoubleVectorVector dvvIntensity;
	DoubleVectorVector dvvArea;
	vector <IsotopePeakStats*> vips;
	if ( o18Flag ) {
		vips.push_back ( new IsotopePeakStats ( formulae [0], charge, 0.0001 ) );	// A lower probability limit has been used for increased
																					// speed. The only effect is likely to be with very high ratios.
	}
	else {		// Correct intensities for purity
		dvvIntensity = intensity1;
		dvvArea = area1;
		for ( int i = 0 ; i < numQuanStates ; i++ ) {
			vips.push_back ( new IsotopePeakStats ( formulae [i], charge, 0.0001 ) );
		}
		for ( int j = 0 ; j < numPeaks ; j++ ) {
			DoubleVector dvIntensity (numQuanStates);
			DoubleVector dvArea (numQuanStates);
			for ( int k = 0 ; k < numQuanStates ; k++ ) {
				dvIntensity [k] = dvvIntensity [k][j];
				dvArea [k] = dvvArea [k][j];
			}
			MSPurityCorrection pc ( vips, j );
			bool singMatrix = false;
			try {
				pc.correction ( dvIntensity );
				pc.correction ( dvArea );
			}
			catch ( lNrecGaussjSingularMatrix1 ) {
				singMatrix = true;
			}
			catch ( lNrecGaussjSingularMatrix2 ) {
				singMatrix = true;
			}
			if ( !singMatrix ) {
				for ( int m = 0 ; m < numQuanStates ; m++ ) {
					dvvIntensity [m][j] = dvIntensity [m];
					dvvArea [m][j] = dvArea [m];
				}
			}
		}
	}
	for ( int k = 1 ; k < numQuanStates ; k++ ) {
		for ( int m = 0 ; m < numPeaks ; m++ ) {
			if ( ok [0][m] && ok [k][m] ) {
				double measuredLightHeavyRatioInt;
				double measuredLightHeavyRatioArea;
				if ( o18Flag ) {
					IsotopePeakStats* ipsLight = vips [0];
					double m4overm0 = ipsLight->getProbability ( m + 4 ) / ipsLight->getProbability ( m );	// Eqn. 1 from Zang et al 2004, J. Proteome Research Vol. 3, No. 3, pp 604-612
					double m2overm0 = ipsLight->getProbability ( m + 2 ) / ipsLight->getProbability ( m );
					double oneMinusm2overm0 = 1.0 - m2overm0;

					double intensityPlus2 = coeff [0][m+2][0];
					double bottomIntensity = intensity1 [k][0];
					bottomIntensity -= m4overm0 * intensity1 [0][m];
					bottomIntensity += oneMinusm2overm0 * intensityPlus2;
					bottomIntensity -= oneMinusm2overm0 * m2overm0 * intensity1 [0][m];
					measuredLightHeavyRatioInt = intensity1 [0][m] / bottomIntensity;

					double areaPlus2 = intensityPlus2 * coeff [0][m+2][2] * sqrtPi;
					double bottomArea = area1 [k][0];
					bottomArea -= m4overm0 * area1 [0][m];
					bottomArea += oneMinusm2overm0 * areaPlus2;
					bottomArea -= oneMinusm2overm0 * m2overm0 * area1 [0][m];
					measuredLightHeavyRatioArea = area1 [0][m] / bottomArea;
				}
				else {
					measuredLightHeavyRatioInt = dvvIntensity [0][m] / dvvIntensity [k][m];
					measuredLightHeavyRatioArea = dvvArea [0][m] / dvvArea [k][m];
				}
				if ( measuredLightHeavyRatioInt < 0.0 ) lightHeavyIntRatioType [k-1][m] = RATIO_HIGH;
				lightHeavyIntRatio [k-1][m] = measuredLightHeavyRatioInt;
				if ( measuredLightHeavyRatioArea < 0.0 ) lightHeavyAreaRatioType [k-1][m] = RATIO_HIGH;
				lightHeavyAreaRatio [k-1][m] = measuredLightHeavyRatioArea;
			}
		}
	}
	for ( int x = 0 ; x < vips.size () ; x++ ) {
		delete vips [x];
	}
	if ( reportPeakIntensity )	intensity = intensity1;
	if ( reportPeakSNR )		snr = snr1;
	if ( reportPeakFWHM )		fwhm = fwhm1;
	if ( reportPeakResolution )	resolution = resolution1;
	if ( reportPeakArea )		area = area1;
}
void QuantitationRatio::printHTML ( ostream& os ) const
{
	string styleID1 = "sc_stripe_1";
	tableStart ( os, true );
	for ( int i = 0 ; i < numPeaks ; i++ ) {
		if ( i == 0 ) {
			tableRowStart ( os );
			printHTMLHeader1 ( os, "m/z", styleID1 );
			printHTMLNoiseHeader ( os, styleID1 );
			tableRowEnd ( os );
		}
		tableRowStart ( os );
		printHTMLMOverZ ( os, i, styleID1 );
		printHTMLNoiseLine ( os, styleID1 );
		tableRowEnd ( os );
	}
	tableEnd ( os );

	tableStart ( os, true );
	for ( int j = 0 ; j < numPeaks ; j++ ) {
		if ( j == 0 ) {
			tableRowStart ( os );
			printHTMLHeader ( os, 0, styleID1 );
			tableRowEnd ( os );
		}
		tableRowStart ( os );
		printHTMLLine ( os, j, styleID1 );
		tableRowEnd ( os );
	}
	tableEnd ( os );
	if ( reportFormulaString && !formulaStrings.empty () ) {
		os << "<p>" << endl;
		MSIsotopeLink isotopeLink;
		startJavascript ( os );
		isotopeLink.printHTML ( os );
		endJavascript ( os );
		for ( StringVectorSizeType k = 0 ; k < formulaStrings.size () ; k++ ) {
			if ( !formulaStrings [k].empty () ) {
				os << "<b>Elemental Composition Pk " << k+1 << ": </b>";
				ElementalFormula efTemp ( formulaStrings [k] ); 
				isotopeLink.write ( os, efTemp, charge );
				os << "<br />" << endl;
				CosSimilarity cs ( formulaStrings [k], charge, intensity [k], area [k] );
				cs.printHTML ( os, k+1 );
			}
		}
		os << "<br />" << endl;
	}
}
int QuantitationRatio::getColspan ( DoubleVectorVectorSizeType peakNumber )
{
	int colspan = 0;
	if ( reportPeakIntensity )				colspan += numQuanStates;
	if ( reportPeakSNR )					colspan += numQuanStates;
	if ( reportPeakResolution )				colspan += numQuanStates;
	if ( reportPeakCSIntensity && peakNumber == 0 )			colspan += numQuanStates;
	if ( reportPeakFWHM )					colspan += numQuanStates;
	if ( reportActualLightHeavyIntensityRatio )		colspan += numRatios;
	if ( reportPeakArea )							colspan += numQuanStates;
	if ( reportPeakCSArea && peakNumber == 0 )					colspan += numQuanStates;
	if ( reportActualLightHeavyAreaRatio )			colspan += numRatios;
	return colspan;
}
void QuantitationRatio::printHTMLHeader ( ostream& os, DoubleVectorVectorSizeType peakNumber, const string& styleID )
{
	if ( reportPeakIntensity )	printHTMLHeader1 ( os, "Int", styleID );
	if ( reportPeakSNR )		printHTMLHeader1 ( os, "SNR", styleID );
	if ( reportPeakResolution )	printHTMLHeader1 ( os, "Resolution", styleID );
	if ( reportPeakCSIntensity && peakNumber == 0 ) printHTMLHeader1 ( os, "CS", styleID );
	if ( reportPeakFWHM )		printHTMLHeader1 ( os, "FWHM", styleID );

	if ( reportActualLightHeavyIntensityRatio )		printHTMLHeader2 ( os, "Intensity", styleID );

	if ( reportPeakArea )		printHTMLHeader1 ( os, "Area", styleID );
	if ( reportPeakCSArea && peakNumber == 0 )		printHTMLHeader1 ( os, "CS", styleID );

	if ( reportActualLightHeavyAreaRatio )			printHTMLHeader2 ( os, "Area", styleID );
}
void QuantitationRatio::printHTMLHeader1 ( ostream& os, const string& str, const string& styleID )
{
	for ( int i = 1 ; i <= numQuanStates ; i++ ) {
		tableHeader ( os, str + " Pk " + gen_itoa ( i ), styleID );
	}
}
void QuantitationRatio::printHTMLHeader2 ( ostream& os, const string& str, const string& styleID )
{
	for ( int i = 0 ; i < numRatios ; i++ ) {
		tableHeader ( os, getRatioString ( i ) + " " + str, styleID );
	}
}
void QuantitationRatio::printHTMLBlankLine ( ostream& os, DoubleVectorVectorSizeType peakNumber, const string& styleID )
{
	if ( reportPeakIntensity )						tableEmptyNCells ( os, numQuanStates, styleID );
	if ( reportPeakSNR )							tableEmptyNCells ( os, numQuanStates, styleID );
	if ( reportPeakResolution )						tableEmptyNCells ( os, numQuanStates, styleID );
	if ( reportPeakCSIntensity && peakNumber == 0 )	tableEmptyNCells ( os, numQuanStates, styleID );
	if ( reportPeakFWHM )							tableEmptyNCells ( os, numQuanStates, styleID );
	if ( reportActualLightHeavyIntensityRatio )		tableEmptyNCells ( os, numRatios, styleID );
	if ( reportPeakArea )							tableEmptyNCells ( os, numQuanStates, styleID );
	if ( reportPeakCSArea && peakNumber == 0 )		tableEmptyNCells ( os, numQuanStates, styleID );
	if ( reportActualLightHeavyAreaRatio )			tableEmptyNCells ( os, numRatios, styleID );
}
void QuantitationRatio::printHTMLMOverZ ( ostream& os, DoubleVectorVectorSizeType peakNumber, const string& styleID ) const
{
	for ( int i = 0 ; i < mOverZ.size () ; i++ ) {
		if ( ok [i][peakNumber] ) tableCell ( os, mOverZ [i][peakNumber], 4, false, styleID );
		else tableEmptyCell ( os, styleID );
	}
}
void QuantitationRatio::printHTMLLine ( ostream& os, DoubleVectorVectorSizeType peakNumber, const string& styleID ) const
{
	CosSimilarityList* csl = 0;
	if ( ( reportPeakCSIntensity || reportPeakCSArea ) && !formulaStrings.empty () && peakNumber == 0 ) {
		csl = new CosSimilarityList ( formulaStrings, charge, intensity, area, numQuanStates );
	}
	if ( reportPeakIntensity )						printHTMLLine1 ( os, intensity, peakNumber, styleID );
	if ( reportPeakSNR )							printHTMLLine1 ( os, snr, peakNumber, styleID );
	if ( reportPeakResolution )						printHTMLLine2 ( os, resolution, peakNumber, styleID );
	if ( reportPeakCSIntensity && peakNumber == 0 ) csl->printHTMLIntensity ( os, styleID );
	if ( reportPeakFWHM )							printHTMLLine1 ( os, snr, peakNumber, styleID );
	if ( reportActualLightHeavyIntensityRatio )		printHTMLLine4 ( os, lightHeavyIntRatio, lightHeavyIntRatioType, peakNumber, styleID );
	if ( reportPeakArea )							printHTMLLine1 ( os, area, peakNumber, styleID );
	if ( reportPeakCSArea && peakNumber == 0 )		csl->printHTMLArea ( os, styleID );
	if ( reportActualLightHeavyAreaRatio )			printHTMLLine4 ( os, lightHeavyAreaRatio, lightHeavyAreaRatioType, peakNumber, styleID );
	if ( csl != 0 ) delete csl;
}
void QuantitationRatio::printHTMLLine1 ( ostream& os, const DoubleVectorVector& dvv, DoubleVectorVectorSizeType peakNumber, const string& styleID ) const
{
	for ( int i = 0 ; i < dvv.size () ; i++ ) {
		if ( ok [i][peakNumber] ) tableCellSigFig ( os, dvv [i][peakNumber], 4, false, styleID );
		else tableEmptyCell ( os, styleID );
	}
}
void QuantitationRatio::printHTMLLine2 ( ostream& os, const DoubleVectorVector& dvv, DoubleVectorVectorSizeType peakNumber, const string& styleID ) const
{
	for ( int i = 0 ; i < dvv.size () ; i++ ) {
		if ( ok [i][peakNumber] ) tableCell ( os, dvv [i][peakNumber], 0, false, styleID );
		else tableEmptyCell ( os, styleID );
	}
}
void QuantitationRatio::printHTMLLine3 ( ostream& os, const DoubleVectorVector& dvv, DoubleVectorVectorSizeType peakNumber, const string& styleID ) const
{
	for ( int i = 0 ; i < dvv.size () ; i++ ) {
		if ( ok [0][peakNumber] && ok [i+1][peakNumber] && quanPep )
			tableCellSigFig ( os, dvv [i][peakNumber], 4, false, styleID );
		else
			tableEmptyCell ( os, styleID );
	}
}
void QuantitationRatio::printHTMLLine4 ( ostream& os, const DoubleVectorVector& ratio, const CharVectorVector& ratioType, DoubleVectorVectorSizeType peakNumber, const string& styleID ) const
{
	for ( int i = 0 ; i < ratio.size () ; i++ ) {
		if ( ok [0][peakNumber] && ok [i+1][peakNumber] && quanPep )
			printHTMLRatio ( os, ratio [i][peakNumber], ratioType [i][peakNumber], styleID );
		else
			tableEmptyCell ( os, styleID );
	}
}
void QuantitationRatio::printDelimitedHeader ( ostream& os, DoubleVectorVectorSizeType peakNumber )
{
	if ( reportPeakIntensity )	printDelimitedHeader1 ( os, "Int" );
	if ( reportPeakSNR )		printDelimitedHeader1 ( os, "SNR" );
	if ( reportPeakResolution )	printDelimitedHeader1 ( os, "Resolution" );
	if ( reportPeakCSIntensity && peakNumber == 0 )printDelimitedHeader1 ( os, "CS" );
	if ( reportPeakFWHM )		printDelimitedHeader1 ( os, "FWHM" );

	if ( reportActualLightHeavyIntensityRatio )		printDelimitedHeader2 ( os, "Intensity" );

	if ( reportPeakArea )		printDelimitedHeader1 ( os, "Area" );
	if ( reportPeakCSArea && peakNumber == 0 )		printDelimitedHeader1 ( os, "CS" );

	if ( reportActualLightHeavyAreaRatio )			printDelimitedHeader2 ( os, "Area" );
}
void QuantitationRatio::printDelimitedHeader1 ( ostream& os, const string& str )
{
	for ( int i = 1 ; i <= numQuanStates ; i++ ) delimitedHeader ( os, str + " Pk " + gen_itoa ( i ) );
}
void QuantitationRatio::printDelimitedHeader2 ( ostream& os, const string& str )
{
	for ( int i = 0 ; i < numRatios ; i++ ) delimitedHeader ( os, getRatioString ( i ) + " " + str );
}
void QuantitationRatio::printDelimitedBlankLine ( ostream& os, DoubleVectorVectorSizeType peakNumber )
{
	if ( reportPeakIntensity )						delimitedEmptyNCells ( os, numQuanStates );
	if ( reportPeakSNR )							delimitedEmptyNCells ( os, numQuanStates );
	if ( reportPeakResolution )						delimitedEmptyNCells ( os, numQuanStates );
	if ( reportPeakCSIntensity && peakNumber == 0 )					delimitedEmptyNCells ( os, numQuanStates );
	if ( reportPeakFWHM )							delimitedEmptyNCells ( os, numQuanStates );
	if ( reportActualLightHeavyIntensityRatio )		delimitedEmptyNCells ( os, numRatios );
	if ( reportPeakArea )							delimitedEmptyNCells ( os, numQuanStates );
	if ( reportPeakCSArea && peakNumber == 0 )							delimitedEmptyNCells ( os, numQuanStates );
	if ( reportActualLightHeavyAreaRatio )			delimitedEmptyNCells ( os, numRatios );
}
void QuantitationRatio::printDelimitedLine ( ostream& os, DoubleVectorVectorSizeType peakNumber ) const
{
	CosSimilarityList* csl = 0;
	if ( ( reportPeakCSIntensity || reportPeakCSArea ) && !formulaStrings.empty () && peakNumber == 0 ) {
		csl = new CosSimilarityList ( formulaStrings, charge, intensity, area, numQuanStates );
	}
	if ( reportPeakIntensity )						printDelimitedLine1 ( os, intensity, peakNumber );
	if ( reportPeakSNR )							printDelimitedLine1 ( os, snr, peakNumber );
	if ( reportPeakResolution )						printDelimitedLine2 ( os, resolution, peakNumber );
	if ( reportPeakCSIntensity && peakNumber == 0 )	csl->printDelimitedIntensity ( os );
	if ( reportPeakFWHM )							printDelimitedLine1 ( os, fwhm, peakNumber );
	if ( reportActualLightHeavyIntensityRatio )		printDelimitedLine4 ( os, lightHeavyIntRatio, lightHeavyIntRatioType, peakNumber );
	if ( reportPeakArea )							printDelimitedLine1 ( os, area, peakNumber );
	if ( reportPeakCSArea && peakNumber == 0 )		csl->printDelimitedArea ( os );
	if ( reportActualLightHeavyAreaRatio )			printDelimitedLine4 ( os, lightHeavyAreaRatio, lightHeavyAreaRatioType, peakNumber );
}
void QuantitationRatio::printDelimitedLine1 ( ostream& os, const DoubleVectorVector& dvv, DoubleVectorVectorSizeType peakNumber ) const
{
	for ( int i = 0 ; i < dvv.size () ; i++ ) {
		if ( ok [i][peakNumber] ) delimitedCellSigFig ( os, dvv [i][peakNumber], 4 );
		else delimitedEmptyCell ( os );
	}
}
void QuantitationRatio::printDelimitedLine2 ( ostream& os, const DoubleVectorVector& dvv, DoubleVectorVectorSizeType peakNumber ) const
{
	for ( int i = 0 ; i < dvv.size () ; i++ ) {
		if ( ok [i][peakNumber] ) delimitedCell ( os, dvv [i][peakNumber], 0 );
		else delimitedEmptyCell ( os );
	}
}
void QuantitationRatio::printDelimitedLine3 ( ostream& os, const DoubleVectorVector& dvv, DoubleVectorVectorSizeType peakNumber ) const
{
	for ( int i = 0 ; i < dvv.size () ; i++ ) {
		if ( ok [0][peakNumber] && ok [i+1][peakNumber] && quanPep )
			delimitedCellSigFig ( os, dvv [i][peakNumber], 4 );
		else
			delimitedEmptyCell ( os );
	}
}
void QuantitationRatio::printDelimitedLine4 ( ostream& os, const DoubleVectorVector& ratio, const CharVectorVector& ratioType, DoubleVectorVectorSizeType peakNumber ) const
{
	for ( int i = 0 ; i < ratio.size () ; i++ ) {
		if ( ok [0][peakNumber] && ok [i+1][peakNumber] && quanPep )
			printDelimitedRatio ( os, ratio [i][peakNumber], ratioType [i][peakNumber] );
		else
			delimitedEmptyCell ( os );
	}
}
void QuantitationRatio::getModInfo ( const string& modLine, string& mName, StringVector& aas )
{
	int end = modLine.find ( " " );
	if ( end != string::npos ) {
		string lab;
		if ( !isPrefix ( modLine, "Label" ) ) {			// The light form is a modified amino acid
			StringSizeType idx = modLine.find_last_of ( ":" );
			if ( idx != string::npos ) {
				lab = modLine.substr ( 0, idx );
			}
		}
		int end2 = modLine.find ( "(", end );
		int end3 = modLine.find ( ")", end2 );
		string str = modLine.substr ( end2+1, end3-end2-1 );
		for ( StringVectorSizeType i = 0 ; i < str.length () ; i++ ) {
			string s = string ( 1, str [i] );
			if ( !lab.empty () ) s += "(" + lab;
			aas.push_back ( s );
		}
		mName = modLine.substr ( 0, end );
	}
	else {
		mName = modLine;
		aas.push_back ( "." );
	}
}
class SortQuanModsByDiffMass {
	public:
		int operator () ( const QuanModInfo& a, const QuanModInfo& b ) const
		{
			return a.getMonoMass () < b.getMonoMass ();
		}
};
void QuantitationRatio::setQuanResidue ( const string& quanType )
{
	static bool init = false;
	if ( !init ) {
		if ( quanType == "Label:15N" ) {
			n15Flag = ( quanType == "Label:15N" );
			numRatios = 1;
			numQuanStates = 2;
			init = true;
		}
		else {
			try {
				StringVector sv = QuantitationInfo::instance ().getQuanInfo ( quanType );
				mod.clear ();
				for ( StringVectorSizeType i = 0 ; i < sv.size () ; i++ ) {
					string modName;
					StringVector aas;
					getModInfo ( sv [i], modName, aas );
					mod.push_back ( modName );
					for ( StringVectorSizeType j = 0 ; j < aas.size () ; j++ ) {
						string lab = aas [j];
						char aa = lab[0];
						quanMods [lab].push_back ( QuanModInfo ( string ( 1, aa ), modName ) );
					}
				}
				numRatios = quanMods.begin ()->second.size ();
				numQuanStates = numRatios + 1;
				for ( QuanModInfoMapIterator k = quanMods.begin () ; k != quanMods.end () ; k++ ) {
					if ( (*k).second.size () + 1 != numQuanStates ) {
						throw runtime_error ( "Invalid entry in the file quan.txt." );
					}
					if ( numQuanStates > 2 ) stable_sort ( (*k).second.begin (), (*k).second.end (), SortQuanModsByDiffMass () );
				}
				o18Flag = ( quanType == "Label:18O(2)" );
				init = true;
			}
			catch ( runtime_error e ) {
				ErrorHandler::genError ()->error ( e );
			}
		}
	}
}
string QuantitationRatio::getUnmodifiedPeptide ( const string& peptide )
{
	string unmodPep = peptide;
	for ( StringVectorSizeType i = 0 ; i < mod.size () ; i++ ) {
		string m = mod [i];
		if ( isPrefix ( m, "Label:" ) ) {
			unmodPep = genReplaceSubstrings ( unmodPep, "(" + m + ")", "" );
			unmodPep = genReplaceSubstrings ( unmodPep, "(" + m + "+", "(" );
			unmodPep = genReplaceSubstrings ( unmodPep, m + "-", "-" );				// N-terminus
			unmodPep = genReplaceSubstrings ( unmodPep, "-" + m, "-" );				// C-terminus
			unmodPep = genReplaceSubstrings ( unmodPep, m + "+", "" );				// Should work for N and C terminus
		}
		else {
			StringSizeType end = m.find_last_of ( ":" );
			if ( end != string::npos ) {
				string u = m.substr ( 0, end );
				unmodPep = genReplaceSubstrings ( unmodPep, "(" + m + ")", "(" + u + ")" );
				unmodPep = genReplaceSubstrings ( unmodPep, "(" + m + "+", "(" + u + "+" );
				unmodPep = genReplaceSubstrings ( unmodPep, m + "-", u + "-" );		// N-terminus
				unmodPep = genReplaceSubstrings ( unmodPep, "-" + m, "-" + u );		// C-terminus
				unmodPep = genReplaceSubstrings ( unmodPep, m + "+", u + u );		// Should work for N and C terminus
			}
		}

	}
	return unmodPep;
}
bool QuantitationRatio::getDataRange ( const QuanPeptide& qPeptide, double mOverZ, int charge, double& startMass, double& endMass )
{
	IntVector totalNumMods ( numQuanStates, 0 );
	DoubleVector massDelta ( numQuanStates, 0.0 );
	massDelta [0] = 0.0;
	for ( QuanModInfoMapConstIterator i = quanMods.begin () ; i != quanMods.end () ; i++ ) {	// Check if the peptide is a heavy one
		string res = (*i).first;
		int nm = ( res == "." ) ? 1 : getNumPotentialMods ( qPeptide, res );	// Number of times the residue appears
		totalNumMods [0] += nm;
		vector <QuanModInfo> vqmi = (*i).second;
		for ( VectorQuanModInfoSizeType j = 0 ; j < vqmi.size () ; j++ ) {	// Different mods on the residue
			int numMods = vqmi [j].getNumMods ( qPeptide );	// Times the modified residue appears
			totalNumMods [j+1] += numMods;
			if ( nm ) massDelta [j+1] += vqmi [j].getMonoMass ( nm );
		}
	}
	if ( totalNumMods [0] == 0 ) return false;	// No modifiable residues
	int idx = 0;
	for ( IntVectorSizeType k = 1 ; k < numQuanStates ; k++ ) {		// Check for heavy labels
		if ( totalNumMods [k] ) {
			if ( totalNumMods [k] < totalNumMods [0] ) return false;		// Not fully labelled
			else
				idx = k;
		}
	}
	if ( idx != 0 ) {	// deltas initially set for light peptide
		massDelta [0] = -massDelta [idx];
		for ( int i = 1 ; i < numQuanStates ; i++ ) {
			if ( i == idx ) massDelta [i] = 0.0;
			else massDelta [i] += massDelta [0];
		}
	}
	double mPlusH = mOverZToMPlusH ( mOverZ, charge, true );
	startMass = mPlusHToMOverZ ( mPlusH + massDelta [0], charge, true ) - 6.0;
	endMass = mPlusHToMOverZ ( mPlusH + massDelta [numQuanStates-1], charge, true ) + 10.0;
	return true;
}
bool QuantitationRatio::getQuanMasses ( const ElementalFormula* ef, vector <ElementalFormula>& efDelta )
{
	ElementalFormula f = *ef;
	int changed = 0;
	for ( f.first () ; f.isDone () ; f.next () ) {
		if ( f.element () == "N" ) {
			changed++;
			if ( changed == 2 ) return false;			// Incomplete labelling - peptide contains N and 15N
			string multi = gen_itoa ( f.multiplier () );
			efDelta [0] = "";
			efDelta [1] = "N-" + multi + " 15N" + multi;
		}
		if ( f.element () == "15N" ) {
			changed++;
			if ( changed == 2 ) return false;			// Incomplete labelling - peptide contains N and 15N
			string multi = gen_itoa ( f.multiplier () );
			efDelta [0] = "15N-" + multi + " N" + multi;
			efDelta [1] = "";
		}
	}
	if ( changed == 0 ) return false;	// No nitrogen content
	else				return true;
}
bool QuantitationRatio::getQuanMasses ( const QuanPeptide& qPeptide, vector <ElementalFormula>& efDelta )
{
	IntVector totalNumMods ( numQuanStates, 0 );

	efDelta [0] = "";
	for ( QuanModInfoMapConstIterator i = quanMods.begin () ; i != quanMods.end () ; i++ ) {	// Check if the peptide is a heavy one
		string res = (*i).first;
		int nm = ( res == "." ) ? 1 : getNumPotentialMods ( qPeptide, res );	// Number of times the residue appears
		totalNumMods [0] += nm;
		vector <QuanModInfo> vqmi = (*i).second;
		for ( VectorQuanModInfoSizeType j = 0 ; j < vqmi.size () ; j++ ) {	// Different mods on the residue
			int numMods = vqmi [j].getNumMods ( qPeptide );					// Times the modified residue appears
			totalNumMods [j+1] += numMods;
			if ( nm ) efDelta [j+1] += vqmi [j].getDiffFormula ( nm );
		}
	}
	if ( totalNumMods [0] == 0 ) return false;	// No modifiable residues
	int idx = 0;
	for ( IntVectorSizeType k = 1 ; k < numQuanStates ; k++ ) {		// Check for heavy labels
		if ( totalNumMods [k] ) {
			if ( totalNumMods [k] < totalNumMods [0] ) return false;		// Not fully labelled
			else
				idx = k;
		}
	}
	if ( idx != 0 ) {	// efDeltas initially set for light peptide
		efDelta [0] = efDelta [idx];
		efDelta [0].multiply ( -1 );
		for ( int i = 1 ; i < numQuanStates ; i++ ) {
			if ( i == idx ) efDelta [i] = "";
			else efDelta [i] += efDelta [0];
		}
	}
	return true;
}
void getResiduesFromFormula ( const string& pep, StringVector& residues )
{
	StringSizeType len = pep.length ();
	for ( StringSizeType i = 0 ; i < len ; i++ ) {
		if ( pep [i] == '(' ) {
			int bracket = 0;
			StringSizeType start = i+1;
			for ( ; i < len ; i++ ) {
				char a = pep [i];
				if ( a == '(' ) bracket++;
				if ( a == ')' ) bracket--;
				if ( bracket == 0 ) break;
			}
			residues.push_back ( pep.substr ( start, i-start ) );
		}
	}
}
bool QuantitationRatio::getQuanMassesN15 ( const QuanPeptide& qPeptide, const ElementalFormula* ef, vector <ElementalFormula>& efDelta )
{
	static MapStringToInt balance = Usermod::getN15Balance ();
	bool heavy = false;
	string p1 = qPeptide.getPeptide1 ();
	string p2 = qPeptide.getPeptide2 ();
	if ( p1.find ( "Label:15N" ) != string::npos ) heavy = true;
	if ( heavy ) {
		ElementalFormula f = *ef;
		for ( f.first () ; f.isDone () ; f.next () ) {
			if ( f.element () == "15N" ) {
				string multi = gen_itoa ( f.multiplier () );
				efDelta [0] = "15N-" + multi + " N" + multi;
				efDelta [1] = "";
			}
		}
	}
	else {
		ElementalFormula f = *ef;
		int multi = 0;
		for ( f.first () ; f.isDone () ; f.next () ) {
			if ( f.element () == "N" ) {
				multi = f.multiplier ();
			}
		}
		StringVector residues;
		getResiduesFromFormula ( p1, residues );
		string nTerm1 = qPeptide.getNTerm1 ();
		if ( !nTerm1.empty () ) residues.push_back ( nTerm1 );
		string cTerm1 = qPeptide.getCTerm1 ();
		if ( !cTerm1.empty () ) residues.push_back ( cTerm1 );
		if ( !p2.empty () ) {
			string nTerm2 = qPeptide.getNTerm2 ();
			if ( !nTerm2.empty () ) residues.push_back ( nTerm2 );
			string cTerm2 = qPeptide.getCTerm2 ();
			if ( !cTerm2.empty () ) residues.push_back ( cTerm2 );
			getResiduesFromFormula ( p2, residues );
		}
		int b = 0;
		for ( StringVectorSizeType i = 0 ; i < residues.size () ; i++ ) {
			MapStringToIntConstIterator cur = balance.find ( residues [i] );
			if ( cur != balance.end () ) {
				b += (*cur).second;
			}
		}
		multi -= b;
		string sMulti = gen_itoa ( multi );
		efDelta [0] = "";
		efDelta [1] = "N-" + sMulti + " 15N" + sMulti;
	}
	return true;
}
int QuantitationRatio::getNumPotentialMods ( const QuanPeptide& qPeptide, const string& aas )
{
	int num = 0;
	if ( aas.empty () ) return 0;									// There aren't any aa's to look for
	string pep;
	if ( aas.length () == 1 )
		pep = gen_strstriptags2 ( qPeptide.getPeptide1 (), '(', ')' );				// Remove modifications
	else {
		if ( isPrefix ( qPeptide.getNTerm1 (), aas.substr ( 2 ) ) ) num++;
		if ( isPrefix ( qPeptide.getCTerm1 (), aas.substr ( 2 ) ) ) num++;
		pep = qPeptide.getPeptide1 ();
	}
	int startInd = 0;
	for ( ; ; ) {
		int ind = pep.find ( aas, startInd );
		if ( ind == string::npos ) break;
		else {
			num++;
			startInd = ind + 1;
		}
	}
	string pep2 = qPeptide.getPeptide2 ();
	if ( !pep2.empty () ) {
		if ( aas.length () == 1 )
			pep2 = gen_strstriptags2 ( pep2, '(', ')' );				// Remove modifications
		else {
			if ( isPrefix ( qPeptide.getNTerm2 (), aas.substr ( 2 ) ) ) num++;
			if ( isPrefix ( qPeptide.getCTerm2 (), aas.substr ( 2 ) ) ) num++;
		}
		int startInd = 0;
		for ( ; ; ) {
			int ind = pep2.find ( aas, startInd );
			if ( ind == string::npos ) break;
			else {
				num++;
				startInd = ind + 1;
			}
		}
		if ( aas [0] == 'X' ) num++;		// One for the crosslinker
	}
	return num;
}
bool QuantitationRatio::outputQuanResults ( ostream& os, const string& searchName, int numRepeats, bool area ) const
{
	if ( !quanPep ) return false;
	bool flag = false;
	for ( int i = 0 ; i < numRatios ; i++ ) {
		string s = searchName + " " + getRatioString ( i );
		if ( area ) {
			if ( outputQuanRatio ( os, s, lightHeavyAreaRatio [i][0], lightHeavyAreaRatioType [i][0], numRepeats ) ) {
				flag = true;
			}
		}
		else {
			if ( outputQuanRatio ( os, s, lightHeavyIntRatio [i][0], lightHeavyIntRatioType [i][0], numRepeats ) ) {
				flag = true;
			}
		}
	}
	return flag;
}
DoubleVector QuantitationRatio::getAreaRatios () const
{
	DoubleVector dv;
	for ( int i = 0 ; i < numRatios ; i++ ) {
		if ( quanPep )
			dv.push_back ( getQuanRatio ( lightHeavyAreaRatio [i][0], lightHeavyAreaRatioType [i][0] ) );
		else
			dv.push_back ( 0.0 );
	}
	return dv;
}
DoubleVector QuantitationRatio::getIntensityRatios () const
{
	DoubleVector dv;
	for ( int i = 0 ; i < numRatios ; i++ ) {
		if ( quanPep )
			dv.push_back ( getQuanRatio ( lightHeavyIntRatio [i][0], lightHeavyIntRatioType [i][0] ) );
		else
			dv.push_back ( 0.0 );
	}
	return dv;
}
string QuantitationRatio::getRatioString ( int i )
{
	if ( numRatios == 1 ) return "L/H";
	else if ( numRatios == 2 ) {
		if ( i == 0 ) return "L/M";
		else return "L/H";
	}
	else
		return "1/" + gen_itoa ( i + 2 );
}
int QuantitationRatio::getNumRatios ( const string& quanType )
{
	setQuanResidue ( quanType );
	return numRatios;
}
int QuantitationRatio::getNumPeaks ( int numPeaks )
{
	return o18Flag ? genMin ( numPeaks, 2 ) : ( reportPeakCSIntensity || reportPeakCSArea ) ? genMax ( numPeaks, 3 ) : numPeaks;
}
void QuantitationInfo::initialise ()
{
	GenIFStream fromFile ( MsparamsDir::instance ().getParamPath ( "quan.txt" ) );
	string line;
	while ( getline ( fromFile, line ) ) {
		if ( line.length () != 0 && line [0] != '#' ) {
			name.push_back ( line );
			string n = name.back ();
			StringVector sv;
			for ( ; ; ) {
				string modificationName;
				getline ( fromFile, modificationName );
				if ( modificationName [0] == '>' ) {
					if ( sv.empty () && n != "Label:15N" ) {
						ErrorHandler::genError ()->error ( "No quantitation modifications specified for quantitation type " + n + " in file quan.txt.\n" );
					}
					break;
				}
				else
					sv.push_back ( modificationName );
			}
			singQuanInfo [n] = sv;
		}
	}
}
QuantitationInfo::QuantitationInfo ()
{
	initialise ();
}
QuantitationInfo& QuantitationInfo::instance ()
{
	static QuantitationInfo quantitationInfo;
	return quantitationInfo;
}
QuanModInfo::QuanModInfo ( const string& aa, const string& name ) :
	name ( name ),
	name1 ( aa + "(" + name + ")" ),
	name2 ( aa + "(" + name + "+" ),
	diffFormula ( convertPSIlabelToFormula ( name ) ),
	monoMass ( formula_to_monoisotopic_mass ( diffFormula.getFormula ().c_str () ) )
{
}
int QuanModInfo::getNumMods ( const QuanPeptide& qPeptide )
{
	string pep1 = qPeptide.getPeptide1 ();
	int numMods = genNumSubstrings ( pep1, name1 ) + genNumSubstrings ( pep1, name2 );
	if ( name == qPeptide.getNTerm1 () ) numMods += 1;
	if ( name == qPeptide.getCTerm1 () ) numMods += 1;
	string pep2 = qPeptide.getPeptide2 ();
	if ( !pep2.empty () ) {
		numMods += genNumSubstrings ( pep2, name1 ) + genNumSubstrings ( pep2, name2 );
		if ( name == qPeptide.getNTerm2 () ) numMods += 1;
		if ( name == qPeptide.getCTerm2 () ) numMods += 1;
		if ( name1 [0] == 'X' && qPeptide.getLinkName () == name ) numMods += 1; // One for the crosslinker
	}
	return numMods;
}
ElementalFormula QuanModInfo::getDiffFormula ( int numMods ) const
{
	ElementalFormula ef = diffFormula;
	ef.multiply ( numMods );
	return ef;
}
string QuanModInfo::convertPSIlabelToFormula ( const string& modStr )
{
	string s = modStr.substr ( modStr.find_last_of ( ":" ) + 1 );
	string s2;
	for ( string::size_type ii = 0 ; ii < s.length () ; ) {		// Remove non-isotopes from the formula
		int isotope = isdigit ( s [ii] );
		do {
			if ( isotope ) s2 += s [ii];
		} while ( ii < s.length () && s [ii++] != ')' );
	}
	string pos = genReplaceSubstrings ( s2, "(", "" );
	pos = genReplaceSubstrings ( pos, ")", " " );
	pos = pos.substr ( 0, pos.length () - 1 );
	string neg = genReplaceSubstrings ( s2, "(", "-" );
	neg = genReplaceSubstrings ( neg, ")", " " );
	neg = neg.substr ( 0, neg.length () - 1 );
	string neg2;
	for ( string::size_type i = 0 ; i < neg.length () ; ) {
		while ( isdigit ( neg [i] ) )	i++;
		while ( isalpha ( neg [i] ) )	neg2 += neg [i++];
		while ( genIsNumberStart ( neg [i] ) ) neg2 += neg [i++];
		while ( neg [i] == ' ' )		neg2 += neg [i++];
	}
	return pos + " " + neg2;
}
#endif
