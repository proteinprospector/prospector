/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_parent.cpp                                                 *
*                                                                             *
*  Created    : December 26th 2003                                            *
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
#include <memory>
#include <iomanip>
#include <lg_string.h>
#include <lu_parent.h>
#include <lu_iso.h>
#include <lu_table.h>
#include <lu_getfil.h>
#include <lu_delim.h>
#include <lu_mass_elem.h>
#include <lu_cgi_val.h>
#include <lu_param_list.h>
#include <lu_spec_id.h>
#include <lu_app_gr.h>
#include <lu_quan_multi.h>
#include <lu_iso_par.h>
#include <lu_html.h>
using std::cout;
using std::endl;
using std::ostream;
using std::setprecision;
using std::string;
using std::ostringstream;
using std::auto_ptr;

bool PeakFit::graphs = false;
bool PeakFit::diagnostics = false;
void PeakFit::getNoise ( const XYData& xyData, double& mean, double& stddev )
{
	XYData noiseRange = xyData;
	for ( ; ; ) {
		XYData newData;
		mean = noiseRange.mean ();
		stddev = noiseRange.stddev ();
		for ( noiseRange.first () ; noiseRange.isDone () ; noiseRange.next () ) {
			if ( noiseRange.y () < mean + ( stddev * 2 ) ) {
				newData.add ( noiseRange.x (), noiseRange.y () );
			}
		}
		if ( newData.size () == noiseRange.size () ) break;
		noiseRange = newData;
	}
	if ( stddev < 0.0001 ) {
		stddev = 0.1;
	}
	if ( mean < 0.0001 ) {
		mean = 0.1;
		stddev = 0.01;
	}
}
void PeakFit::getStats ( const XYData& xyData, double& mean, double& stddev )
{
	mean = xyData.mean ();
	stddev = xyData.stddev ();
	if ( stddev < 0.0001 ) {
		stddev = 0.1;
	}
	if ( mean < 0.0001 ) {
		mean = 0.1;
		stddev = 0.01;
	}
}
DoubleVectorVector PeakFit::getCoefficients ( const XYData& xyData, ColouredGraphData* graphData, double monoMass, int numPeaks, int charge, double inputResolution, double noiseMean, double noiseStDev, double minSNR )
{
	double extent = 2000 / inputResolution;
	double graphInc = 100 / inputResolution;
	using nr::sqrtTwo;

	DoubleVectorVector ret;
	if ( !xyData.empty () ) {
		for ( int i = 0 ; i < numPeaks ; i++ ) {
			DoubleVector gues (3);
			double mass = monoMass + ( (double) i / (double) charge );
			gues [2] = ( sqrtTwo * mass ) / ( 2.35 * inputResolution );	// pg 23, Sivia, D. S. Data Analysis, A Bayesian Tutorial
			DoubleVector defaultA (3);
			defaultA [0] = 0.0;
			defaultA [1] = mass;
			defaultA [2] = gues [2];
			try {
				gues [0] = xyData.getMaxYInTolRange ( mass, 2.0 * gues [2] ) - noiseMean;
				if ( gues [0] / noiseStDev < minSNR ) {
					ret.push_back ( defaultA );
					continue;
				}
			}
			catch ( lNrecEmptyXYDataRange ) {
				ret.push_back ( defaultA );
				continue;
			}
			try {
				gues [1] = xyData.getXAtMaxYInTolRange ( mass, 2.0 * gues [2] );
			}
			catch ( lNrecEmptyXYDataRange ) {
				ret.push_back ( defaultA );
				continue;
			}
			int numCoeff = 3;
			double peakStartMass = mass - extent;
			double peakEndMass = mass + extent;
			if ( peakEndMass < xyData.maxX () ) { 
				XYData peakRange = xyData.getXRange ( peakStartMass, peakEndMass );
				try {
					DoubleVector a = minimize ( peakRange, noiseStDev, noiseMean, gues, 1 );
					a [2] = fabs ( a [2] );	// Deal with an occasional negative width.
					if ( a [1] > peakStartMass && a [1] < peakEndMass ) {	// Is mass within data range
						if ( graphData ) {
							for ( double x = peakStartMass ; x < peakEndMass ; x += graphInc ) {
								double y = 0.0;
								for ( int b = 0 ; b < numCoeff ; b += 3 ) {
									double sd = a [b+2];
									double twoSDSquared = sd * sd;
									double top = x - a [b+1];
									top *= top;
									y += a [b] * exp ( - top / twoSDSquared );
								}
								y += noiseMean;
								graphData->add ( x, y, 8 );
							}
						}
						ret.push_back ( a );
					}
					else {
						ret.push_back ( defaultA );
					}
				}
				catch ( lNrecGaussjSingularMatrix1 ) {
					break;
				}
				catch ( lNrecGaussjSingularMatrix2 ) {
					ret.push_back ( defaultA );
					continue;
				}
			}
		}
	}
	return ret;
}
DoubleVector PeakFit::minimize ( XYData& xyData, double noiseWidth, double noiseValue, const DoubleVector& gues, int numDist )
{
	int numCoeff = numDist * 3;
	DoubleVector a ( numCoeff );
	int* lista = inrvector(1,numCoeff);
	int numPts = xyData.size ();
	DoubleVector x ( numPts );
	DoubleVector y ( numPts );
	DoubleVector sig ( numPts );
	double** covar=nrmatrix(1,numCoeff,1,numCoeff);
	double** alpha=nrmatrix(1,numCoeff,1,numCoeff);
	int i = 0;
	for ( xyData.first () ; xyData.isDone () ; xyData.next () ) {
		x[i]=xyData.x ();
		y[i]=xyData.y () - noiseValue;
		sig[i]=noiseWidth;
		i++;
	}
	int mfit = numCoeff;
	for ( i = 1 ; i <= mfit ; i++ ) {
		lista [i] = i;
	}
	double alamda = -1;
	for ( i = 0 ; i < numCoeff ; i++ ) a[i] = gues[i];
	double chisq;
	double ochisq;
	mrqmin ( &x[0]-1, &y[0]-1, &sig[0]-1, numPts, &a[0]-1, numCoeff, lista, mfit, covar, alpha, &chisq, fgauss, &alamda );
	int k=1;
	int itst=0;
	while (itst < 2) {
		if ( diagnostics ) {
			cout << "Iteration: " << k;
			cout << "  ";
			cout << "chi-squared: " << setprecision ( 4 ) << chisq;
			cout << "  ";
			cout << "alamda: " << setprecision ( 2 ) << alamda;
			cout << "<br />";
			cout << "m/z: " << setprecision ( 5 ) << a[1] << " ";
			cout << "intensity: " << setprecision ( 1 ) << a[0] << " ";
			cout << "width: " << setprecision ( 4 ) << a[2];
			cout << "<br />" << endl;
		}
		k++;
		ochisq=chisq;
		mrqmin ( &x[0]-1, &y[0]-1, &sig[0]-1, numPts, &a[0]-1, numCoeff, lista, mfit, covar, alpha, &chisq, fgauss, &alamda );
		if (chisq > ochisq)
			itst=0;
		else if (fabs(ochisq-chisq) < 0.1)
			itst++;
	}
	alamda=0.0;
	mrqmin ( &x[0]-1, &y[0]-1, &sig[0]-1, numPts, &a[0]-1, numCoeff, lista, mfit, covar, alpha, &chisq, fgauss, &alamda );
	if ( diagnostics ) {
		cout << "<br />Uncertainties:<br />";
		for ( i = 1 ; i <= numCoeff ; i += 3 ) {
			cout << "m/z: ";
			cout << setprecision ( 8 ) << a[i] << " +/- ";
			cout << setprecision ( 3 ) << sqrt(covar[i+1][i+1]) << "<br />";

			cout << "intensity: ";
			cout << setprecision ( 4 ) << a[i-1] << " +/- ";
			cout << setprecision ( 3 ) << sqrt(covar[i][i]) << "<br />";

			cout << "width: ";
			cout << setprecision ( 4 ) << a[i+1] << " +/- ";
			cout << setprecision ( 3 ) << sqrt(covar[i+2][i+2]) << "<br />";

			cout << "resolution: ";
			using nr::sqrtTwo;
			double factor = sqrtTwo * a[i] / 2.35;
			double resolution = factor / a[i+1];
			cout << setprecision ( 4 ) << resolution << " +/- ";
			cout << setprecision ( 3 ) << sqrt(covar[i+2][i+2]) * resolution / a[i+1] << "<br />";

			cout << "area: ";
			using nr::sqrtPi;
			cout << setprecision ( 4 ) << a[i-1] * a[i+1] * sqrtPi;
			cout << "<br />";
			cout << "<br />";
		}
		cout << "<br />";
	}
	free_nrmatrix(alpha,1,numCoeff,1,numCoeff);
	free_nrmatrix(covar,1,numCoeff,1,numCoeff);
	free_inrvector(lista,1,numCoeff);
	return a;
}
void PeakFit::drawGraph ( GraphData& graphData, bool sorted )
{
	SpectrumGraph s ( "pr_graph.par.txt" );
	s.drawGraph ( cout, graphData, sorted );
	cout << "<br />" << endl;
	cout << "<br />" << endl;
}

bool PeakFitData::reportPeakIntensity = true;
bool PeakFitData::reportPeakSNR = true;
bool PeakFitData::reportPeakResolution = true;
bool PeakFitData::reportPeakCSIntensity = false;
bool PeakFitData::reportPeakFWHM = true;
bool PeakFitData::reportPeakArea = true;
bool PeakFitData::reportPeakCSArea = false;
bool PeakFitData::reportNoiseMean = true;
bool PeakFitData::reportStdDev = true;
bool PeakFitData::reportFormulaString = false;
double PeakFitData::snrThreshold = 0.0;
double PeakFitData::areaThreshold = 0.0;
double PeakFitData::intensityThreshold = 0.0;

PeakFitData::~PeakFitData () {}
int PeakFitData::getNoiseColspan ()
{
	int colspan = 0;
	if ( reportNoiseMean )	colspan++;
	if ( reportStdDev )		colspan++;
	return colspan;
}
void PeakFitData::printHTMLNoiseHeader ( ostream& os, const string& styleID )
{
	if ( reportNoiseMean )	tableHeader ( os, "Noise Mean", styleID );
	if ( reportStdDev )		tableHeader ( os, "Noise St Dev", styleID );
}
void PeakFitData::printHTMLNoiseBlankLine ( ostream& os, const string& styleID )
{
	if ( reportNoiseMean )	tableEmptyCell ( os, styleID );
	if ( reportStdDev )		tableEmptyCell ( os, styleID );
}
void PeakFitData::printHTMLNoiseLine ( ostream& os, const string& styleID ) const
{
	if ( reportNoiseMean )	tableCellSigFig ( os, noiseMean, 3, false, styleID );
	if ( reportStdDev )		tableCellSigFig ( os, noiseStDev, 3, false, styleID );
}
void PeakFitData::printDelimitedNoiseHeader ( ostream& os )
{
	if ( reportNoiseMean )	delimitedHeader ( os, "Noise Mean" );
	if ( reportStdDev )		delimitedHeader ( os, "Noise St Dev" );
}
void PeakFitData::printDelimitedNoiseBlankLine ( ostream& os )
{
	if ( reportNoiseMean )	delimitedEmptyCell ( os );
	if ( reportStdDev )		delimitedEmptyCell ( os );
}
void PeakFitData::printDelimitedNoiseLine ( ostream& os ) const
{
	if ( reportNoiseMean )	delimitedCellSigFig	( os, noiseMean, 3 );
	if ( reportStdDev )		delimitedCellSigFig	( os, noiseStDev, 3 );
}
ParentData::ParentData ( const XYData& xyData, double monoMass, int charge, ElementalFormula* ef, bool efFlag, double inputResolution, int numPeaks )
{
	ColouredGraphData* graphData = 0;
	bool graphFlag = PeakFit::getGraphs ();
	if ( graphFlag ) graphData = new ColouredGraphData ( xyData );
	PeakFit::getNoise ( xyData, noiseMean, noiseStDev );
	DoubleVectorVector coeff = PeakFit::getCoefficients ( xyData, graphData, monoMass, genMax ( numPeaks, 3 ), charge, inputResolution, noiseMean, noiseStDev, snrThreshold );
	for ( DoubleVectorVectorSizeType i = 0 ; i < coeff.size () ; i++ ) {
		double width = coeff [i][2];
		mOverZ.push_back ( coeff [i][1] );
		intensity.push_back ( coeff [i][0] );
		fwhm.push_back ( 2.35 * width / nr::sqrtTwo );
		resolution.push_back ( mOverZ [i] / fwhm [i] );
		area.push_back ( intensity [i] * width * nr::sqrtPi );
		snr.push_back ( intensity [i] / noiseStDev );
	}
	auto_ptr <IsotopePeakStats> ips;
	if ( ef && efFlag )	{
		ips = auto_ptr <IsotopePeakStats> ( new IsotopePeakStats ( ef->getFormula (), charge ) );
		fString = ef->getFormula ();
	}
	else {
		ips = auto_ptr <IsotopePeakStats> ( new IsotopePeakStats ( monoMass, charge ) );
		fString = "";
	}
	ch = charge;
	double maximumProbability = ips->getGroupMaximumProbability ();
	for ( DoubleVectorVectorSizeType j = 0 ; j < coeff.size () ; j++ ) {
		theoreticalPercentMax.push_back ( ips->getProbability ( j ) * 100.0 / maximumProbability );
	}
	if ( graphFlag ) {
		PeakFit::drawGraph ( *graphData, false );
		delete graphData;
	}
}
void ParentData::printHTML ( ostream& os ) const
{
	if ( !mOverZ.empty () ) {
		tableStart ( os, true );
		double maxArea = *(std::max_element ( area.begin (), area.end () ));
		double maxIntensity = *(std::max_element ( intensity.begin (), intensity.end () ));
		for ( DoubleVectorSizeType i = 0 ; i < mOverZ.size () ; i++ ) {
			if ( i == 0 ) {
				tableRowStart ( os );
					tableHeader ( os, "m/z" );
					tableHeader ( os, "Int" );
					tableHeader ( os, "%" );
					tableHeader ( os, "Resolution" );
					tableHeader ( os, "FWHM" );
					tableHeader ( os, "Area" );
					tableHeader ( os, "%" );
					tableHeader ( os, "Theor %" );
					tableHeader ( os, "Noise Mean" );
					tableHeader ( os, "Noise St Dev" );
				tableRowEnd ( os );
			}
			tableRowStart ( os );
				tableCell ( os, mOverZ [i], 4 );
				tableCellSigFig ( os, intensity [i], 4 );
				if ( maxIntensity == 0.0 ) tableCell ( os, "---" );
				else tableCell ( os, intensity [i] * 100.0 / maxIntensity, 1 );
				tableCell ( os, resolution [i], 0 );
				tableCellSigFig ( os, fwhm [i], 4 );
				tableCellSigFig ( os, area [i], 4 );
				if ( maxArea == 0.0 ) tableCell ( os, "---" );
				else tableCell ( os, area [i] * 100.0 / maxArea, 1 );
				tableCell ( os, theoreticalPercentMax [i], 1 );
				tableCellSigFig ( os, noiseMean, 3 );
				tableCellSigFig ( os, noiseStDev, 3 );
			tableRowEnd ( os );
		}
		tableEnd ( os );
		if ( !fString.empty () ) {
			os << "<p>" << endl;
			MSIsotopeLink isotopeLink;
			startJavascript ( os );
			isotopeLink.printHTML ( os );
			endJavascript ( os );
			os << "<b>Elemental Composition: </b>";
			ElementalFormula efTemp ( fString ); 
			isotopeLink.write ( os, efTemp, ch );
			os << "<br />" << endl;
			CosSimilarity cs ( fString, ch, intensity, area );
			cs.printHTML ( os, 1 );
			os << "</p>" << endl;
		}
	}
}
void ParentData::printHTMLMOverZ ( ostream& os, DoubleVectorVectorSizeType peakNumber, const string& styleID ) const
{
}
void ParentData::printHTMLLine ( ostream& os, DoubleVectorVectorSizeType peakNumber, const string& styleID ) const
{
	if ( peakNumber < mOverZ.size () ) {
		double cosSimilarityIntensity = -100.0;
		double cosSimilarityArea = -100.0;
		if ( ( reportPeakCSIntensity || reportPeakCSArea ) && !fString.empty () ) {
			CosSimilarity cs ( fString, ch, intensity, area );
			cosSimilarityIntensity = cs.getIntensity ();
			cosSimilarityArea = cs.getArea ();
		}
		if ( reportPeakIntensity )	tableCellSigFig ( os, intensity [peakNumber], 4, false, styleID );
		if ( reportPeakSNR )		tableCellSigFig ( os, snr [peakNumber], 4, false, styleID );
		if ( reportPeakResolution )	tableCell		( os, resolution [peakNumber], 0, false, styleID );
		if ( reportPeakCSIntensity ) {
			if ( cosSimilarityIntensity == -100.0 )	tableEmptyCell ( os, styleID );
			else									tableCellSigFig ( os, cosSimilarityIntensity, 3, false, styleID );
		}
		if ( reportPeakFWHM )		tableCellSigFig ( os, fwhm [peakNumber], 4, false, styleID );
		if ( reportPeakArea )		tableCellSigFig ( os, area [peakNumber], 4, false, styleID );
		if ( reportPeakCSArea ) {
			if ( cosSimilarityArea == -100.0 )	tableEmptyCell ( os, styleID );
			else								tableCellSigFig ( os, cosSimilarityArea, 3, false, styleID );
		}
	}
	else printHTMLBlankLine ( os, styleID );
}
void ParentData::printDelimitedLine ( ostream& os, DoubleVectorVectorSizeType peakNumber ) const
{
	if ( peakNumber < mOverZ.size () ) {
		double cosSimilarityIntensity = -100.0;
		double cosSimilarityArea = -100.0;
		if ( ( reportPeakCSIntensity || reportPeakCSArea ) && !fString.empty () ) {
			CosSimilarity cs ( fString, ch, intensity, area );
			cosSimilarityIntensity = cs.getIntensity ();
			cosSimilarityArea = cs.getArea ();
		}
		if ( reportPeakIntensity )	delimitedCellSigFig ( os, intensity [peakNumber], 4 );
		if ( reportPeakSNR )		delimitedCellSigFig ( os, snr [peakNumber], 4 );
		if ( reportPeakResolution )	delimitedCellSigFig ( os, resolution [peakNumber], 0 );
		if ( reportPeakCSIntensity ) {
			if ( cosSimilarityIntensity == -100.0 )	delimitedEmptyCell ( os );
			else									delimitedCellSigFig ( os, cosSimilarityIntensity, 3 );
		}
		if ( reportPeakFWHM )		delimitedCellSigFig ( os, fwhm [peakNumber], 4 );
		if ( reportPeakArea )		delimitedCellSigFig ( os, area [peakNumber], 4 );
		if ( reportPeakCSArea ) {
			if ( cosSimilarityArea == -100.0 )	delimitedEmptyCell ( os );
			else								delimitedCellSigFig ( os, cosSimilarityArea, 3 );
		}
	}
	else printDelimitedBlankLine ( os );
}
void ParentData::printHTMLBlankLine ( ostream& os, const string& styleID )
{
	if ( reportPeakIntensity )	tableEmptyCell ( os, styleID );
	if ( reportPeakSNR )		tableEmptyCell ( os, styleID );
	if ( reportPeakResolution )	tableEmptyCell ( os, styleID );
	if ( reportPeakCSIntensity )tableEmptyCell ( os, styleID );
	if ( reportPeakFWHM )		tableEmptyCell ( os, styleID );
	if ( reportPeakArea )		tableEmptyCell ( os, styleID );
	if ( reportPeakCSArea )		tableEmptyCell ( os, styleID );
}
void ParentData::printDelimitedBlankLine ( ostream& os )
{
	if ( reportPeakIntensity ) 	delimitedEmptyCell ( os );
	if ( reportPeakSNR )		delimitedEmptyCell ( os );
	if ( reportPeakResolution )	delimitedEmptyCell ( os );
	if ( reportPeakCSIntensity )delimitedEmptyCell ( os );
	if ( reportPeakFWHM )		delimitedEmptyCell ( os );
	if ( reportPeakArea )		delimitedEmptyCell ( os );
	if ( reportPeakCSArea )		delimitedEmptyCell ( os );
}
int ParentData::getColspan ()
{
	int colspan = 0;
	if ( reportPeakIntensity )	colspan++;
	if ( reportPeakSNR )		colspan++;
	if ( reportPeakResolution )	colspan++;
	if ( reportPeakCSIntensity )colspan++;
	if ( reportPeakFWHM )		colspan++;
	if ( reportPeakArea )		colspan++;
	if ( reportPeakCSArea )		colspan++;
	return colspan;
}
void ParentData::printHTMLHeader ( ostream& os, const string& styleID )
{
	if ( reportPeakIntensity )	tableHeader ( os, "Intensity", styleID );
	if ( reportPeakSNR )		tableHeader ( os, "SNR", styleID );
	if ( reportPeakResolution )	tableHeader ( os, "Resolution", styleID );
	if ( reportPeakCSIntensity )tableHeader ( os, "CS", styleID );
	if ( reportPeakFWHM )		tableHeader ( os, "FWHM", styleID );
	if ( reportPeakArea )		tableHeader ( os, "Area", styleID );
	if ( reportPeakCSArea )		tableHeader ( os, "CS", styleID );
}
void ParentData::printDelimitedHeader ( ostream& os )
{
	if ( reportPeakIntensity )	delimitedHeader ( os, "Intensity" );
	if ( reportPeakSNR )		delimitedHeader ( os, "SNR" );
	if ( reportPeakResolution )	delimitedHeader ( os, "Resolution" );
	if ( reportPeakCSIntensity )delimitedHeader ( os, "CS" );
	if ( reportPeakFWHM )		delimitedHeader ( os, "FWHM" );
	if ( reportPeakArea )		delimitedHeader ( os, "Area" );
	if ( reportPeakCSArea )		delimitedHeader ( os, "CS" );
}
bool QuantitationData::reportActualLightHeavyIntensityRatio = false;
bool QuantitationData::reportActualLightHeavyAreaRatio = false;
const char QuantitationData::RATIO_GREATER_THAN = '>';
const char QuantitationData::RATIO_LESS_THAN = '<';
const char QuantitationData::RATIO_NOT_CALCULATED = '-';
const char QuantitationData::RATIO_OK = 'K';
const char QuantitationData::RATIO_HIGH = 'H';
const char QuantitationData::RATIO_LOW = 'L';
QuantitationData::QuantitationData ( int numQuanPeaks ) :
	lightHeavyIntRatioType ( 1 ),
	lightHeavyAreaRatioType ( 1 ),
	lightHeavyIntRatio ( 1 ),
	lightHeavyAreaRatio ( 1 )
{
	for ( int i = 0 ; i < 1 ; i++ ) {
		for ( int j = 0 ; j < numQuanPeaks ; j++ ) {
			lightHeavyIntRatioType [i].push_back ( RATIO_NOT_CALCULATED );
			lightHeavyAreaRatioType [i].push_back ( RATIO_NOT_CALCULATED );
			lightHeavyIntRatio [i].push_back ( 0.0 );
			lightHeavyAreaRatio [i].push_back ( 0.0 );
		}
	}
}
double QuantitationData::getQuanRatio ( double ratio, char ratioType )
{
	if ( ratioType != RATIO_NOT_CALCULATED ) {
		if ( ratio >= 100.0 || ratioType == RATIO_HIGH ) {
			return 100.0;
		}
		else if ( ratio <= 0.01 || ratioType == RATIO_LOW ) {
			return 0.01;
		}
		return ratio;
	}
	return 0.0;
}
void QuantitationData::setRatio ( char& ratioType, double& ratio, double threshold, double refValue, double value )
{
	if ( threshold > 0 ) {
		if ( refValue < threshold && value < threshold ) ratioType = RATIO_NOT_CALCULATED;
		else if ( refValue < threshold )setLTRatio ( ratioType, ratio, threshold, value );
		else if ( value < threshold )	setGTRatio ( ratioType, ratio, refValue, threshold );
		else setRatio ( ratioType, ratio, refValue, value );
	}
	else setRatio ( ratioType, ratio, refValue, value );
}
void QuantitationData::setRatio ( char& ratioType, double& ratio, double refValue, double value )
{
	if ( refValue <= 0.0 && value <= 0.0 )	ratioType = RATIO_NOT_CALCULATED;
	else if ( refValue <= 0.0 )				ratioType = RATIO_LOW;
	else if ( value <= 0.0 )				ratioType = RATIO_HIGH;
	else {
		ratioType = RATIO_OK;
		ratio = refValue / value;
	}
}
void QuantitationData::setGTRatio ( char& ratioType, double& ratio, double refValue, double value )
{
	ratioType = RATIO_GREATER_THAN;
	ratio = refValue / value;
}
void QuantitationData::setLTRatio ( char& ratioType, double& ratio, double refValue, double value )
{
	ratioType = RATIO_LESS_THAN;
	ratio = refValue / value;
}
bool QuantitationData::outputQuanRatio ( ostream& os, const string& searchName, double ratio, char ratioType, int numRepeats )
{
	if ( ratioType != RATIO_NOT_CALCULATED ) {
		os << searchName;
		os << '\t';
		if ( ratio >= 100.0 || ratioType == RATIO_HIGH ) {
			os << 100.0;
		}
		else if ( ratio <= 0.01 || ratioType == RATIO_LOW ) {
			os << 0.01;
		}
		else os << ratio;
		os << '\t';
		os << ratioType;
		os << '\t';
		os << numRepeats;
		os << endl;
		return true;
	}
	return false;
}
void QuantitationData::printHTMLRatio ( ostream& os, double ratio, char ratioType, const string& styleID )
{
	if ( ratioType == RATIO_NOT_CALCULATED ) {
		tableCell ( os, "---", false, true, styleID );
	}
	else if ( ratioType == RATIO_HIGH ) {
		tableCell ( os, "high", false, true, styleID );
	}
	else if ( ratioType == RATIO_LOW ) {
		tableCell ( os, "low", false, true, styleID );
	}
	else {
		ostringstream ostr;
		if ( ratioType == RATIO_GREATER_THAN ) {
			ostr << '>';
		}
		if ( ratioType == RATIO_LESS_THAN ) {
			ostr << '<';
		}
		genPrintSigFig ( ostr, ratio, 3 );
		tableCell ( os, ostr.str (), false, true, styleID );
	}
}
void QuantitationData::printDelimitedRatio ( ostream& os, double ratio, char ratioType )
{
	if ( ratioType == RATIO_NOT_CALCULATED ) {
		delimitedCell ( os, "---" );
	}
	else if ( ratioType == RATIO_HIGH ) {
		delimitedCell ( os, "high" );
	}
	else if ( ratioType == RATIO_LOW ) {
		delimitedCell ( os, "low" );
	}
	else {
		ostringstream ostr;
		if ( ratioType == RATIO_GREATER_THAN ) {
			ostr << '>';
		}
		if ( ratioType == RATIO_LESS_THAN ) {
			ostr << '<';
		}
		genPrintSigFig ( ostr, ratio, 3 );
		delimitedCell ( os, ostr.str () );
	}
}
void MSParentLink::putCGI ( ostream& os ) const
{
	const ParameterList* params = ProgramLink::getParams ();
	printCGIString ( os, "version", Version::instance ().getVersion () );
	params->copyToCGI ( os, "raw_type" );
	params->copyToCGI ( os, "quan_type" );
	IsotopePurity::copyToCGI ( os, params );
	if ( isQuanMSMS ( params->getStringValue ( "quan_type", "" ) ) ) {
		params->copyToCGI ( os, "purity_correction" );
		params->copyToCGI ( os, "reporter_ion_window" );
	}
	params->copyToCGI ( os, "area_threshold" );
	params->copyToCGI ( os, "intensity_threshold" );
	params->copyToCGI ( os, "resolution" );
	printCGI ( os, "num_peaks", 3 );
}
void MSParentLink::write ( ostream& os, const SpecID& specID, double mOverZ, int charge, double rtIntervalStart, double rtIntervalEnd, double snrThreshold, const string& formula, const string& nterm, const string& cterm, const string& nloss, const string& searchKey, const string& systematicError, const string& toleranceUnits ) const
{
	ProgramLink::openLink ( os, "msdisplayLink", -1 );
	specID.putCGI ( os );
	printCGIString ( os, "m_over_z", gen_ftoa ( mOverZ, "%.4f" ) );
	printCGI ( os, "charge", charge );
	printCGI ( os, "rt_int_start", rtIntervalStart, 4 );
	printCGI ( os, "rt_int_end", rtIntervalEnd, 4 );
	printCGI ( os, "snr_threshold", snrThreshold, 4 );
	printCGIString ( os, "search_key", searchKey );
	printCGIString ( os, "formula", formula );
	if ( !nterm.empty () ) printCGIString ( os, "nterm", nterm );
	if ( !cterm.empty () ) printCGIString ( os, "cterm", cterm );
	if ( !nloss.empty () ) printCGIString ( os, "nloss", nloss );
	printCGI ( os, "msms_parent_mass_systematic_error", systematicError );
	printCGIString ( os, "msms_parent_mass_tolerance_units", toleranceUnits );
	os << "\\\">";
	genPrint ( os, mOverZ, 4 );
	ProgramLink::closeLink ( os );
}
void MSParentLink::write ( ostream& os, const SpecID& specID, double mOverZ, int charge, double rtIntervalStart, double rtIntervalEnd, double snrThreshold, const string& formula, const string& nterm, const string& cterm, const string& nloss, const string& formula2, const string& nterm2, const string& cterm2, const string& linkSearchType, const string& searchKey, const string& systematicError, const string& toleranceUnits ) const
{
	ProgramLink::openLink ( os, "msdisplayLink", -1 );
	specID.putCGI ( os );
	printCGIString ( os, "m_over_z", gen_ftoa ( mOverZ, "%.4f" ) );
	printCGI ( os, "charge", charge );
	printCGI ( os, "rt_int_start", rtIntervalStart, 4 );
	printCGI ( os, "rt_int_end", rtIntervalEnd, 4 );
	printCGI ( os, "snr_threshold", snrThreshold, 4 );
	printCGIString ( os, "search_key", searchKey );
	printCGIString ( os, "formula", formula );
	if ( !nterm.empty () ) printCGIString ( os, "nterm", nterm );
	if ( !cterm.empty () ) printCGIString ( os, "cterm", cterm );
	printCGIString ( os, "formula2", formula2 );
	if ( !nterm2.empty () ) printCGIString ( os, "nterm", nterm2 );
	if ( !cterm2.empty () ) printCGIString ( os, "cterm", cterm2 );
	if ( !nloss.empty () ) printCGIString ( os, "nloss", nloss );
	printCGIString ( os, "link_search_type", linkSearchType );
	printCGI ( os, "msms_parent_mass_systematic_error", systematicError );
	printCGIString ( os, "msms_parent_mass_tolerance_units", toleranceUnits );
	os << "\\\">";
	genPrint ( os, mOverZ, 4 );
	ProgramLink::closeLink ( os );
}
void MSParentLink::printHTML ( ostream& os ) const
{
	os << "msdisplayLink" << "=\"";
	os << ProgramLink::getURLStart ( "msdisplay" );
	os << "?";
	putCGI ( os );
	os << "\";\n";
}
#endif
