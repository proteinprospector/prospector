/******************************************************************************
*                                                                             *
*  Program    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_histogram.cpp                                              *
*                                                                             *
*  Created    : February 26th 2004                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2004-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lu_app_gr.h>
#include <lu_histogram.h>
#include <lu_parent.h>
#include <lu_delim.h>
#include <lu_inst.h>
#include <lu_r_plot.h>
#include <lu_param_list.h>
#include <lu_getfil.h>
using std::sort;
using std::istream;
using std::ostream;
using std::endl;
using std::string;
using nr::pi;

class Gaussian {
	double intensity;
	double mean;
	double width;
	double stdev;
	double twoStDevSquared;
public:
	Gaussian ( double intensity, double mean, double width ) :
		intensity ( intensity ), mean ( mean ), width ( width ),
		stdev ( width / nr::sqrtTwo ), twoStDevSquared ( width * width ) {}
	double getInterval () const { return ( getStdev3 () - getStdevm3 () ) / 30.0; }
	double getArea () const { return intensity * width * nr::sqrtPi; }
	double getMean () const { return mean; }
	double getStdev () const { return stdev; }
	double getStdev1 () const { return mean + stdev; }
	double getStdev2 () const { return mean + ( 2.0 * stdev ); }
	double getStdev3 () const { return mean + ( 3.0 * stdev ); }
	double getStdevm3 () const { return mean - ( 3.0 * stdev ); }
	void printStatsDelimited ( ostream& os ) const;
	void printStatsHTML ( ostream& os ) const;
	void printStatsHTML2 ( ostream& os ) const;
	double getY ( double x ) const
	{
		double top = x - mean;
		top *= top;
		return intensity * exp ( - top / twoStDevSquared );
	};
	double prob ( double x ) const
	{
		double factor = 1.0 / ( stdev * nr::sqrtTwoPi );
		double top = x - mean;
		top *= top;
		return factor * exp ( - top / twoStDevSquared );
	}
};
void Gaussian::printStatsDelimited ( ostream& os ) const
{
	delimitedRowStart ( os );
		delimitedCell ( os, "Mean" );
		delimitedCellSigFig ( os, mean, 4 );
	delimitedRowEnd ( os );
	delimitedRowStart ( os );
		delimitedCell ( os, "St Dev" );
		delimitedCellSigFig ( os, stdev, 4 );
	delimitedRowEnd ( os );
	delimitedRowStart ( os );
		delimitedCell ( os, "2 St Dev" );
		delimitedCellSigFig ( os, getStdev2 (), 4 );
	delimitedRowEnd ( os );
	delimitedRowStart ( os );
		delimitedCell ( os, "3 St Dev" );
		delimitedCellSigFig ( os, getStdev3 (), 4 );
	delimitedRowEnd ( os );
}
void Gaussian::printStatsHTML ( ostream& os ) const
{
	ParameterList::printDoubleHTMLSigFig ( os, "Mean", mean, 4 );
	ParameterList::printDoubleHTMLSigFig ( os, "St Dev", stdev, 4 );
	ParameterList::printDoubleHTMLSigFig ( os, "2 St Dev", getStdev2 (), 4 );
	ParameterList::printDoubleHTMLSigFig ( os, "3 St Dev", getStdev3 (), 4 );
}
void Gaussian::printStatsHTML2 ( ostream& os ) const
{
	ParameterList::printDoubleHTMLSigFig ( os, "Mean", mean, 4 );
	ParameterList::printDoubleHTMLSigFig ( os, "St Dev", stdev, 4 );
	ParameterList::printDoubleHTMLSigFig ( os, "St Dev * 3", stdev * 3, 4 );
}
class MixtureGaussian {
	Gaussian negDistribution;
	Gaussian posDistribution;
	double noiseMean;
	double priorNeg;
	double priorPos;
public:
	MixtureGaussian ( const Gaussian& negDistribution, const Gaussian& posDistribution, double noiseMean );
	double prob ( double x ) const
	{
		double a = posDistribution.prob ( x ) * priorPos;
		double b = negDistribution.prob ( x ) * priorNeg;
		return a / ( a + b );
	}
	double getInterval () const { return negDistribution.getInterval (); };
	double getY ( double x ) const { return negDistribution.getY ( x ) + posDistribution.getY ( x ) + noiseMean; };
	void printStatsDelimited ( ostream& os ) const;
	void printStatsHTML ( ostream& os ) const;
};
MixtureGaussian::MixtureGaussian ( const Gaussian& negDistribution, const Gaussian& posDistribution, double noiseMean ) :
	negDistribution ( negDistribution ), posDistribution ( posDistribution ), noiseMean ( noiseMean )
{
	priorNeg = negDistribution.getArea () / ( negDistribution.getArea () + posDistribution.getArea () );
	priorPos = 1.0 - priorNeg;
}
void MixtureGaussian::printStatsDelimited ( ostream& os ) const
{
	negDistribution.printStatsDelimited ( os );
	posDistribution.printStatsDelimited ( os );

	delimitedRowStart ( os );
		delimitedCell ( os, "Prior -" );
		delimitedCellSigFig ( os, priorNeg, 4 );
	delimitedRowEnd ( os );
	delimitedRowStart ( os );
		delimitedCell ( os, "Prior +" );
		delimitedCellSigFig ( os, priorPos, 4 );
	delimitedRowEnd ( os );
}
void MixtureGaussian::printStatsHTML ( ostream& os ) const
{
	negDistribution.printStatsHTML ( os );
	os << "<p />" << endl;
	posDistribution.printStatsHTML ( os );
	os << "<p />" << endl;

	ParameterList::printDoubleHTMLSigFig ( os, "Prior -", priorNeg, 4 );
	ParameterList::printDoubleHTMLSigFig ( os, "Prior +", priorPos, 4 );
	os << "<p />" << endl;
}
Histogram::Histogram () :
	negDistribution ( 0 )
{
}
Histogram::~Histogram ()
{
	delete negDistribution;
}
void Histogram::makeHistogram () const
{
	sort ( val.begin (), val.end () );
	double min = val.front ();
	double max = val.back ();
	if ( max == std::numeric_limits<double>::infinity() ) return;
	int nbins = 60;
	if ( val.size () > 60 ) {
		double mean, sdev, var;
		moment ( &val[0]-1, val.size (), &mean, &sdev, &var );
		double binWidth = 3.49 * sdev * pow ( val.size (), - ( 1.0 / 3.0 ) );
		nbins = ( max - min ) / binWidth;
		nbins = genMax ( 60, nbins );
	}
	double interval = ( max - min ) / nbins;
	if ( interval == 0.0 ) interval = 0.000000001;
	double halfInterval = interval / 2;
	int n = 0;
	int num = 0;
	double end = min + interval;
	for ( DoubleVectorSizeType i = 0 ; i < val.size () ; ) {
		if ( val [i] < end ) {
			num++;
			i++;
		}
		else {
			xyData.add ( end - halfInterval, num );
			end += interval;
			num = 0;
			n++;
		}
	}
	xyData.add ( end - halfInterval, num );
}
void Histogram::compute () const
{
	if ( !val.empty () ) {
		makeHistogram ();
	}
}
void Histogram::drawGraph ( ostream& os ) const
{
	if ( negDistribution == 0 ) compute ();

	GraphData graphData ( xyData );

	if ( graphData.size () ) {
		SpectrumGraph s1 ( "hist.par.txt" );
		s1.drawGraph ( os, graphData, false );
		os << "<p />" << endl;
	}
}
istream& operator>> ( istream& is, Histogram& h )
{
	is >> h.xyData;
	return is;
}
ostream& operator<< ( ostream& os, const Histogram& h )
{
	os << h.xyData;
	return os;
}
ErrorHistogram::ErrorHistogram () :
	Histogram ()
{
}
ErrorHistogram::~ErrorHistogram () {}
void ErrorHistogram::compute () const
{
	if ( !val.empty () ) {
		makeHistogram ();
	}
}
void ErrorHistogram::drawGraph ( ostream& os ) const
{
	compute ();
	if ( xyData.size () ) {
		os << "<p>" << endl;
		SpectrumGraph s1 ( "error_hist.par.txt" );
		ColouredGraphData gd ( xyData, 1 );
		s1.drawGraph ( os, gd, false );
		os << "</p>" << endl;
	}
}
string SurvivalHistogram::expectationMethod = "None";
int SurvivalHistogram::minUsedPeptides = ExpectationParameters::instance ().getMinUsedPeptides ();
double SurvivalHistogram::tailPercent = ExpectationParameters::instance ().getTailPercent ();

SurvivalHistogram::SurvivalHistogram ( size_t sizeLimit ) :
	Histogram (),
	a ( 0.0 ),
	b ( 0.0 ),
	size ( 0 ),
	sizeLimit ( sizeLimit )
{
}
SurvivalHistogram::~SurvivalHistogram () {}
void SurvivalHistogram::compute ( int numSavedPeptides ) const
{
	int numUsedPeptides = val.size () - numSavedPeptides;
	if ( expectationMethod == "Linear Tail Fit" ) numUsedPeptides -= 20;
	if ( numUsedPeptides > 0 ) {
		makeHistogram ();
		if ( expectationMethod == "Linear Tail Fit" ) {
			if ( numUsedPeptides >= minUsedPeptides ) {
				double tot = 0.0;
				double actTot = 0.0;
				XYData xyTemp;
				for ( int i = xyData.size () ; i-- ; ) {
					tot += xyData.y ( i );
					if ( actTot > numUsedPeptides / tailPercent ) break;
					actTot += xyData.y ( i );
					if ( tot > (numSavedPeptides + 20) ) {
						xyTemp.add ( xyData.x ( i ), actTot / numUsedPeptides );
					}
				}
				for ( int j = 0 ; j < xyTemp.size () ; j++ ) {
					if ( xyTemp.y ( j ) ) {
						double x = xyTemp.x ( j );			// Method 1
						double y = log10 ( xyTemp.y ( j ) );
						xySurv.add ( x, y );
					}
				}
				xySurv.linearRegression ( &b, &a );
			}
		}
		else if ( expectationMethod == "Method of Moments" ) {	//	Method of moments
			if ( numUsedPeptides >= 500 ) {
				double mean, sdevError, var;
				moment ( &val[0]-1, numUsedPeptides, &mean, &sdevError, &var );
				b = ( sqrt (6.0 * var) ) / pi;
				a = mean - ( 0.57721566490153286060651209008240243104 * b );
			}
		}
		else if ( expectationMethod == "Closed Form Max Likelihood" ) {	//	Closed Form Maximum Likelihood
			if ( numUsedPeptides >= 500 ) {
				double mean, sdevError, var;
				sort ( val.begin (), val.begin () + numSavedPeptides );
				moment ( &val[0]-1, numUsedPeptides, &mean, &sdevError, &var );
				double sum1 = 0.0;
				for ( int i = 0 ; i < numUsedPeptides ; i++ ) {
					sum1 += val [i] * log ( ( (i+1) - 0.5 ) / ( numUsedPeptides + 0.5 )  );
				}
				sum1 /= numUsedPeptides;
				b = mean + sum1;
				double sum2 = 0.0;
				for ( int j = 0 ; j < numUsedPeptides ; j++ ) {
					sum2 += exp ( - val [j] / b  );
				}
				a = - b * log ( sum2 / numUsedPeptides );
			}
		}
	}
}
void SurvivalHistogram::printHTML ( ostream& os, double score, int numSavedPeptides ) const
{
	if ( expectationMethod != "None" ) {
		score /= SCORE_TYPE_MULTIPLIER;
		if ( a == 0.0 ) {
			compute ( numSavedPeptides );
		}
		ColouredGraphData graphData ( xyData, 1 );

		if ( graphData.size () ) {
			SpectrumGraph s1 ( "hist.par.txt" );
			s1.drawGraph ( os, graphData, false );
			os << "<p />" << endl;
		}
		if ( a != 0.0 ) {
			if ( expectationMethod == "Linear Tail Fit" ) {
				if ( RPlot::getRFlag () ) {
					RPlot rplot ( "scatterSurvival.R" );
					GenOFStream ofs ( rplot.getDataFileFullPath () );
					for ( int ii = 0 ; ii < xySurv.size () ; ii++ ) {
						genPrintSigFig ( ofs, xySurv.x ( ii ), 4 );
						ofs << " ";
						genPrintSigFig ( ofs, xySurv.y ( ii ), 4 );
						ofs << endl;
					}
					ofs.close ();
					rplot.printImageAndLink ( os );
				}
				double y = (a * score) + b;
				os << "expectation value = ";
				genPrintSigFig ( os, pow ( 10.0, y ) * size, 3 );
				os << "<br />" << endl;
			}
			else {
				os << "expectation value = ";
				genPrintSigFig ( os, ( 1.0 - exp ( -exp ( (- 1 / b) * (score - a) ) ) ) * size, 3 );
				os << "<br />" << endl;
			}
			os << "num peptides considered = " << size << "<br />" << endl;
		}
		else {
			os << "num peptides considered = " << size << "<br />" << endl;
		}
	}
}
void SurvivalHistogram::printXML ( ostream& os, int numSavedPeptides ) const
{
	if ( expectationMethod != "None" ) {
		if ( a == 0.0 ) {
			compute ( numSavedPeptides );
		}
		delimitedRowStart ( os );
			delimitedCellSigFig ( os, a, 14 );
			delimitedCellSigFig ( os, b, 14 );
			delimitedCell ( os, (int)size );
		delimitedRowEnd ( os );
	}
}
void SurvivalHistogram::init ( int numSavedPeptides ) const
{
	if ( expectationMethod != "None" ) {
		if ( a == 0.0 ) {
			compute ( numSavedPeptides );
		}
	}
}
double SurvivalHistogram::getEValue ( double score ) const
{
	if ( a != 0.0 ) {
		score /= SCORE_TYPE_MULTIPLIER;
		if ( expectationMethod == "Linear Tail Fit" ) {
			double y = (a * score) + b;
			return pow ( 10.0, y ) * size;
		}
		else
			return ( 1.0 - exp ( -exp ( (- 1 / b) * (score - a) ) ) ) * size;
	}
	else
		return -1.0;
}
ExpectationParameters::ExpectationParameters ()
{
	GenNameValueStream nvs ( MsparamsDir::instance ().getParamPath ( "expectation.txt" ) );
	nvs.getValue ( "min_used_peptides", minUsedPeptides );
	nvs.getValue ( "max_used_peptides", maxUsedPeptides );
	nvs.getValue ( "tail_percent", tailPercent );
}
ExpectationParameters& ExpectationParameters::instance ()
{
	static ExpectationParameters ep;
	return ep;
}
