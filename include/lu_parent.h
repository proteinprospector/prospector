/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_parent.h                                                   *
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

#ifndef __lu_parent_h
#define __lu_parent_h

#include <nr.h>
#include <lu_prog.h>
#include <lu_formula.h>

class SpecID;
class GraphData;
class ColouredGraphData;

class MSParentLink {
	void putCGI ( std::ostream& os ) const;
public:
	MSParentLink () {}
	void write ( std::ostream& os, const SpecID& specNumber, double mOverZ, int charge, double rtIntervalStart, double rtIntervalEnd, double snrThreshold, const std::string& formula, const std::string& nterm, const std::string& cterm, const std::string& nloss, const std::string& searchKey, const std::string& systematicError, const std::string& toleranceUnits ) const;
	void write ( std::ostream& os, const SpecID& specID, double mOverZ, int charge, double rtIntervalStart, double rtIntervalEnd, double snrThreshold, const std::string& formula, const std::string& nterm, const std::string& cterm, const std::string& nloss, const std::string& formula2, const std::string& nterm2, const std::string& cterm2, const std::string& linkSearchType, const std::string& searchKey, const std::string& systematicError, const std::string& toleranceUnits ) const;
	void printHTML ( std::ostream& os ) const;
};

class PeakFit {
	static bool diagnostics;
	static bool graphs;
public:
	static void getNoise ( const XYData& xyData, double& mean, double& stddev );
	static void getStats ( const XYData& xyData, double& mean, double& stddev );
	static DoubleVectorVector getCoefficients ( const XYData& xyData, ColouredGraphData* graphData, double monoMass, int numPeaks, int charge, double inputResolution, double noiseMean, double noiseStDev, double minSNR );
	static DoubleVector minimize ( XYData& xyData, double noiseWidth, double noiseValue, const DoubleVector& gues, int numDist = 1 );

	static void drawGraph ( GraphData& graphData, bool sorted );

	static void setDiagnostics ( bool d ) { diagnostics = d; }
	static void setGraphs ( bool g ) { graphs = g; }
	static bool getGraphs () { return graphs; }
};

class PeakFitData {
protected:
	double noiseMean;
	double noiseStDev;

	static bool reportPeakIntensity;
	static bool reportPeakSNR;
	static bool reportPeakResolution;
	static bool reportPeakCSIntensity;
	static bool reportPeakFWHM;
	static bool reportPeakArea;
	static bool reportPeakCSArea;
	static bool reportNoiseMean;
	static bool reportStdDev;
	static bool reportFormulaString;
	static double snrThreshold;
	static double areaThreshold;
	static double intensityThreshold;
	virtual void printHTMLMOverZ ( std::ostream& os, DoubleVectorVectorSizeType peakNumber, const std::string& styleID ) const = 0;
public:
	virtual ~PeakFitData () = 0;
	static int getNoiseColspan ();
	static void printHTMLNoiseHeader ( std::ostream& os, const std::string& styleID );
	void printHTMLNoiseLine ( std::ostream& os, const std::string& styleID ) const;
	static void printHTMLNoiseBlankLine ( std::ostream& os, const std::string& styleID );

	static void printDelimitedNoiseHeader ( std::ostream& os );
	static void printDelimitedNoiseBlankLine ( std::ostream& os );
	void printDelimitedNoiseLine ( std::ostream& os ) const;

	static void setReportNoiseMean	( bool f )	{ reportNoiseMean = f; }
	static void setReportStdDev		( bool f )	{ reportStdDev = f; }
	static void setReportPeakIntensity	( bool f )	{ reportPeakIntensity = f; }
	static void setReportPeakSNR		( bool f )	{ reportPeakSNR = f; }
	static void setReportPeakResolution	( bool f )	{ reportPeakResolution = f; }
	static void setReportPeakCSIntensity( bool f )	{ reportPeakCSIntensity = f; }
	static void setReportPeakFWHM		( bool f )	{ reportPeakFWHM = f; }
	static void setReportPeakArea		( bool f )	{ reportPeakArea = f; }
	static void setReportPeakCSArea		( bool f )	{ reportPeakCSArea = f; }
	static void setReportFormulaString	( bool f )	{ reportFormulaString = f; }

	static void setSNRThreshold ( double t ){ snrThreshold = t; }
	static void setAreaThreshold ( double t ){ areaThreshold = t; }
	static void setIntensityThreshold ( double t ){ intensityThreshold = t; }
	static double getSNRThreshold () { return snrThreshold; }
	static double getAreaThreshold () { return areaThreshold; }
	static double getIntensityThreshold () { return intensityThreshold; }

	virtual void printHTML ( std::ostream& os ) const = 0;

	virtual void printHTMLLine ( std::ostream& os, DoubleVectorVectorSizeType peakNumber, const std::string& styleID ) const = 0;
	virtual void printDelimitedLine ( std::ostream& os, DoubleVectorVectorSizeType peakNumber ) const = 0;
	virtual bool outputQuanResults ( std::ostream& os, const std::string& searchName, int numRepeats, bool area ) const = 0;
	virtual DoubleVector getAreaRatios () const = 0;
	virtual DoubleVector getIntensityRatios () const = 0;
};

class QuantitationData : public PeakFitData {
protected:
	static bool reportActualLightHeavyIntensityRatio;
	static bool reportActualLightHeavyAreaRatio;
	static const char RATIO_GREATER_THAN;
	static const char RATIO_LESS_THAN;
	static const char RATIO_NOT_CALCULATED;
	static const char RATIO_OK;
	static const char RATIO_HIGH;
	static const char RATIO_LOW;
	static double getQuanRatio ( double ratio, char ratioType );
	static bool outputQuanRatio ( std::ostream& os, const std::string& searchName, double ratio, char ratioType, int numRepeats );
	static void printHTMLRatio ( std::ostream& os, double ratio, char ratioType, const std::string& styleID );
	static void printDelimitedRatio ( std::ostream& os, double ratio, char ratioType );

	CharVectorVector lightHeavyIntRatioType;
	CharVectorVector lightHeavyAreaRatioType;
	DoubleVectorVector lightHeavyIntRatio;
	DoubleVectorVector lightHeavyAreaRatio;
	static void setRatio ( char& ratioType, double& ratio, double threshold, double refValue, double value );
	static void setRatio ( char& ratioType, double& ratio, double refValue, double value );
	static void setGTRatio ( char& ratioType, double& ratio, double refValue, double value );
	static void setLTRatio ( char& ratioType, double& ratio, double refValue, double value );
public:
	QuantitationData ( int numQuanPeaks );
	static void setReportActualLightHeavyIntensityRatio		( bool f )	{ reportActualLightHeavyIntensityRatio = f; }
	static void setReportActualLightHeavyAreaRatio			( bool f )	{ reportActualLightHeavyAreaRatio = f; }
	static bool getQuanReport ()
	{
		return	reportPeakIntensity ||
				reportPeakResolution ||
				reportPeakCSIntensity ||
				reportPeakSNR ||
				reportPeakFWHM ||
				reportActualLightHeavyIntensityRatio ||
				reportPeakArea ||
				reportPeakCSArea ||
				reportActualLightHeavyAreaRatio ||
				reportNoiseMean ||
				reportStdDev;
	}
	static bool getIntRatioReport () { return reportActualLightHeavyIntensityRatio; }
	static bool getAreaRatioReport () { return reportActualLightHeavyAreaRatio; }
};

class AACalculator;

class CosSimilarity {
	double cosSimilarityIntensity;
	double cosSimilarityArea;
public:
	CosSimilarity ( const std::string& formulaString, int charge, const DoubleVector& intensity, const DoubleVector& area );
	void printHTML ( std::ostream& os );
	void printHTML ( std::ostream& os, int num );
	double getIntensity () const { return cosSimilarityIntensity; }
	double getArea () const { return cosSimilarityArea; }
};

class CosSimilarityList {
	DoubleVector csi;
	DoubleVector csa;
	void printHTML ( std::ostream& os, const DoubleVector& c, const std::string& styleID ) const;
	void printDelimited ( std::ostream& os, const DoubleVector& c ) const;
public:
	CosSimilarityList ( const StringVector& formulaString, int charge, const DoubleVectorVector& intensity, const DoubleVectorVector& area, int numQuanStates );
	void printHTMLIntensity ( std::ostream& os, const std::string& styleID ) const;
	void printHTMLArea ( std::ostream& os, const std::string& styleID ) const;
	void printDelimitedIntensity ( std::ostream& os ) const;
	void printDelimitedArea ( std::ostream& os ) const;
};

class ParentData : public PeakFitData {
	DoubleVector mOverZ;
	DoubleVector intensity;
	DoubleVector fwhm;
	DoubleVector resolution;
	DoubleVector area;
	DoubleVector theoreticalPercentMax;
	DoubleVector snr;
	std::string fString;
	int ch;
	void printHTMLMOverZ ( std::ostream& os, DoubleVectorVectorSizeType peakNumber, const std::string& styleID ) const;
public:
	ParentData ( const XYData& xyData, double monoMass, int charge, ElementalFormula* ef, bool efFlag, double inputResolution, int numPeaks );
	void printHTML ( std::ostream& os ) const;

	static int getColspan ();
	static void printHTMLHeader ( std::ostream& os, const std::string& styleID );
	static void printDelimitedHeader ( std::ostream& os );
	void printHTMLLine ( std::ostream& os, DoubleVectorVectorSizeType peakNumber, const std::string& styleID ) const;
	void printDelimitedLine ( std::ostream& os, DoubleVectorVectorSizeType peakNumber ) const;
	static void printHTMLBlankLine ( std::ostream& os, const std::string& styleID );
	static void printDelimitedBlankLine ( std::ostream& os );
	bool outputQuanResults ( std::ostream& os, const std::string& searchName, int numRepeats, bool area ) const { return false; }
	DoubleVector getAreaRatios () const { return DoubleVector (); }
	DoubleVector getIntensityRatios () const { return DoubleVector (); }
};

#endif /* ! __lu_parent_h */
