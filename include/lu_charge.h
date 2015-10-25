/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_charge.h                                                   *
*                                                                             *
*  Created    : January 28th 1997                                             *
*                                                                             *
*  Purpose    : Functions for dealing with multiply charged data.             *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1997-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_charge_h
#define __lu_charge_h

#include <vector>
#include <iostream>
#include <string>

#include <lu_data.h>
#include <lu_tol.h>
#include <lu_inst.h>

class MSPeakFilterOptions;
class MSMSPeakFilterOptions;
class MSDataSetInfo;

class PeakContainerInfo {
	ToleranceInfo toleranceInfo;
	std::string massType;
	bool monoisotopicFlag;
	double systematicError;
	IntVector averageToMonoConvert;
	DoubleVector parentContaminantMasses;
public:
	PeakContainerInfo ( const ParameterList* params );
	Tolerance* getTolerance () const { return toleranceInfo.getTolerance (); }
	std::string getMassType () const { return massType; }
	bool getMonoisotopicFlag () const { return monoisotopicFlag; }
	double getSystematicError () const { return systematicError; }
	IntVector getAverageToMonoConvert () const { return averageToMonoConvert; }
	DoubleVector getParentContaminantMasses () const { return parentContaminantMasses; }
};

class Peak {
protected:
	double mOverZ;
	double tolerance;
	int charge;
	double intensity;
	mutable double mass;
	double adductMass;
	double averageMass;
	std::string specNumber;
	friend class sortPeaksByMass;
	friend class sortPeaksByMOverZ;
	friend class sortPeaksByDescendingIntensity;
	friend class PeakContainer;
public:
	Peak ( double mOverZ, double tolerance, int charge, double intensity, double adductMass, int averageToMonoConvertFlag = false );
	Peak ( double mOverZ, double tolerance, int charge, double intensity, double mass, double adductMass, double averageMass, const std::string& specNumber );
	double getMOverZ () const { return mOverZ; }
	int getCharge () const { return charge; }
	bool getNonUnitChargeData () const { return charge != 1; }
	double getTolerance () const { return tolerance; }
	double getIntensity () const { return intensity; }
	double getMass () const { return mass; }
	double getMassMinusTol () const { return mass - tolerance; }
	double getMassPlusTol () const { return mass + tolerance; }
	double getAverageMass () const { return averageMass; }
	void setTolerance ( double t ) { tolerance = t; }
	void setMassAverage () { mass = averageMass; }
	void setSpecNumber ( std::string sNum ) { specNumber = sNum; }
	std::string getSpecNumber () const { return specNumber; }
	bool isMatch ( double modelMass ) const { return genAbsDiff ( modelMass, mass ) < tolerance; }
	bool isLowerMatch ( double modelMass ) const { return modelMass > mass - tolerance; }
	bool isUpperMatch ( double modelMass ) const { return modelMass < mass + tolerance; }
	bool isInRange ( double modelMass, double range ) const { return genAbsDiff ( modelMass, mass ) < tolerance + range; }
	bool isAtLeastRatio ( double inten, double ratio ) const { return inten >= intensity * ratio; }
	bool getAverageFlag () const { return averageMass != 0.0; }

	void calibrate ( const Tolerance* t, double gradient, double offset );
	void putCGI ( std::ostream& os ) const;
	void putParentCGI ( std::ostream& os ) const;
	void putHiddenFormEntry ( std::ostream& os ) const;
	void putHiddenFormJavascriptEntry ( std::ostream& os ) const;
	void printMOverZHTML ( std::ostream& os, int precision ) const;
	void printMOverZDelimited ( std::ostream& os, int precision ) const;
	void printHTML ( std::ostream& os, const std::string& label, const PeakPrecision& pp ) const;
	void printXML ( std::ostream& os, const PeakPrecision& pp ) const;
	void printXML ( std::ostream& os, const std::string& label, const PeakPrecision& pp ) const;
	void printDelimited ( std::ostream& os, const PeakPrecision& pp ) const;
	static void printDelimited ( std::ostream& os, double mOverZ, int charge, double intensity, int mPrecision, int iPrecision );
	static void printXML ( std::ostream& os, double mOverZ, int charge, double intensity, int mPrecision, int iPrecision );
};

typedef std::vector <Peak> PeakVector;
typedef PeakVector::size_type PeakVectorSizeType;
typedef PeakVector::iterator PeakVectorIterator;
typedef PeakVector::const_iterator PeakVectorConstIterator;

typedef std::vector <Peak*> PeakPtrVector;
typedef PeakPtrVector::size_type PeakPtrVectorSizeType;

class PeakContainer {
	typedef Peak* T;
	std::vector <T> peaks;
	const Tolerance* tolerance;
	bool monoisotopicFlag;
	bool spectrumRetained;
	static bool multiChargeAssign;

	void initMS ( MSDataPoint* dataPoint, const MSPeakFilterOptions& msPeakFilterOptions, const PeakContainerInfo& peakContainerInfo, const std::string& specNum );
	void makePeaks ( PeakVector& pks );
	static void getPeaks ( PeakVector& pks, const DataFilePeakVector& dataPeaks, const Tolerance* tolerance, bool monoisotopicFlag, const IntVector& average_to_mono_convert_array, double systematic_error = 0.0, const std::string& specNumber = "" );
	static void deleteContaminantPeaks ( PeakVector& pks, const DoubleVector& contaminantPeaks );
	static DataFilePeakVector joinSplitPeaks ( const DataFilePeakVector& pks );
	PeakVector deisotope ( const PeakVector& pks );
	PeakVector deisotopeHighResolution ( const PeakVector& pks );
	PeakVector setSingleCharge ( const PeakVector& pks );
	static void removeAboveMOverZ ( DataFilePeakVector& pks, double maxMOverZ );
	static void removeMOverZRange ( DataFilePeakVector& pks, double minMOverZ, double maxMOverZ );
	static void filterMSMSPeaks ( PeakVector& pks, double minMass, double precursorExclusion, const Peak* parentPeak );
	static void filterPeaks ( PeakVector& pks, double minMass, double maxMass );
	static void filterQuantitationPeaks ( DataFilePeakVector& pks );
	static void filterIntensity ( DataFilePeakVector& pks, double minIntensity );
	static void retainPeaks ( PeakVector& pks, PeakVectorSizeType maxPeaks, PeakVectorSizeType minPeaks );
	static void retainPeaksMSMS ( PeakVector& pks, PeakVectorSizeType maxPeaks, PeakVectorSizeType minPeaks );
	static void retainNPeaksPerDaRangeMSMS ( PeakVector& pks, int n, double daRange, PeakVectorSizeType minPeaks );
	static void removeMatrix ( PeakVector& pks, double maxMass );
	static void removeFTPeaks ( DataFilePeakVector& pks, const Peak* parentPeak, const Tolerance* tol );
	static void removeECDorETDSideChainPeaks ( DataFilePeakVector& pks, const Peak* parentPeak, const Tolerance* parTol, const Tolerance* tol );
	static void removeIsotopeDistribution ( DataFilePeakVector& pks, double mz, int z, double tol );
	static double splitInterval ( double mass )
	{
		if ( mass <= 500.0 ) return 0.2;
		if ( mass <= 650.0 ) return 0.25;
		if ( mass <= 800.0 ) return 0.3;
		return 0.35;
	}
public:
	PeakContainer () {}
	PeakContainer ( MSDataSetInfo* dsi, const MSPeakFilterOptions& msPeakFilterOptions, const PeakContainerInfo& peakContainerInfo );
	PeakContainer ( MSDataPoint* dataPoint, const MSPeakFilterOptions& msPeakFilterOptions, const PeakContainerInfo& peakContainerInfo );
	PeakContainer ( MSMSDataPoint* dataPoint, const MSMSPeakFilterOptions* msmsPeakFilterOptions, const Peak* parentPeak, const Tolerance* parTol, const Tolerance* tolerance, bool monoisotopicFlag, bool averageToMonoConvertFlag );
	~PeakContainer ();
	typedef std::vector <T>::size_type size_type;
	bool getNonUnitChargeData () const;
	int getMaxCharge () const;
	double getMaxIntensity () const;
	double getMinMassMinusTol () const;
	double getMaxMassPlusTol () const;
	size_type size () const { return peaks.size (); }
	bool empty () const { return peaks.empty (); }
	void truncate ( int size );
	void sortPeaks ();
	T& operator [] ( int index ) { return peaks [index]; }
	const T& operator [] ( int index ) const { return peaks [index]; }

	const Tolerance* getTolerance () const { return tolerance; }
	const std::string getToleranceUnits () const { return tolerance->getUnitsString (); }
	bool getSpectrumRetained () const { return spectrumRetained; }
	bool getPrecursorLessThanFragments ( double parentMass ) const;
	bool getMonoisotopicFlag () const { return monoisotopicFlag; }
	void setToleranceValue ( double val );
	void calibrate ( double gradient, double offset );
	void printPeaksHTML ( std::ostream& os, const PeakPrecision& pp ) const;
	void printPeaksXML ( std::ostream& os, const std::string& label, const PeakPrecision& pp ) const;
	void putCGI ( std::ostream& os, const std::string& name, const Peak* parentPeak ) const;
	void putCGI ( std::ostream& os, const std::string& name, const BoolDeque& peakUsed ) const;
	void putCGI ( std::ostream& os, const std::string& name ) const;
	void putHiddenFormEntry ( std::ostream& os, const std::string& name, const BoolDeque& peakUsed ) const;
	void putHiddenFormEntry ( std::ostream& os, const std::string& name ) const;
	void putHiddenFormJavascriptEntry ( std::ostream& os, const std::string& name, const BoolDeque& peakUsed, const std::string& str ) const;

	void setMassAverage () const;
	static void setMultiChargeAssign ( bool flag ) { multiChargeAssign = flag; }
};

typedef PeakContainer::size_type PeakContainerSizeType;

class PeakMatchContext {
	const Tolerance* tolerance;
	PeakPrecision peakPrecision;
	bool nonUnitChargeData;
public:
	PeakMatchContext ( const Tolerance* tolerance, const PeakPrecision& pp, bool nonUnitChargeData );
	std::string getToleranceUnits () const { return tolerance->getUnitsString (); }
	PeakPrecision getPeakPrecision () const { return peakPrecision; }
	int getPrecision () const { return peakPrecision.getMassDecimalPlaces (); }
	int getErrorSigFig () const { return peakPrecision.getErrorSigFig (); }
	double getError ( double measuredMass, double matchedMass, int charge ) const
		{ return tolerance->getError ( measuredMass, matchedMass, charge ); }
	bool getNonUnitChargeData () const { return nonUnitChargeData; }
};

class PeakMatch {
	const Peak* dataPeak;
	double matchedMass;
public:
	PeakMatch ( const Peak* dataPeak, double matchedMass );
	double getMassDiff () const { return dataPeak->getMass () - matchedMass; }
	double getError ( const PeakMatchContext& pmc ) const;
	double getDataMass () const { return dataPeak->getMass (); }
	double getDataMOverZ () const { return dataPeak->getMOverZ (); }
	double getMatchedMass () const { return matchedMass; }
	int getCharge () const { return dataPeak->getCharge (); }
	std::string getSpecNumber () const { return dataPeak->getSpecNumber (); }
	void printHeaderHTML ( std::ostream& os, const PeakMatchContext& peakMatchContext );
	static void printHeaderDelimited ( std::ostream& os, const std::string& errUnits, bool multi );
	void printHTML ( std::ostream& os, const PeakMatchContext& peakMatchContext, bool printLine = true ) const;
	void printDelimited ( std::ostream& os, const PeakMatchContext& pmc, bool printLine ) const;
	void printXML ( std::ostream& os, const PeakMatchContext& peakMatchContext ) const;
	void printTagHitHTML ( std::ostream& os, const PeakMatchContext& pmc, bool zeroErr ) const;
};

#endif /* ! __lu_charge_h */
