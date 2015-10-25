/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_pk_filter.h                                                *
*                                                                             *
*  Created    : June 26th 2003                                                *
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

#ifndef __lu_pk_filter_h
#define __lu_pk_filter_h

#include <string>
#include <iostream>

class ParameterList;

class MSPeakFilterOptions {
	static bool peakExclusionDefault;
	static unsigned int maxPeaksDefault;
	static unsigned int minPeaksDefault;
	static double minIntensityDefault;
	static bool massExclusionDefault;
	static double minMassDefault;
	static double maxMassDefault;
	static bool matrixExclusionDefault;
	static double maxMatrixMassDefault;
	bool peakExclusion;
	unsigned int maxPeaks;
	unsigned int minPeaks;
	double minIntensity;
	bool massExclusion;
	double minMass;
	double maxMass;
	bool matrixExclusion;
	double maxMatrixMass;
public:
	MSPeakFilterOptions ( const ParameterList* params );
	bool getPeakExclusion () const { return peakExclusion; }
	unsigned int getMaxPeaks () const { return maxPeaks; }
	unsigned int getMinPeaks () const { return minPeaks; }
	double getMinIntensity () const { return minIntensity; }
	bool getMassExclusion () const { return massExclusion; }
	double getMinMass () const { return minMass; }
	double getMaxMass () const { return maxMass; }
	bool getMatrixExclusion () const { return matrixExclusion; }
	double getMaxMatrixMass () const { return maxMatrixMass; }
	static void copyToCGI ( std::ostream& os, const ParameterList* params );
	static void copyToHiddenFormEntry ( std::ostream& os, const ParameterList* params );
	static void setDefault ( const MSPeakFilterOptions& ms );
};

class MSMSPeakFilterOptions {
	static std::string peakFilterDefault;
	static bool highResETDDeisotopeDefault;
	static bool ftPeakExclusionDefault;
	static bool ECDorETDSideChainExclusionDefault;
	static bool peakExclusionDefault;
	static unsigned int maxPeaksDefault;
	static unsigned int minPeaksDefault;
	static double minIntensityDefault;
	static bool massExclusionDefault;
	static double minMassDefault;
	static double precursorExclusionDefault;
	static double minPrecursorMassDefault;
	static bool matrixExclusionDefault;
	static double maxMatrixMassDefault;
	static bool joinPeaksDefault;
	static bool deisotopeDefault;
	static bool deisotopeHiResDefault;

	std::string peakFilter;
	bool highResETDDeisotope;
	bool ftPeakExclusion;
	bool ECDorETDSideChainExclusion;
	bool peakExclusion;
	unsigned int maxPeaks;
	unsigned int minPeaks;
	double minIntensity;
	bool massExclusion;
	double minMass;
	double precursorExclusion;
	double minPrecursorMass;
	bool matrixExclusion;
	double maxMatrixMass;
	bool joinPeaks;
	bool deisotope;
	bool deisotopeHiRes;
public:
	MSMSPeakFilterOptions ( const ParameterList* params );
	std::string getPeakFilter () const { return peakFilter; }
	bool getRawSpectrum () const { return ( peakFilter == "Unprocessed MSMS" ); }
	bool getPer100Da () const { return peakFilter == "Max MSMS Pks / 100 Da"; }
	bool getHighResETDDeisotope () const { return highResETDDeisotope; }
	bool getFTPeakExclusion () const { return ftPeakExclusion; }
	bool getECDorETDSideChainExclusion () const { return ECDorETDSideChainExclusion; }
	bool getPeakExclusion () const { return peakExclusion; }
	unsigned int getMaxPeaks () const { return maxPeaks; }
	void setPeakFilter ( const std::string& pf ) { peakFilter = pf; }
	void setMaxPeaks ( unsigned int mp ) { maxPeaks = mp; }
	unsigned int getMinPeaks () const { return minPeaks; }
	double getMinIntensity () const { return minIntensity; }
	bool getMassExclusion () const { return massExclusion; }
	double getMinMass () const { return minMass; }
	double getPrecursorExclusion () const { return precursorExclusion; }
	double getMinPrecursorMass () const { return minPrecursorMass; }
	bool getMatrixExclusion () const { return matrixExclusion; }
	double getMaxMatrixMass () const { return maxMatrixMass; }
	bool getJoinPeaks () const { return joinPeaks; }
	bool getDeisotope () const { return deisotope; }
	bool getDeisotopeHiRes () const { return deisotopeHiRes; }

	static void copyToCGI ( std::ostream& os, const ParameterList* params );
	static void copyToHiddenFormEntry ( std::ostream& os, const ParameterList* params );
	static std::string getCommandLineNVPair ( const ParameterList* params );
	static void setDefault ( const MSMSPeakFilterOptions& msms );
	static int getMaxPeaksDefault () { return maxPeaksDefault; }
};

#endif /* ! __lu_pk_filter_h */
