/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_filter_par.h                                               *
*                                                                             *
*  Created    : September 3rd 2012                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2012-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_filter_par_h
#define __lu_filter_par_h

#include <lu_mass.h>
#include <lu_prog_par.h>
#include <lu_tol.h>

class MSMSPeakFilterOptions;
class MassInfo;

class MSFilterParameters : public MSProgramParameters {
	AAInitInfo aaInitInfo;
	double lowMPlusH;
	double highMPlusH;
	bool fullMPlusHRange;
	StringVector chargeFilter;
	bool allCharges;
	std::string lossFormula;
	double fragmentMZ;
	DoubleVector fragmentMZs;
	int minMatches;
	MSMSPeakFilterOptions* msmsPeakFilterOptions;
	ToleranceInfo parentMassTolerance;
	ToleranceInfo productMassTolerance;
	double systematicError;
	MassInfo* massInfo;
	std::string instrumentName;
	std::string keepOrRemove;

	std::string archiveName;
	std::string inPeakListFpath;
	std::string outPeakListFpath;
	std::string outPeakListURL;

	std::string initKeepOrRemove ( const std::string& s );
	void initPeakList ( const ParameterList* params );
public:
	MSFilterParameters ( const ParameterList* params );
	~MSFilterParameters ();

	double getLowMPlusH () const { return lowMPlusH; }
	double getHighMPlusH () const { return highMPlusH; }
	bool getFullMPlusHRange () const { return fullMPlusHRange; }
	StringVector getChargeFilter () const { return chargeFilter; }
	bool getAllCharges () const { return allCharges; }
	std::string getLossFormula () const { return lossFormula; }
	double getFragmentMZ () const { return fragmentMZ; }
	DoubleVector getFragmentMZs () const { return fragmentMZs; }
	int getMinMatches () const { return minMatches; }
	MSMSPeakFilterOptions* getMSMSPeakFilterOptions () const { return msmsPeakFilterOptions; }
	Tolerance* getParentMassTolerance () const { return parentMassTolerance.getTolerance (); }
	Tolerance* getProductMassTolerance () const { return productMassTolerance.getTolerance (); }
	double getSystematicError () const { return systematicError; }
	bool getMonoisotopicFlag () const { return massInfo->getMonoisotopicFlag (); }
	bool getMonoParentAverageFragments () const { return massInfo->getMonoParentAverageFragments (); }
	bool getAverageParentMonoFragments () const { return massInfo->getAverageParentMonoFragments (); }
	std::string getKeepOrRemove () const { return keepOrRemove; }

	std::string getInPeakListFpath () const { return inPeakListFpath; }
	std::string getOutPeakListFpath () const { return outPeakListFpath; }
	std::string getOutPeakListURL () const { return outPeakListURL; }
	std::string getArchiveName () const { return archiveName; }
};

#endif /* ! __lu_filter_par_h */
