/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_fit_par.h                                                  *
*                                                                             *
*  Created    : June 16th 2001                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_fit_par_h
#define __lu_fit_par_h

#include <lu_srch_par.h>
#include <lu_mut_mtrx.h>
#include <lu_spec_id.h>
#include <lu_prog.h>
#include <lu_charge.h>
#include <lu_pk_filter.h>
#include <lu_df_info.h>
#include <lu_usermod.h>

class MowseInfo {
	bool mowseOn;
	double mowsePFactor;
public:
	MowseInfo ( const ParameterList* params );
	bool getMowseOn () const { return mowseOn; }
	double getMowsePFactor () const { return mowsePFactor; }
	void printHTML ( std::ostream& os ) const;
	static void copyToCGI ( std::ostream& os, const ParameterList* params );
	static void copyToHiddenFormEntry ( std::ostream& os, const ParameterList* params );
	static void copyToHiddenFormJavascriptEntry ( std::ostream& os, const ParameterList* params );
};

class MSFitParameters : public MSSearchParameters {

	MSPeakFilterOptions msPeakFilterOptions;

	MSDataSetInfo dataSetInfo;
	PeakContainerInfo peakContainerInfo;
	SpecID specID;

	int minMatches;
	MowseInfo mowseInfo;
	std::string sortType;
	MultipleModification multipleModification;
	MultipleModification2 multipleModification2;
	int minParentIonMatches;
	std::string reportHomologousProteins;

	ModificationParameters modificationParameters;
	static int maxReportedHitsLimit;
public:
	MSFitParameters ( const ParameterList* params );

	MSPeakFilterOptions getMSPeakFilterOptions () const { return msPeakFilterOptions; }

	MSDataSetInfo* getDataSetInfo () { return &dataSetInfo; }
	const PeakContainerInfo& getPeakContainerInfo () const { return peakContainerInfo; }
	bool getMonoisotopicFlag () const { return peakContainerInfo.getMonoisotopicFlag (); }

	int getMinMatches () const { return minMatches; }
	MowseInfo getMowseInfo () const { return mowseInfo; }
	bool getMowseFlag () const { return !getPreSearchInfoFromFile () && mowseInfo.getMowseOn (); }
	std::string getSortType () const { return sortType; }
	MultipleModification getMultipleModification () const { return multipleModification; }
	const MultipleModification2& getMultipleModification2 () const { return multipleModification2; }
	int getMinParentIonMatches () const { return minParentIonMatches; }
	std::string getReportHomologousProteins () const { return reportHomologousProteins; }

	ModificationParameters getModificationParameters () const { return modificationParameters; }

	void printHTML ( std::ostream& os ) const;
};

class MSFitLink {
	static int num;
	int ind;
public:
	MSFitLink ();
	static void write ( std::ostream& os, int minMatches, const PeakContainer& peaks, const BoolDeque& peakUsed, const std::string& str );
	void putHidden ( std::ostream& os ) const;
};

#endif /* ! __lu_fit_par_h */
