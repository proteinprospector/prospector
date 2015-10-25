/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_nspec_par.h                                                *
*                                                                             *
*  Created    : February 12th 2002                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2002-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_nspec_par_h
#define __lu_nspec_par_h

#include <lu_dig_par.h>
#include <lu_charge.h>

class MSNonSpecificParameters : public MSDigestParameters {

	MSPeakFilterOptions msPeakFilterOptions;
	MSDataSetInfo dataSetInfo;
	PeakContainerInfo peakContainerInfo;
	int maxHits;
public:
	MSNonSpecificParameters ( const ParameterList* params );

	int getMaxHits () const { return maxHits; }
	MSPeakFilterOptions getMSPeakFilterOptions () const { return msPeakFilterOptions; }
	MSDataSetInfo* getDataSetInfo () { return &dataSetInfo; }
	const PeakContainerInfo& getPeakContainerInfo () const { return peakContainerInfo; }

	void printHTML ( std::ostream& os ) const;
};

class MSNonSpecificLink {
	static StringVector database;
	static int num;
	int ind;
public:
	MSNonSpecificLink ( const std::string& db );
	static void write ( std::ostream& os, int indexNumber, int dnaReadingFrame, int openReadingFrame, const std::string& str, int index );
	void putHidden ( std::ostream& os ) const;
};

#endif /* ! __lu_nspec_par_h */
