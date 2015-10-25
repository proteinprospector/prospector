/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_bdg_par.h                                                  *
*                                                                             *
*  Created    : June 14th 2001                                                *
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

#ifndef __lu_bdg_par_h
#define __lu_bdg_par_h

#include <lu_get_link.h>
#include <lu_dig_par.h>
#include <lu_charge.h>

class MSBridgeParameters : public MSDigestParameters {

	MSPeakFilterOptions msPeakFilterOptions;
	StringVector knownSequences;
	MSDataSetInfo dataSetInfo;
	PeakContainerInfo peakContainerInfo;
	SpecID specID;
	LinkInfo* linkInfo;
	StringVector initKnownSequences ( const ParameterList* params );
public:
	MSBridgeParameters ( const ParameterList* params );

	MSPeakFilterOptions getMSPeakFilterOptions () const { return msPeakFilterOptions; }
	MSDataSetInfo* getDataSetInfo () { return &dataSetInfo; }
	const PeakContainerInfo& getPeakContainerInfo () const { return peakContainerInfo; }

	StringVector getAccessionNumbers () { return getSingleEntryParameters ().getAccessionNums (); }
	LinkInfo* getLinkInfo () const { return linkInfo; }

	StringVector getKnownSequences () const { return knownSequences; }
	void printHTML ( std::ostream& os ) const;
};

class MSBridgeLink {
	static StringVector database;
	static int num;
	int ind;
public:
	MSBridgeLink ( const std::string& db );
	static void write ( std::ostream& os, int indexNumber, int dnaReadingFrame, int openReadingFrame, const MultipleModification2& mm2, const std::string& str, int index );
	void putHidden ( std::ostream& os ) const;
};

#endif /* ! __lu_bdg_par_h */
