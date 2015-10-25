/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_dbst_srch.h                                                *
*                                                                             *
*  Created    : September 5th 2001                                            *
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

#ifndef __lu_dbst_srch_h
#define __lu_dbst_srch_h

#include <lu_hit.h>
#include <lu_srch_par.h>
#include <lu_db_srch.h>

class DBStatParameters : public MSSearchParameters {
	bool showAAStatistics;
	double minHistogramMass;
	double maxHistogramMass;
	double bandwidth;
public:
	DBStatParameters ( const ParameterList* params );
	bool getShowAAStatistics () const { return showAAStatistics; }
	double getMinHistogramMass () const { return minHistogramMass; }
	double getMaxHistogramMass () const { return maxHistogramMass; }
	double getBandwidth () const { return bandwidth; }
};

class AminoAcidStats;

class DBStatSearch : public DBSearch {
	AminoAcidStats* aaStats;
	double bandwidth;
	int maxMissedCleavages;
	bool outputFlag;

	int longestProteinDBIndex;
	int longestProteinIndex;
	double longestProteinMW;
	int largestNumEnzymeFragmentsDBIndex;
	int largestNumEnzymeFragmentsIndex;
	int longestProteinOrfNumber;
	int largestNumEnzymeFragmentsOrfNumber;

	int maxNumAA;
	int largestNumEnzymeFragments;
	double totalNumFragments;
	double averageProteinMW;
	double totalAA;

	void printParamsBodyHTML ( std::ostream& os ) const { params.printHTML ( os ); }
	void calculateDigestHistogram ( const char* frame, const IntVector& cleavageIndex, double minHistMass, double maxHistMass );
	static int getNumFragments ( int numEnzymeFragments, int maxMissedCleavages );
public:
	DBStatSearch ( const DBStatParameters& params );
	~DBStatSearch ();
	void printHTMLHits ( std::ostream& os );
	void printXMLHits ( std::ostream& os ) const;
};

#endif /* ! __lu_dbst_srch_h */
