/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_filter_srch.h                                              *
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

#ifndef __lu_filter_srch_h
#define __lu_filter_srch_h

#include <lgen_define.h>
#include <lu_filter_par.h>
#include <lu_program.h>

class MSFilterSearch : public MSProgram {
	const MSFilterParameters& filterParams;

	std::string inPeakListFPath;
	std::string outPeakListFPath;
	std::string outPeakListURL;
	std::string archiveName;
	StringVector peakListFractionNames;
	StringVector peakListCentroidFiles;

	void processPeakListFile ();
	void filterPeakLists ( bool keepFlag, const std::string& aSuffix = "" );
	bool checkMPlusHRange ( MSMSDataPoint& mmdp );
	bool checkCharge ( MSMSDataPoint& mmdp );
	bool checkNeutralLoss ( MSMSDataPoint& mmdp );
	bool checkNeutralLoss ( MSMSDataPoint* mmdp, double testMass );
	bool checkFragmentMZ ( MSMSDataPoint& mmdp );
	bool checkFragmentMZs ( MSMSDataPoint& mmdp );
public:
	MSFilterSearch ( const MSFilterParameters& params );
	~MSFilterSearch ();
	void printBodyHTML ( std::ostream& os );
};

#endif /* ! __lu_filter_srch_h */
