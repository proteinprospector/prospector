/******************************************************************************
*                                                                             *
*  Library    : libanalyst                                                    *
*                                                                             *
*  Filename   : la_wiff.h                                                     *
*                                                                             *
*  Created    : June 18th 2004                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2004-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __la_wiff_h
#define __la_wiff_h

#include <string>
#include <vector>
#include <nr.h>
#include <lgen_define.h>

class WiffExploreInstance;

class WiffFile {
protected:
	WiffExploreInstance* inst;
	DoubleVector times;
	double timeWindowStart;
	double timeWindowEnd;
	void initCache (  const std::string& file );
public:
	WiffFile ( const std::string& file, double timeWindowStart = 0.0, double timeWindowEnd = 0.0 );
	~WiffFile ();
	void getDataFromRT ( XYData& xyData, double spot, double startMass, double endMass );
	void getDataFromRT ( XYData& xyData, const std::string& spot, double startMass, double endMass );
	void getDataFromSpecList ( XYData& xyData, const std::string& msmsInfo, double startMass, double endMass, bool msData = false );
	DoubleVector getTimes () const { return times; }
};

class WiffCentroidFile : public WiffFile {
	bool smoothFlag;
	bool useIntensityNotArea;
	bool mergeDistInPPM;
	double heightPercentage;
	double mergeDistance;
public:
	WiffCentroidFile ( const std::string& file, bool smoothFlag, bool useIntensityNotArea, bool mergeDistInPPM, double heightPercentage, double mergeDistance );
	void getData ( XYData& xyData, int fraction, double spot, int msType, double startMass, double endMass ) const;
	void getData ( XYData& xyData, int fraction, const std::string& spot, int msType, double startMass, double endMass ) const;
};

#endif /* ! __la_wiff_h */
