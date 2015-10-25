/******************************************************************************
*                                                                             *
*  Library    : libraw                                                        *
*                                                                             *
*  Filename   : lr_main.h                                                     *
*                                                                             *
*  Created    : April 23rd 2007                                               *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2007-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lr_main_h
#define __lr_main_h

#include <string>
#include <nr.h>

class ParameterList;
class SpecID;
class RawInstance;

class RawFile {
	RawInstance* rawInstance;
public:
	RawFile ( const ParameterList* params, int fraction, double rtIntervalStart, double rtIntervalEnd );
	~RawFile ();
	void getXYData ( std::vector <XYData>& vXYData, const SpecID& specID, double mOverZ, double startMass = 0.0, double endMass = 0.0 );
	void getXYData ( std::vector <XYData>& vXYData, PairStringVectorString& sTimes, bool msFullScan, const SpecID& specID, double mOverZ, double startMass = 0.0, double endMass = 0.0 );
	void getQuantitationXYData ( std::vector <XYData>& vXYData, PairStringVectorString& sTimes, bool ITRAQ, const SpecID& specID, double mOverZ, const std::string& version, double startMass = 0.0, double endMass = 0.0 );
	std::string getRawType () const;
};

#endif /* ! __lr_main_h */
