/******************************************************************************
*                                                                             *
*  Library    : liboracle                                                     *
*                                                                             *
*  Filename   : lo_calibration.h                                              *
*                                                                             *
*  Created    :                                                               *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2002-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifndef __lo_calibration_h
#define __lo_calibration_h

#include <lgen_define.h>
#include <lo_init.h>

class CalibrationCoeffs : public OracleStatement {
	bool okFlag;
	unsigned int peakListID;
	DoubleVector curCoeff;
	DoubleVector defCoeff;
	std::string serialNumber; 
	int eqnID;
	int ionMode;
	int calState;
	int numCoeff;
public:
	CalibrationCoeffs ( double dJobRunItem, double dAcqSourceJobRunItemID, OracleConnection& oc );
	virtual ~CalibrationCoeffs ();
	void writeFile ( const std::string& filename ) const;
	bool getOKFlag () const { return okFlag; };
	bool getMS () const { return ionMode == 1 || ionMode == 2; };
	bool getMSMS () const { return ionMode == 4 || ionMode == 8; };
};


#endif /* ! __lo_calibration_h */
