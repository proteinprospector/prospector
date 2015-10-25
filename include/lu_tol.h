/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_tol.h                                                      *
*                                                                             *
*  Created    : June 7th 2001                                                 *
*                                                                             *
*  Purpose    : Functions concerned with tolerance and tolerance units.       *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2011) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_tol_h
#define __lu_tol_h

#include <string>
#include <iostream>

class ParameterList;

class Tolerance {
protected:
	double value;
	std::string valueStr;
	double offset;
	double lowVal;
	double highVal;
	std::string unitsString;
public:
	Tolerance ( double value, const std::string& valueStr, double offset, const std::string& unitsString );
	virtual ~Tolerance ();
	double getValue () const { return value; }
	std::string getValueStr () const { return valueStr; }
	double getOffset () const { return offset; }
	void setValue ( double v )
	{
		value = v;
		lowVal = -value + offset;
		highVal = value + offset;
	}
	std::string getUnitsString () const { return unitsString; }
	void print ( std::ostream& os ) const;
	virtual double getTolerance ( double mass, int charge = 1 ) const = 0;
	virtual double getLowerTolerance ( double mass, int charge = 1 ) const = 0;
	virtual double getUpperTolerance ( double mass, int charge = 1 ) const = 0;
	virtual double getCorrection ( double mOverZ, double gradient, double offset ) const = 0;
	virtual double getCorrection ( double mass, double systematicError ) const = 0;
	virtual double getError ( double measuredMass, double actualMass, int charge ) const = 0;
	virtual double getActualMass ( double measuredMass, double error, int charge ) const = 0;

	void printHTML ( std::ostream& os, const std::string& label ) const;
	void putCGI ( std::ostream& os, const std::string& prefix ) const;
	void putHiddenFormEntry ( std::ostream& os, const std::string& prefix ) const;
	void putHiddenFormJavascriptEntry ( std::ostream& os, const std::string& prefix ) const;
	std::string getCommandLineNVPair ( const std::string& prefix ) const;
};

class ToleranceInfo {
	Tolerance* tolerance;
	std::string prefix;
	static std::string TOLERANCE_VALUE;
	static std::string TOLERANCE_OFFSET;
	static std::string TOLERANCE_UNITS;
	Tolerance* initTolerance ( const std::string& prefix, const ParameterList* params );
	Tolerance* getNewTolerance ( double value, const std::string& valueStr, double offset, const std::string& units );
public:
	ToleranceInfo ( const std::string& prefix, const ParameterList* params );
	ToleranceInfo ( const std::string& prefix, double offset, const std::string& units );
	~ToleranceInfo ();
	ToleranceInfo ( const ToleranceInfo& rhs );
	ToleranceInfo ( const ToleranceInfo* rhs );
	ToleranceInfo& operator= ( ToleranceInfo& rhs );
	ToleranceInfo& operator= ( ToleranceInfo* rhs );
	Tolerance* getTolerance () const { return tolerance; }
	void setValue ( double value ) { tolerance->setValue ( value ); }
	void printHTML ( std::ostream& os ) const;
	void putCGI ( std::ostream& os ) const;
	friend class Tolerance;
};

#endif /* ! __lu_tol_h */
