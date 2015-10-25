/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_tol.cpp                                                    *
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
#include <lu_tol.h>
#include <lu_html_form.h>
#include <lu_cgi_val.h>
#include <lu_param_list.h>
using std::string;
using std::ostream;

class DaltonTolerance : public Tolerance {
public:
	DaltonTolerance ( double value, const string& valueStr, double offset ) :
		Tolerance ( value, valueStr, offset, "Da" ) {}
	~DaltonTolerance ();
	double getTolerance ( double mass, int charge = 1 ) const;
	double getLowerTolerance ( double mass, int charge = 1 ) const;
	double getUpperTolerance ( double mass, int charge = 1 ) const;
	double getCorrection ( double mOverZ, double gradient, double offset ) const;
	double getCorrection ( double mass, double systematicError ) const;
	double getError ( double measuredMass, double actualMass, int charge ) const;
	double getActualMass ( double measuredMass, double error, int charge ) const;
};
class PPMTolerance : public Tolerance {
public:
	PPMTolerance ( double value, const string& valueStr, double offset ) :
		Tolerance ( value, valueStr, offset, "ppm" ) {}
	~PPMTolerance ();
	double getTolerance ( double mass, int charge = 1 ) const;
	double getLowerTolerance ( double mass, int charge = 1 ) const;
	double getUpperTolerance ( double mass, int charge = 1 ) const;
	double getCorrection ( double mOverZ, double gradient, double offset ) const;
	double getCorrection ( double mass, double systematicError ) const;
	double getError ( double measuredMass, double actualMass, int charge ) const;
	double getActualMass ( double measuredMass, double error, int charge ) const;
};
class PercentTolerance : public Tolerance {
public:
	PercentTolerance ( double value, const string& valueStr, double offset ) :
		Tolerance ( value, valueStr, offset, "%" ) {}
	~PercentTolerance ();
	double getTolerance ( double mass, int charge = 1 ) const;
	double getLowerTolerance ( double mass, int charge = 1 ) const;
	double getUpperTolerance ( double mass, int charge = 1 ) const;
	double getCorrection ( double mOverZ, double gradient, double offset ) const;
	double getCorrection ( double mass, double systematicError ) const;
	double getError ( double measuredMass, double actualMass, int charge ) const;
	double getActualMass ( double measuredMass, double error, int charge ) const;
};
class MMUTolerance : public Tolerance {
public:
	MMUTolerance ( double value, const string& valueStr, double offset ) :
		Tolerance ( value, valueStr, offset, "mmu" ) {}
	~MMUTolerance ();
	double getTolerance ( double mass, int charge = 1 ) const;
	double getLowerTolerance ( double mass, int charge = 1 ) const;
	double getUpperTolerance ( double mass, int charge = 1 ) const;
	double getCorrection ( double mOverZ, double gradient, double offset ) const;
	double getCorrection ( double mass, double systematicError ) const;
	double getError ( double measuredMass, double actualMass, int charge ) const;
	double getActualMass ( double measuredMass, double error, int charge ) const;
};
Tolerance::Tolerance ( double value, const string& valueStr, double offset, const string& unitsString ) :
	value ( value ),
	valueStr ( valueStr ),
	offset ( offset ),
	lowVal ( -value + offset ),
	highVal ( value + offset ),
	unitsString ( unitsString )
{
}
Tolerance::~Tolerance () {}
void Tolerance::print ( ostream& os ) const
{
	os << " (+/- ";
	genPrint ( os, value, 2 );
	os << " ";
	os << unitsString;
	os << ")";
}
void Tolerance::printHTML ( ostream& os, const string& label ) const
{
	os << label << ": <b>" << value << " " << unitsString << "</b><br />\n";
}
void Tolerance::putCGI ( ostream& os, const string& prefix ) const
{
	printCGI ( os, prefix + ToleranceInfo::TOLERANCE_VALUE, valueStr );
	printCGIString ( os, prefix + ToleranceInfo::TOLERANCE_UNITS, unitsString );
}
void Tolerance::putHiddenFormEntry ( ostream& os, const string& prefix ) const
{
	printHTMLFORMHidden ( os, prefix + ToleranceInfo::TOLERANCE_VALUE, valueStr );
	printHTMLFORMHidden ( os, prefix + ToleranceInfo::TOLERANCE_UNITS, unitsString );
}
void Tolerance::putHiddenFormJavascriptEntry ( ostream& os, const string& prefix ) const
{
	printHTMLFORMJavascriptHidden ( os, prefix + ToleranceInfo::TOLERANCE_VALUE, valueStr );
	printHTMLFORMJavascriptHidden ( os, prefix + ToleranceInfo::TOLERANCE_UNITS, unitsString );
}
string Tolerance::getCommandLineNVPair ( const string& prefix ) const
{
	string s;
	s += ::getCommandLineNVPair ( prefix + ToleranceInfo::TOLERANCE_VALUE, valueStr );
	s += " ";
	s += ::getCommandLineNVPair ( prefix + ToleranceInfo::TOLERANCE_UNITS, unitsString );
	s += " ";
	return s;
}
DaltonTolerance::~DaltonTolerance () {}
double DaltonTolerance::getTolerance ( double mass, int charge ) const
{
	return ( value * charge );
}
double DaltonTolerance::getLowerTolerance ( double mass, int charge ) const
{
	return ( lowVal * charge );
}
double DaltonTolerance::getUpperTolerance ( double mass, int charge ) const
{
	return highVal * charge;
}
double DaltonTolerance::getCorrection ( double mass, double systematicError ) const
{
	return systematicError;
}
double DaltonTolerance::getCorrection ( double mOverZ, double gradient, double offset ) const
{
	return offset + ( gradient * mOverZ );
}
double DaltonTolerance::getError ( double measuredMass, double actualMass, int charge ) const
{
	return ( measuredMass - actualMass ) / charge;
}
double DaltonTolerance::getActualMass ( double measuredMass, double error, int charge ) const
{
	return measuredMass - ( error * charge );
}
PPMTolerance::~PPMTolerance () {}
double PPMTolerance::getTolerance ( double mass, int charge ) const
{
	return ( mass * value * charge / 1000000.0 );
}
double PPMTolerance::getLowerTolerance ( double mass, int charge ) const
{
	return ( mass * lowVal * charge / 1000000.0 );
}
double PPMTolerance::getUpperTolerance ( double mass, int charge ) const
{
	return ( mass * highVal * charge / 1000000.0 );
}
double PPMTolerance::getCorrection ( double mass, double systematicError ) const
{
	return ( mass * systematicError / 1000000.0 );
}
double PPMTolerance::getCorrection ( double mOverZ, double gradient, double offset ) const
{
	return ( offset + ( gradient * mOverZ ) ) * mOverZ / 1000000.0;
}
double PPMTolerance::getError ( double measuredMass, double actualMass, int charge ) const
{
	return ( ( measuredMass - actualMass ) * 1000000.0 / actualMass );
}
double PPMTolerance::getActualMass ( double measuredMass, double error, int charge ) const
{
	return ( ( measuredMass * 1000000.0 ) / ( error + 1000000.0 ) );
}
PercentTolerance::~PercentTolerance () {}
double PercentTolerance::getTolerance ( double mass, int charge ) const
{
	return ( mass * value * charge / 100.0 );
}
double PercentTolerance::getLowerTolerance ( double mass, int charge ) const
{
	return ( mass * lowVal * charge / 100.0 );
}
double PercentTolerance::getUpperTolerance ( double mass, int charge ) const
{
	return ( mass * highVal * charge / 100.0 );
}
double PercentTolerance::getCorrection ( double mass, double systematicError ) const
{
	return ( mass * systematicError / 100.0 );
}
double PercentTolerance::getCorrection ( double mOverZ, double gradient, double offset ) const
{
	return ( offset + ( gradient * mOverZ ) ) * mOverZ / 100.0;
}
double PercentTolerance::getError ( double measuredMass, double actualMass, int charge ) const
{
	return ( ( measuredMass - actualMass ) * 100.0 / actualMass );
}
double PercentTolerance::getActualMass ( double measuredMass, double error, int charge ) const
{
	return ( ( measuredMass * 100.0 ) / ( error + 100.0 ) );
}
MMUTolerance::~MMUTolerance () {}
double MMUTolerance::getTolerance ( double mass, int charge ) const
{
	return ( value * charge / 1000.0 );
}
double MMUTolerance::getLowerTolerance ( double mass, int charge ) const
{
	return ( lowVal * charge / 1000.0 );
}
double MMUTolerance::getUpperTolerance ( double mass, int charge ) const
{
	return ( highVal * charge / 1000.0 );
}
double MMUTolerance::getCorrection ( double mass, double systematicError ) const
{
	return ( systematicError / 1000.0 );
}
double MMUTolerance::getCorrection ( double mOverZ, double gradient, double offset ) const
{
	return ( offset + ( gradient * mOverZ ) ) / 1000.0;
}
double MMUTolerance::getError ( double measuredMass, double actualMass, int charge ) const
{
	return ( ( measuredMass - actualMass ) * 1000.0 / charge );
}
double MMUTolerance::getActualMass ( double measuredMass, double error, int charge ) const
{
	return measuredMass - ( error * 1000.0 * charge );
}

string ToleranceInfo::TOLERANCE_VALUE = "_tolerance";
string ToleranceInfo::TOLERANCE_OFFSET = "_offset";
string ToleranceInfo::TOLERANCE_UNITS = "_tolerance_units";

ToleranceInfo::ToleranceInfo ( const string& prefix, const ParameterList* params ) :
	tolerance ( initTolerance ( prefix, params ) ),
	prefix ( prefix )
{
}
ToleranceInfo::ToleranceInfo ( const string& prefix, double offset, const string& units ) :
	tolerance ( getNewTolerance ( 0.0, "0.0", offset, units ) ),
	prefix ( prefix )
{
}
ToleranceInfo::~ToleranceInfo ()
{
	delete tolerance;
}
ToleranceInfo::ToleranceInfo ( const ToleranceInfo& rhs ) :
	tolerance ( getNewTolerance ( rhs.tolerance->getValue (), rhs.tolerance->getValueStr (), rhs.tolerance->getOffset (), rhs.tolerance->getUnitsString () ) ),
	prefix ( rhs.prefix )
{
}
ToleranceInfo::ToleranceInfo ( const ToleranceInfo* rhs ) :
	tolerance ( getNewTolerance ( rhs->tolerance->getValue (), rhs->tolerance->getValueStr (), rhs->tolerance->getOffset (), rhs->tolerance->getUnitsString () ) ),
	prefix ( rhs->prefix )
{
}
ToleranceInfo& ToleranceInfo::operator= ( ToleranceInfo& rhs )
{
	if ( this != &rhs ) {
		delete tolerance;
		tolerance = getNewTolerance ( rhs.tolerance->getValue (), rhs.tolerance->getValueStr (), rhs.tolerance->getOffset (), rhs.tolerance->getUnitsString () );
		prefix = rhs.prefix;
	}
	return *this;
}
ToleranceInfo& ToleranceInfo::operator= ( ToleranceInfo* rhs )
{
	if ( this != rhs ) {
		delete tolerance;
		tolerance = getNewTolerance ( rhs->tolerance->getValue (), rhs->tolerance->getValueStr (), rhs->tolerance->getOffset (), rhs->tolerance->getUnitsString () );
		prefix = rhs->prefix;
	}
	return *this;
}
Tolerance* ToleranceInfo::initTolerance ( const string& prefix, const ParameterList* params )
{
	double value = params->getDoubleValue ( prefix + TOLERANCE_VALUE, 1.0 );
	string valueStr = params->getStringValue ( prefix + TOLERANCE_VALUE, "1.0" );
	double offset = params->getDoubleValue ( prefix + TOLERANCE_OFFSET, 0.0 );
	string units = params->getStringValue ( prefix + TOLERANCE_UNITS, "Da" );
	return getNewTolerance ( value, valueStr, offset, units );
}
Tolerance* ToleranceInfo::getNewTolerance ( double value, const string& valueStr, double offset, const string& units )
{
	if ( units == "Da" )	return new DaltonTolerance ( value, valueStr, offset );
	if ( units == "ppm" )	return new PPMTolerance ( value, valueStr, offset );
	if ( units == "%" )		return new PercentTolerance ( value, valueStr, offset );
	if ( units == "mmu" )	return new MMUTolerance ( value, valueStr, offset );
	return 0;
}
void ToleranceInfo::printHTML ( ostream& os ) const
{
	tolerance->printHTML ( os, prefix );
}
void ToleranceInfo::putCGI ( ostream& os ) const
{
	tolerance->putCGI ( os, prefix );
}
