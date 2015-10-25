/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_histogram.h                                                *
*                                                                             *
*  Created    : February 25th 2004                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2004-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_histogram_h
#define __lu_histogram_h

#include <nr.h>
#include <lgen_define.h>
#include <limits>

class Gaussian;
class MixtureGaussian;

class Histogram {
protected:
	mutable DoubleVector val;
	mutable XYData xyData;
	mutable Gaussian* negDistribution;
	mutable bool valid;
	virtual void compute () const;
	void makeHistogram () const;
	virtual DoubleVector getGuess () const { return DoubleVector (); }
public:
	Histogram ();
	virtual ~Histogram ();
	virtual void add ( const double x ) { val.push_back ( x ); }
	virtual void drawGraph ( std::ostream& os ) const;
	int size () const { return val.size (); }
	bool getValid () const { return valid; }
	friend std::ostream& operator<< ( std::ostream& os, const Histogram& h );
	friend std::istream& operator>> ( std::istream& is, Histogram& h );
};
std::ostream& operator<< ( std::ostream& os, const Histogram& h );
std::istream& operator>> ( std::istream& is, Histogram& h );

class ErrorHistogram : public Histogram {
	virtual void compute () const;
public:
	ErrorHistogram ();
	virtual ~ErrorHistogram ();
	void drawGraph ( std::ostream& os ) const;
};

class ExpectationParameters {
	int minUsedPeptides;
	int maxUsedPeptides;
	double tailPercent;
	ExpectationParameters ();
public:
	static ExpectationParameters& instance ();
	int getMinUsedPeptides () const { return minUsedPeptides; }
	int getMaxUsedPeptides () const { return maxUsedPeptides; }
	double getTailPercent () const { return tailPercent; }
};

class SurvivalHistogram : public Histogram {
	mutable XYData xySurv;
	mutable int numSpectra;
	mutable double b;
	mutable double a;
	void compute ( int numSavedPeptides ) const;
	static std::string expectationMethod;
	static int minUsedPeptides;
	static double tailPercent;
	size_t size;
	size_t sizeLimit;
public:
	SurvivalHistogram ( size_t sizeLimit );
	virtual ~SurvivalHistogram ();
	void printHTML ( std::ostream& os, double score, int numSavedPeptides ) const;
	void printXML ( std::ostream& os, int numSavedPeptides ) const;
	bool addValue ( const double x )
	{
		if ( sizeLimit == 0 ) {		// Just add up the number of spectra
			size++;
			return true;
		}
		if ( size < sizeLimit ) {
			val.push_back ( x );
			size++;
			return true;
		}
		return false;
	}
	static void setExpectationMethod ( const std::string& method )
	{
		expectationMethod = method;
	}
	int getSize () const { return size; }
	void init ( int numSavedPeptides ) const;
	double getEValue ( double score ) const;
	bool getEValueFlag () const { return a != 0.0; }
};

#endif /* ! __lu_histogram_h */
