/******************************************************************************
*                                                                             *
*  Library    : libnrec                                                       *
*                                                                             *
*  Filename   : nr.h                                                          *
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
*  Copyright (1997-2010) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __nr_h
#define __nr_h

#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <lgen_define.h>

#ifdef NR_MAIN
#define NR_EXTERN
#else
#define NR_EXTERN extern
#endif

NR_EXTERN int nrerrno;
#define ZBRENT_ROOT_NOT_BRACKETED 1
#define ZBRENT_MAX_ITERATIONS_EXCEEDED 2

struct lNrecMomentLessThanOneDataValue {};
struct lNrecMomentLessThanTwoDataValues {};
struct lNrecMomentZeroVariance {};

struct lNrecGaussjSingularMatrix1 {};
struct lNrecGaussjSingularMatrix2 {};

// This is a list of the numerical recipies currently used.
double bico ( int n, int k );
void covsrt ( double** covar, int ma, int* lista, int mfit );
double factln ( int n );
double factrl ( int n );
void fit ( const double* x, const double* y, int ndata, double* sig, int mwt, double* a, double* b, double* siga, double* sigb, double* chi2, double* q );
void fgauss ( double x, double* a, double* y, double* dyda, int na );
double gammln ( double xx );
double gammq ( double a, double x );
void gaussj ( double** a, int n, double** b, int m );
void gcf ( double* gammcf, double a, double x, double* gln );
void gser ( double* gamser, double a, double x, double* gln );
void moment ( double* data, int n, double* ave );
void moment_weighted ( double* data, double* wt, double s_wt, int n, double* ave );
void moment ( double* data, int n, double* ave, double* sdev, double* svar );
void moment ( double* data, int n, double* ave, double* adev, double* sdev, double* svar, double* skew, double* curt);
double median ( const DoubleVector& sortedData );
double lowerQ ( const DoubleVector& sortedData );
double upperQ ( const DoubleVector& sortedData );
void mrqcof ( double* x, double* y, double* sig, int ndata, double* a, int ma, int* lista, int mfit, double** alpha, double* beta, double* chisq, void (*funcs)(double, double*, double*, double*, int ) );
void mrqmin ( double* x, double* y, double* sig, int ndata, double* a, int ma, int* lista, int mfit, double** covar, double** alpha, double* chisq, void (*funcs)(double, double*, double*, double*, int ), double* alamda );
double zbrent ( double (*func)(double), double x1, double x2, double tol );

void nrerror ( const char* errorText );
double* nrvector ( int nl, int nh );
void free_nrvector ( double* v, int nl, int nh );
int* inrvector ( int nl, int nh );
void free_inrvector ( int* v, int nl, int nh );
double** nrmatrix ( int nrl, int nrh, int ncl, int nch );
void free_nrmatrix ( double** m, int nrl, int nrh, int ncl, int nch );
// End of list.

typedef std::pair <double,double> XYPair;

struct lNrecEmptyXYDataRange {};

class XYData {
	DoubleVector xList;
	DoubleVector yList;
	mutable DoubleVectorSizeType cur;
public:
	XYData () {};
	void add ( double x, double y )
	{
		xList.push_back ( x );
		yList.push_back ( y );
	}
	void resize ( int size ) { xList.resize ( size ); yList.resize ( size ); }
	DoubleVectorSizeType size () const { return xList.size (); }
	bool empty () const { return xList.size () == 0; }
	void first () const { cur = 0; }
	void next () const { cur++; }
	bool isDone () const { return cur < xList.size (); }
	bool isDone ( int num ) const { return cur < num; }
	const DoubleVectorSizeType& index () const { return cur; }
	const double& x () const { return xList [cur]; }
	const double& y () const { return yList [cur]; }
	const double& x ( const DoubleVectorSizeType& ind ) const { return xList [ind]; }
	const double& y ( const DoubleVectorSizeType& ind ) const { return yList [ind]; }
	std::vector <double> getXList () const { return xList; }
	std::vector <double> getYList () const { return yList; }
	double minX () const { return *(std::min_element ( xList.begin (), xList.end () )); }
	double maxX () const { return *(std::max_element ( xList.begin (), xList.end () )); }
	double minY () const { return *(std::min_element ( yList.begin (), yList.end () )); }
	double maxY () const { return *(std::max_element ( yList.begin (), yList.end () )); }
	double xMaxY () const { return xList [std::max_element ( yList.begin (), yList.end () ) - yList.begin ()]; }
	double xInc () const { return xList [1] - xList [0]; }
	XYData getXRange ( double startX, double endX ) const;
	double getY ( double x ) const;
	double getXAtMaxYInTolRange ( double x, double tol ) const;
	double getMaxYInTolRange ( double x, double tol ) const;
	double mean () const;
	double stddev () const;
	void linearRegression ( double* offset, double* gradient ) const;
	void linearRegression ( double* offset, double* gradient, double* siga, double* sigb, double* chi2, double* q ) const;
	double getSumOfYInTolRange ( double x, double tol ) const;
	void offsetCalibration ( double offset, const std::string& units );
};

std::istream& operator>> ( std::istream& input, XYData& xyData );
std::ostream& operator<< ( std::ostream& os, const XYData& xyData );
void getXYDataRange ( std::istream& input, XYData& xyData, double lowX, double highX );

class SplineInterpolator {
	const DoubleVector xa;
	const DoubleVector ya;
	DoubleVector y2;
	int na;
public:
	SplineInterpolator ( const XYData& xyData, double yp1 = 1.0e30, double yp2 = 1.0e30 );
	double getY ( double x ) const;
};
XYData averageScans ( const std::vector <XYData>& vxyd );

namespace nr {
	const double pi = 3.141592653589793;
	const double sqrtPi = sqrt ( pi );
	const double sqrtTwo = sqrt ( 2.0 );
	const double sqrtTwoPi = sqrt ( 2.0 * pi );
}

#endif /* ! __nr_h */
