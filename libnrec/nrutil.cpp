/******************************************************************************
*                                                                             *
*  Library    : libnrec                                                       *
*                                                                             *
*  Filename   : nrutil.cpp                                                    *
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
*  Copyright (1998-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <cstdio>
#include <numeric>
#define NR_MAIN
#include <nr.h>
using std::istream;
using std::ostream;
using std::string;
using std::lower_bound;
using std::upper_bound;
using std::max_element;
using std::accumulate;
using std::endl;

void nrerror ( const char* errorText )
{
	fprintf ( stderr, "Error: %s\n", errorText );
	exit ( 1 );
}
double* nrvector ( int nl, int nh )
{
	double *v;

	v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl;
}
void free_nrvector ( double* v, int nl, int nh )
{
	free((char*) (v+nl));
}
int* inrvector ( int nl, int nh )
{
	int *v;

	v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl;
}
void free_inrvector ( int* v, int nl, int nh )
{
	free((char*) (v+nl));
}
double** nrmatrix ( int nrl, int nrh, int ncl, int nch )
{
	int i;
	double **m;

	m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
		if (!m[i]) nrerror("allocation failure 2 in matrix()");
		m[i] -= ncl;
	}
	return m;
}
void free_nrmatrix ( double** m, int nrl, int nrh, int ncl, int nch )
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}
istream& operator>> ( istream& is, XYData& xyData )
{
	xyData.resize ( 0 );
	double x;
	double y;

	while ( is >> x ) {
		is >> y;
		xyData.add ( x, y );
	}
	if ( xyData.size () == 0 && is.fail () ) throw std::ios_base::failure ( "Illegal Data Format." );
	return is;
}
ostream& operator<< ( ostream& os, const XYData& xyData )
{
	for ( xyData.first () ; xyData.isDone () ; xyData.next () ) {
		os << xyData.x () << " " << xyData.y () << endl;
	}
	return os;
}
void getXYDataRange ( istream& is, XYData& xyData, double lowX, double highX )
{
	xyData.resize ( 0 );
	double x;
	double y;

	while ( is >> x ) {
		is >> y;
		if ( x > highX ) return;
		if ( x >= lowX ) xyData.add ( x, y );
	}
	if ( xyData.size () == 0 && is.fail () ) throw std::ios_base::failure ( "Illegal Data Format." );
}
XYData XYData::getXRange ( double startX, double endX ) const
{
	XYData xyData;
	int start = lower_bound ( xList.begin (), xList.end (), startX ) - xList.begin ();
	int end = upper_bound ( xList.begin (), xList.end (), endX ) - xList.begin ();
	for ( int i = start ; i != end ; i++ ) {
		xyData.add ( xList [i], yList [i] );
	}
	return xyData;
}
double XYData::getY ( double x ) const
{
	int start = lower_bound ( xList.begin (), xList.end (), x ) - xList.begin ();
	int end = start + 1;
	double gradient = ( yList [end] - yList [start] ) / ( xList [end] - xList [start] );
	double offset = yList [end] - ( gradient * xList [end] );
	return ( ( gradient * x ) + offset );
}
double XYData::getXAtMaxYInTolRange ( double x, double tol ) const
{
	int start = lower_bound ( xList.begin (), xList.end (), x - tol ) - xList.begin ();
	int end = upper_bound ( xList.begin (), xList.end (), x + tol ) - xList.begin ();
	if ( start != end )
		return xList [max_element ( yList.begin () + start, yList.begin () + end ) - yList.begin ()];
	else
		throw lNrecEmptyXYDataRange ();
}
double XYData::getMaxYInTolRange ( double x, double tol ) const
{
	int start = lower_bound ( xList.begin (), xList.end (), x - tol ) - xList.begin ();
	int end = upper_bound ( xList.begin (), xList.end (), x + tol ) - xList.begin ();
	if ( start != end )
		return yList [max_element ( yList.begin () + start, yList.begin () + end ) - yList.begin ()];
	else
		throw lNrecEmptyXYDataRange ();
}
double XYData::mean () const
{
	double s = 0.0;
	int n = yList.size ();
	for ( int i = 0 ; i < n ; i++ ) {
		s += yList [i];
	}
	return ( s / n );
}
double XYData::stddev () const
{
	double sd = 0.0;
	double m = mean ();
	int n = yList.size ();
	for ( int i = 0 ; i < n ; i++ ) {
		double s = yList [i] - m;
		sd += s * s;
	}
	sd /= ( n - 1 );
	sd = sqrt ( sd );
	return sd;
}
void XYData::linearRegression ( double* offset, double* gradient ) const
{
	double* sig = 0;
	double siga;
	double sigb;
	double chi2;
	double q;
	int mwt = 0; // sig (standard deviations) unavailable
	fit ( &xList[0]-1, &yList[0]-1, size (), sig, mwt, offset, gradient, &siga, &sigb, &chi2, &q );
}
void XYData::linearRegression ( double* offset, double* gradient, double* siga, double* sigb, double* chi2, double* q ) const
{
	double* sig = 0;
	int mwt = 0; // sig (standard deviations) unavailable
	fit ( &xList[0]-1, &yList[0]-1, size (), sig, mwt, offset, gradient, siga, sigb, chi2, q );
}
double XYData::getSumOfYInTolRange ( double x, double tol ) const
{
	int start = lower_bound ( xList.begin (), xList.end (), x - tol ) - xList.begin ();
	int end = upper_bound ( xList.begin (), xList.end (), x + tol ) - xList.begin ();
	return accumulate ( yList.begin () + start, yList.begin () + end, 0.0 );
}
void XYData::offsetCalibration ( double offset, const string& units )
{
	double denominator;
	if ( units == "ppm" || units == "%" ) {
		if ( units == "%" )		denominator = 100.0;
		if ( units == "ppm" )	denominator = 1000000.0;
		for ( DoubleVectorSizeType i = 0 ; i < xList.size () ; i++ ) {
			xList [i] -= ( xList [i] * offset / denominator );
		}
	}
	if ( units == "mmu" || units == "Da" ) {
		if ( units == "mmu" )	denominator = 1000.0;
		if ( units == "Da" )	denominator = 1.0;
		for ( DoubleVectorSizeType i = 0 ; i < xList.size () ; i++ ) {
			xList [i] -= ( offset / denominator );
		}
	}
}
