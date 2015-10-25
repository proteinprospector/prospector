/******************************************************************************
*                                                                             *
*  Library    : libnrec                                                       *
*                                                                             *
*  Filename   : moment.cpp                                                    *
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
*  Copyright (1998-2010) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <nr.h>

void moment ( double* data, int n, double* ave )	// This version just calculates the mean value
{
	if ( n <= 0 ) throw lNrecMomentLessThanOneDataValue ();
	if ( n <= 1 ) *ave = data [1];
	double s = 0.0;
	for ( int j = 1 ; j <= n ; j++ ) {
		s += data [j];
	}
	*ave = s / n;
}
void moment_weighted ( double* data, double* wt, int n, double s_wt, double* ave )	// This version just calculates the weighted mean value
{
	if ( n <= 0 ) throw lNrecMomentLessThanOneDataValue ();
	if ( n <= 1 ) *ave = data [1];
	if ( s_wt == 0.0 ) {					// sum of weights may already have been calculated
		for ( int i = 1 ; i <= n ; i++ ) {
			s_wt += wt [i];
		}
	}
	double s = 0.0;
	for ( int j = 1 ; j <= n ; j++ ) {
		s += data [j] * wt [j];
	}
	*ave = s / ( n * s_wt );
}
void moment ( double* data, int n, double* ave, double* sdev, double* svar )
{
	int j;

	if ( n <= 0 ) throw lNrecMomentLessThanOneDataValue ();
	if ( n <= 1 ) {
		*ave = data [1];
		throw lNrecMomentLessThanTwoDataValues ();
	}
	double s = 0.0;
	for ( j = 1 ; j <= n ; j++ ) {
		s += data [j];
	}
	*ave = s / n;
	*svar = 0.0;
	for ( j = 1 ; j <= n ; j++ ) {
		s = data [j] - (*ave);
		*svar += s * s;
	}
	*svar /= ( n - 1 );
	*sdev = sqrt ( *svar );
}
void moment ( double* data, int n, double* ave, double* adev, double* sdev, double* svar, double* skew, double* curt )
{
	int j;
	double s;
	double p;

	if ( n <= 0 ) throw lNrecMomentLessThanOneDataValue ();
	if ( n <= 1 ) {
		*ave = data [1];
		throw lNrecMomentLessThanTwoDataValues ();
	}
	s = 0.0;
	for ( j = 1 ; j <= n ; j++ ) {
		s += data [j];
	}
	*ave = s / n;
	*adev = (*svar) = (*skew) = (*curt) = 0.0;
	for ( j = 1 ; j <= n ; j++ ) {
		*adev += fabs (s = data [j] - (*ave) );
		*svar += (p = s * s );
		*skew += (p *= s);
		*curt += (p *= s);
	}
	*adev /= n;
	*svar /= ( n - 1 );
	*sdev = sqrt ( *svar );
	if ( *svar ) {
		*skew /= (n*(*svar)*(*sdev));
		*curt = (*curt) / ( n * (*svar) * (*svar) ) - 3.0;
	}
	else throw lNrecMomentZeroVariance ();
}
double median ( const DoubleVector& sortedData )
{
	int siz = sortedData.size ();
	int ind = (int) (0.5 * ( siz - 1 ));
	if ( genIsOdd ( siz ) )
		return sortedData [ind];
	else
		return 0.5 * (sortedData[ind+1] + sortedData[ind]);
}
double lowerQ ( const DoubleVector& sortedData )
{
	return median ( DoubleVector ( sortedData.begin (), sortedData.begin () + 1 + (int) (0.5 * ( sortedData.size () - 1 )) ) );
}
double upperQ ( const DoubleVector& sortedData )
{
	return median ( DoubleVector ( sortedData.begin () + (int) (0.5 * ( sortedData.size () )), sortedData.end () ) );
}
