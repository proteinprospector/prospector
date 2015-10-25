/******************************************************************************
*                                                                             *
*  Library    : libnrec                                                       *
*                                                                             *
*  Filename   : ln_dot.cpp                                                    *
*                                                                             *
*  Created    : July 13th 2010                                                *
*                                                                             *
*  Purpose    : Dot product functions.                                        *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2010-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <cmath>
#include <lgen_define.h>

double dotProduct ( const DoubleVector& a, const DoubleVector& b )
{
	int siz = genMin ( a.size (), b.size () );
	double dot = 0.0;
	for ( int i = 0 ; i < siz ; i++ ) {
		dot += a[i] * b[i];
	}
	return dot;
}
double magnitude ( const DoubleVector& a )
{
	double mag = 0.0;
	int siz = a.size ();
	for ( int i = 0 ; i < siz ; i++ ) {
		mag += a[i] * a[i];
	}
	return sqrt ( mag );
}
double cosSimilarity ( const DoubleVector& a, const DoubleVector& b )
{
	int siz = genMin ( a.size (), b.size () );
	double dot = 0.0;
	double magA = 0.0;
	double magB = 0.0;
	for ( int i = 0 ; i < siz ; i++ ) {
		const double& a1 = a[i];
		const double& b1 = b[i];
		dot += a1 * b1;
		magA += a1 * a1;
		magB += b1 * b1;
	}
	return dot / sqrt ( magA * magB );
}
double correlation ( const DoubleVector& a, const DoubleVector& b )
{
	int siz = genMin ( a.size (), b.size () );
	double sumA = 0.0;
	double sumB = 0.0;
	double sumASquared = 0.0;
	double sumBSquared = 0.0;
	double sumAB = 0.0;
	for ( int i = 0 ; i < siz ; i++ ) {
		const double& a1 = a[i];
		const double& b1 = b[i];
		sumA += a1;
		sumB += b1;
		sumASquared += a1 * a1;
		sumBSquared += b1 * b1;
		sumAB += a1 * b1;
	}
	return ( ( siz * sumAB ) - ( sumA * sumB ) ) / sqrt ( ( ( siz * sumASquared ) - ( sumA * sumA ) ) * ( ( siz * sumBSquared ) - ( sumB * sumB ) ) );
}
