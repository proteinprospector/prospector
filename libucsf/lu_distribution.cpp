/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_distribution.cpp                                           *
*                                                                             *
*  Created    : July 14th 2010                                                *
*                                                                             *
*  Purpose    :                                                               *
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
#include <algorithm>
#include <lg_io.h>
#include <lu_getfil.h>
#include <lu_distribution.h>
#include <lu_mass_conv.h>

using std::string;
using std::getline;
using std::istringstream;
using std::lower_bound;

//1.0029 - mono to first
//1.0025 - first to second
//1.002 - second to third

TheoreticalDistribution::TheoreticalDistribution ( const string& type )
{
	GenIFStream fromFile ( MsparamsDir::instance ().getParamPath ( "distribution.txt" ) );
	string line;
	bool flag = false;
	while ( getline ( fromFile, line ) ) {
		if ( line.length () != 0 && line [0] != '#' ) {
			if ( line == type ) flag = true;
			while ( getline ( fromFile, line ) ) {
				if ( line [0] == '>' ) break;
				if ( flag ) {
					istringstream ist ( line );
					double d;
					DoubleVector dv;
					DoubleVector normCoeff;
					double normFactor = 0.0;
					if ( ist >> d ) {
						m.push_back ( d );
						while ( ist >> d ) {
							if ( d ) {
								dv.push_back ( d );
								normFactor += d * d;
							}
						}
					}
					normFactor = sqrt ( normFactor );
					for ( int i = 0 ; i < dv.size () ; i++ ) {
						normCoeff.push_back ( dv [i] / normFactor );
					}
					id.push_back ( dv );
					normCoeffVector.push_back ( normCoeff );
				}
			}
		}
	}
}
DoubleVector& TheoreticalDistribution::getDistribution ( double mOverZ, int charge )
{
	double mph = mOverZToMPlusH ( mOverZ, charge, true );
	int ind = lower_bound ( m.begin (), m.end (), mph ) - m.begin ();
	return ind < id.size () ? id [ind] : id [id.size ()-1];
}
const DoubleVector& TheoreticalDistribution::getNormDistribution ( double mOverZ, int charge ) const
{
	double mph = mOverZToMPlusH ( mOverZ, charge, true );
	int ind = lower_bound ( m.begin (), m.end (), mph ) - m.begin ();
	return ind < normCoeffVector.size () ? normCoeffVector [ind] : normCoeffVector [normCoeffVector.size ()-1];
}
DoubleVector TheoreticalDistribution::getOverlapNormDistribution ( double mOverZ, int charge, double int1, double int2, int offset ) const
{
	double mph = mOverZToMPlusH ( mOverZ, charge, true );
	int ind = lower_bound ( m.begin (), m.end (), mph ) - m.begin ();
	DoubleVector dv = id [ind];
	double ratio = ( int2 - ( ( int1 * dv [offset] )/dv [0]) ) / int1;
	DoubleVector dv2;
	double normFactor = 0.0;
	for ( DoubleVectorSizeType i = 0 ; i < dv.size () ; i++ ) {
		double d = dv [i];
		if ( i >= offset )	d += ratio * dv [i-offset];
		dv2.push_back ( d );
		normFactor += d * d;
	}
	for ( DoubleVectorSizeType j = offset ; j < dv.size () ; j++ ) {
		double d = ratio * dv [j];
		dv2.push_back ( d );
		normFactor += d * d;
	}
	normFactor = sqrt ( normFactor );
	for ( DoubleVectorSizeType k = 0 ; k < dv2.size () ; k++ ) {
		dv2 [k] /= normFactor;
	}
	return dv2;
}
// This function calculates if there is a shift between a measured distribtion and a theoretical one based
// on their maxima
int TheoreticalDistribution::maxShift ( double mOverZ, int charge, const DoubleVector& intensity )
{
	double mph = mOverZToMPlusH ( mOverZ, charge, true );
	int ind = lower_bound ( m.begin (), m.end (), mph ) - m.begin ();
	DoubleVector& av = id [ind];
	int idx1 = std::max_element ( av.begin (), av.end () ) - av.begin ();
	int idx2 = std::max_element ( intensity.begin (), intensity.end () ) - intensity.begin ();
std::cout << mOverZ << " " << idx2-idx1 << "<br />" << std::endl;
	return idx2-idx1;
}
