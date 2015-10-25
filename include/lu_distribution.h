/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_distribution.h                                             *
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

#ifndef __lu_distribution_h
#define __lu_distribution_h

#include <lgen_define.h>

class TheoreticalDistribution {
	DoubleVector m;
	DoubleVectorVector id;
	DoubleVectorVector normCoeffVector;
public:
	TheoreticalDistribution ( const std::string& type );
	DoubleVector& getDistribution ( double mass, int charge );
	const DoubleVector& getNormDistribution ( double mOverZ, int charge ) const;
	DoubleVector getOverlapNormDistribution ( double mOverZ, int charge, double int1, double int2, int offset ) const;
	int maxShift ( double mOverZ, int charge, const DoubleVector& intensity );
};

#endif /* ! __lu_distribution_h */
