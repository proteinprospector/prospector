/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_mass_conv.cpp                                              *
*                                                                             *
*  Created    : August 9th 2005                                               *
*                                                                             *
*  Purpose    : Functions for dealing with multiply charged data.             *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2005-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lu_mass.h>
#include <lu_mass_conv.h>
#include <lu_mass_elem.h>

static const char* adduct = "H";

double getAdductMass ( bool monoisotopicFlag )
{
	return monoisotopicFlag ? formula_to_monoisotopic_mass ( adduct ) : formula_to_average_mass ( adduct );
}
double mOverZToMPlusH ( double mOverZ, int charge, bool monoisotopicFlag )
{
	return mOverZToMPlusH ( mOverZ, charge, getAdductMass ( monoisotopicFlag ) );
}
// X = M + qH - qe   ---- (1)
//     -----------
//         q
//
// Y = M + H - e     ---- (2)
//
// From (1)
//
// M = qX - qH + qe
//
// Combining (2) and (3)
//
// Y = qX -qH + qe + H - e 
//
// Y = qX - H ( q - 1 ) + e ( q - 1)
//
// Y = qX + ( e - H ) ( q - 1 )
//
double mOverZToMPlusH ( double mOverZ, int charge, double adductMass )
{
	adductMass -= ELECTRON_REST_MASS;
	if ( charge < 0 )
		return ( mOverZ * -charge ) + ( ( -charge + 1 ) * adductMass );
	else if ( charge == 0 )
		return mOverZ + adductMass;
	else
		if ( charge != 1 )
			return ( mOverZ * charge ) - ( ( charge - 1 ) * adductMass );
		else
			return mOverZ;
}
double mOverZToM ( double mOverZ, int charge, bool monoisotopicFlag )
{
	return mOverZToM ( mOverZ, charge, getAdductMass ( monoisotopicFlag ) );
}
double mOverZToM ( double mOverZ, int charge, double adductMass )
{
	return charge * ( mOverZ - adductMass + ELECTRON_REST_MASS );
}
double mPlusHToMOverZ ( double mass, int precursorCharge, int charge, bool monoisotopicFlag )
{
	double adductMass = getAdductMass ( monoisotopicFlag ) - ELECTRON_REST_MASS;
	return ( ( mass + (( precursorCharge - 1 ) * adductMass ) ) / charge );
}
double mPlusHToMOverZ ( double mass, int charge, bool monoisotopicFlag )
{
	return mPlusHToMOverZ ( mass, charge, getAdductMass ( monoisotopicFlag ) );
}
double mPlusHToMOverZ ( double mass, int charge, double adductMass )
{
	adductMass -= ELECTRON_REST_MASS;
	if ( charge < 0 )
		return ( ( mass - ( ( -charge + 1 ) * adductMass ) ) / -charge );
	else if ( charge == 0 )
		return mass - adductMass;
	else
		if ( charge != 1 )
			return ( ( mass + (( charge - 1 ) * adductMass ) ) / charge );
		else
			return mass;
}
double mToMOverZ ( double mass, int charge, bool monoisotopicFlag )
{
	return ( mass / charge ) + getAdductMass ( monoisotopicFlag ) - ELECTRON_REST_MASS;
}
double mToMOverZ ( double mass, int charge, double adductMass )
{
	return ( mass / charge ) + adductMass - ELECTRON_REST_MASS;
}
