/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_mass_conv.h                                                *
*                                                                             *
*  Created    : August 9th 2005                                               *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
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

#ifndef __lu_mass_conv_h
#define __lu_mass_conv_h

double getAdductMass ( bool monoisotopicFlag );
double mOverZToMPlusH ( double mOverZ, int charge, bool monoisotopicFlag );
double mOverZToMPlusH ( double mOverZ, int charge, double adductMass );
double mOverZToM ( double mOverZ, int charge, bool monoisotopicFlag );
double mOverZToM ( double mOverZ, int charge, double adductMass );
double mPlusHToMOverZ ( double mass, int precursorCharge, int charge, bool monoisotopicFlag );
double mPlusHToMOverZ ( double mass, int charge, bool monoisotopicFlag );
double mPlusHToMOverZ ( double mass, int charge, double adductMass );
double mToMOverZ ( double mass, int charge, bool monoisotopicFlag );
double mToMOverZ ( double mass, int charge, double adductMass );

#endif /* ! __lu_mass_conv_h */
