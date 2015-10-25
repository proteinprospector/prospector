/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_mask.cpp                                                   *
*                                                                             *
*  Created    : July 26th 2001                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2009) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lu_mask.h>
#include <lu_mass.h>
using std::string;

CompositionMask::CompositionMask ( const string& str )
{
	mask = 0;
	for ( StringSizeType i = 0 ; i < str.length () ; i++ ) {
		mask |= aa_composition_mask [str [i]];
	}
}
