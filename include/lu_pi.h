/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_pi.h                                                       *
*                                                                             *
*  Created    : May 30th 2001                                                 *
*                                                                             *
*  Purpose    : Functions to calculate the pI of a protein.                   *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_pi_h
#define __lu_pi_h

#include <ostream>
#include <string>

class ProteinPI {
	double pi;
public:
	ProteinPI ( const std::string& protein_string );
	double getProteinPI () const {return pi;};
	static int piPrecision;
};
std::ostream& operator<< ( std::ostream& os, const ProteinPI& ppi );

const double NO_PK_VALUE = -98.0;
const double PI_NOT_CALCULATED = -99.0;

#endif /* ! __lu_pi_h */
