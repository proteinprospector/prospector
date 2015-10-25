/******************************************************************************
*                                                                             *
*  Library    : libnrec                                                       *
*                                                                             *
*  Filename   : ln_dot.h                                                      *
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
*  Copyright (2010-2010) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __ln_dot_h
#define __ln_dot_h

#include <lgen_define.h>

double dotProduct ( const DoubleVector& a, const DoubleVector& b );
double magnitude ( const DoubleVector& a );
double cosSimilarity ( const DoubleVector& a, const DoubleVector& b );
double correlation ( const DoubleVector& a, const DoubleVector& b );

#endif /* ! __ln_dot_h */
