/******************************************************************************
*                                                                             *
*  Library    : libnrec                                                       *
*                                                                             *
*  Filename   : factln.cpp                                                    *
*                                                                             *
*  Created    :                                                               *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  :                                                               *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1998-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <nr.h>

double factln ( int n )
{
	static double a [101];

	if ( n < 0 ) nrerror ( "Negative factorial in routine FACTLN" );
	if ( n <= 1 ) return 0.0;
	if ( n <= 100 ) return a[n] ? a[n] : (a[n]=gammln(n+1.0));
	else return gammln ( n + 1.0 );
}
