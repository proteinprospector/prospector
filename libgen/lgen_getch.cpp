/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lgen_getch.cpp                                                *
*                                                                             *
*  Created    : September 25th 1996                                           *
*                                                                             *
*  Purpose    : getyes function.                                              *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1996-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <cctype>
#ifdef VIS_C
#include <conio.h>
#else
#include <lg_stdio.h>
#endif

int gen_getyes ()
{
	int ch;

#ifdef VIS_C
	do {
		ch = _getch ();
		ch = toupper ( ch );
	} while ( ch != 'Y' && ch != 'N' );

	_putch ( ch );
	_putch ( '\r' );
	_putch ( '\n' );
	if ( ch == 'Y' ) return ( 1 );
	else return ( 0 );
#else
	scanf( "%c", &ch );
	if ( toupper ( ch ) == 'Y' ) return ( 1 );
	else return ( 0 );
#endif
}
