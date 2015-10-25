/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lg_memory.cpp                                                 *
*                                                                             *
*  Created    : August 1st 1996                                               *
*                                                                             *
*  Purpose    : Machine independent memory manipulation functions.            *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1996-2009) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <cstring>
#include <lg_stdlib.h>

void gen_bzero ( char* block, int size )
{
	if ( size == 0 ) return;

	memset ( (void*) block, 0, (size_t)size );
}
void gen_bcopy ( char* source, char* dest, int size )
{
	if ( size == 0 ) return;

	memcpy ( (void*) dest, (void*) source, (size_t)size );
}
