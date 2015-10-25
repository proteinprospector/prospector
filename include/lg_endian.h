/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lg_endian.h                                                   *
*                                                                             *
*  Created    : July 29th 2005                                                *
*                                                                             *
*  Purpose    : Function to swap "endainness".                                *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2005-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lg_endian_h
#define __lg_endian_h

#include <lgen_define.h>

inline void genEndianConvert ( unsigned short& x )
{
	x = (x>>8) | (x<<8);
}
inline void genEndianConvert ( unsigned int& x )
{
	x = (x>>24) | ((x<<8) & 0x00FF0000) | ((x>>8) & 0x0000FF00) | (x<<24);
}
inline void genEndianConvert ( GENUINT64& x )
{
#ifdef VIS_C
	x = (x>>56) | ((x<<40) & 0x00FF000000000000) | ((x<<24) & 0x0000FF0000000000) | ((x<<8) & 0x000000FF00000000) | ((x>>8) & 0x00000000FF000000) | ((x>>24) & 0x0000000000FF0000) | ((x>>40) & 0x000000000000FF00) | (x<<56);
#else
	x = (x>>56) | ((x<<40) & 0x00FF000000000000LL) | ((x<<24) & 0x0000FF0000000000LL) | ((x<<8) & 0x000000FF00000000LL) | ((x>>8) & 0x00000000FF000000LL) | ((x>>24) & 0x0000000000FF0000LL) | ((x>>40) & 0x000000000000FF00LL) | (x<<56);
#endif
}

#endif /* ! __lg_endian_h */
