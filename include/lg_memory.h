/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lg_memory.h                                                   *
*                                                                             *
*  Created    : August 1st 1996                                               *
*                                                                             *
*  Purpose    : Machine independent interface to memory.h.                    *
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

#ifndef __lg_memory_h
#define __lg_memory_h

void gen_bzero ( char* block, int size );
void gen_bcopy ( char* source, char* dest, int size );

#endif /* ! __lg_memory_h */
