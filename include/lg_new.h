/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lg_new.h                                                      *
*                                                                             *
*  Created    : July 10th 2000                                                *
*                                                                             *
*  Purpose    : Machine independent memory allocation functions based on new. *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2000-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lg_new_h
#define __lg_new_h

char* gen_new_string ( const char* string );
char* genNewString ( const char* start, const char* end );
void gen_set_new_handler ();

#endif /* ! __lg_new_h */
