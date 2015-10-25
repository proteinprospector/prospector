/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lg_stdlib.h                                                   *
*                                                                             *
*  Created    : June 23rd 1996                                                *
*                                                                             *
*  Purpose    : Machine independent interface to stdlib.h.                    *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1996-2010) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lg_stdlib_h
#define __lg_stdlib_h

#include <string>
#include <cstdlib>

#ifndef VIS_C
#include <unistd.h>			/* Declaration of unlink */
#endif

int genSystem ( const std::string& command, const std::string& workingDirectory = "", bool cont = false );
char* gen_getenv ( char* name, const std::string& error_location_string );
void gen_qsort ( char* base, int num, int width, int ( *compare ) ( const void*, const void* ) );

#endif /* ! __lg_stdlib_h */
