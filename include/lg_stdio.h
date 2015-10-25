/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lg_stdio.h                                                    *
*                                                                             *
*  Created    : June 23rd 1996                                                *
*                                                                             *
*  Purpose    : Machine independent interface to stdio.h.                     *
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

#ifndef __lg_stdio_h
#define __lg_stdio_h

#include <string>

#ifdef VIS_C
#include <io.h>
#endif

#include <cstdio>

int gen_open ( const std::string& path, int flags, const std::string& errorLocation );
int gen_close ( int fd, const std::string& errorLocation );
FILE* gen_fopen ( const std::string& filename, const char* type, const std::string& errorLocation );
FILE* gen_fopen_text ( const std::string& filename, const char* type, const std::string& errorLocation );
FILE* gen_fopen_binary ( const std::string& filename, const char* type, const std::string& errorLocation );
int gen_fclose ( FILE* fp, const std::string& errorLocation );
int gen_fread ( char* ptr, size_t size, size_t nitems, FILE* fp, const std::string& errorLocation );
int gen_fwrite ( char* ptr, size_t size, size_t nitems, FILE* fp, const std::string& errorLocation );
void gen_stdin_buffering_off ();
void gen_stdout_buffering_off ();
void gen_stderr_buffering_off ();
void gen_stdin_binary_mode ();

#endif /* ! __lg_stdio_h */
