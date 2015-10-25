/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lgen_process.h                                                *
*                                                                             *
*  Created    : June 23rd 2006                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2006-2009) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lgen_process_h
#define __lgen_process_h

#include <lgen_define.h>

bool killProcess ( int pid );
void genReduceProcessPriority ();
IntVector getProcessNumberList ();
bool isProcessRunning ( int pid );
bool isProcessRunning ( const std::string& processName );
BoolDeque areProcessesRunning ( const IntVector& pid );
void genSleep ( int millisecs );
int genCreateProcess ( const std::string& fullExePath, const std::string& params = "", const std::string& dir = "" );

#ifndef VIS_C
void genInitSigchld ( void ( *handler ) ( int sigNum ) );
void genInitSigterm ( void ( *handler ) ( int sigNum ) );
void genInitSighup ( void ( *handler ) ( int sigNum ) );
#endif

#endif /* ! __lgen_process_h */
