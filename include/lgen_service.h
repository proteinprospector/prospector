/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lgen_service.h                                                *
*                                                                             *
*  Created    : October 22nd 2007                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2007-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lgen_service_h
#define __lgen_service_h

#ifdef VIS_C
#include <string>
void genStartDaemon ( void ( *passedDaemon ) () );
bool genIsServiceRunning ( const std::string& serviceName );
bool genStartService ( const std::string& serviceName, int argc = 0, const char** argv = 0 );
bool genStopService ( const std::string& serviceName );
void genInstallService ( const std::string& serviceName, const std::string& serviceDisplayName, const std::string& serviceBinaryPath, const std::string& user, const std::string& password );
void genUninstallService ( const std::string& serviceName );
#endif

#endif /* ! __lgen_service_h */
