/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_version.h                                                  *
*                                                                             *
*  Created    : February 18th 1997                                            *
*                                                                             *
*  Purpose    : Function to define and print version numbers on program       *
*               outputs.                                                      *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1997-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_version_h
#define __lu_version_h

#include <string>
#include <ostream>

std::string getVersionFromPPXMLFile ( const std::string& filename );
void printProgramInformation ( const std::string& programName );
void printProgramInformationHTML ( std::ostream& os, const std::string& programName );

#endif /* ! __lu_version_h */
