 /******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_script.h                                                   *
*                                                                             *
*  Created    : January 28th 2003                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2003-2008) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_script_h
#define __lu_script_h

#include <string>

class ParameterList;

void writeScript ( ParameterList* paramList );
void writePerlScript ( const ParameterList* paramList, const std::string& filename );
void writeBatchFile ( const ParameterList* paramList, const std::string& filename );
void writeParamsXML ( ParameterList* paramList, const std::string& programName );

#endif /* ! __lu_script_h */
