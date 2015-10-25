/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_pp_param.h                                                 *
*                                                                             *
*  Created    : July 17th 2003                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2003-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_pp_param_h
#define __lu_pp_param_h

#include <string>

namespace PPParamParser {
int getIntValue ( const std::string& str, const std::string& tag );
double getDoubleValue ( const std::string& str, const std::string& tag );
std::string getStringValue ( const std::string& str, const std::string& tag );
}

#endif /* ! __lu_pp_param_h */
