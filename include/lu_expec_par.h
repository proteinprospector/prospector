/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_expec_par.h                                                *
*                                                                             *
*  Created    : July 20st 2006                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2006-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_expec_par_h
#define __lu_expec_par_h

class ParameterList;

ParameterList* getExpectationParams ( const ParameterList* searchParamList, std::string& outputName, bool force );

#endif /* ! __lu_expec_par_h */
