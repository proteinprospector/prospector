/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_del_proj.h                                                 *
*                                                                             *
*  Created    : September 25th 2007                                           *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2007-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_del_proj_h
#define __lu_del_proj_h

#ifdef MYSQL_DATABASE
void deleteResults ( const std::string& user, const std::string& project, const std::string& results );
bool deleteProject ( const std::string& user, const std::string& project );
#endif

#endif /* ! __lu_del_proj_h */
