/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_export_proj.h                                              *
*                                                                             *
*  Created    : June 29th 2012                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2012-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_export_proj_h
#define __lu_export_proj_h

#ifdef MYSQL_DATABASE
class ParameterList;

void exportData ( const ParameterList& paramList );
void compressProjects ( const ParameterList& paramList );
void uncompressProjects ( const ParameterList& paramList );
bool uncompressProjects2 ( const ParameterList& paramList );
void checkProjects ( const ParameterList& paramList );
#endif

#endif /* ! __lu_export_proj_h */
