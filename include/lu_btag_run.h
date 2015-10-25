/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_btag_run.h                                                 *
*                                                                             *
*  Created    : October 2nd 2007                                              *
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

#ifndef __lu_btag_run_h
#define __lu_btag_run_h

class ParameterList;

void initialiseProspectorBTag ( const std::string& searchKey = "" );
int getMSMSMaxSpectra ( ParameterList* paramList );
void runBatchTag ( ParameterList* paramList, int numSearches, int searchNumber, const std::string& searchJobID, int startSerial );
void joinResultsFiles ( const ParameterList* paramList, int numSearches );
std::string getBatchTagOutputPath ( const ParameterList* params );

#endif /* ! __lu_btag_run_h */
