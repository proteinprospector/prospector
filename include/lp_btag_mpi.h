/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_btag_mpi.h                                                 *
*                                                                             *
*  Created    : September 20th 2004                                           *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2004-2011) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_btag_mpi_h
#define __lu_btag_mpi_h

#include <lgen_define.h>

class ParameterList;
class FileSplit;

class BatchTagMPI {
	int argc;
	char** argv;
	int numSearches;
	int rank;
	int numRuns;
	int searchStage;
	int startSerial;
#ifdef MYSQL_DATABASE
	std::string searchJobID;
#endif
	std::vector <ParameterList*> paramList;
	void rankZeroSearchLoop ( const FileSplit* fs, const StringVector& filesToDelete, int startSerial );
	void otherRankSearchLoop ( const char* message, int startSerial );
public:
	BatchTagMPI ( int argc, char** argv );
	~BatchTagMPI ();
	void run ();
};

#endif /* ! __lu_btag_mpi_h */
