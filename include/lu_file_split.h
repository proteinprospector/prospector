/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_file_split.h                                               *
*                                                                             *
*  Created    : September 17th 2004                                           *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2004-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_file_split_h
#define __lu_file_split_h

#include <lgen_define.h>

class ParameterList;

class FileSplit {
	IntVector startFraction;
	IntVector startSpec;
	IntVector endFraction;
	IntVector endSpec;
	int numSerial;
	int totSpec;
	void init ( const IntVector& nFileSpec, const IntVector& nSearchSpec );
public:
	FileSplit ( const IntVector& nFileSpec, int numProcesses, int maxSpectra );
	FileSplit ( const IntVector& nFileSpec, const DoubleVector& prop, int maxSpectra );
	IntVector getSendData ( int searchNumber, int index ) const;
	int getNumSerial () const { return numSerial; }
	int getTotalSpectra () const { return totSpec; }
	static void setParams ( ParameterList* paramList, const IntVector& iv );
	static const int NUM_PARAMS;
};

#endif /* ! __lu_file_split_h */
