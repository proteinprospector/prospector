/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_coverage.h                                                 *
*                                                                             *
*  Created    : July 4th 2003                                                 *
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

#ifndef __lu_coverage_h
#define __lu_coverage_h

#include <lgen_define.h>

class CoverageMap {
	CharVector aaCovered;
	bool error;
	mutable bool entrySet;
	mutable int coverCount;
	mutable double percentCoverage;
	void calculateStats () const;
public:
	CoverageMap ();
	CoverageMap ( int protLen );
	void resize ( int protLen );
	void setCoverage ( const int start, const int end, const unsigned char val = 1 );
	int getCoverCount () const;
	double getPercentCoverage () const;
	CharVector getAACovered () const { return aaCovered; }
	void putCGI ( std::ostream& os ) const;
	void printCoverCountHTML ( std::ostream& os ) const;
};

#endif /* ! __lu_coverage_h */
