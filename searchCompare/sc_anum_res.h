/******************************************************************************
*                                                                             *
*  Program    : searchCompare                                                 *
*                                                                             *
*  Filename   : sc_anum_res.h                                                 *
*                                                                             *
*  Created    : April 2nd 2003                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2003-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __sc_anum_res_h
#define __sc_anum_res_h

#include <map>
#include <ostream>
#include <vector>
#include <string>
#include <sc_search_res.h>

class ANumResults {
	StringVector aNumList;
public:
	ANumResults ( const MapStringToStringVector& aNumMap );
	StringVector getAccessionNumbers () { return aNumList; }
};

class CheckAccessionNumber {
	const StringVector& aNumFilter;
	bool checkEquality; 
public:
	CheckAccessionNumber ( const StringVector& aNumFilter, bool checkEquality ) :
		aNumFilter ( aNumFilter ),
		checkEquality ( checkEquality ) {}
	bool operator () ( const SearchResultsPeptideHit* pp )
	{
		bool flag = false;
		const std::string& aNum = pp->getAccessionNumber ();
		for ( StringVectorSizeType i = 0 ; i < aNumFilter.size () ; i++ ) {
			if ( aNumFilter [i] == aNum ) {
				flag = true;
				break;
			}
		}
		if ( checkEquality ) return flag;
		else return !flag;
	}
};

#endif /* ! __sc_anum_res_h */
