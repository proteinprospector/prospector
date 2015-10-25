/******************************************************************************
*                                                                             *
*  Program    : searchCompare                                                 *
*                                                                             *
*  Filename   : sc_anum_res.cpp                                               *
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
*  Copyright (2003-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef VIS_C
#pragma warning( disable : 4503 )	// Disable warnings about decorated name length being too long
#endif
#include <lgen_define.h>
#include <sc_anum_res.h>
using std::string;
using std::stable_sort;
using std::unique;

ANumResults::ANumResults ( const MapStringToStringVector& aNumMap )
{
	for ( MapStringToStringVectorConstIterator i = aNumMap.begin () ; i != aNumMap.end () ; i++ ) {
		string db = (*i).first;
		StringVector sv = (*i).second;
		int dbIdx = db.empty () ? 1 : ProteinInfo::getDBIndex ( db );
		for ( StringVectorSizeType j = 0 ; j < sv.size () ; j++ ) {
			string a = sv [j];
			int n = gen_strcharcount ( a, '$' );
			if ( n == 0 || n == 2 ) {
				a = gen_itoa ( dbIdx ) + "$" + sv [j];
				n++;
			}
			ProteinInfo pi ( a );
			string ac = pi.getFirstAccessionNumber ();
			string::size_type start = 0;
			string::size_type end = 0;
			string s = genNextString ( a, "$", start, end );
			if ( n == 1 ) {	// Protein
				a = s + "$" + ac;
			}
			else {	// n == 3 ---- DNA
				string s2 = genNextString ( a, "$", start, end );
				a = s + "$" + ac + "$" + a.substr ( start );
			}
			aNumList.push_back ( a );
		}
	}
	stable_sort ( aNumList.begin (), aNumList.end () );
	aNumList.erase ( unique ( aNumList.begin (), aNumList.end () ), aNumList.end () );
	ProteinInfo::resetANumMap ();
}
