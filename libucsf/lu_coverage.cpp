/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_coverage.cpp                                               *
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
*  Copyright (2003-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <stdexcept>
#include <iostream>
#include <lg_string.h>
#include <lu_coverage.h>
using std::fill;
using std::ostream;
using std::endl;
using std::string;
using std::runtime_error;

CoverageMap::CoverageMap () :
	entrySet ( false ),
	error ( false )
{
}
CoverageMap::CoverageMap ( int protLen ) :
	entrySet ( false )
{
	resize ( protLen );
}
void CoverageMap::resize ( int protLen )
{
	error = ( protLen == 0 );
	if ( !error ) {
		aaCovered.resize ( protLen );
		fill ( aaCovered.begin (), aaCovered.end (), false );
	}
}
void CoverageMap::setCoverage ( const int start, const int end, const unsigned char val )
{
	if ( !error ) {
		int len = aaCovered.size ();
		if ( start > len || end > len ) {
			string err ( "A problem has been encountered setting a coverage map.\n Start AA = " );
			err += gen_itoa ( start );
			err += ". End AA = ";
			err += gen_itoa ( end );
			err += ". Protein Length = ";
			err += gen_itoa ( len );
			err += ".\n";
			throw runtime_error ( err );
		}
		else
			fill ( aaCovered.begin () + start - 1, aaCovered.begin () + end, val );
	}
}
int CoverageMap::getCoverCount () const
{
	if ( !entrySet ) calculateStats ();
	return coverCount;
}
double CoverageMap::getPercentCoverage () const
{
	if ( !entrySet ) calculateStats ();
	return percentCoverage;
}
void CoverageMap::calculateStats () const
{
	coverCount = 0;
	if ( !error ) {
		for ( CharVectorSizeType i = 0 ; i < aaCovered.size () ; i++ ) {
			if ( aaCovered [i] > 0 ) coverCount++;
		}
		percentCoverage = coverCount * 100.0 / aaCovered.size ();
	}
	else {
		percentCoverage = 0.0;
	}
	entrySet = true;
}
void CoverageMap::putCGI ( ostream& os ) const
{
	if ( !entrySet ) calculateStats ();
	if ( !aaCovered.empty () ) {
		CharVectorSizeType i, j;
		unsigned char val = aaCovered [0];
		os << "coverage_map=" << ( val ? 1 : 0 );
		for ( i = 1, j = 1 ; i < aaCovered.size () ; i++ ) {
			if ( aaCovered [i] == val ) j++;
			else {
				val = aaCovered [i];
				os << "+" << j;
				j = 1;
			}
		}
		os << "+" << j << "&";
	}
}
void CoverageMap::printCoverCountHTML ( ostream& os ) const
{
	if ( !entrySet ) calculateStats ();
	os << "The&nbsp;matched&nbsp;peptides&nbsp;cover&nbsp;<b>";
	os << percentCoverage;
	os << "%</b>&nbsp;";
	os << "(" << coverCount << "/" << aaCovered.size ();
	os << "AA's)&nbsp;of&nbsp;the&nbsp;protein.";
	os << "<br />" << endl ;
}
