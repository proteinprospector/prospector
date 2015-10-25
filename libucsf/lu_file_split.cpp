/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_file_split.cpp                                             *
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
#include <numeric>
#include <lg_string.h>
#include <lu_file_split.h>
#include <lu_param_list.h>
using std::accumulate;

const int FileSplit::NUM_PARAMS = 4;

FileSplit::FileSplit ( const IntVector& nFileSpec, int numProcesses, int maxSpectra )
{
	totSpec = accumulate ( nFileSpec.begin (), nFileSpec.end (), 0 );
	int nSpecPerProcess = totSpec / numProcesses;			// Number of spectra per process
	numSerial = ( nSpecPerProcess / maxSpectra ) + 1;		// Number of serial searches
	int numSearches = numProcesses * numSerial;				// Number of searches
	int nSpecPerSearch = totSpec / numSearches;				// Number of spectra per search
	int rem = totSpec % numSearches;
	IntVector nSearchSpec;
	for ( int i = 0 ; i < numSearches ; i++ ) {
		int n = nSpecPerSearch;
		if ( rem != 0 ) {
			n += 1;
			rem--;
		}
		nSearchSpec.push_back ( n );
	}
	init ( nFileSpec, nSearchSpec );
}
FileSplit::FileSplit ( const IntVector& nFileSpec, const DoubleVector& prop, int maxSpectra )
{
	totSpec = accumulate ( nFileSpec.begin (), nFileSpec.end (), 0 );
	int maxSpecPerProcess = static_cast<int> ( totSpec * prop [0] );
	numSerial = ( maxSpecPerProcess / maxSpectra ) + 1;		// Number of serial searches
	IntVector nSearchSpec;
	for ( DoubleVectorSizeType i = 0 ; i < prop.size () ; i++ ) {
		for ( int j = 0 ; j < numSerial ; j++ ) {
			nSearchSpec.push_back ( static_cast<int> ( ( totSpec / numSerial ) * prop [i] ) );
		}
	}
	int rem = accumulate ( nSearchSpec.begin (), nSearchSpec.end (), 0 ) - totSpec;
	while ( rem != 0 ) {
		int inc = ( rem < 0 ) ? 1 : -1;
		for ( DoubleVectorSizeType i = 0 ; i < prop.size () ; i++ ) {
			nSearchSpec [i] += inc;
			rem += inc;
			if ( rem == 0 ) break;
		}
	}
	init ( nFileSpec, nSearchSpec );
}
void FileSplit::init ( const IntVector& nFileSpec, const IntVector& nSearchSpec )
{
	int spec = 0;
	for ( IntVectorSizeType j = 0, k = 0 ; j < nSearchSpec.size () ; j++ ) {
		startFraction.push_back ( k+1 );
		startSpec.push_back ( spec+1 );
		int r = nFileSpec [k] - spec;		// Number of spectra remaining in this fraction
		int s = nSearchSpec [j];			// Number of spectra for this search
		if ( r >= s ) {						// All spectra available from this fraction
			s += spec;
		}
		else {								// Get some spectra from other fractions
			while ( r < s ) {
				k++;
				s -= r;
				r = nFileSpec [k];
			}
		}
		endFraction.push_back ( k+1 );
		endSpec.push_back ( s );
		if ( nFileSpec [k] == s ) {
			if ( k+1 < nFileSpec.size () ) k++;
			spec = 0;
		}
		else spec = s;
	}
}
IntVector FileSplit::getSendData ( int searchNumber, int index ) const
{
	int numProcesses = startFraction.size () / numSerial;
	int i = searchNumber + ( index * numProcesses );
	IntVector sendData;
	sendData.push_back ( startFraction [i] );
	sendData.push_back ( startSpec [i] );
	sendData.push_back ( endFraction [i] );
	sendData.push_back ( endSpec [i] );
	return sendData;
}
void FileSplit::setParams ( ParameterList* paramList, const IntVector& iv )
{
	paramList->setValue ( "start_fraction",	gen_itoa ( iv [0] ) );
	paramList->setValue ( "start_spectrum",	gen_itoa ( iv [1] ) );
	paramList->setValue ( "end_fraction",	gen_itoa ( iv [2] ) );
	paramList->setValue ( "end_spectrum",	gen_itoa ( iv [3] ) );
}
