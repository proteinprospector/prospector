/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_indicies.cpp                                               *
*                                                                             *
*  Created    : April 9th 2000                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2000-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lu_getfil.h>
using std::string;
using std::fill;

static DoubleVectorVector indexV;
static DoubleVector chosenIndex;
static StringVector names;
static bool initialised = false;

static void initialiseIndexScore ();
static inline void initialise ();

static inline void initialise ()
{
	if ( initialised == false ) initialiseIndexScore ();
}
void setIndexType ( const string& indexType )
{
	initialise ();
	for ( StringVectorSizeType i = 0 ; i < names.size () ; i++ ) {
		if ( indexType == names [i] ) {
			chosenIndex = indexV [i];
			return;
		}
	}
}
double indexScore ( const string& qseq )
{
	double score = 0.0;
	for ( StringSizeType i = 0 ; i < qseq.length () ; i++ ) {
		score += chosenIndex [qseq [i] - 'A'];
	}
	return ( score );
}
static void initialiseIndexScore ()
{
	int numEntries;
	char* info = getFileInfo ( MsparamsDir::instance ().getParamPath ( "indicies.txt" ), '>', 1, true, &numEntries );

	for ( int i = 0 ; i < numEntries ; i++ ) {
		names.push_back ( ( i == 0 ) ? strtok ( info, "\n" ) : strtok ( NULL, "\n" ) );
		DoubleVector dv (52);
		indexV.push_back ( dv );
		fill ( dv.begin (), dv.end (), 0.0 );
		for ( ; ; ) {
			char aa;
			double value;
			char* line = strtok ( NULL, "\n" );

			if ( !strcmp ( line, ">" ) ) break;
			sscanf ( line, "%c %lf", &aa, &value );
			indexV [i][aa-'A'] = value;
		}
	}
	initialised = true;
}
