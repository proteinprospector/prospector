/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_get_imm.cpp                                                *
*                                                                             *
*  Created    : March 12th 1997                                               *
*                                                                             *
*  Purpose    : Gets the immonium ion masses and corresponding compositional  *
*               information from the params/imm.txt file.                     *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1997-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lu_immonium.h>
#include <lu_getfil.h>
#include <lu_mass.h>
#include <lu_mass_elem.h>
#include <lu_charge.h>
using std::string;

const CharVector& ImmoniumInfo::getMajorIons ()
{
	initialise ();
	return immoniumMajorIons;
}
const DoubleVector& ImmoniumInfo::getMajorIonMasses ()
{
	initialise ();
	return immoniumMajorIonMasses;
}
const StringVector& ImmoniumInfo::getMajorIonCompositions ()
{
	initialise ();
	return immoniumMajorIonCompositions;
}
double ImmoniumInfo::getTolerance ()
{
	initialise ();
	return immoniumTolerance;
}
double ImmoniumInfo::getMinFragmentMass ()
{
	initialise ();
	return minimumFragmentMass;
}

bool ImmoniumInfo::initialised = false;

DoubleVector ImmoniumInfo::immoniumMasses;
StringVector ImmoniumInfo::immoniumCompositions;
UIntVector ImmoniumInfo::immoniumMasks;
CharVector ImmoniumInfo::immoniumMajorIons;

DoubleVector ImmoniumInfo::immoniumMajorIonMasses;
StringVector ImmoniumInfo::immoniumMajorIonCompositions;

DoubleVector ImmoniumInfo::immoniumExcludeAAMasses;
StringVector ImmoniumInfo::immoniumExcludeAACompositions;

double ImmoniumInfo::immoniumTolerance;
double ImmoniumInfo::minimumFragmentMass;

void ImmoniumInfo::initialise ()
{
	if ( initialised == false ) initialiseImmonium ();
}
void ImmoniumInfo::initialiseImmonium ()
{
	static const int NUM_HEADER_LINES = 2;
	int numLines;
	char* info = getParamsFileInfo ( "imm.txt", &numLines );
	int numImmoniumIons = numLines - NUM_HEADER_LINES;

	immoniumTolerance = atof ( strtok ( info, "\n" ) );
	minimumFragmentMass = atof ( strtok ( NULL, "\n" ) );
	for ( int i = 0 ; i < numImmoniumIons ; i++ ) {
		ElementalFormula ef = strtok ( NULL, "|" );
		immoniumMasses.push_back ( formula_to_monoisotopic_mass ( ef ) - ELECTRON_REST_MASS );
		immoniumCompositions.push_back ( strtok ( NULL, "|" ) );
		immoniumMasks.push_back ( string_to_mask ( immoniumCompositions [i] ) );
		immoniumMajorIons.push_back ( ( strtok ( NULL, "|" ) ) [0] );
		if ( immoniumMajorIons [i] == 'M' ) {
			immoniumMajorIonMasses.push_back ( immoniumMasses [i] );
			immoniumMajorIonCompositions.push_back ( immoniumCompositions [i] );
		}
		strtok ( NULL, "|" );	// Immonium flag
		char* allowExclude = strtok ( NULL, "\n" );
		if ( allowExclude [0] != '-' ) {
			immoniumExcludeAAMasses.push_back ( immoniumMasses [i] );
			immoniumExcludeAACompositions.push_back ( allowExclude );
		}
	}
	initialised = true;
}
string ImmoniumInfo::getBracketedString ( const string& str )
{
	int len = str.length ();
	if ( len != 1 ) return str.substr ( 1, len-2 );
	else return str;
}
bool ImmoniumInfo::inImmoniumRegion ( double mass )
{
	initialise ();
	double maxImmoniumMass = immoniumMasses [immoniumMasses.size ()-1] + immoniumTolerance;
	if ( mass <= maxImmoniumMass && mass <= minimumFragmentMass ) {
		return true;
	}
	return false;
}
unsigned int ImmoniumInfo::massToMask ( const Peak* peak )
{
	initialise ();
	unsigned int mask = 0;
	for ( DoubleVectorSizeType i = 0 ; i < immoniumMasses.size () ; i++ ) {
		if ( peak->isMatch ( immoniumMasses [i] ) ) {
			mask |= immoniumMasks [i];
		}
	}
	return mask;
}
string ImmoniumInfo::getExcludeListFromPeaks ( const PeakContainer& peaks, const string& excludeAA )
{
	initialise ();
	string excludeList = excludeAA;
	for ( DoubleVectorSizeType i = 0 ; i < immoniumExcludeAAMasses.size () ; i++ ) {
		bool flag = true;
		for ( PeakContainerSizeType j = 0 ; j < peaks.size () ; j++ ) {
			if ( genAbsDiff ( immoniumExcludeAAMasses [i], peaks [j]->getMass () ) <= immoniumTolerance ) {
				flag = false;
				break;
			}
		}
		if ( flag ) excludeList += immoniumExcludeAACompositions [i];
	}
	string excludeList2;
	for ( StringSizeType a = 0 ; a < excludeList.size () ; a++ ) {			// Remove duplicates.
		bool flag = false;
		for ( StringSizeType b = a + 1 ; b < excludeList.size () ; b++ ) {
			if ( excludeList [b] == excludeList [a] ) flag = true;
		}
		if ( flag == false ) excludeList2 += excludeList [a];
	}
	return excludeList2;
}
DoubleVector ImmoniumInfo::getRelatedIons ( const char aa )
{
	initialise ();
	DoubleVector ret;
	unsigned int mask = string_to_mask ( string ( 1, aa ) );
	for ( DoubleVectorSizeType i = 0 ; i < immoniumMasses.size () ; i++ ) {
		if ( mask & immoniumMasks [i] ) ret.push_back ( immoniumMasses [i] );
	}
	return ret;
}
