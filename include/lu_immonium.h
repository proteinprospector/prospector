/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_immonium.h                                                 *
*                                                                             *
*  Created    : May 30th 2001                                                 *
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
*  Copyright (2001-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_immonium_h
#define __lu_immonium_h

#include <string>
#include <lgen_define.h>

class PeakContainer;
class Peak;

class ImmoniumInfo {
	static bool initialised;
	static DoubleVector immoniumMasses;
	static StringVector immoniumCompositions;
	static UIntVector immoniumMasks;
	static CharVector immoniumMajorIons;

	static DoubleVector immoniumMajorIonMasses;
	static StringVector immoniumMajorIonCompositions;

	static DoubleVector immoniumExcludeAAMasses;
	static StringVector immoniumExcludeAACompositions;

	static double immoniumTolerance;
	static double minimumFragmentMass;

	static void initialiseImmonium ();
	static void initialise ();
	static std::string getBracketedString ( const std::string& str );
public:
	static const CharVector& getMajorIons ();

	static const DoubleVector& getMajorIonMasses ();
	static const StringVector& getMajorIonCompositions ();

	static double getTolerance ();
	static double getMinFragmentMass ();

	static bool inImmoniumRegion ( double mass );
	static unsigned int massToMask ( const Peak* peak );
	static std::string getExcludeListFromPeaks ( const PeakContainer& peaks, const std::string& excludeAA = "" );
	static DoubleVector getRelatedIons ( const char aa );
};

#endif /* ! __lu_immonium_h */
