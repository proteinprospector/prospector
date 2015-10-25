/******************************************************************************
*                                                                             *
*  Program    : searchCompare                                                 *
*                                                                             *
*  Filename   : sc_quan.h                                                     *
*                                                                             *
*  Created    : July 9th 2012                                                 *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2012-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __sc_quan_h
#define __sc_quan_h

#include <lgen_define.h>
#include <lu_const_mod.h>
#include <lu_pros_form.h>

class RawFile;
class AACalculator;
class Tolerance;
class PeakFitData;
class QuanPeptide;
class SpecID;
class SearchResults;
typedef std::vector <SearchResults*> SearchResultsPtrVector;
typedef SearchResultsPtrVector::size_type SearchResultsPtrVectorSizeType;

class PeptidePositionQuan {
	static RawFile* rf;
	static bool quanMSMSFlag;
	static bool n15Flag;
	static std::string qType;
	static bool qFlag;
	static double defaultResolution;
	static int numIsotopePeaks;
	static StringVector instrument;
	static StringVectorVector rawTypes;
	static StringVector version;
	static VectorConstParameterListPtr params;
	static std::vector <Tolerance*> parentTolerances;
	static DoubleVectorVector offsets;
	static AACalculator* aaCalc;
	static bool quanMultiNormalFlag;
	static double rtIntervalStart;
	static double rtIntervalEnd;
public:
	static void initialiseQuan ( int i, int fraction );
	static void deleteQuan ( int i );
	static PeakFitData* getQuanRatio ( double mOverZ, int charge, const QuanPeptide& qp, const SpecID& specID, int searchIndex );
	static int getColspan ( int searchNumber );
	static void printHeaderHTML ( std::ostream& os, int searchNumber, const std::string& styleID );
	static void printHTML ( std::ostream& os, const PeakFitData* quanRatio, const std::string& styleID, bool joint );
	static void printQuanBlankHTML ( std::ostream& os, int searchNumber, const std::string& styleID );
	static void printHeaderDelimited ( std::ostream& os, int searchNumber );
	static void printDelimited ( std::ostream& os, const PeakFitData* quanRatio, bool joint );
	static void printQuanBlankDelimited ( std::ostream& os, int searchNumber );
	static bool outputQuanResults ( std::ostream& os, const PeakFitData* quanRatio, const std::string& searchName, int numRepeats, bool area );
	static DoubleVector getIntensityRatios ( const PeakFitData* quanRatio );
	static DoubleVector getAreaRatios ( const PeakFitData* quanRatio );
	static void initialiseParams ( const ParameterList* p );
	static void initialiseDataSetInfo ( double timeWindowStart, double timeWindowEnd );
	static void initialise ( const SearchResultsPtrVector& searchResults );
	static void initialiseAACalculator ( const MapStringConstModPtr& constMods, const std::string& rawType, const std::string& quanType, double resolution );

	static int getNumQuanRatioColumns ();
	static bool getQuanMSMSFlag () { return quanMSMSFlag; }
};

#endif /* ! __sc_quan_h */
