/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_composit.h                                                 *
*                                                                             *
*  Created    : January 7th 1997                                              *
*                                                                             *
*  Purpose    : Composition searches based on bit maps.                       *
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

#ifndef __lu_composit_h
#define __lu_composit_h

#include <string>
#include <vector>

class PeakContainer;

class CompositionSearchParameters {
	StringVector compIons;
	std::string compositionExclude;
	bool excludeFlag;
	bool searchFlag;
public:
	CompositionSearchParameters ( const ParameterList* params );

	StringVector getCompIons () const { return compIons; }
	std::string getCompositionExclude () const { return compositionExclude; }
	bool getExcludeFlag () const { return excludeFlag; }
	bool getSearchFlag () const { return searchFlag; }
	void printHTML ( std::ostream& os ) const;
};

class CompositionSearch {
	std::string specifiedComposition;
	std::string excludedComposition;
	std::string knownCompositionSequence;
	unsigned int includeMask;
	unsigned int excludeMask;
	UIntVector eitherOrMasks;
	bool compositionSearchFlag;

	StringVector getIncludeList ( PeakContainer& peaks, const StringVector& compAA, bool removeDuplicates ) const;
	std::string calculateKnownCompositionSequence ( const StringVector& includeList ) const;
	std::string calculateCompositionSequence ( const StringVector& compList ) const;
	UIntVector calculateEitherOrMasks ( const StringVector& includeList, unsigned int includeMask ) const;
public:
	CompositionSearch () {};
	CompositionSearch ( const CompositionSearchParameters& compSearchParams, PeakContainer& peaks );
	bool doCompositionSearch ( const std::string& peptide ) const;
	const char* getKnownCompositionSequence () const { return knownCompositionSequence.c_str (); }
	bool isCompositionSearch () const { return compositionSearchFlag; }
	double subtractKnownComposition ( double mass ) const;
	void printHTML ( std::ostream& os ) const;
	void printXML ( std::ostream& os ) const;
	static unsigned int calculateIncludeMask ( const StringVector& includeList );
};

#endif /* ! __lu_composit_h */
