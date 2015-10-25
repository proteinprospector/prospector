/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_species.h                                                  *
*                                                                             *
*  Created    : July 10th 2008                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2008-2009) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_species_h
#define __lu_species_h

#include <lgen_define.h>

class SpeciesMap {
protected:
	MapStringToInt specMap;
public:
	SpeciesMap () {}
	int getNode ( const std::string& s )
	{
		MapStringToIntConstIterator cur = specMap.find ( s );
		if ( cur != specMap.end () )
			return (*cur).second;
		else
			return 0;
	}
};

class UniProtSpeciesMap : public SpeciesMap {
public:
	UniProtSpeciesMap ();
};

class NCBISpeciesMap : public SpeciesMap {
	MapIntToString nodeNameMap;
public:
	NCBISpeciesMap ( bool nodeNames = false );
	std::string getScientificName ( int node ) const;
};

class NCBIMergedNodes;

class NCBINodes {
	MapIntToIntVector childNodes;
	MapIntToInt parent;
	MapIntToInt nodeRank;
	void parseRank ( const std::string& r, int node );
	int getMergedNode ( int node ) const;
	mutable NCBIMergedNodes* mNodes;
public:
	NCBINodes ( bool parentFlag = false, bool childFlag = false, bool rankFlag = false );
	~NCBINodes ();
	void getChildNodes ( int node, IntVector& childNodes ) const;
	int getRank ( int node ) const;
	std::string getRankAsString ( int node ) const;
	IntVector getTaxonomy ( int node ) const;
};

class TaxonomyMatch {
	static int count;
	static NCBINodes* nodes;
	static UniProtSpeciesMap* usm;
	static NCBISpeciesMap* nsm;
	std::string preferredSpeciesStr;
	MapStringToInt matchLengths;
	IntVector preferredMatchTaxonomy;
	int preferredNode;
	size_t prefLen;
public:
	TaxonomyMatch ( const std::string& preferred, bool parentFlag = true, bool childFlag = false, bool rankFlag = false );
	~TaxonomyMatch ();
	int getPreferredMatchLength () const { return prefLen; }
	int getSpeciesMatchLength ( const std::string& sp );
	int getBestSpeciesIndex ( const StringVector& sv );
	bool equal ( const std::string& sp1, const std::string& sp2 );
	bool greaterThan ( const std::string& sp1, const std::string& sp2 );
	void getNodes ( IntVector& n, bool children ) const;
};

#endif /* ! __lu_species_h */
