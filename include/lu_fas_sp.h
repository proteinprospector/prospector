/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_fas_sp.h                                                   *
*                                                                             *
*  Created    : July 13th 2000                                                *
*                                                                             *
*  Purpose    : Functions to do species pre-search.                           *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2000-2009) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_fas_sp_h
#define __lu_fas_sp_h

#include <lgen_define.h>

class TaxonomyMatch;

class TaxonomyNodeList {
	std::vector <TaxonomyMatch*> tm;
public:
	TaxonomyNodeList ();
	~TaxonomyNodeList ();
	void getTaxonomyNodes ( const std::string& species, IntVector& nodes, bool sort );
};

class TaxonomyNames {
	StringVector names;
	TaxonomyNames ();
public:
	static TaxonomyNames& instance ();
	StringVector getList () const { return names; }
};

class TaxonomyGroupNames {
	StringVector names;
	MapStringToStringVector nameAliases;
	TaxonomyGroupNames ();
public:
	static TaxonomyGroupNames& instance ();
	StringVector getList () const { return names; }
	bool check ( const std::string& s ) const;
	StringVector getAliases ( const std::string& s ) const;
};

class TaxonomySearch {
	int numDatabaseEntries;
	StringVector speciesList;
	IntVector indicies;
	IntVector taxList;
	IntVector taxListSize;
public:
	TaxonomySearch ( const std::string& databaseName, int numDatabaseEntries, const StringVector& speciesList, bool useCache );
	IntVector getIndicies () const { return indicies; }
	IntVector getTaxList () const { return taxList; }
	IntVector getTaxListSize () const { return taxListSize; }
};

#endif /* ! __lu_fas_sp_h */
