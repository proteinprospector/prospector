/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_pre_srch.h                                                 *
*                                                                             *
*  Created    : July 18th 2000                                                *
*                                                                             *
*  Purpose    : Functions to do all pre-searches.                             *
*                                                                             *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2000-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_pre_srch_h
#define __lu_pre_srch_h

#include <vector>
#include <lgen_define.h>

class ParameterList;
class PreSearch;
class FastaServer;
typedef std::vector <FastaServer*> FastaServerPtrVector;

class PreSearchInfo {
	std::vector <PreSearch*> preSearch;
	IntVector preSearchIndex;
	bool fromFile;
	static std::string RESULTS_FROM_FILE;
	static bool RESULTS_FROM_FILE_DEFAULT;
	PreSearchInfo& operator= ( PreSearchInfo& rhs );
	PreSearchInfo ( const PreSearchInfo& rhs );
	static PairBoolStringVector getAccNumSearch ( const std::string& database, const MapStringToStringVector& accessionNumbers );
	static StringVector getAddAccessionNumbers ( const std::string& database, const MapStringToStringVector& addAccessionNumbers );
public:
	PreSearchInfo ( const ParameterList* params, const std::string& prefix = "" );
	~PreSearchInfo ();
	void doSearch ( const FastaServerPtrVector& fs ) const;
	const IntVector& getIndicies ( int num ) const;
	int getNumIndicies ( int num ) const;
	void printHTML ( std::ostream& os ) const;
	void printBodyXML ( std::ostream& os ) const;
	bool getFromFile () const { return fromFile; }
	void putCGI ( std::ostream& os ) const;
	void putHiddenFormEntryJavascript ( std::ostream& os ) const;
	static void setReportTaxonomy ();
	static void setAccessionNumbers ( MapStringToStringVector& accessionNumbers, const StringVector& aNum );
};

#endif /* ! __lu_pre_srch_h */
