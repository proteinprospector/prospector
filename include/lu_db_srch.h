/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_db_srch.h                                                  *
*                                                                             *
*  Created    : September 5th 2001                                            *
*                                                                             *
*  Purpose    : Database search and hit base classes.                         *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_db_srch_h
#define __lu_db_srch_h

#include <ostream>
#include <string>
#include <lu_dig_par.h>
#include <lu_program.h>

class MSSearchParameters;
class PreSearch;
class FastaServer;
typedef std::vector <FastaServer*> FastaServerPtrVector;

class DatabaseHits {
public:
	DatabaseHits () {}
	virtual ~DatabaseHits ();
	virtual int getIndex ( int num, int dataSet = 0 ) const = 0;
	virtual int getDBIndex ( int num, int dataSet = 0 ) const = 0;
	virtual std::string getAccessionNumber ( int num, int dataSet = 0 ) const = 0;
	virtual int getOpenReadingFrame ( int num, int dataSet = 0 ) const = 0;
	virtual int getDNAReadingFrame ( int num, int dataSet = 0 ) const = 0;
	virtual int size ( int dataSet = 0 ) const = 0;
	virtual void printHTMLHeader ( std::ostream& os, int num, int dataSet = 0 ) const = 0;
	virtual void printHTMLHit ( std::ostream& os, int num, int dataSet = 0 ) const = 0;
	virtual void printXMLHit ( std::ostream& os, int num, int dataSet = 0 ) const = 0;
};

class DBSearch : public MSProgram {
	void printBodyHTML ( std::ostream& os );
	void printBodyTabDelimitedText ( std::ostream& os );
	void printBodyXML ( std::ostream& os );
	static StringVector tempDirs;
	static FastaServerPtrVector userFS;
protected:
	FastaServerPtrVector fs;
	const MSSearchParameters& params;
	VectorPairIntInt dnaFrameTranslationPairVector;
	DatabaseHits* databaseHits;
	int numHits;
	virtual void printHTMLHits ( std::ostream& os ) = 0;
	virtual void printDelimitedHits ( std::ostream& os ) const {}
	virtual void printXMLHits ( std::ostream& os ) const = 0;
public:
	DBSearch ( const MSSearchParameters& params );
	virtual ~DBSearch ();
	static void deleteTempDirs ();
};

class DatabaseSearch : public DBSearch {
protected:
	virtual void printHTMLHits ( std::ostream& os );
	virtual void printHTMLHitsReport ( std::ostream& os );
	virtual void printHTMLHitsTable ( std::ostream& os ) const;
	virtual void printXMLHits ( std::ostream& os ) const;
public:
	DatabaseSearch ( const MSSearchParameters& params );
	virtual ~DatabaseSearch ();
	virtual void printHTMLHitsJavascript ( std::ostream& os ) const;
	virtual void saveHits () const;
};

#endif /* ! __lu_db_srch_h */
