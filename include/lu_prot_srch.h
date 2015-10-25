/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_prot_srch.h                                                *
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
*  Copyright (2001-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_prot_srch_h
#define __lu_prot_srch_h

#include <lu_hit.h>
#include <lu_db_srch.h>

class MSSearchParameters;

class ProteinHits : public DatabaseHits {
	std::vector <ProteinHit> hits;
public:
	ProteinHits () : DatabaseHits () {}
	int getIndex ( int num, int dataSet = 0 ) const { return hits [num].getIndex (); }
	int getDBIndex ( int num, int dataSet = 0 ) const { return hits [num].getDBIndex (); }
	std::string getAccessionNumber ( int num, int dataSet = 0 ) const { return hits [num].getAccessionNumber (); }
	int getOpenReadingFrame ( int num, int dataSet = 0 ) const { return hits [num].getOpenReadingFrame (); }
	int getDNAReadingFrame ( int num, int dataSet = 0 ) const { return hits [num].getDNAReadingFrame (); }
	int size ( int dataSet = 0 ) const { return hits.size (); }
	void printHTMLHeader ( std::ostream& os, int num, int dataSet = 0 ) const { hits [num].printHTMLHeader ( os ); }
	void printHTMLHit ( std::ostream& os, int num, int dataSet = 0 ) const
		{ hits [num].printHTMLHit2 ( os ); }
	void printXMLHit ( std::ostream& os, int num, int dataSet = 0 ) const { hits [num].printXML ( os ); }

	void addHit ( const ProteinHit& hit ) { hits.push_back ( hit ); }
};

class ProteinSearch : public DatabaseSearch {
	void printParamsBodyHTML ( std::ostream& os ) const;
public:
	ProteinSearch ( const MSSearchParameters& params );
};

#endif /* ! __lu_prot_srch_h */
