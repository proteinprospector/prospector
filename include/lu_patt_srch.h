/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_patt_srch.h                                                *
*                                                                             *
*  Created    : September 5th 2001                                            *
*                                                                             *
*  Purpose    : MS-Pattern search and hit classes.                            *
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

#ifndef __lu_patt_srch_h
#define __lu_patt_srch_h

#include <lu_hit.h>
#include <lu_db_srch.h>

class RegularExpression;
class RegularExpressionWithErrors;
class MSPatternParameters;

class PatternHit : public ProteinHit {
	RegularExpression* regExp;
	std::string sequence;
	char previousAA;
	char nextAA;
	short startAA;
	double matchedPepMW;
	int numSubstitutions;
	static void printSequenceWithErrors ( std::ostream& os, const std::string& sequenceHit, RegularExpressionWithErrors* regExp );
public:
	PatternHit ( FastaServer* fs, int index, int dnaReadingFrame, int openReadingFrame, RegularExpression* regExp );
	void printHTMLHeader ( std::ostream& os, bool reportSubstitutions ) const;
	void printHTML ( std::ostream& os, bool reportSubstitutions ) const;
	void printXML ( std::ostream& os ) const;
	friend class sort_hits;
};

class sort_hits {
	public:
		int operator () ( const PatternHit& a, const PatternHit& b ) const
		{
			if ( a.numSubstitutions == b.numSubstitutions ) {
				if ( a.matchedPepMW == b.matchedPepMW ) {
					return ( a.getIndex () < b.getIndex () );
				}
				return ( a.matchedPepMW < b.matchedPepMW );
			}
			return ( a.numSubstitutions < b.numSubstitutions );
		}
};

class PatternHits : public DatabaseHits {
	std::vector <PatternHit> hits;
	bool reportSubstitutions;
public:
	PatternHits ( bool reportSubstitutions ) :
		DatabaseHits (),
		reportSubstitutions ( reportSubstitutions ) {}

	int getIndex ( int num, int dataSet = 0 ) const { return hits [num].getIndex (); }
	int getDBIndex ( int num, int dataSet = 0 ) const { return hits [num].getDBIndex (); }
	std::string getAccessionNumber ( int num, int dataSet = 0 ) const { return hits [num].getAccessionNumber (); }
	int getOpenReadingFrame ( int num, int dataSet = 0 ) const { return hits [num].getOpenReadingFrame (); }
	int getDNAReadingFrame ( int num, int dataSet = 0 ) const { return hits [num].getDNAReadingFrame (); }
	int size ( int dataSet = 0 ) const { return hits.size (); }
	void printHTMLHeader ( std::ostream& os, int num, int dataSet = 0 ) const { hits [num].printHTMLHeader ( os, reportSubstitutions ); }
	void printHTMLHit ( std::ostream& os, int num, int dataSet = 0 ) const
		{ hits [num].printHTML ( os, reportSubstitutions ); }
	void printXMLHit ( std::ostream& os, int num, int dataSet = 0 ) const { hits [num].printXML ( os ); }

	void sortHits () { std::sort ( hits.begin (), hits.end (), sort_hits () ); }
	void addHit ( const PatternHit& hit ) { hits.push_back ( hit ); }
};

class PatternSearch : public DatabaseSearch {
	const MSPatternParameters& patternParams;
	int maxPeptideHits;
	void printParamsBodyHTML ( std::ostream& os ) const;
	void doSearch ( PatternHits* patternHits, FastaServer* fsPtr, int i );
public:
	PatternSearch ( const MSPatternParameters& params );
};

#endif /* ! __lu_patt_srch_h */
