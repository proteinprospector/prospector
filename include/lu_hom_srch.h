/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_hom_srch.h                                                 *
*                                                                             *
*  Created    : September 5th 2001                                            *
*                                                                             *
*  Purpose    : MS-Homology search and hit classes.                           *
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

#ifndef __lu_hom_srch_h
#define __lu_hom_srch_h

#include <lu_hit.h>
#include <lu_db_srch.h>

class MSHomologyParameters;

class HomologyHit : public ProteinHit {
	std::string sequence;
	std::string previousAA;
	std::string nextAA;
	short startAA;
	std::string measuredSequence;
	int sequenceSet;
	int numSubstitutions;
	mutable double score;
	mutable bool scoreSet;
	static int numPreviousAA;
	static int numNextAA;
	std::string getPreviousAA ( const char* protein ) const;
	std::string getNextAA () const;
public:
	double proteinScore;
	HomologyHit ( FastaServer* fs, int index, int dnaReadingFrame, int openReadingFrame, const char* protein, const std::string& measuredSequence, int sequenceSet, int numSubstitutions );
	void printHTMLHeader ( std::ostream& os ) const;
	void printDelimitedHeader ( std::ostream& os ) const;
	void printHTML ( std::ostream& os, bool newProtein ) const;
	void printDelimited ( std::ostream& os, bool newProtein ) const;
	void printXML ( std::ostream& os, bool newProtein ) const;
	double getScore () const;
	double getSequenceSet () const { return sequenceSet; }
	friend class sort_multi_sequence_list_hits;
	friend class sort_hits_by_score;
	static void setNumPreviousAA ( int num );
	static void setNumNextAA ( int num );
};

class sort_multi_sequence_list_hits {
	public:
		bool operator () ( const HomologyHit& a, const HomologyHit& b ) const
		{
			if ( a.getIndex () == b.getIndex () ) {
				if ( a.sequenceSet == b.sequenceSet ) {
					return ( strcmp ( a.measuredSequence.c_str (), b.measuredSequence.c_str () ) < 0 );
				}
				return ( a.sequenceSet < b.sequenceSet );
			}
			return ( a.getIndex () < b.getIndex () );
		}
};

class sort_hits_by_score {
	public:
		bool operator () ( const HomologyHit& a, const HomologyHit& b ) const
		{
			if ( a.proteinScore == b.proteinScore ) {
				if ( a.getIndex () == b.getIndex () ) {
					return ( strcmp ( a.measuredSequence.c_str (), b.measuredSequence.c_str () ) < 0 );
				}
				return ( a.getIndex () < b.getIndex () );
			}
			return ( b.proteinScore < a.proteinScore );
		}
};

class ProteinHomologyHit {
	int dbIndex;
	int index;
	std::string accessionNumber;
	int dnaReadingFrame;
	int openReadingFrame;
	double score;
public:
	ProteinHomologyHit ( int dbIndex, int index, const std::string& accessionNumber, int dnaReadingFrame, int openReadingFrame, double score ) :
		dbIndex ( dbIndex ),
		index ( index ),
		accessionNumber ( accessionNumber ),
		dnaReadingFrame ( dnaReadingFrame ),
		openReadingFrame ( openReadingFrame ),
		score ( score ) {}
	int getDBIndex () const { return dbIndex; }
	std::string getAccessionNumber () const { return accessionNumber; }
	int getIndex () const { return index; }
	int getDNAReadingFrame () const { return dnaReadingFrame; }
	int getOpenReadingFrame () const { return openReadingFrame; }
};

class HomologyHits : public DatabaseHits {
	std::vector <HomologyHit> peptideHits;
	std::vector <ProteinHomologyHit> proteinHits;
public:
	HomologyHits () : DatabaseHits () {}

	int getIndex ( int num, int dataSet = 0 ) const { return proteinHits [num].getIndex (); }
	int getDBIndex ( int num, int dataSet = 0 ) const { return proteinHits [num].getDBIndex (); }
	std::string getAccessionNumber ( int num, int dataSet = 0 ) const { return proteinHits [num].getAccessionNumber (); }
	int getOpenReadingFrame ( int num, int dataSet = 0 ) const { return proteinHits [num].getOpenReadingFrame (); }
	int getDNAReadingFrame ( int num, int dataSet = 0 ) const { return proteinHits [num].getDNAReadingFrame (); }
	int size ( int dataSet = 0 ) const { return proteinHits.size (); }
	void printHTMLHeader ( std::ostream& os, int num, int dataSet = 0 ) const { peptideHits [num].printHTMLHeader ( os ); }
	void printHTMLHit ( std::ostream& os, int num, int dataSet = 0 ) const {}
	void printHTMLHit2 ( std::ostream& os, bool newProtein, int num ) const
		{ peptideHits [num].printHTML ( os, newProtein ); }
	void printDelimitedHeader ( std::ostream& os, int num, int dataSet = 0 ) const { peptideHits [num].printDelimitedHeader ( os ); }
	void printDelimitedHit2 ( std::ostream& os, bool newProtein, int num ) const
		{ peptideHits [num].printDelimited ( os, newProtein ); }
	void printXMLHit ( std::ostream& os, int num, int dataSet = 0 ) const {}
	void printXMLHit ( std::ostream& os, bool newProtein, int num ) const
		{ peptideHits [num].printXML ( os, newProtein ); }

	void calculateScores ( int minMatches, const std::string& scoreMatrix );
	void addHit ( const HomologyHit& hit ) { peptideHits.push_back ( hit ); }

	bool differentProteins ( int num1, int num2 ) const
		{ return peptideHits [num1].getIndex () != peptideHits [num2].getIndex (); }
	int peptideHitsSize () const { return peptideHits.size (); }
};

class HomologySearch : public DatabaseSearch {
	HomologyHits* homologyHits;
	const MSHomologyParameters& homologyParams;
	int maxPeptideHits;
	void sequenceSearch ();
	void doSearch ( FastaServer* fsPtr, int i );
	void add_hit ( FastaServer* fs, int index, int dnaReadingFrame, int openReadingFrame, char* protein );
	void printParamsBodyHTML ( std::ostream& os ) const;
	void printHTMLHitsTable ( std::ostream& os ) const;
	void printDelimitedHits ( std::ostream& os ) const;
	void printXMLHits ( std::ostream& os ) const;
public:
	HomologySearch ( const MSHomologyParameters& params );
};

#endif /* ! __lu_hom_srch_h */
