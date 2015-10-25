/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_sim_ent.h                                                  *
*                                                                             *
*  Created    : May 17th 2001                                                 *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_sim_ent_h
#define __lu_sim_ent_h

#include <cstdio>
#include <string>
#include <vector>
#include <lu_pi.h>
#include <lu_mass_seq.h>

class ParameterList;

class SingleEntryParameters {
	std::string database;
	std::string accessMethod;

	StringVector accessionNum;
	IntVector indexNum;
	IntVector dnaReadingFrame;
	IntVector openReadingFrame;
	int maxNTermAA;
	StringVector userProteinSequence;
	StringVector names;
	mutable std::string tempDBFullPath;
	void addUserProteinSequence ( const std::string& value );
	void addDBUserProteinSequence ( const std::string& value );
	std::string checkSequence ( const std::string& s ) const;
	std::string convertSequence ( const std::string& s ) const;
	std::string convertDBSequence ( const std::string& s ) const;
	void parseForUserProtein ( const ParameterList* params );
	void parseForEntryData ( const ParameterList* params );
	static void writeProteinEntry ( FILE* fp, const std::string& n, const std::string& protein );
public:
	SingleEntryParameters ( const ParameterList* params, const std::string& tempDatabase = "" );

	std::string getDatabase () const { return database; }
	std::string getAccessMethod () const { return accessMethod; }
	std::string getAccessionNum ( int i = 0 ) const { return accessionNum [i]; }
	StringVector getAccessionNums () const { return accessionNum; }
	bool getDNADatabase () const { return !dnaReadingFrame.empty (); }
	int getIndexNum ( int i = 0 ) const { return indexNum [i]; }
	int getDNAReadingFrame ( int i = 0 ) const { return !dnaReadingFrame.empty () ?  dnaReadingFrame [i] : 0; }
	int getOpenReadingFrame ( int i = 0 ) const { return !openReadingFrame.empty () ? openReadingFrame [i] : 0; }
	int getMaxNTermAA () const { return maxNTermAA; }
	std::string getUserProteinSequence ( int i = 0 ) const { return userProteinSequence [i]; }
	std::string getName ( int i = 0 ) const { return i < names.size () ? names [i] : ""; }
	int getNumEntries () const;

	void printHTML ( std::ostream& os ) const;
	void putCGI ( std::ostream& os ) const;
	static void putIndexCGI ( std::ostream& os, bool dnaFlag, int indexNum, int dnaReadingFrame, int openReadingFrame );
	static void putIndexHiddenFormEntry ( std::ostream& os, bool dnaFlag, int indexNum, int dnaReadingFrame, int openReadingFrame );
	void putHiddenFormEntry ( std::ostream& os ) const;
	static void copyToHiddenFormEntry ( std::ostream& os, const ParameterList* params );
	PairStringBool createTemporaryDatabase ( bool random, bool reverse, const std::string& tempDir = "" ) const;
	std::string getTempDBFullPath () const { return tempDBFullPath; }
};

class SingleEntry {
	std::string prot;
	std::string singleEntryName;
	ProteinPI ppi;
	ProteinMW pmw;
	std::string aaComposition;
protected:
	bool dnaDatabase;
	int indexNum;
	int drf;
	int orf;
public:
	SingleEntry ( const std::string& protein, const std::string& singleEntryName, bool dnaDatabase, int indexNumber, int dnaReadingFrame, int openReadingFrame );
	const char* getProtein () const { return prot.c_str (); }
	int getProteinLength () const { return prot.length (); }
	virtual int getIndex () const { return indexNum; }
	virtual void printHTML ( std::ostream& os ) const;
	void printXML ( std::ostream& os ) const;
	virtual void printXMLBody ( std::ostream& os ) const;
};

typedef std::vector <SingleEntry*> SingleEntryPtrVector;
typedef SingleEntryPtrVector::size_type SingleEntryPtrVectorSizeType;

std::vector <SingleEntry*> getSingleEntry ( const SingleEntryParameters& p );

#endif /* ! __lu_sim_ent_h */
