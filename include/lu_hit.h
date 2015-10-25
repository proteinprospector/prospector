/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_hit.h                                                      *
*                                                                             *
*  Created    : May 21st 2001                                                 *
*                                                                             *
*  Purpose    :                                                               *
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

#ifndef __lu_hit_h
#define __lu_hit_h

#include <map>
#include <iostream>
#include <lu_db_entry.h>
#include <lu_prog.h>

class FastaServer;

class ProteinHit {
	static bool uniprot;
protected:
	static std::map <const FastaServer*, int> idxMap;
	FastaServer* fs;
	static bool dnaDatabase;
	DatabaseEntry databaseEntry;
	mutable double proteinMW;
	mutable double proteinPI;
	mutable std::string accessionNumber;
	mutable std::string species;
	mutable std::string name;
	mutable bool entrySet;
	void printAccessionNumber ( std::ostream& os, const std::string& an, bool multi, int num ) const;
	void printDelimitedAccessionNumber ( std::ostream& os, const std::string& an, bool multi ) const;
	void printProteinPI ( std::ostream& os ) const;
	void printProteinMW ( std::ostream& os ) const;
	void getEntry () const;
public:
	ProteinHit () {};
	ProteinHit ( FastaServer* fs, int ind, int drf, int orf );
	void printDelimitedAccNum ( std::ostream& os ) const;
	void printXMLHit ( std::ostream& os ) const;
	void printXML ( std::ostream& os ) const;
	void printHTMLHit ( std::ostream& os ) const;
	void printHTMLHit2 ( std::ostream& os ) const;
	void printHTMLMultiHit ( std::ostream& os ) const;
	void printDelimitedHit ( std::ostream& os ) const;
	void printDelimitedHit2 ( std::ostream& os ) const;
	void printHit ( std::ostream& os ) const;
	virtual void printHTMLHeader ( std::ostream& os ) const;
	void printDelimitedHeader ( std::ostream& os ) const;
	void printDelimitedHeader2 ( std::ostream& os ) const;
	void printHTMLHeader2 ( std::ostream& os ) const;
	void printHTMLHeader3 ( std::ostream& os ) const;
	DatabaseEntry getDatabaseEntry () const { return databaseEntry; }
	int getIndex () const { return databaseEntry.getIndex (); }
	int getDNAReadingFrame () const { return databaseEntry.getDNAReadingFrame (); }
	int getOpenReadingFrame () const { return databaseEntry.getOpenReadingFrame (); }
	double getProteinMW () const;
	double getProteinPI () const;
	std::string getSpecies () const;
	std::string getName () const;
	int getDBIndex () const
	{
		std::map <const FastaServer*, int>::const_iterator cur = idxMap.find ( fs );
		return (*cur).second;
	}
	std::string getAccessionNumber () const;
	bool isDecoy () const;
	static void addFS ( const FastaServer* fs, int num );
	static void reset ();
};

#endif /* ! __lu_hit_h */
