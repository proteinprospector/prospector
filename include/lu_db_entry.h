/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_db_entry.h                                                 *
*                                                                             *
*  Created    :                                                               *
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

#ifndef __lu_db_entry_h
#define __lu_db_entry_h

#include <ostream>

class DatabaseEntry {
	int index;
	int dnaReadingFrame;
	int openReadingFrame;
public:
	DatabaseEntry () {};
	DatabaseEntry ( int index, int dnaReadingFrame, int openReadingFrame );
	bool operator== ( const DatabaseEntry& rhs ) const;
	int getIndex () const { return index; }
	int getDNAReadingFrame () const { return dnaReadingFrame; }
	int getOpenReadingFrame () const { return openReadingFrame; }
	void printHTMLHeader ( std::ostream& os, bool dnaDatabase ) const;
	void printDelimitedHeader ( std::ostream& os, bool dnaDatabase ) const;
	void printHTMLHit ( std::ostream& os, int num, bool dnaEntry ) const;
	void printHTML ( std::ostream& os, int num, bool dnaDatabase, bool dnaEntry ) const;
	void printDelimited ( std::ostream& os, int num, bool dnaDatabase, bool dnaEntry ) const;
	void printGapHTML ( std::ostream& os, bool dnaDatabase ) const;
	void printXML ( std::ostream& os, bool dnaEntry ) const;
};

#endif /* ! __lu_db_entry_h */
