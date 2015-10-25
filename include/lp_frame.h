/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_frame.h                                                    *
*                                                                             *
*  Created    : June 7th 2000                                                 *
*                                                                             *
*  Purpose    :                                                               *
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

#ifndef __lp_frame_h
#define __lp_frame_h

#include <lg_time.h>
#include <lgen_define.h>

class FastaServer;

class FrameIterator {
	FastaServer* fs;
	IntVector indicies;
	int numIndicies;
	int frameTranslationStart;
	int frameTranslationEnd;
	bool tempOverride;
	int indexNum;
	int transNum;
	int orfNum;
	int index;
	bool output;
	CharPtrVector readingFrames;
	int numReadingFrames;
	bool getNextEntry ();
	bool getNextFrameTranslation ();
	static int getIncrement ();
	static int numDatabaseEntries;
	static int numSearches;
	static char outputCharacter;
	static bool mpi;
	static GenElapsedTime elapsedTime;
	static int num;
	static int timeUpdate;
	static int timeout;
	static std::string searchJobID;
public:
	FrameIterator ( FastaServer* fs, int index, int dnaReadingFrame, int openReadingFrame );
	FrameIterator ( FastaServer* fs, const IntVector& indicies, const PairIntInt& frameTransPair, bool temp_override );
	~FrameIterator ();
	char* getNextFrame ();
	int getEntry () const { return index; }
	int getFrameTranslation () const { return ( transNum - 1 ); }
	int getFrame () const { return ( orfNum - 1 ); }
	static bool updateProgress ( bool reset );
	FastaServer* getFS () const { return fs; }
	static void setOutputCharacter ( char ch ) { outputCharacter = ch; }
	static void setMPI () { mpi = true; }
	static void setNumDatabaseEntries ( int n ) { numDatabaseEntries = n; }
	static void setNumSearches ( int n ) { numSearches = n; }
	static void resetElapsedTime ( int startFraction )
	{
		GenElapsedTime et;
		elapsedTime = et;
		num = startFraction;
	}
	static void setSearchJobID ( const std::string& s ) { searchJobID = s; }
	static void searchFinished ();
	static void decrementNum () { num--; }
	static void incrementNum () { num++; }
	static int getNum () { return num; }
	static int getNumSearches () { return numSearches; }
};

#endif /* ! __lp_frame_h */
