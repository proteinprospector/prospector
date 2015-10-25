/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_frame.cpp                                                  *
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
*  Copyright (2000-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef PP_MULTI
#ifdef VIS_C
#define MPICH_SKIP_MPICXX
// This is required for MPICH2 version 1.0.2, it may be fixed by the next version also the mpicxx.h needs replacing with the one on homer
#endif
#include <mpi.h>
#endif
#include <lgen_error.h>
#include <lu_getfil.h>
#include <lu_fasta.h>
#include <lu_html.h>
#include <lp_frame.h>
#ifdef MYSQL_DATABASE
#include <ld_init.h>
#endif
using std::cout;
using std::endl;
using std::string;
using std::ostringstream;

FrameIterator::FrameIterator ( FastaServer* fs, int index, int dnaReadingFrame, int openReadingFrame ) :
	fs ( fs ),
	index ( index ),
	transNum ( dnaReadingFrame + 1 ),
	orfNum ( openReadingFrame ),
	output ( false )
{
}
FrameIterator::FrameIterator ( FastaServer* fs, const IntVector& indicies, const PairIntInt& frameTransPair, bool tempOverride ) :
	fs ( fs ),
	indicies ( indicies ),
	numIndicies ( indicies.size () ),
	frameTranslationStart ( frameTransPair.first ),
	frameTranslationEnd ( frameTransPair.second ),
	tempOverride ( tempOverride ),
	output ( false )
{
#ifdef PP_MULTI
	if ( mpi ) {
		char c = 'n';
		MPI_Bsend ( &c, 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD );
		MPI_Bsend ( &numIndicies, 1, MPI_INT, 0, 0, MPI_COMM_WORLD );
	}
#endif
	numDatabaseEntries = numIndicies;
	indexNum = 0;
	transNum = frameTranslationStart;
	orfNum = 0;
	if ( getNextEntry () ) {
		getNextFrameTranslation ();
	}
}
FrameIterator::~FrameIterator ()
{
#ifdef PP_MULTI
	if ( mpi ) {
		char c = 'x';
		MPI_Bsend ( &c, 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD );
	}
	else
#endif
	if ( output ) updateProgress ( true );
}
bool FrameIterator::getNextEntry ()
{
	if ( indexNum >= numIndicies ) {
		return false;
	}
	if ( (indexNum+1) % getIncrement () == 0 ) {
		if ( tempOverride == false ) {
#ifdef PP_MULTI
			if ( mpi ) {
				MPI_Bsend ( &outputCharacter, 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD );
			}
			else {
#endif
				if ( outputCharacter != ' ' ) {
					output = true;
					if ( updateProgress ( false ) ) {
						ErrorHandler::genError ()->error ( "The search has been abandoned as it is likely the timeout will be exceeded.\n" );
					}
					cout.flush ();
					if ( cout.fail () )
						exit ( 0 );
				}
#ifdef PP_MULTI
			}
#endif
		}
	}
	index = indicies [indexNum];
	indexNum++;
	return true;
}
bool FrameIterator::getNextFrameTranslation ()
{
	if ( transNum > frameTranslationEnd ) {
		if ( getNextEntry () == false ) return false;
		transNum = frameTranslationStart;
	}
	char* protein = fs->get_fasta_protein ( index, transNum++ );
	readingFrames = fs->split_protein_to_reading_frames ( protein );
	numReadingFrames = readingFrames.size ();
	return true;
}
char* FrameIterator::getNextFrame ()
{
	if ( numIndicies == 0 ) return NULL;
	if ( orfNum >= numReadingFrames ) {
		if ( getNextFrameTranslation () == false ) {
			return NULL;
		}
		orfNum = 0;
	}
	return ( readingFrames [orfNum++] );
}
char FrameIterator::outputCharacter = '.';
int FrameIterator::numDatabaseEntries = 0;
int FrameIterator::numSearches = 1;
bool FrameIterator::mpi = false;
GenElapsedTime FrameIterator::elapsedTime;
int FrameIterator::num = 1;
int FrameIterator::timeUpdate = 2;
int FrameIterator::timeout = InfoParams::instance ().getIntValue ( "timeout" );
string FrameIterator::searchJobID;
bool FrameIterator::updateProgress ( bool reset )
{
	bool ret = false;
	static int numSoFar = 0;
	numSoFar += getIncrement ();
	static int previousOutputPercent = 0;
	static int previousTimeSec = 0;
	int outputPercent;
	if ( numDatabaseEntries )
		outputPercent = int( (double)numSoFar * (double)100 / (double)numDatabaseEntries );
	else
		outputPercent = 100;
	bool percentChanged = outputPercent != previousOutputPercent;
	int timeSec = elapsedTime.getElapsedTime ();
	bool timeChanged = ( timeSec >= previousTimeSec + timeUpdate );
	if ( timeChanged ) {
		previousTimeSec = timeSec;
		previousOutputPercent = outputPercent;
	}
	if ( reset ) {
		numSoFar = 0;
		previousOutputPercent = 0;
		previousTimeSec = 0;
		//resetElapsedTime ( 1 );
	}
	if ( reset || ( percentChanged && timeChanged ) ) {
		static UpdatingJavascriptMessage ujm;
		ostringstream ostr;
		ostr << "Search ";
		if ( numSearches > 1 ) ostr << num << "/" << numSearches << " ";
		if ( reset ) {
			num++;
		}
		else {
			ostr << "Processing " << outputPercent << "%";
#ifdef MYSQL_DATABASE
			if ( !searchJobID.empty () ) {
				MySQLPPSDDBase::instance ().updatePercentComplete ( searchJobID, outputPercent );
			}
#endif
		}
		ostr << " completed.";
		double ratioComplete;
		if ( reset ) ratioComplete = (double)( num - 1 ) / (double)numSearches;
		else ratioComplete = ( (double)( num - 1 ) + ( outputPercent / 100.0 ) ) / (double)numSearches;
		int elapsedTimeSec = elapsedTime.getElapsedTime ();
		int remainingTime = (elapsedTimeSec / ratioComplete) - elapsedTimeSec;
		ostr << " " << elapsedTime.getElapsedTimeString () << " elapsed." << " " << elapsedTime.getTimeString ( remainingTime ) << " remaining.";
		ujm.writeMessage ( cout, ostr.str () );
		if ( timeout && remainingTime > timeout ) ret = true;
	}
	return ret;
}
int FrameIterator::getIncrement ()
{
	if ( numDatabaseEntries > 5000 )	return 1000;
	else if ( numDatabaseEntries > 500 )return 100;
	else if ( numDatabaseEntries > 50 )	return 10;
	else return 1;
}
void FrameIterator::searchFinished ()
{
#ifdef PP_MULTI
	if ( mpi ) {
		char c = 'y';
		MPI_Bsend ( &c, 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD );
	}
#endif
}
