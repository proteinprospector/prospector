/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_btag_mpi.cpp                                               *
*                                                                             *
*  Created    : September 20th 2004                                           *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2004-2012) The Regents of the University of California.         *
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
#include <lg_new.h>
#include <lgen_file.h>
#include <lgen_process.h>
#include <lp_frame.h>
#include <lu_param_list.h>
#include <lu_tag_srch.h>
#include <lu_html.h>
#include <lu_file_split.h>
#include <lu_proj_file.h>
#include <lu_expec_par.h>
#include <lu_btag_run.h>
#include <lu_getfil.h>
#include <lp_btag_mpi.h>
#ifdef MYSQL_DATABASE
#include <ld_init.h>
#endif
using std::ostringstream;
using std::string;
using std::cout;
using std::generate;
using std::remove;

static void abortProgram ( const string& message )
{
	ErrorHandler::genError ()->message ( message );
	MPI_Abort ( MPI_COMM_WORLD, -1 );
}
class MPIErrorHandler : public ErrorHandler {
	void messageDisplay ( const string& messageString );
	void errorDisplay ( const string& errorString, bool endProgram );
public:
	MPIErrorHandler ();
};
MPIErrorHandler::MPIErrorHandler ()
{
	ErrorHandler::registration ( this );
}
void MPIErrorHandler::errorDisplay ( const string& errorString, bool endProgram )
{
	char c = 'e';
	MPI_Bsend ( &c, 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD );
	char* message = gen_new_string ( errorString.c_str () );
	int len = errorString.length ();
	MPI_Bsend ( &len, 1, MPI_INT, 0, 0, MPI_COMM_WORLD );
	MPI_Bsend ( message, len + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD );	// Add 1 for the null terminator

	MPI_Ssend ( &len, 1, MPI_INT, 0, 0, MPI_COMM_WORLD ); // Dummy message to wait for abort
}
void MPIErrorHandler::messageDisplay ( const string& messageString )
{
	/* Don't print messages with MPI */
}

BatchTagMPI::BatchTagMPI ( int argc, char** argv ) :
	argc ( argc ),
	argv ( argv )
{
	MPI_Init ( &argc, &argv );

	int ntasks = 0;
	MPI_Comm_size ( MPI_COMM_WORLD, &ntasks);
	numSearches = ntasks - 1;

	MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
	genReduceProcessPriority ();
	if ( rank == 0 ) {
		if ( string ( argv [1] ) == "-k" )	initialiseProspectorBTag ( argv [2] );
		else								initialiseProspectorBTag ();	
	}
	else gen_set_new_handler ();
}
BatchTagMPI::~BatchTagMPI ()
{
}
void BatchTagMPI::run ()
{
	bool expectationSearchFirst = InfoParams::instance ().getBoolValue ( "expectation_search_first" );
	if ( rank == 0 ) {
		bool expectationSearchDone = false;
		ParameterList* searchParamList = new ParameterList ( argc, argv );
		startSerial = 1;
#ifdef MYSQL_DATABASE
		string searchKey = searchParamList->getStringValue ( "search_key" );
		if ( !searchKey.empty () ) {
			JobItem* jobItem = MySQLPPSDDBase::instance ().getSearchJobByKey ( searchKey );
			searchJobID = jobItem->getSearchJobID ();
			FrameIterator::setSearchJobID ( searchJobID );
			searchStage = jobItem->getSearchStage ();
			startSerial = jobItem->getStartSerial ();
			expectationSearchDone = expectationSearchFirst && ( searchStage == 2 );
		}
#endif
		if ( searchParamList->getStringValue ( "expect_calc_method" ) != "None" && !expectationSearchDone ) {
			string outputFilename;
			ParameterList* expParamList = getExpectationParams ( searchParamList, outputFilename, startSerial != 1 );
			if ( expParamList ) {
				paramList.push_back ( expParamList );
			}
			searchParamList->addOrReplaceName ( "expect_coeff_file", outputFilename );
		}
		paramList.push_back ( searchParamList );
		numRuns = paramList.size ();
	}
	MPI_Barrier ( MPI_COMM_WORLD );
	MPI_Bcast ( &numRuns, 1, MPI_INT, 0, MPI_COMM_WORLD );
	MPI_Bcast ( &searchStage, 1, MPI_INT, 0, MPI_COMM_WORLD );
	MPI_Bcast ( &startSerial, 1, MPI_INT, 0, MPI_COMM_WORLD );
	for ( int s = 0 ; s < numRuns ; s++ ) {
		int index = expectationSearchFirst ? s : numRuns - s - 1;
		char* message;
		int len;
		StringVector filesToDelete;
		if ( rank == 0 ) {
			string filename = getBatchTagOutputPath ( paramList [index] );
			for ( int i = 0 ; i < numSearches ; i++ ) {
				filesToDelete.push_back ( filename + string ( "_" ) + gen_itoa ( i ) );
			}
			init_html_premature_stop ( paramList [index]->getStringValue ( "report_title" ), true, filesToDelete );
			ostringstream ost;
			paramList [index]->copyToCGI ( ost );
			message = gen_new_string ( ost.str ().c_str () );
			len = ost.str ().length () + 1;
		}
		else {
			FrameIterator::setMPI ();
			static MPIErrorHandler e;
		}
		MPI_Barrier ( MPI_COMM_WORLD );
		MPI_Bcast ( &len, 1, MPI_INT, 0, MPI_COMM_WORLD );
		if ( rank != 0 ) {
			message = new char [len];
		}
		MPI_Bcast ( message, len, MPI_CHAR, 0, MPI_COMM_WORLD );	// Send parameters
		FileSplit* fs;
		if ( rank == 0 ) {
			ProjectFile pf ( paramList [index] );
			fs = new FileSplit ( pf.getNumMSMSSpectra (), numSearches, getMSMSMaxSpectra ( paramList [index] ) );
#ifdef MYSQL_DATABASE
			if ( !searchJobID.empty () ) MySQLPPSDDBase::instance ().setNumSerial ( searchJobID, fs->getNumSerial () );
#endif
			if ( fs->getTotalSpectra () == 0 ) {
				abortProgram ( "No spectra in the data file.\n" );
			}
		}
		MPI_Barrier ( MPI_COMM_WORLD );
		int startSerialActual = ( index == 0 ? startSerial : 1 );
		if ( rank == 0 ) {
#ifdef MYSQL_DATABASE
			if ( !searchJobID.empty () ) MySQLPPSDDBase::instance ().updateSearchStage ( searchJobID, numRuns == 1 ? 2 : index+1 );
#endif
			rankZeroSearchLoop ( fs, filesToDelete, startSerialActual );
		}
		else {
			otherRankSearchLoop ( message, startSerialActual );
		}
		MPI_Barrier ( MPI_COMM_WORLD );
		if ( rank == 0 ) joinResultsFiles ( paramList [index], numSearches );
	}
#ifdef MYSQL_DATABASE
	if ( rank == 0 && !searchJobID.empty () ) MySQLPPSDDBase::instance ().setJobDone ( searchJobID );
#endif
	MPI_Finalize ();
}
void BatchTagMPI::rankZeroSearchLoop ( const FileSplit* fs, const StringVector& filesToDelete, int startSerial )
{
	MPI_Status status;
	int tag = 0;
	int numDatabaseEntries;
	FrameIterator::setNumSearches ( fs->getNumSerial () );
	FrameIterator::resetElapsedTime ( startSerial );
	for ( int i = startSerial ; i <= fs->getNumSerial () ; i++ ) {
#ifdef MYSQL_DATABASE
		if ( !searchJobID.empty () ) MySQLPPSDDBase::instance ().updateSearchNumber ( searchJobID, i );
#endif
		int numActualSearches = genMin ( fs->getTotalSpectra (), numSearches );	// Deal with the case where there are less spectra than processes
		IntVector s ( numActualSearches );
		generate ( s.begin (), s.end (), GenUnique ( 1 ) );
		for ( int j = 0 ; ; j++ ) {		// This deals with the situation of multiple passes through the database
			for ( ; ; ) {
				char c;
				IntVector copyS = s;
				for ( int k = 0 ; k < s.size () ; k++ ) {
					int source = s [k];
					MPI_Recv ( &c, 1, MPI_CHAR, source, tag, MPI_COMM_WORLD, &status );
					if ( c == 'e' ) {
						int len;
						MPI_Recv ( &len, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status );
						char* message = new char [len+1];
						MPI_Recv ( message, len+1, MPI_CHAR, source, tag, MPI_COMM_WORLD, &status );	// Add 1 for the null terminator
						genUnlink ( filesToDelete );
						abortProgram ( message );
					}
					else if ( c == 'n' ) {
						MPI_Recv ( &numDatabaseEntries, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status );
						FrameIterator::setNumDatabaseEntries ( numDatabaseEntries );
					}
					if ( c == 'y' ) copyS.erase ( remove ( copyS.begin (), copyS.end (), source ), copyS.end () );
				}
				s = copyS;
				if ( s.empty () ) goto lab;
				if ( c != 'n' && c != 'x' ) {
					if ( FrameIterator::updateProgress ( false ) ) {
						genUnlink ( filesToDelete );
						abortProgram ( "The search has been abandoned as it is likely the timeout will be exceeded.\n" );
					}
					cout.flush ();
				}
				if ( c == 'x' ) {
					FrameIterator::updateProgress ( true );
					FrameIterator::decrementNum ();
					break;
				}
			}
		}
		lab:;
		FrameIterator::incrementNum ();
	}
}
void BatchTagMPI::otherRankSearchLoop ( const char* message, int startSerial )
{
	ParameterList pList ( message, false, false );
	int bufferSize = 10000;
	char* buffer = new char [bufferSize];
	MPI_Buffer_attach ( buffer, bufferSize );
	runBatchTag ( &pList, numSearches, rank - 1, "", startSerial );
	MPI_Buffer_detach ( &buffer, &bufferSize );
	delete [] buffer;
}
#endif
