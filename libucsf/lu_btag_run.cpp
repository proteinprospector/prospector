/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_btag_run.cpp                                               *
*                                                                             *
*  Created    : October 2nd 2007                                              *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2007-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef BATCHTAG
#include <lg_new.h>
#include <lg_string.h>
#include <lgen_error.h>
#include <lgen_file.h>
#include <lu_getfil.h>
#include <lu_param_list.h>
#include <lu_prog.h>
#include <lu_program.h>
#include <lu_proj_file.h>
#include <lu_file_split.h>
#include <lp_frame.h>
#ifdef MYSQL_DATABASE
#include <ld_init.h>
#endif
#include <lu_tag_par.h>
#include <lu_tag_srch.h>

using std::string;
using std::cout;
using std::endl;

class BTagRootErrorHandler : public ErrorHandler {
	string searchKey;
	void messageDisplay ( const string& messageString );
	void errorDisplay ( const string& errorString, bool endProgram );
public:
	BTagRootErrorHandler ( const string& searchKey );
};
BTagRootErrorHandler::BTagRootErrorHandler ( const string& searchKey ) :
	searchKey ( searchKey )
{
	ErrorHandler::registration ( this );
}
void BTagRootErrorHandler::errorDisplay ( const string& errorString, bool endProgram )
{
#ifdef MYSQL_DATABASE
	if ( !searchKey.empty () ) MySQLPPSDDBase::instance ().submitAbortedJob ( searchKey, errorString );
	else
#endif
		cout << "Error message: " << errorString << endl;
}
void BTagRootErrorHandler::messageDisplay ( const string& messageString )
{
#ifdef MYSQL_DATABASE
	if ( !searchKey.empty () ) MySQLPPSDDBase::instance ().submitAbortedJob ( searchKey, messageString );
	else
#endif
		cout << "Message: " << messageString << endl;
}
void initialiseProspectorBTag ( const string& searchKey )
{
	gen_set_new_handler ();

	static BTagRootErrorHandler h ( searchKey );
}
string getBatchTagOutputPath ( const ParameterList* params )
{
	string path;
	string outputFilename = params->getStringValue ( "output_filename" );
#ifdef MYSQL_DATABASE
	string searchKey = params->getStringValue ( "search_key" );
	if ( !searchKey.empty () ) {
		static BatchJobItem* bji = MySQLPPSDDBase::instance ().getBatchJobByKey ( searchKey );
		path = bji->getResultsPath ();
		if ( outputFilename.find ( ".exp." ) != string::npos ) {	// Random search
			path = bji->getProjectPath ();
		}
		else {														// Normal search
			outputFilename = searchKey;
		}
	}
	else {
#endif
		if ( outputFilename.find ( ".exp." ) != string::npos ) {	// Random search
			path = params->getStringValue ( "project_filepath" );
		}
		else {														// Normal search
			path = params->getStringValue ( "output_filepath" );
		}
#ifdef MYSQL_DATABASE
	}
#endif
	return adjustPPOutputPath ( path ) + outputFilename + ".xml";
}
int getMSMSMaxSpectra ( ParameterList* paramList )
{
	int maxSpectra = InfoParams::instance ().getIntValue ( "msms_max_spectra", 500 );
	string linkSearchType = paramList->getStringValue ( "link_search_type" );
	if ( linkSearchType != "No Link" && linkSearchType != "" ) {
		int maxSavedTagHits = paramList->getIntValue ( "max_saved_tag_hits" );
		int msmsMaxReportedHits = paramList->getIntValue ( "msms_max_reported_hits" );
		int ratio = maxSavedTagHits / msmsMaxReportedHits;
		maxSpectra /= ratio;
		maxSpectra = genMax ( maxSpectra, 4 );
	}
	return maxSpectra;
}
void runBatchTag ( ParameterList* paramList, int numSearches, int searchNumber, const string& searchJobID, int startSerial )
{
	string outputFilename = getBatchTagOutputPath ( paramList ) + string ( "_" ) + gen_itoa ( searchNumber );

	MSProgramParameters mpp ( paramList );
	MSProgram* ds = new MSProgram ( mpp );
	MSProgram::setParams ( paramList );
	ProjectFile pf ( paramList );
	FileSplit fs ( pf.getNumMSMSSpectra (), numSearches, getMSMSMaxSpectra ( paramList ) );
#ifdef MYSQL_DATABASE
	if ( searchJobID.empty () ) MySQLPPSDDBase::instance ( false, true );
#endif

	std::ios_base::openmode mode = std::ios_base::out;
	if ( startSerial != 1 ) mode |= std::ios_base::app;

	if ( searchNumber == 0 && startSerial == 1 ) {
		GenOFStream os ( outputFilename, std::ios_base::out );
		ds->printXMLTop ( os );
	}
	FrameIterator::resetElapsedTime ( startSerial );
	FrameIterator::setNumSearches ( fs.getNumSerial () );

#ifdef MYSQL_DATABASE
	if ( !searchJobID.empty () ) MySQLPPSDDBase::instance ().setNumSerial ( searchJobID, fs.getNumSerial () );
#endif
	int startSerialIndex = startSerial - 1;
	MSTagParameters* params = 0;
	for ( int i = startSerialIndex ; i < fs.getNumSerial () ; i++ ) {
#ifdef MYSQL_DATABASE
		if ( !searchJobID.empty () ) MySQLPPSDDBase::instance ().updateSearchNumber ( searchJobID, i+1 );
#endif
		IntVector iv = fs.getSendData ( searchNumber, i );
		if ( searchNumber == 0 || ( iv [2] != iv [0] ) || iv [3] >= iv [1] ) {	// This is to catch the case where the number of spectra is less than the number of processors
			FileSplit::setParams ( paramList, iv );
			if ( i == startSerialIndex )	params = new MSTagParameters ( paramList );
			else							params->setDataSetInfo ( paramList );
			TagSearch* ts = getTagSearch ( *params );
			GenOFStream os ( outputFilename, std::ios_base::out | std::ios_base::app );
			ts->printBodyXML ( os, searchNumber == 0 && i == 0 );	// Only show pre search results for one search
			delete ts;
		}
	}
	delete params;
	if ( searchNumber == genMin ( fs.getTotalSpectra (), numSearches ) - 1 ) {
		GenOFStream os ( outputFilename, std::ios_base::out | std::ios_base::app );
		ds->printXMLBottom ( os );
	}
}
void joinResultsFiles ( const ParameterList* paramList, int numSearches )
{
	bool joinResFiles = InfoParams::instance ().getBoolValue ( "join_results_files", true );
	string filename = getBatchTagOutputPath ( paramList );
	if ( numSearches > 1 && joinResFiles ) {
		GenOFStream os ( filename + string ( "_0" ), std::ios_base::out | std::ios_base::app | std::ios_base::binary );
		for ( int i = 1 ; i < numSearches ; i++ ) {
			string f = filename + string ( "_" ) + gen_itoa ( i );
			if ( genFileExists ( f ) ) {
				int size = (int) genFileSize ( f );
				GenIFStream is ( f, std::ios_base::out | std::ios_base::binary );
				char* info = new char [size];
				is.read ( (char*) info, size );
				os.write ( info, size );
				delete [] info;
				is.close ();
				genUnlink ( f );
			}
		}
	}
	genUnlink ( filename );
	if ( numSearches == 1 || joinResFiles ) {
		genRename ( ( filename + string ( "_0" ) ), filename );
	}
}
#endif
