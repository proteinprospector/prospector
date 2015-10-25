/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_db_srch.cpp                                                *
*                                                                             *
*  Created    : Septmber 5th 2001                                             *
*                                                                             *
*  Purpose    : Database search and hit base classes.                         *
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
#include <lg_stdio.h>
#include <lgen_file.h>
#include <lu_acc_link.h>
#include <lu_check_db.h>
#include <lu_db_srch.h>
#include <lu_cgi_val.h>
#include <lu_srch_par.h>
#include <lu_param_list.h>
#include <lu_html.h>
#include <lu_table.h>
#include <lu_getfil.h>
#include <lu_fas_ind.h>
#include <lu_hit.h>
using std::ostream;
using std::string;
using std::stringstream;
using std::endl;
using std::make_pair;

DatabaseHits::~DatabaseHits () {}

StringVector DBSearch::tempDirs;
FastaServerPtrVector DBSearch::userFS;

DBSearch::DBSearch ( const MSSearchParameters& params ) :
	MSProgram ( params ),
	params ( params )
{
	StringVector db = params.getDatabase ();
	DBSearchFlags dbFlags ( db );
	if ( dbFlags.getRandomAndReverseFlag () ) {
		ErrorHandler::genError ()->error ( "Random and reverse database cannot be simultaneously selected.\n" );
	}
	for ( int i = 0 ; i < db.size () ; i++ ) {
		const string& d = db [i];
		if ( d == "User Protein" ) {
			PairStringBool psb = params.createTemporaryDatabase ( dbFlags.getRandomFlag (), dbFlags.getReverseFlag () );
			tempDirs.push_back ( genDirectoryFromPath ( psb.first ) );
			string fName = psb.first;
			fs.push_back ( new FastaServer ( fName, psb.second ) );
			userFS.push_back ( fs.back () );
			if ( dbFlags.getRandomFlag () ) {
				fs.push_back ( new FastaServer ( fName.substr ( 0, fName.length () - 6 ) + ".random.fasta", psb.second ) );
				userFS.push_back ( fs.back () );
			}
			if ( dbFlags.getReverseFlag () ) {
				fs.push_back ( new FastaServer ( fName.substr ( 0, fName.length () - 6 ) + ".reverse.fasta", psb.second ) );
				userFS.push_back ( fs.back () );
			}
		}
		else {
			if ( dbFlags.getConcatFlag () && !isSuffix ( d, ".concat" ) ) {
				ErrorHandler::genError ()->error ( "If you select a concat database then all selected database must be concat databases.\n" );
			}
			PairStringString pss;
			if ( getConcatDBPair ( d, pss ) ) {
				fs.push_back ( new FastaServer ( pss.first ) );
				fs.push_back ( new FastaServer ( pss.second ) );
			}
			else {
				fs.push_back ( new FastaServer ( d ) );
			}
		}
	}
	int maxNTermAA = params.getMaxNTermAA ();
	for ( int j = 0 ; j < fs.size () ; j++ ) {
		fs [j]->setMaxNTermAA ( maxNTermAA );
		dnaFrameTranslationPairVector.push_back ( getDNAFrameTranslationPair ( fs [j]->getFileName (), params.getDNAFrameTranslation () ) );
	}
	params.doSearch ( fs );
}
DBSearch::~DBSearch ()
{
	for ( int i = 0 ; i < fs.size () ; i++ ) {
		delete fs [i];
	}
	userFS.resize ( 0 );
	ProteinHit::reset ();
}
void DBSearch::deleteTempDirs ()
{
	for ( int i = 0 ; i < userFS.size () ; i++ ) {
		delete userFS [i];
	}
	for ( StringVectorSizeType j = 0 ; j < tempDirs.size () ; j++ ) {
		genUnlinkDirectory ( tempDirs [j] );
	}
}
void DBSearch::printBodyHTML ( ostream& os )
{
	params.printPreSearchHTML ( os );
	printHTMLHits ( os );
}
void DBSearch::printBodyTabDelimitedText ( ostream& os )
{
	printDelimitedHits ( os );
}
void DBSearch::printBodyXML ( ostream& os )
{
	params.printPreSearchBodyXML ( os );
	printXMLHits ( os );
}
DatabaseSearch::DatabaseSearch ( const MSSearchParameters& params ) :
	DBSearch ( params )
{
}
DatabaseSearch::~DatabaseSearch ()
{
}
void DatabaseSearch::saveHits () const
{
	if ( numHits ) {
		GenOFStream os ( params.getOutputPath ( ".res" ) );

		printCGIString ( os, "output_filename", params.getOutputFilename () );
		printCGI ( os, "num_indicies", numHits );

		for ( int i = 0 ; i < fs.size () ; i++ ) {
			string db = fs [i]->getFileName ();
			if ( db == "UserProtein.fasta" )			db = "User Protein";
			if ( db == "UserProtein.random.fasta" )		db = "User Protein Random";
			if ( db == "UserProtein.reverse.fasta" )	db = "User Protein Reverse";
			bool hit = false;
			int dbIdx = i+1;
			string dbIdxStr = ( dbIdx == 1 ) ? "" : gen_itoa ( dbIdx );
			int j;
			printCGIString ( os, "database" + dbIdxStr, db );
			if ( !isPrefix ( db, "User Protein" ) ) {
				printCGI ( os, "database_date" + dbIdxStr, ParameterList::spaceToPlus ( getIndexFileLastModifiedTime ( db ) ) );
			}
			os << "indicies" << dbIdxStr << "=";
			for ( j = 0 ; j < numHits ; j++ ) {
				if ( databaseHits->getDBIndex ( j ) == i ) {
					if ( hit ) os << "+";
					os << databaseHits->getIndex ( j );
					hit = true;
				}
			}
			os << "&";
			os << "accession_numbers" << dbIdxStr << "=";
			for ( j = 0 ; j < numHits ; j++ ) {
				if ( databaseHits->getDBIndex ( j ) == i ) {
					os << databaseHits->getAccessionNumber ( j );
					os << "%0D%0A";
				}
			}
			os << "&";
			hit = false;
			os << "hit_open_reading_frames" << dbIdxStr << "=";
			for ( j = 0 ; j < numHits ; j++ ) {
				if ( databaseHits->getDBIndex ( j ) == i ) {
					if ( hit ) os << "+";
					os << databaseHits->getOpenReadingFrame ( j );
					hit = true;
				}
			}
			os << "&";
			hit = false;
			os << "hit_dna_reading_frames" << dbIdxStr << "=";
			for ( j = 0 ; j < numHits ; j++ ) {
				if ( databaseHits->getDBIndex ( j ) == i ) {
					if ( hit ) os << "+";
					os << databaseHits->getDNAReadingFrame ( j );
					hit = true;
				}
			}
			if ( i != db.size () - 1 ) os << "&";
		}
	}
}
void DatabaseSearch::printHTMLHits ( ostream& os )
{
	printHTMLHitsJavascript ( os );

	printHTMLHitsReport ( os );
}
void DatabaseSearch::printHTMLHitsJavascript ( ostream& os ) const
{
	startJavascript ( os );
	StringVector dbPaths;
	for ( int i = 0 ; i < fs.size () ; i++ ) {
		string fPath = fs [i]->getFilePath ();
		if ( fPath.empty () ) fPath = fs [i]->getFileName ();
		dbPaths.push_back ( fPath );
		AccessionNumberLinkInfo anli;
		anli.printHTML ( os, fs [i]->getFileName () );
		MSDigestLink digLink ( params.getSearchName (), dbPaths.back () );
		digLink.printHTML ( os );
	}
	endJavascript ( os );
}
void DatabaseSearch::printHTMLHitsReport ( ostream& os )
{
	int maxReportedHits = params.getMaxReportedHits ();
	printNumHits ( os, params.getReportTitle (), numHits, maxReportedHits );
	numHits = genMin( numHits, maxReportedHits );

	if ( numHits == 0 ) {
		os << "<p />" << endl;
		os << params.getReportTitle () << " could not match any entries to the input data.<br />" << endl;
	}
	else
		printHTMLHitsTable ( os );
}
void DatabaseSearch::printHTMLHitsTable ( ostream& os ) const
{
	os << "<table cellspacing=\"3\">" << endl;
	for ( int i = 0 ; i < numHits ; i++ ) {
		if ( i == 0 ) {
			tableRowStart ( os );
				databaseHits->printHTMLHeader ( os, i );
			tableRowEnd ( os );
		}
		tableRowStart ( os );
			databaseHits->printHTMLHit ( os, i );
		tableRowEnd ( os );
	}
	os << "</table>" << endl;
}
void DatabaseSearch::printXMLHits ( ostream& os ) const
{
	os << "<" << params.getSearchName () << "_hits>" << endl;
	for ( int i = 0 ; i < numHits ; i++ ) {
		databaseHits->printXMLHit ( os, i );
	}
	os << "</" << params.getSearchName () << "_hits>" << endl;
}
