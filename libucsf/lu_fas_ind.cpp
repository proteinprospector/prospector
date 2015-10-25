/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_fas_ind.cpp                                                *
*                                                                             *
*  Created    : July 20th 1996                                                *
*                                                                             *
*  Purpose    : Functions to create the index file for fasta formatted        *
*               databases.                                                    *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1996-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lg_string.h>
#include <lg_time.h>
#include <lgen_error.h>
#include <lgen_file.h>
#include <lu_fas_ind.h>
#include <lu_html.h>
using std::string;
using std::ios;
using std::ios_base;
using std::cout;
using std::endl;

DatabaseIndicies::DatabaseIndicies ( const string& filename, bool createFlag )
{
	string databasePath;
	string fullDatabasePath;
	if ( filename.find ( "UserProtein" ) != string::npos ) {
		databasePath = filename.substr ( 0, filename.length () - 6 );	// Cut off .fasta
		fullDatabasePath = filename;
	}
	else {
		databasePath = SeqdbDir::instance ().getSeqdbDir () + filename;
		fullDatabasePath = SeqdbDir::instance ().getDatabasePath ( filename );
	}
	fullIDCPath = databasePath + ".idc";
	fullIDPPath = databasePath + ".idp";
	fullIDIPath = databasePath + ".idi";

	if ( !genFileExists ( fullDatabasePath ) ) {
		ErrorHandler::genError ()->error ( "The database is not present.\n" );
	}
	databaseMap = new MMapFile <char> ( fullDatabasePath );
	if ( createFlag ) {
		writeIndexFile ();
	}
	readIndexFile ();
}
DatabaseIndicies::~DatabaseIndicies ()
{
	delete databaseMap;
	delete commentIndexMap;
	delete proteinIndexMap;
}
char* DatabaseIndicies::getCommentPointer ( unsigned int serialNumber, int* length )
{
	if ( serialNumber > numEntries ) {
		ErrorHandler::genError ()->error ( "The index number you have entered is greater than the number of entries in the database.\n" );
	}
	GENINT64 thisEntry = commentIndexMap->subscript(serialNumber - 1);
	GENINT64 nextEntry = proteinIndexMap->subscript(serialNumber - 1);
	*length = static_cast <int> ( nextEntry - thisEntry );

	return ( databaseMap->getRange ( thisEntry, nextEntry - 1 ) );
}
char* DatabaseIndicies::getProteinPointer ( unsigned int serialNumber, int* length )
{
	if ( serialNumber > numEntries ) {
		ErrorHandler::genError ()->error ( "The index number you have entered is greater than the number of entries in the database.\n" );
	}
	GENINT64 thisEntry = proteinIndexMap->subscript(serialNumber - 1);
	GENINT64 nextEntry = ( serialNumber == numEntries ) ? databaseMap->getFileSize () : commentIndexMap->subscript(serialNumber);
	*length = static_cast <int> ( nextEntry - thisEntry );

	return ( databaseMap->getRange ( thisEntry, nextEntry - 1 ) );
}
char* DatabaseIndicies::getEntryPointer ( unsigned int serialNumber, int* length )
{
	if ( serialNumber > numEntries ) {
		ErrorHandler::genError ()->error ( "The index number you have entered is greater than the number of entries in the database.\n" );
	}
	GENINT64 thisEntry = commentIndexMap->subscript(serialNumber - 1);
	GENINT64 nextEntry = ( serialNumber == numEntries ) ? databaseMap->getFileSize () : commentIndexMap->subscript(serialNumber);
	*length = static_cast <int> ( nextEntry - thisEntry );

	return ( databaseMap->getRange ( thisEntry, nextEntry - 1 ) );
}
void DatabaseIndicies::readIndexFile ()
{
	GenIFStream idiFile ( fullIDIPath, ios::binary );
	idiFile.read ( (char*) &numEntries, sizeof (unsigned int) );
	idiFile.read ( (char*) &maxCommentLength, sizeof (unsigned int) );
	idiFile.read ( (char*) &maxProteinLength, sizeof (unsigned int) );

	commentIndexMap = new MMapFile <GENINT64> ( fullIDCPath, 0, 0x80000 );
	proteinIndexMap = new MMapFile <GENINT64> ( fullIDPPath, 0, 0x80000 );
}
void DatabaseIndicies::writeIndexFile ()
{
	GenOFStream idcFile ( fullIDCPath, ios_base::binary );
	GenOFStream idpFile ( fullIDPPath, ios_base::binary );
	maxProteinLength = 0;
	maxCommentLength = 0;
	char previous = '\n';
	char* pointer = databaseMap->getStartPointer ();
	char* endPointer = databaseMap->getEndPointer ();
	numEntries = 0;
	GENINT64 commentOffset = 0;
	GENINT64 proteinOffset = 0;
	bool comment = false;
	UpdatingJavascriptMessage ujm;
	for ( ; ; ) {
		while ( pointer <= endPointer && previous != '\n' ) previous = *pointer++; 
		if ( pointer > endPointer ) {
			GENINT64 offset2 = databaseMap->getMapOffset ( pointer );
			if ( offset2 == databaseMap->getFileSize () ) break;
			pointer = databaseMap->nextRange ();
			endPointer = databaseMap->getEndPointer ();
		}
		if ( previous == '\n' ) {
			previous = *pointer;
			if ( previous == '>' ) {
				if ( comment == true ) {	// Blank entry no protein defined
					proteinOffset = databaseMap->getMapOffset ( pointer );
 					idpFile.write ( (char*) &proteinOffset, sizeof (GENINT64) );
					maxCommentLength = genMax( maxCommentLength, (unsigned int)(proteinOffset - commentOffset) );
					comment = false;
				}
				commentOffset = databaseMap->getMapOffset ( pointer );
				numEntries++;
				if ( numEntries % 100000 == 0 ) ujm.writeMessage ( cout, numEntries );
				idcFile.write ( (char*) &commentOffset, sizeof (GENINT64) );
				maxProteinLength = genMax( maxProteinLength, (unsigned int)(commentOffset - proteinOffset) );
				comment = true;
			}
			else {
				if ( comment == true ) {
					proteinOffset = databaseMap->getMapOffset ( pointer );
					idpFile.write ( (char*) &proteinOffset, sizeof (GENINT64) );
					maxCommentLength = genMax( maxCommentLength, (unsigned int)(proteinOffset - commentOffset) );
					comment = false;
				}
			}
			pointer++;
		}
	}
	ujm.deletePreviousMessage ( cout );

	maxProteinLength = genMax( maxProteinLength, (unsigned int)(databaseMap->getFileSize () - proteinOffset) );

	GenOFStream idiFile ( fullIDIPath, ios_base::binary );
	idiFile.write ( (char*) &numEntries, sizeof (unsigned int) );
	idiFile.write ( (char*) &maxCommentLength, sizeof (unsigned int) );
	idiFile.write ( (char*) &maxProteinLength, sizeof (unsigned int) );
}
string getIndexFileLastModifiedTime ( const string& filename )	
{
	return ( genTimeAndDateString ( genLastModifyTime ( SeqdbDir::instance ().getSeqdbDir () + filename + ".idi" ) ) );
}
