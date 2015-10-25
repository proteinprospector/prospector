/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_blib.cpp                                                   *
*                                                                             *
*  Created    : July 18th 2013                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2013-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <cmath>
#include <iostream>
#include <string>
#include <sqlite3.h>
#include <zlib.h>
#include <lg_io.h>
#include <lg_string.h>
#include <lg_time.h>
#include <lgen_define.h>
#include <lgen_error.h>
#include <lgen_file.h>
#include <lu_data.h>
#include <lu_blib.h>
using std::ostream;
using std::ostringstream;
using std::string;
using std::cout;
using std::endl;
using std::make_pair;
using std::vector;
using std::multimap;
using std::pair;

namespace {
unsigned long compressBound ( unsigned long sourceLen )
{
	return ( ( (int)ceil ( sourceLen * 1.001 ) ) + 12 );
}
}

Blib::Blib () :
	PPSQLite (),
	pM ( 0 ),
	pI ( 0 ),
	pSiz ( 0 ),
	compMZ ( 0 ),
	compMZLen ( 0 ),
	compInt ( 0 ),
	compIntLen ( 0 )
{
}
Blib::~Blib ()
{
	if ( pSiz ) {
		delete [] pM;
		delete [] pI;
	}
	if ( compMZLen ) delete [] compMZ;
	if ( compIntLen ) delete [] compInt;
}
void Blib::updateQuery ( const string& table, const PairStringString& setClause, const PairStringString& whereClause )
{
	string sql = "UPDATE " + table + " ";
	sql += "SET " + setClause.first + " = ";
	const string& sec = setClause.second;
	if ( sec == "NULL" || isSuffix ( sec, "()" ) )
		sql += sec;
	else
		sql += "'" + sec + "'";
	sql += " WHERE " + whereClause.first + " = '" + whereClause.second + "'";
	rc = sqlite3_exec ( db, sql.c_str (), callback, 0, &zErrMsg );
}
void Blib::resizePMPI ( int numPeaks )
{
	if ( numPeaks > pSiz ) {
		delete [] pM;
		delete [] pI;
		pM = new double [numPeaks];
		pI = new float [numPeaks];
		pSiz = numPeaks;
	}
}
void Blib::resizeCompMZ ( unsigned long compMZSiz )
{
	if ( compMZSiz > compMZLen ) {
		delete [] compMZ;
		compMZ = new unsigned char [compMZSiz];
		compMZLen = compMZSiz;
	}
}
void Blib::resizeCompInt ( unsigned long compIntSiz )
{
	if ( compIntSiz > compIntLen ) {
		delete [] compInt;
		compInt = new unsigned char [compIntSiz];
		compIntLen = compIntSiz;
	}
}
BlibWrite::BlibWrite ( const string& name, bool append ) :
	Blib ()
{
	if ( !append ) {		// Get rid of the old database if creating a new one
		genUnlink ( name );
	}
	rc = sqlite3_open_v2 ( name.c_str (), &db, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, NULL );
}
BlibWrite::~BlibWrite ()
{
}

void BlibWrite::create ()
{
	StringVector sql;
	sql.push_back ( "CREATE TABLE LibInfo(libLSID TEXT, createTime TEXT, numSpecs INTEGER, majorVersion INTEGER, minorVersion INTEGER)" );
	sql.push_back ( "CREATE TABLE RefSpectra (id INTEGER primary key autoincrement not null, peptideSeq VARCHAR(150), precursorMZ REAL, precursorCharge INTEGER, peptideModSeq VARCHAR(200), prevAA CHAR(1), nextAA CHAR(1), copies INTEGER, numPeaks INTEGER, retentionTime REAL, fileID INTEGER, SpecIDinFile VARCHAR(256), score REAL, scoreType TINYINT)" );
	sql.push_back ( "CREATE TABLE sqlite_sequence(name,seq)" );
	sql.push_back ( "CREATE TABLE Modifications (id INTEGER primary key autoincrement not null,RefSpectraID INTEGER, position INTEGER, mass REAL)" );
	sql.push_back ( "CREATE TABLE RefSpectraPeaks(RefSpectraID INTEGER, peakMZ BLOB, peakIntensity BLOB)" );
	sql.push_back ( "CREATE TABLE SpectrumSourceFiles (id INTEGER PRIMARY KEY autoincrement not null,fileName VARCHAR(512) )" );
	sql.push_back ( "CREATE TABLE ScoreTypes (id INTEGER PRIMARY KEY, scoreType VARCHAR(128) )" );
	sql.push_back ( "CREATE TABLE RetentionTimes (RefSpectraID INTEGER, RedundantRefSpectraID INTEGER, SpectrumSourceID INTEGER, retentionTime REAL, bestSpectrum INTEGER, FOREIGN KEY(RefSpectraID) REFERENCES RefSpectra(id) )" );
	sql.push_back ( "CREATE INDEX idxPeptide ON RefSpectra (peptideSeq, precursorCharge)" );
	sql.push_back ( "CREATE INDEX idxPeptideMod ON RefSpectra (peptideModSeq, precursorCharge)" );
	sql.push_back ( "CREATE INDEX idxRefIdPeaks ON RefSpectraPeaks (RefSpectraID)" );
	for ( StringVectorSizeType i = 0 ; i < sql.size () ; i++ ) {
		rc = sqlite3_exec ( db, sql [i].c_str (), callback, 0, &zErrMsg );
	}
}
void BlibWrite::insertLibInfo ( int numSpecs, const string& server, bool redundant, const string& name )
{
	StringVector names;
	StringVector values;

	string sRedundant = "nr";
	if ( redundant ) sRedundant = "redundant";

	names.push_back ( "libLSID" );

	string libLSID;
	libLSID += "urn";				// Stands for Uniform Resource Name
	libLSID += ":";
	libLSID += "lsid";				// Stands for Life Science Identifier
	libLSID += ":";
	libLSID += server;				// This is the lab authority
	libLSID += ":";
	libLSID += "spectral_library";	// This is the namespace
	libLSID += ":";
	libLSID += "bibliospec";
	libLSID += ":";
	libLSID += sRedundant;
	libLSID += ":";
	libLSID += name;

	values.push_back ( libLSID );

	names.push_back ( "createTime" );
	string date = genCurrentTimeAndDateString ();
	values.push_back ( date );

	names.push_back ( "numSpecs" );
	string sNumSpecs = gen_itoa ( numSpecs );
	values.push_back ( sNumSpecs );

	names.push_back ( "majorVersion" );
	values.push_back ( "1" );

	names.push_back ( "minorVersion" );
	values.push_back ( "1" );

	insertQuery ( "LibInfo", names, values );
}
void BlibWrite::insertModifications ( int refSpectraID, const BlibRefSpectraEntryValue& val )
{
	const VectorPairIntDouble& vpid = val.getVpid ();
	if ( !vpid.empty () ) insertModifications ( refSpectraID, vpid );
}
void BlibWrite::insertModifications ( int refSpectraID, const VectorPairIntDouble& vpid )
{
	for ( int i = 0 ; i < vpid.size () ; i++ ) {
		StringVector names;
		StringVector values;

		names.push_back ( "RefSpectraID" );
		string sRefSpectraID = gen_itoa ( refSpectraID );
		values.push_back ( sRefSpectraID );

		names.push_back ( "position" );
		string sPosition = gen_itoa ( vpid [i].first );
		values.push_back ( sPosition );

		names.push_back ( "mass" );
		string sMass = gen_ftoa ( vpid [i].second, "%.4f" );
		values.push_back ( sMass );

		insertQuery ( "Modifications", names, values );
	}
}
void BlibWrite::insertRefSpectra ( const PairStringInt& key, const BlibRefSpectraEntryValue& val )
{
	StringVector names;
	StringVector values;

	names.push_back ( "peptideSeq" );
	values.push_back ( val.getPeptideSeq () );

	names.push_back ( "precursorMZ" );
	string sPrecursorMZ = gen_ftoa ( val.getPrecursorMZ (), "%.4f" );
	values.push_back ( sPrecursorMZ );

	names.push_back ( "precursorCharge" );
	string sPrecursorCharge = gen_itoa ( key.second );
	values.push_back ( sPrecursorCharge );

	names.push_back ( "peptideModSeq" );
	values.push_back ( key.first );

	names.push_back ( "prevAA" );
	values.push_back ( val.getPrevAA () );

	names.push_back ( "nextAA" );
	values.push_back ( val.getNextAA () );

	names.push_back ( "copies" );
	string sCopies = gen_itoa ( val.getCopies () );
	values.push_back ( sCopies );

	names.push_back ( "numPeaks" );
	string sNumPeaks = gen_itoa ( val.getNumPeaks () );
	values.push_back ( sNumPeaks );

	names.push_back ( "retentionTime" );
	string sRetentionTime = gen_ftoa ( val.getRetentionTime (), "%.4f" );
	values.push_back ( sRetentionTime );

	names.push_back ( "fileID" );
	string sFileID = gen_itoa ( val.getFileID () );
	values.push_back ( sFileID );

	names.push_back ( "SpecIDinFile" );
	values.push_back ( val.getSpecIDinFile () );

	names.push_back ( "score" );
	ostringstream ostr;
	genPrintSigFig ( ostr, val.getEValue (), 2 );
	values.push_back ( ostr.str () );

	names.push_back ( "scoreType" );
	string sScoreType = gen_itoa ( val.getScoreType () );
	values.push_back ( sScoreType );

	insertQuery ( "RefSpectra", names, values );
}
void BlibWrite::populatePMPI ( const DataFilePeakVector& dfpv )
{
	int numPeaks = dfpv.size ();
	resizePMPI ( numPeaks );
	for ( int i = 0 ; i < numPeaks ; i++ ) {
		pM [i] = dfpv [i].getMOverZ ();
		pI [i] = dfpv [i].getIntensity ();
	}
}
void BlibWrite::insertRefSpectraPeaks ( int spectraID, const DataFilePeakVector& dfpv )
{
	populatePMPI ( dfpv );
	int numPeaks = dfpv.size ();
	unsigned long sizeMZ = (unsigned long) numPeaks * sizeof(double);			// compress mz
	unsigned long compMZSiz = compressBound ( sizeMZ );
	resizeCompMZ ( compMZSiz );
	unsigned char* compMZLocal = compMZ;
	unsigned long compMZLocalLen = compMZSiz;
	compress ( compMZLocal, &compMZLocalLen, (const unsigned char*)pM, sizeMZ );
	if ( compMZLocalLen >= sizeMZ ) {
		compMZLocal = (unsigned char*)pM;
		compMZLocalLen = sizeMZ;
	}

	unsigned long sizeI = (unsigned long) numPeaks * sizeof(float);				// compress intensity
	unsigned long compIntSiz = compressBound ( sizeI );
	resizeCompInt ( compIntSiz );
	unsigned char* compIntLocal = compInt;
	unsigned long compIntLocalLen = compIntSiz;
	compress ( compIntLocal, &compIntLocalLen, (const unsigned char*)pI, sizeI );
	if ( compIntLocalLen >= sizeI ) {
		compIntLocal = (unsigned char*)pI;
		compIntLocalLen = sizeI;
	}
	updateRefSpectraNumPeaks ( spectraID, numPeaks );
	insertRefSpectraPeaks ( spectraID, compMZLocal, (int)compMZLocalLen, compIntLocal, (int)compIntLocalLen );
}
void BlibWrite::updateRefSpectraNumPeaks ( int spectraID, int numPeaks )
{
	PairStringString setClause = make_pair ( string("numPeaks"), gen_itoa ( numPeaks ) );
	PairStringString whereClause = make_pair ( string("id"), gen_itoa ( spectraID ) );
	updateQuery ( "RefSpectra", setClause, whereClause );
}
void BlibWrite::insertRefSpectraPeaks ( int spectraID, unsigned char* mzComp, int mzCompSiz, unsigned char* intComp, int intCompSiz )
{
	string sql;
	sql += "INSERT INTO RefSpectraPeaks VALUES(";
	sql += gen_itoa ( spectraID );
	sql += ",?,?)";
	
	sqlite3_stmt* pStmt;
	rc = sqlite3_prepare_v2 ( db, sql.c_str (), -1, &pStmt, 0 );
/*
To execute an SQL query, it must first be compiled into a byte-code program.

The first argument, "db", is a database connection obtained from a prior successful call to sqlite3_open().
The database connection must not have been closed.

The second argument, "sql", is the statement to be compiled, encoded as UTF-8.

If the nByte argument is less than zero, then sql is read up to the first zero terminator. If nByte
is non-negative, then it is the maximum number of bytes read from zSql. When nByte is non-negative,
the zSql string ends at either the first '\000' or the nByte-th byte, whichever comes first.
If the caller knows that the supplied string is nul-terminated, then there is a small performance
advantage to be gained by passing an nByte parameter that is equal to the number of bytes in the input
string including the nul-terminator bytes as this saves SQLite from having to make a copy of the input
string.

If pzTail is not NULL then *pzTail is made to point to the first byte past the end of the first SQL
statement in zSql. These routines only compile the first statement in zSql, so *pzTail is left
pointing to what remains uncompiled.

*pStmt is left pointing to a compiled prepared statement that can be executed using sqlite3_step().
If there is an error, *pStmt is set to NULL. If the input text contains no SQL (if the input is an
empty string or a comment) then *pStmt is set to NULL. The calling procedure is responsible for
deleting the compiled SQL statement using sqlite3_finalize() after it has finished with it. pStmt
may not be NULL.

On success, the sqlite3_prepare() family of routines return SQLITE_OK; otherwise an error code is
returned.
*/
	if ( rc != SQLITE_OK ) {
		error ( "sqlite3_prepare_v2", sql, rc );
	}
	rc = sqlite3_bind_blob ( pStmt, 1, mzComp, mzCompSiz, SQLITE_STATIC );
	if ( rc != SQLITE_OK ) {
		error ( "sqlite3_bind_blob", sql, rc );
	}
	rc = sqlite3_bind_blob ( pStmt, 2, intComp, intCompSiz, SQLITE_STATIC );
	if ( rc != SQLITE_OK ) {
		error ( "sqlite3_bind_blob", sql, rc );
	}
/*
The sqlite3_bind_* routines return SQLITE_OK on success or an error code if anything goes wrong.
SQLITE_RANGE is returned if the parameter index is out of range. SQLITE_NOMEM is returned if
malloc() fails.
*/
	for ( int i = 0 ; ; i++ ) {
		rc = sqlite3_step ( pStmt );
		if ( rc == SQLITE_DONE ) break;
		if ( i > 10 ) {
			exit ( 0 );
		}
	}
/*
After a prepared statement has been prepared using either sqlite3_prepare_v2(), this function
must be called one or more times to evaluate the statement.

In the legacy interface, the return value will be either SQLITE_BUSY, SQLITE_DONE, SQLITE_ROW,
SQLITE_ERROR, or SQLITE_MISUSE. With the "v2" interface, any of the other result codes or
extended result codes might be returned as well.

SQLITE_BUSY means that the database engine was unable to acquire the database locks it needs
to do its job. If the statement is a COMMIT or occurs outside of an explicit transaction,
then you can retry the statement. If the statement is not a COMMIT and occurs within an
explicit transaction then you should rollback the transaction before continuing.

SQLITE_DONE means that the statement has finished executing successfully. sqlite3_step()
should not be called again on this virtual machine without first calling sqlite3_reset()
to reset the virtual machine back to its initial state.

*/
	if ( rc != SQLITE_OK && rc != SQLITE_DONE ) {
		error ( "sqlite3_step", sql, rc );
	}
	sqlite3_finalize ( pStmt );
/*
The sqlite3_finalize() function is called to delete a prepared statement. If the most
recent evaluation of the statement encountered no errors or if the statement is never
been evaluated, then sqlite3_finalize() returns SQLITE_OK. If the most recent evaluation
of statement S failed, then sqlite3_finalize(S) returns the appropriate error code or
extended error code.

The sqlite3_finalize(S) routine can be called at any point during the life cycle of
prepared statement S: before statement S is ever evaluated, after one or more calls
to sqlite3_reset(), or after any call to sqlite3_step() regardless of whether or not
the statement has completed execution.

Invoking sqlite3_finalize() on a NULL pointer is a harmless no-op.

The application must finalize every prepared statement in order to avoid resource leaks.
It is a grievous error for the application to try to use a prepared statement after it
has been finalized. Any use of a prepared statement after it has been finalized can
result in undefined and undesirable behavior such as segfaults and heap corruption.
*/
}
void BlibWrite::insertRetentionTimes ( int refSpectraID, const BlibRefSpectraEntryValue& val )
{
	IntVector fileIDs = val.getFileIDs ();
	DoubleVector rts = val.getRTs ();
	IntVector indicies = val.getIndicies ();
	int bestSpectrumIdx = val.getBestSpectrumIdx ();
	for ( DoubleVectorSizeType i = 0 ; i < rts.size () ; i++ ) {
		insertRetentionTimes ( refSpectraID, indicies [i], fileIDs [i], rts [i], ( i == bestSpectrumIdx ) ? 1 : 0 );
	}
}
void BlibWrite::insertRetentionTimes ( int refSpectraID, int redundantRefSpectraID, int spectrumSourceID, double retentionTime, int bestSpectrum )
{
	StringVector names;
	StringVector values;

	names.push_back ( "RefSpectraID" );
	string sRefSpectraID = gen_itoa ( refSpectraID );
	values.push_back ( sRefSpectraID );

	names.push_back ( "RedundantRefSpectraID" );
	string sRedundantRefSpectraID = gen_itoa ( redundantRefSpectraID );
	values.push_back ( sRedundantRefSpectraID );

	names.push_back ( "SpectrumSourceID" );
	string sSpectrumSourceID = gen_itoa ( spectrumSourceID );
	values.push_back ( sSpectrumSourceID );

	names.push_back ( "retentionTime" );
	string sRetentionTime = gen_ftoa ( retentionTime, "%.4f" );
	values.push_back ( sRetentionTime );

	names.push_back ( "bestSpectrum" );
	string sBestSpectrum = gen_itoa ( bestSpectrum );
	values.push_back ( sBestSpectrum );

	insertQuery ( "RetentionTimes", names, values );
}
void BlibWrite::insertScoreTypes ()
{
	StringVector names;
	names.push_back ( "scoreType" );

	StringVector v;
	v.push_back ( "UNKNOWN" );
	v.push_back ( "PERCOLATOR QVALUE" );
	v.push_back ( "PEPTIDE PROPHET SOMETHING" );
	v.push_back ( "SPECTRUM MILL" );
	v.push_back ( "IDPICKER FDR" );
	v.push_back ( "MASCOT IONS SCORE" );
	v.push_back ( "TANDEM EXPECTATION VALUE" );
	v.push_back ( "PROTEIN PILOT CONFIDENCE" );
	v.push_back ( "SCAFFOLD SOMETHING" );
	v.push_back ( "WATERS MSE PEPTIDE SCORE" );
	v.push_back ( "OMSSA EXPECTATION SCORE" );
	v.push_back ( "PROTEIN PROSPECTOR EXPECTATION SCORE" );
	v.push_back ( "SEQUEST XCORR" );
	v.push_back ( "MAXQUANT SCORE" );

	for ( StringVectorSizeType i = 0 ; i < v.size () ; i++ ) {
		StringVector values;
		values.push_back ( v [i] );
		insertQuery ( "ScoreTypes", names, values );
	}
}
void BlibWrite::insertSpectrumSourceFiles ( const string& file )
{
	StringVector names;
	names.push_back ( "fileName" );
	StringVector values;
	values.push_back ( file );

	insertQuery ( "SpectrumSourceFiles", names, values );
}
BlibRefSpectraEntryValue::BlibRefSpectraEntryValue ( const string& peptideSeq, const VectorPairIntDouble& vpid, double precursorMZ, const string& prevAA, const string& nextAA, int copies, int numPeaks, int fileID, double retentionTime, int index, const string& specIDinFile, const string& specID, double pValue, double eValue ) :
	peptideSeq ( peptideSeq ),
	vpid ( vpid ),
	precursorMZ ( precursorMZ ),
	prevAA ( prevAA ),
	nextAA ( nextAA ),
	copies ( copies ),
	numPeaks ( numPeaks ),
	fileID ( fileID ),
	retentionTime ( retentionTime ),
	specIDinFile ( specIDinFile ),
	specID ( specID ),
	pValue ( pValue ),
	eValue ( eValue ),
	scoreType ( 12 )
{
	fileIDs.push_back ( fileID );
	rts.push_back ( retentionTime );
	indicies.push_back ( index );
	bestSpectrumIdx = 0;
}
void BlibRefSpectraEntryValue::update ( double precursorMZ, const string& prevAA, const string& nextAA, int numPeaks, int fileID, double retentionTime, int index, const string& specIDinFile, const string& specID, double pValue, double eValue )
{
	this->precursorMZ = precursorMZ;
	this->prevAA = prevAA;
	this->nextAA = nextAA;
	this->numPeaks = numPeaks;
	this->fileID = fileID;
	this->retentionTime = retentionTime;
	this->specIDinFile = specIDinFile;
	this->specID = specID;
	this->pValue = pValue;
	this->eValue = eValue;

	bestSpectrumIdx = rts.size ();
	fileIDs.push_back ( fileID );
	rts.push_back ( retentionTime );
	indicies.push_back ( index );
}
void BlibRefSpectraEntryValue::update ( int fileID, double retentionTime, int index )
{
	fileIDs.push_back ( fileID );
	rts.push_back ( retentionTime );
	indicies.push_back ( index );
}
namespace {

StringVector peptideSeq;
StringVector modifications;
DoubleVector precursorMZ;
IntVector precursorCharge;
StringVector prevAA;
StringVector nextAA;
IntVector copies;
IntVector numPeaks;
StringVector retentionTime;
IntVector fileID;
StringVector specIDInFile;
DoubleVector score;

int refSpectraCallback ( void* NotUsed, int argc, char** argv, char** azColName )
{
	peptideSeq.push_back		( argv [0] );
	precursorMZ.push_back		( atof ( argv [1] ) );
	precursorCharge.push_back	( atoi ( argv [2] ) );
	prevAA.push_back			( argv [3] );
	nextAA.push_back			( argv [4] );
	copies.push_back			( atoi ( argv [5] ) );
	numPeaks.push_back			( atoi ( argv [6] ) );
	retentionTime.push_back		( argv [7] );
	fileID.push_back			( atoi ( argv [8] ) );
	specIDInFile.push_back		( argv [9] == 0 ? "" : argv [9] );
	score.push_back				( atof ( argv [10] ) );
	return 0;
}

IntVector localRefSpectraID;
IntVector localPosition;
DoubleVector localMass;

int modificationsCallback ( void* NotUsed, int argc, char** argv, char** azColName )
{
	localRefSpectraID.push_back	( atoi ( argv [0] ) );
	localPosition.push_back		( atoi ( argv [1] ) );
	localMass.push_back			( atof ( argv [2] ) );
	return 0;
}

StringVector fileName;

int spectrumSourceFilesCallback ( void* NotUsed, int argc, char** argv, char** azColName )
{
	fileName.push_back	( argv [0] );
	return 0;
}

}


BlibRead::BlibRead ( const string& name ) :
	Blib ()
{
	rc = sqlite3_open_v2 ( name.c_str (), &db, SQLITE_OPEN_READONLY, NULL );
}
BlibRead::~BlibRead ()
{
}
void BlibRead::read ( const string& peakListPath, StringVectorVector& header, StringVectorVector& rows )
{
	genCreateDirectory ( peakListPath );
	readResults ( header, rows );
	readSpectra ( peakListPath );
}
void BlibRead::readResults ( StringVectorVector& header, StringVectorVector& rows )
{
	string sql = "SELECT \
		peptideSeq, precursorMZ, precursorCharge, prevAA, nextAA, copies, numPeaks, retentionTime, fileID, SpecIDInFile, score \
		from RefSpectra";
	executeSQL ( sql, refSpectraCallback );

	sql = "SELECT \
		RefSpectraID, position, mass \
		from Modifications";
	executeSQL ( sql, modificationsCallback );

	bool modsInString = false;
	int lastRsidx = -1;
	if ( modsInString ) {
		int offset = 0;
		for ( IntVectorSizeType i = 0 ; i < localRefSpectraID.size () ; i++ ) {
			int rsidx = localRefSpectraID [i]-1;
			if ( rsidx != lastRsidx ) {
				lastRsidx = rsidx;
				offset = 0;
			}
			int pos = localPosition [i] + offset;
			string massStr = '(' + gen_ftoa ( localMass [i], "%.4f" ) + ')';
			int lenMassStr = massStr.length ();
			offset += lenMassStr;
			string& ps = peptideSeq [rsidx];
			ps = ps.substr ( 0, pos ) + massStr + ps.substr ( pos );
		}
	}
	else {
		modifications.resize ( peptideSeq.size () );
		for ( IntVectorSizeType i = 0 ; i < localRefSpectraID.size () ; i++ ) {
			int rsidx = localRefSpectraID [i]-1;
			if ( rsidx != lastRsidx ) {
				lastRsidx = rsidx;
			}
			string& m = modifications [rsidx];
			m += gen_ftoa ( localMass [i], "%.4f" ) + "@" + gen_itoa ( localPosition [i] ) + ";";
		}
	}

	sql = "SELECT fileName from SpectrumSourceFiles";
	executeSQL ( sql, spectrumSourceFilesCallback );
	for ( StringVectorSizeType j = 0 ; j < fileName.size () ; j++ ) {
		fileName [j] = genShortFilenameFromPath ( fileName [j] );
	}
	StringVector hCols;
	if ( modsInString ) {
		hCols.push_back ( "Peptide Sequence" );
	}
	else {
		hCols.push_back ( "DB Peptide" );
		hCols.push_back ( "Modifications" );
	}
	hCols.push_back ( "Precursor MZ" );
	hCols.push_back ( "Precursor Charge" );
	hCols.push_back ( "Prev AA" );
	hCols.push_back ( "Next AA" );
	hCols.push_back ( "Copies" );
	hCols.push_back ( "Num Peaks" );
	hCols.push_back ( "RT" );
	hCols.push_back ( "Fraction" );
	hCols.push_back ( "Spec ID in File" );
	hCols.push_back ( "Score" );
	header.push_back ( hCols );
	for ( StringVectorSizeType k = 0 ; k < peptideSeq.size () ; k++ ) {
		StringVector cols;
		cols.push_back ( peptideSeq [k] );
		if ( !modsInString ) {
			cols.push_back ( genStrtrimSemiColon ( modifications [k] ) );
		}
		cols.push_back ( gen_ftoa ( precursorMZ [k], "%.4f" ) );
		cols.push_back ( gen_itoa ( precursorCharge [k] ) );
		cols.push_back ( prevAA [k] );
		cols.push_back ( nextAA [k] );
		cols.push_back ( gen_itoa ( copies [k] ) );
		cols.push_back ( gen_itoa ( numPeaks [k] ) );
		cols.push_back ( retentionTime [k] );
		cols.push_back ( fileName [fileID [k]-1] );
		cols.push_back ( specIDInFile [k] );
		ostringstream ostr;
		genPrintSigFig ( ostr, score [k], 2 );
		cols.push_back ( ostr.str () );
		rows.push_back ( cols );
	}
}
/*
void BlibRead::readSpectra ( const string& path )
{
	int curFileID = 1;
	peakListFractionNames.push_back ( fileName [curFileID-1] );
	peakListCentroidFileNames.push_back ( fileName [curFileID-1] + ".mgf" );
	peakListCentroidFilePaths.push_back ( path + SLASH + peakListCentroidFileNames.back () );
	GenOFStream os ( peakListCentroidFilePaths.back () );

	string sql = "SELECT RefSpectraID, peakMZ, peakIntensity FROM RefSpectraPeaks";	// LIMIT 50 OFFSET 400
	
	sqlite3_stmt* pStmt;
	rc = sqlite3_prepare_v2 ( db, sql.c_str (), -1, &pStmt, 0 );

	if ( rc != SQLITE_OK ) {
		error ( "sqlite3_prepare_v2", sql, rc );
	}
	while ( sqlite3_step ( pStmt ) == SQLITE_ROW ) {
		int rsi = readSpectrumRow ( pStmt );
		int curNumPeaks = numPeaks [rsi];
		if ( fileID [rsi] != curFileID ) {
			os.close ();
			curFileID = fileID [rsi];
			peakListFractionNames.push_back ( fileName [curFileID-1] );
			peakListCentroidFileNames.push_back ( fileName [curFileID-1] + ".mgf" );
			peakListCentroidFilePaths.push_back ( path + SLASH + peakListCentroidFileNames.back () );
			os.open ( peakListCentroidFilePaths.back () );
		}
		writeSpectrum ( os, rsi, curNumPeaks );
	}
	sqlite3_finalize ( pStmt );
	os.close ();
}
*/
void BlibRead::readSpectra ( const string& path )
{
	vector <MultiMapNumberStringToInt> vmsi;
	vmsi.resize ( fileName.size () );
	for ( int a = 0 ; a < retentionTime.size () ; a++ ) {
		vmsi [fileID [a]-1].insert ( PairStringInt ( retentionTime[a], a+1 ) );
	}
	for ( int i = 0 ; i < fileName.size () ; i++ ) {
		peakListFractionNames.push_back ( fileName [i] );
		peakListCentroidFileNames.push_back ( fileName [i] + ".mgf" );
		peakListCentroidFilePaths.push_back ( path + SLASH + peakListCentroidFileNames.back () );
		ErrorHandler::genError ()->message ( "Processing fraction " + fileName [i] + ".\n" );
		GenOFStream os ( peakListCentroidFilePaths.back () );
		MultiMapNumberStringToInt& msi = vmsi [i];
		for ( MultiMapNumberStringToIntConstIterator j = msi.begin () ; j != msi.end () ; j++ ) {
			int rsi = (*j).second;
			string sql = "SELECT RefSpectraID, peakMZ, peakIntensity FROM RefSpectraPeaks";
			sql += " WHERE RefSpectraID = '" + gen_itoa ( rsi ) + "'";

			sqlite3_stmt* pStmt;
			rc = sqlite3_prepare_v2 ( db, sql.c_str (), -1, &pStmt, 0 );

			if ( rc != SQLITE_OK ) {
				error ( "sqlite3_prepare_v2", sql, rc );
			}
			while ( sqlite3_step ( pStmt ) == SQLITE_ROW ) {
				int rsi = readSpectrumRow ( pStmt );
				int curNumPeaks = numPeaks [rsi];
				writeSpectrum ( os, rsi, curNumPeaks );
			}
			sqlite3_finalize ( pStmt );
		}
		os.close ();
	}
}
int BlibRead::readSpectrumRow ( sqlite3_stmt* pStmt )
{
	int rsi = sqlite3_column_int ( pStmt, 0 )-1;
	int curNumPeaks = numPeaks [rsi];
	resizePMPI ( curNumPeaks );

	unsigned char* compMZ = (unsigned char*) sqlite3_column_blob ( pStmt, 1 );
	int len = sqlite3_column_bytes ( pStmt, 1 );
	unsigned long mzBufferLen = curNumPeaks * sizeof (double);
	if ( len >= mzBufferLen )
		memcpy ( pM, compMZ, mzBufferLen );
	else
		uncompress ( (unsigned char*)&pM[0], &mzBufferLen, compMZ, len );

	unsigned char* compInt = (unsigned char*) sqlite3_column_blob ( pStmt, 2 );
	int len2 = sqlite3_column_bytes ( pStmt, 2 );
	unsigned long intBufferLen = curNumPeaks * sizeof (float);
	if ( len2 >= intBufferLen )
		memcpy ( pI, compInt, intBufferLen );
	else
		uncompress ( (unsigned char*)&pI[0], &intBufferLen, compInt, len2 );
	return rsi;
}
void BlibRead::writeSpectrum ( ostream& os, int rsi, int numPeaks ) const
{
	os << "BEGIN IONS" << endl;
	os << "TITLE=Scan ";
	os << specIDInFile [rsi];
	os << " ";
	os << "(rt=";
	os << retentionTime [rsi];
	os << ")";
	os << " ";
	os << "[Prospector Created]";
	os << endl;

	os << "PEPMASS=";
	genPrint ( os, precursorMZ [rsi], 4 );
	os << endl;
	os << "CHARGE=" << precursorCharge [rsi] << "+" << endl;

	for ( DoubleVectorSizeType i = 0 ; i < numPeaks ; i++ ) {
		genPrint ( os, pM [i], 4 );
		os << " ";
		genPrintSigFig ( os, pI [i], 3 );
		os << endl;
	}
	os << "END IONS" << endl;
}
