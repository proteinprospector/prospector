/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_sqlite.cpp                                                 *
*                                                                             *
*  Created    : October 24th 2013                                             *
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
#ifndef VIS_C
#include <cstdio>
#endif

#include <iostream>
#include <string>
#include <sqlite3.h>
#include <lgen_define.h>
#include <lu_sqlite.h>
using std::string;
using std::cout;
using std::endl;

int PPSQLite::rowID;

PPSQLite::PPSQLite () :
	zErrMsg ( 0 )
{
}
PPSQLite::~PPSQLite ()
{
#ifdef VIS_C
	sqlite3_close_v2 ( db );
#else
	sqlite3_close ( db );
#endif
}
void PPSQLite::error ( const string& sql, int rc )
{
	cout << "function = sqlite3_exec<br />" << endl;
	cout << "sql = " << sql << "<br />" << endl;
	cout << "rc = " << rc << "<br />" << endl;
	cout << "err = " << zErrMsg << "<br />" << endl;
	sqlite3_free ( zErrMsg );
	exit ( 0 );
}
void PPSQLite::error ( const string& function, const string& sql, int rc )
{
	cout << "function = " << function << "<br />" << endl;
	cout << "sql = " << sql << "<br />" << endl;
	cout << "rc = " << rc << "<br />" << endl;
	cout << "error = " << sqlite3_errmsg ( db ) << "<br />" << endl;
	exit ( 0 );
}
void PPSQLite::executeSQL ( const string& sql )
{
	rc = sqlite3_exec ( db, sql.c_str (), callback, 0, &zErrMsg );
	if ( rc != SQLITE_OK ) {
		error ( sql, rc );
	}
}
void PPSQLite::executeSQL ( const string& sql, int (*cback)(void*,int,char**,char**) )
{
	rc = sqlite3_exec ( db, sql.c_str (), cback, 0, &zErrMsg );
	if ( rc != SQLITE_OK ) {
		error ( sql, rc );
	}
}
void PPSQLite::beginTransaction ()
{
	sqlite3_exec ( db, "BEGIN TRANSACTION", NULL, NULL, &zErrMsg );
}
void PPSQLite::endTransaction ()
{
	sqlite3_exec ( db, "END TRANSACTION", NULL, NULL, &zErrMsg );
}
int PPSQLite::insertQueryGetIndex ( const string& table, const StringVector& names, const StringVector& values )
{
	insertQuery ( table, names, values );
	string sql = "SELECT last_insert_rowid()";
	rc = sqlite3_exec ( db, sql.c_str (), rowidx_callback, 0, &zErrMsg );
	if ( rc != SQLITE_OK ) {
		error ( sql, rc );
	}
	return rowID;
}
void PPSQLite::insertQuery ( const string& table, const string& name, const string& value )
{
	string sql = "INSERT INTO " + table + " (";
	sql += name;
	sql += ") ";
	sql += "VALUES ( ";
	if ( value.empty () ) {
		sql += "NULL";
	}
	else {
		sql += "'";
		sql += value;
		sql += "'";
	}
	sql += " )";
	rc = sqlite3_exec ( db, sql.c_str (), callback, 0, &zErrMsg );
	if ( rc != SQLITE_OK ) {
		error ( sql, rc );
	}
}
void PPSQLite::insertQuery ( const string& table, const StringVector& names, const StringVector& values )
{

	string sql = "INSERT INTO " + table + " (";
	for ( StringVectorSizeType i = 0 ; i < names.size () ; i++ ) {
		sql += names [i];
		if ( i != names.size () - 1 ) sql += ",";
	}
	sql += ") ";
	sql += "VALUES ( ";
	for ( StringVectorSizeType j = 0 ; j < values.size () ; j++ ) {
		if ( values [j].empty () ) {
			sql += "NULL";
		}
		else {
			sql += "'";
			sql += values [j];
			sql += "'";
		}
		if ( j != values.size () - 1 ) sql += ",";
	}
	sql += " )";
	rc = sqlite3_exec ( db, sql.c_str (), callback, 0, &zErrMsg );
	if ( rc != SQLITE_OK ) {
		error ( sql, rc );
	}
/*
The sqlite3_exec() interface is a convenience wrapper around sqlite3_prepare_v2(),
sqlite3_step(), and sqlite3_finalize(), that allows an application to run multiple
statements of SQL without having to use a lot of C code.

The sqlite3_exec() interface runs zero or more UTF-8 encoded, semicolon-separate
SQL statements passed into its 2nd argument, in the context of the database connection
passed in as its 1st argument. If the callback function of the 3rd argument to
sqlite3_exec() is not NULL, then it is invoked for each result row coming out of
the evaluated SQL statements. The 4th argument to sqlite3_exec() is relayed through
to the 1st argument of each callback invocation. If the callback pointer to
sqlite3_exec() is NULL, then no callback is ever invoked and result rows are ignored.

If an error occurs while evaluating the SQL statements passed into sqlite3_exec(),
then execution of the current statement stops and subsequent statements are skipped.
If the 5th parameter to sqlite3_exec() is not NULL then any error message is written
into memory obtained from sqlite3_malloc() and passed back through the 5th parameter.
To avoid memory leaks, the application should invoke sqlite3_free() on error message
strings returned through the 5th parameter of of sqlite3_exec() after the error message
string is no longer needed. If the 5th parameter to sqlite3_exec() is not NULL and no
errors occur, then sqlite3_exec() sets the pointer in its 5th parameter to NULL before
returning.

If an sqlite3_exec() callback returns non-zero, the sqlite3_exec() routine returns
SQLITE_ABORT without invoking the callback again and without running any subsequent
SQL statements.

The 2nd argument to the sqlite3_exec() callback function is the number of columns in
the result. The 3rd argument to the sqlite3_exec() callback is an array of pointers to
strings obtained as if from sqlite3_column_text(), one for each column. If an element
of a result row is NULL then the corresponding string pointer for the sqlite3_exec()
callback is a NULL pointer. The 4th argument to the sqlite3_exec() callback is an array
of pointers to strings where each entry represents the name of corresponding result
column as obtained from sqlite3_column_name().

If the 2nd parameter to sqlite3_exec() is a NULL pointer, a pointer to an empty string,
or a pointer that contains only whitespace and/or SQL comments, then no SQL statements
are evaluated and the database is not changed.

Restrictions:

The application must insure that the 1st parameter to sqlite3_exec() is a valid and open
database connection.

The application must not close database connection specified by the 1st parameter to
sqlite3_exec() while sqlite3_exec() is running.

The application must not modify the SQL statement text passed into the 2nd parameter of
sqlite3_exec() while sqlite3_exec() is running.

The return value of sqlite3_exec() shall be SQLITE_OK if all SQL statements run
successfully and to completion.

The return value of sqlite3_exec() shall be an appropriate non-zero error code if
any SQL statement fails.
*/
}
int PPSQLite::callback ( void* NotUsed, int argc, char** argv, char** azColName )
{
	for ( int i = 0 ; i < argc ; i++ ) {
		printf ( "%s = %s\n", azColName [i], argv [i] ? argv[i] : "NULL" );
	}
	printf ( "\n" );
	return 0;
}
int PPSQLite::rowidx_callback ( void* NotUsed, int argc, char** argv, char** azColName )
{
	rowID = atoi ( argv [0] );
	return 0;
}

namespace {

StringVector test;

int libInfoCallback ( void* NotUsed, int argc, char** argv, char** azColName )
{
	test.push_back ( argv [0] );
	return 0;
}

}

SQLiteIdentify::SQLiteIdentify ( const string& name )
{
	rc = sqlite3_open_v2 ( name.c_str (), &db, SQLITE_OPEN_READONLY, NULL );
}
SQLiteIdentify::~SQLiteIdentify ()
{
}
void SQLiteIdentify::testID ()
{
	string sql = "SELECT \
		libLSID from LibInfo";
	rc = sqlite3_exec ( db, sql.c_str (), libInfoCallback, 0, &zErrMsg );
}
bool SQLiteIdentify::isBLib () const
{
	if ( !test.empty () ) {
		string s = test [0];
		return s.find ( "bibliospec" ) != string::npos;
	}
	return false;
}
bool SQLiteIdentify::isMSF () const
{
	return test.empty ();
}
