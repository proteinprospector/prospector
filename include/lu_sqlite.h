 /*****************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_sqlite.h                                                   *
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

#ifndef __lu_sqlite_h
#define __lu_sqlite_h

#include <string>

struct sqlite3;

class PPSQLite {
protected:
	sqlite3* db;
	int rc;
	char* zErrMsg;
	static int rowID;
	void executeSQL ( const std::string& sql );
	void executeSQL ( const std::string& sql, int (*callback)(void*,int,char**,char**) );
	void error ( const std::string& sql, int rc );
	void error ( const std::string& function, const std::string& sql, int rc );
	int insertQueryGetIndex ( const std::string& table, const StringVector& names, const StringVector& values );
	void insertQuery ( const std::string& table, const std::string& name, const std::string& value );
	void insertQuery ( const std::string& table, const StringVector& name, const StringVector& value );
	static int callback ( void* NotUsed, int argc, char** argv, char** azColName );
	static int rowidx_callback ( void* NotUsed, int argc, char** argv, char** azColName );
public:
	PPSQLite ();
	virtual ~PPSQLite ();
	void beginTransaction ();
	void endTransaction ();
};

class SQLiteIdentify : public PPSQLite {
public:
	SQLiteIdentify ( const std::string& name );
	~SQLiteIdentify ();
	void testID ();
	bool isBLib () const;
	bool isMSF () const;
};


#endif /* ! __lu_sqlite_h */
