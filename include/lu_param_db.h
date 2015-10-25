 /*****************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_param_db.h                                                 *
*                                                                             *
*  Created    : Dec 18th 2014                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2014-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_param_db_h
#define __lu_param_db_h

#include <string>
#include <lgen_define.h>
#include <lu_sqlite.h>

class ParamDB : public PPSQLite {
public:
	ParamDB ();
	virtual ~ParamDB ();
};

class ParamDBWrite : public ParamDB {
	friend class SQLiteTable;
public:
	ParamDBWrite ( const std::string& name, bool append = false );
	~ParamDBWrite ();
	void create ();
};

class ParamDBRead : public ParamDB {
public:
	ParamDBRead ( const std::string& name );
	~ParamDBRead ();
};

#endif /* ! __lu_param_db_h */
