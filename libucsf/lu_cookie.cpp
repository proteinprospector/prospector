/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_cookie.cpp                                                 *
*                                                                             *
*  Created    : September 24th 2007                                           *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2007-2011) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef MYSQL_DATABASE
#include <lg_string.h>
#include <lg_string.h>
#include <ld_init.h>
#include <lu_html.h>

using std::ostream;
using std::string;

bool setUserCookie ( ostream& os, const string& user )
{
	string key = genRandomString ( 32 );
	string userID = MySQLPPSDDBase::instance ().getUserID ( user );
	if ( userID.empty () ) return false;
	if ( MySQLPPSDDBase::instance ().submitSession ( key, userID ) ) { // Checks it is unique
		setCookie ( os, "user", user, true );	// User is remember for 1 year
		setCookie ( os, "key", key, false );	// Key is remembered only for the current session
		return true;
	}
	return false;
}
#endif
