/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lg_email.cpp                                                  *
*                                                                             *
*  Created    : June 29th 2007                                                *
*                                                                             *
*  Purpose    : Checks email format.                   .                      *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2007-2009) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <cstring>
#include <string>
using std::string;

bool isValidEmailAddress ( const char* address )
{
	int count = 0;
	const char* c;
	const char* domain;
	static const char* rfc822_specials = "()<>@,;:\\\"[]";

	/* first we validate the name portion (name@domain) */
	for ( c = address ;  *c ;  c++ ) {
		if ( *c == '\"' && ( c == address || *(c - 1) == '.' || *(c - 1) == '\"' ) ) {
			while ( *++c ) {
				if ( *c == '\"' ) break;
				if ( *c == '\\' && ( *++c == ' ' ) ) continue;
				if ( *c <= ' ' || *c >= 127 ) return 0;
			}
			if ( !*c++ ) return 0;
			if ( *c == '@' ) break;
			if ( *c != '.' ) return 0;
			continue;
		}
		if ( *c == '@' ) break;
		if ( *c <= ' ' || *c >= 127 ) return 0;
		if ( strchr ( rfc822_specials, *c ) ) return 0;
	}
	if ( c == address || *(c - 1) == '.' ) return 0;

	/* next we validate the domain portion (name@domain) */
	if ( !*(domain = ++c) ) return 0;
	do {
		if ( *c == '.' ) {
			if ( c == domain || *(c - 1) == '.' ) return 0;
			count++;
		}
		if ( *c <= ' ' || *c >= 127 ) return 0;
		if ( strchr ( rfc822_specials, *c ) ) return 0;
	} while ( *++c );

	return ( count >= 1 );
}
bool isValidName ( const string& name )
{
	string::size_type len = name.length ();
	if ( len == 0 ) return false;
	for ( string::size_type i = 0 ; i < len ; i++ ) {
		char c = name [i];
		bool hyphen = ( c == '-' ) && ( i != 0 ) && ( i != len - 1 );
		if ( !isalpha ( c ) && !hyphen ) return false;
	}
	return true;
}
bool isValidUser ( const string& user, string::size_type min, string::size_type max )
{
	string::size_type len = user.length ();
	if ( len < min ) return false;
	if ( len > max ) return false;
	for ( string::size_type i = 0 ; i < len ; i++ ) {
		char c = user [i];
		if ( i == 0 ) {
			if ( !islower ( c ) ) return false;
		}
		else {
			if ( !islower ( c ) && ! isdigit ( c ) ) return false;
		}
	}
	return true;
}
