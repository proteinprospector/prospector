/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lg_string.cpp                                                 *
*                                                                             *
*  Created    : August 1st 1996                                               *
*                                                                             *
*  Purpose    : String manipulation functions.                                *
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
#include <cstdio>
#ifdef VIS_C
#include <process.h>
#else
#include <unistd.h>
#include <errno.h>
#endif
#include <ctime>
#include <lg_string.h>
using std::string;

bool genStringIsFloat ( const string& str )
{
	errno = 0;
	char* ptr;
	double val = strtod ( str.c_str (), &ptr );
	if ( errno != 0 ) return false;
	if ( ptr == str.c_str () ) return false;	// conversion failed (no characters consumed)
	if ( *ptr != 0 ) return false;		        // conversion failed (trailing data)
	return true;
}
bool genStringIsPositiveFloat ( const string& str )
{
	errno = 0;
	char* ptr;
	double val = strtod ( str.c_str (), &ptr );
	if ( errno != 0 ) return false;
	if ( ptr == str.c_str () ) return false;	// conversion failed (no characters consumed)
	if ( *ptr != 0 ) return false;		        // conversion failed (trailing data)
	if ( val < 0.0 ) return false;
	return true;
}
bool genStringIsPositiveNonZeroFloat ( const string& str )
{
	errno = 0;
	char* ptr;
	double val = strtod ( str.c_str (), &ptr );
	if ( errno != 0 ) return false;
	if ( ptr == str.c_str () ) return false;	// conversion failed (no characters consumed)
	if ( *ptr != 0 ) return false;		        // conversion failed (trailing data)
	if ( val <= 0.0 ) return false;
	return true;
}
bool genStringIsInteger ( const string& str )
{
	errno = 0;
	char* ptr;
	int val = strtol ( str.c_str (), &ptr, 10 );
	if ( errno != 0 ) return false;
	if ( ptr == str.c_str () ) return false;	// conversion failed (no characters consumed)
	if ( *ptr != 0 ) return false;		        // conversion failed (trailing data)
	return true;
}
bool genIsLower ( const string& str )
{
	if ( str.empty () ) return false;
	for ( StringSizeType i = 0 ; i < str.length () ; i++ ) {
		if ( !islower ( str [i] ) ) return false;
	}
	return true;
}
bool genIsUpper ( const string& str )
{
	if ( str.empty () ) return false;
	for ( StringSizeType i = 0 ; i < str.length () ; i++ ) {
		if ( !isupper ( str [i] ) ) return false;
	}
	return true;
}
bool genAnyUpper ( const string& str )
{
	if ( str.empty () ) return false;
	for ( StringSizeType i = 0 ; i < str.length () ; i++ ) {
		if ( isupper ( str [i] ) ) return true;
	}
	return false;
}
string genToUpper ( const string& str )
{
	string newStr;
	for ( StringSizeType i = 0 ; i < str.length () ; i++ ) {
		newStr += toupper ( str [i] );
	}
	return newStr;
}
string genToLower ( const string& str )
{
	string newStr;
	for ( StringSizeType i = 0 ; i < str.length () ; i++ ) {
		newStr += tolower ( str [i] );
	}
	return newStr;
}
char* gen_strupr ( char* str )
{
#ifdef VIS_C
	return strupr ( str );
#else
	char* s;
	
	if ( str ) {
		for ( s = str ; *s ; ++s ) *s = toupper(*s);
	}
	return str;
#endif
}
string gen_strstrip ( const string& str )
{
	string newStr;

	for ( StringSizeType i = 0 ; i < str.length () ; i++ ) {
		if ( str [i] != ' ' && str [i] != '\n' && str [i] != '\r' && str [i] != '\t' ) {
			newStr += str [i];
		}
	}
	return newStr;
}
string gen_strstripnum ( const string& str )
{
	string newStr;

	for ( StringSizeType i = 0 ; i < str.length () ; i++ ) {
		if ( !isdigit ( str [i] ) ) {
			newStr += str [i];
		}
	}
	return newStr;
}
string gen_strstripchar ( const string& str, char c )
{
	string newStr;

	for ( StringSizeType i = 0 ; i < str.length () ; i++ ) {
		if ( str [i] != c ) {
			newStr += str [i];
		}
	}
	return newStr;
}
string gen_strstripchars ( const string& str, const string& s2 )
{
	string newStr;

	for ( StringSizeType i = 0 ; i < str.length () ; i++ ) {
		if ( s2.find ( str [i] ) == string::npos ) {
			newStr += str [i];
		}
	}
	return newStr;
}
string gen_strtrim ( const string& str )	// Trim white space off start and end of string
{
	if ( str.empty () ) return str;
	int start = 0;
	while ( str [start] == ' ' || str [start] == '\t' ) start++;
	int end = str.length () - 1;
	while ( str [end] == ' ' || str [end] == '\t' ) end--;
	end += 1;
	return str.substr ( start, end - start );
}
string gen_strtrim2 ( const string& str )	// Trim white space and return characters off start and end of string
{
	if ( str.empty () ) return str;
	int start = 0;
	while ( str [start] == ' ' || str [start] == '\t' || str [start] == '\r' || str [start] == '\n' ) start++;
	int end = str.length () - 1;
	while ( str [end] == ' ' || str [end] == '\t' || str [end] == '\r' || str [end] == '\n' ) end--;
	end += 1;
	return str.substr ( start, end - start );
}
string gen_strstriptags ( const string& str, char ltag, char rtag )
{
	string newStr;
	bool tag = false;
	for ( StringSizeType i = 0 ; i < str.length () ; i++ ) {
		char s = str [i];
		if ( s == ltag ) tag = true;
		if ( !tag ) newStr += s;
		if ( s == rtag ) tag = false;
	}
	return newStr;
}
string gen_strstriptags2 ( const string& str, char ltag, char rtag )
{
	string newStr;
	for ( StringSizeType i = 0 ; i < str.length () ; i++ ) {
		char s = str [i];
		if ( s == ltag ) {
			int ntag = 1;
			i++;
			for ( ; i < str.length () ; i++ ) {
				char a = str [i];
				if ( a == ltag ) ntag++;
				if ( a == rtag ) ntag--;
				if ( ntag == 0 ) break;
			}
		}
		else newStr += s;
	}
	return newStr;
}
string genStrtrimSemiColon ( const string& str )
{
	int len = str.length ();
	if ( len > 2 && str [len-1] == ';' ) {
		return str.substr ( 0, len - 1 ); // Delete trailing ';'
	}
	return str;
}
int gen_strerrcmp ( char* str1, char* str2, int maxerrs )
{
	int nerrs = 0;

	if ( strlen ( str1 ) != strlen ( str2 ) ) return -1;

	while ( *str1 && *str2 ) {
		if ( *str1++ != *str2++ ) {
			nerrs++;
			if ( nerrs > maxerrs ) return -1;
		}
	}
	return nerrs;
}
bool genEmptyString ( const string& str )
{
	for ( int i = 0 ; i < str.length () ; i++ ) {
		char c = str [i];
		if ( !isspace ( c ) && !iscntrl ( c ) ) return false;
	}
	return true;
}
bool gen_strasciiprint ( char* str )
{
	int len = strlen ( str );
	int i;

	for ( i = 0 ; i < len ; i++ ) {
		if ( isprint ( str [i] ) == false ) return false;
	}
	return true;
}
bool gen_strascii ( char* str )
{
	int len = strlen ( str );
	int i;

	for ( i = 0 ; i < len ; i++ ) {
		if ( isascii ( str [i] ) == false ) return false;
	}
	return true;
}
int gen_strcharcount ( const string& str, char character )
{
	int num = 0;
	for ( StringSizeType i = 0 ; i < str.length () ; i++ ) {
		if ( str [i] == character ) num++;
	}
	return num;
}
bool gen_strcontains ( const string& str, char character )
{
	for ( StringSizeType i = 0 ; i < str.length () ; i++ ) {
		if ( str [i] == character ) {
			return true;
		}
	}
	return false;
}
/* Returns a string with the characters in reverse order */
string gen_strrev ( const string& str )
{
	string newStr ( str );
	int len = str.length ();

	for ( int i = 0 ; i < len ; i++ ) {
		newStr [i] = str [len - i - 1];
	}
	return newStr;
}
int gen_longest_matching_prefix ( const char* str, char** prefix_list, int num_prefixes )
{
	int maxlen;
	int len;
	int i;
	int index;

	for ( i = 0, maxlen = 0 ; i < num_prefixes ; i++ ) {
		if ( strstr ( str, prefix_list [i] ) == str ) {
			len = strlen ( prefix_list [i] );
			if ( len > maxlen ) {
				maxlen = len;
				index = i;
			}
		}
	}
	if ( maxlen ) return index;
	else return -1;
}
string gen_itoa ( int number, const char* format )
{
	char temp [100];

	sprintf ( temp, format, number );
	return string ( temp );
}
string gen_ftoa ( double number, const char* format )
{
	char temp [100];

	sprintf ( temp, format, number );
	return string ( temp );
}
string gen_ftoa ( double number, int precision )
{
	char temp [100];
	char format [10];

	sprintf ( format, "%%.%df", precision );
	sprintf ( temp, format, number );
	return string ( temp );
}
bool isPrefix ( const string& str, const string& prefix )
{
	if ( prefix.length () > str.length () ) return false;
	return !str.compare ( 0, prefix.length (), prefix );
}
bool isNoCasePrefix ( const string& str, const string& prefix )
{
	if ( prefix.length () > str.length () ) return false;
	return !genStrncasecmp ( str, prefix, prefix.length () );
}
bool isSuffix ( const string& str, const string& suffix )
{
	if ( suffix.length () > str.length () ) return false;
	return !str.compare ( str.length () - suffix.length (), suffix.length (), suffix );
}
bool isNoCaseSuffix ( const string& str, const string& suffix )
{
	if ( suffix.length () > str.length () ) return false;
	return !genStrncasecmp ( str.substr (str.length () - suffix.length ()), suffix, suffix.length () );
}
int genNumSubstrings ( const string& str, const string& substr )
{
	int num = 0;
	int startInd = 0;
	int len = substr.length ();
	if ( len == 0 ) return 0;
	for ( ; ; ) {
		int ind = str.find ( substr, startInd );
		if ( ind == string::npos ) break;
		else {
			num++;
			startInd = ind + len;
		}
	}
	return num;
}
StringVector genGetSubstrings ( const string& str, char delim )
{
	string::size_type startInd = 0;
	StringVector sv;
	for ( ; ; ) {
		string::size_type ind = str.find ( delim, startInd );
		if ( ind == string::npos ) {
			string s = str.substr ( startInd );
			if ( !s.empty () ) sv.push_back ( s );
			break;
		}
		else {
			string s = str.substr ( startInd, ind-startInd );
			if ( !s.empty () ) sv.push_back ( s );
			startInd = ind + 1;
		}
	}
	return sv;
}
int genNumChars ( const string& str, const string& chars )
{
	if ( chars.empty () ) return 0;
	int num = 0;
	int startInd = 0;
	for ( ; ; ) {
		int ind = str.find_first_of ( chars, startInd );
		if ( ind == string::npos ) break;
		else {
			num++;
			startInd = ind + 1;
		}
	}
	return num;
}
string genStripSubstrings ( const string& str, const StringVector& subStr )
{
	string s = str;
	for ( StringVectorSizeType i = 0 ; i < subStr.size () ; i++ ) {
		s = genReplaceSubstrings ( s, subStr [i], "" );
	}
	return s;
}
string genReplaceSubstrings ( const string& str, const string& subStr, const string& newStr )
{
	string ret;
	int startInd = 0;
	int strLen = str.length ();
	int subStrLen = subStr.length ();
	for ( ; ; ) {
		int ind = str.find ( subStr, startInd );
		if ( ind == string::npos ) {
			ret += str.substr ( startInd, strLen - startInd );
			break;
		}
		else {
			ret += str.substr ( startInd, ind - startInd );
			ret += newStr;
			startInd = ind + subStrLen;
		}
	}
	return ret;
}
string genEscapePath ( const string& path )
{
	string ret;

	for ( StringSizeType i = 0 ; i < path.length () ; i++ ) {
		ret += path [i];
		if ( path [i] == '\\' ) ret += '\\';
	}
	return ret;
}
string genEscapeQuote ( const string& path )
{
	string ret;

	for ( StringSizeType i = 0 ; i < path.length () ; i++ ) {
		ret += path [i];
		if ( path [i] == '\'' ) ret += '\'';
	}
	return ret;
}
string genTranslateCharacter ( const string& s, char c1, char c2 )
{
	string ret;

	for ( StringSizeType i = 0 ; i < s.length () ; i++ ) {
		if ( s [i] == c1 ) ret += c2;
		else ret += s [i];
	}
	return ret;
}
string genNextString ( const string& s, const string& delim, string::size_type& start, string::size_type& end )
{
	end = s.find_first_of ( delim, start );
	if ( end == string::npos ) return "";
	string str = s.substr ( start, end-start );
	start = end + 1;
	return str;
}
int genNextInt ( const string& s, const string& delim, string::size_type& start, string::size_type& end )
{
	end = s.find_first_of ( delim, start );
	if ( end == string::npos ) return 0;
	int i = atoi ( s.substr ( start, end-start ).c_str () );
	start = end + 1;
	return i;
}
double genNextDouble ( const string& s, const string& delim, string::size_type& start, string::size_type& end )
{
	end = s.find_first_of ( delim, start );
	if ( end == string::npos ) return 0.0;
	double d = atof ( s.substr ( start, end-start ).c_str () );
	start = end + 1;
	return d;
}
StringVector genGetSeparatedValues ( const string& s, const string& sep )
{
	StringVector sv;
	if ( s.empty () ) return sv;
	string::size_type start = 0;
	for ( ; ; ) {
		string::size_type end = s.find ( sep, start );
		sv.push_back ( s.substr ( start, end-start ) );
		if ( end == string::npos ) break;
		start = end + sep.length ();
	}
	return sv;
}
string genRandomString ( int len )
{
	static int seeded = false;
	if ( !seeded ) {
		srand ( static_cast<unsigned int> ( time ( 0 ) + getpid () ) );
		seeded = true;
	}
	string s;
	s.reserve ( len );
	for ( int i = 0 ; i < len ; i++ ) {
		char asc;
		char ind = static_cast <char> ( 62 * ( rand () / ( RAND_MAX + 1.0 ) ) );
		if ( ind < 10 ) asc = ind + 48;
		else if ( ind < 36 ) asc = ind + 55;
		else asc = ind + 61;
		s += asc;
	}
	return s;
}
