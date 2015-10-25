/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lgen_sort.h                                                   *
*                                                                             *
*  Created    : June 20th 1996                                                *
*                                                                             *
*  Purpose    : Functions for the gen_qsort function argument.                *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1996-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lgen_sort_h
#define __lgen_sort_h

#include <cstring>
#include <vector>
#include <string>
// The following notice applies to the strcasecmp function below
/*
 * Copyright (c) 1987 Regents of the University of California.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms are permitted
 * provided that this notice is preserved and that due credit is given
 * to the University of California at Berkeley. The name of the University
 * may not be used to endorse or promote products derived from this
 * software without specific written prior permission. This software
 * is provided ``as is'' without express or implied warranty.
 */
static unsigned char charmap[] = {
	'\000', '\001', '\002', '\003', '\004', '\005', '\006', '\007',
	'\010', '\011', '\012', '\013', '\014', '\015', '\016', '\017',
	'\020', '\021', '\022', '\023', '\024', '\025', '\026', '\027',
	'\030', '\031', '\032', '\033', '\034', '\035', '\036', '\037',
	'\040', '\041', '\042', '\043', '\044', '\045', '\046', '\047',
	'\050', '\051', '\052', '\053', '\054', '\055', '\056', '\057',
	'\060', '\061', '\062', '\063', '\064', '\065', '\066', '\067',
	'\070', '\071', '\072', '\073', '\074', '\075', '\076', '\077',
	'\100', '\141', '\142', '\143', '\144', '\145', '\146', '\147',
	'\150', '\151', '\152', '\153', '\154', '\155', '\156', '\157',
	'\160', '\161', '\162', '\163', '\164', '\165', '\166', '\167',
	'\170', '\171', '\172', '\133', '\134', '\135', '\136', '\137',
	'\140', '\141', '\142', '\143', '\144', '\145', '\146', '\147',
	'\150', '\151', '\152', '\153', '\154', '\155', '\156', '\157',
	'\160', '\161', '\162', '\163', '\164', '\165', '\166', '\167',
	'\170', '\171', '\172', '\173', '\174', '\175', '\176', '\177',
	'\200', '\201', '\202', '\203', '\204', '\205', '\206', '\207',
	'\210', '\211', '\212', '\213', '\214', '\215', '\216', '\217',
	'\220', '\221', '\222', '\223', '\224', '\225', '\226', '\227',
	'\230', '\231', '\232', '\233', '\234', '\235', '\236', '\237',
	'\240', '\241', '\242', '\243', '\244', '\245', '\246', '\247',
	'\250', '\251', '\252', '\253', '\254', '\255', '\256', '\257',
	'\260', '\261', '\262', '\263', '\264', '\265', '\266', '\267',
	'\270', '\271', '\272', '\273', '\274', '\275', '\276', '\277',
	'\300', '\341', '\342', '\343', '\344', '\345', '\346', '\347',
	'\350', '\351', '\352', '\353', '\354', '\355', '\356', '\357',
	'\360', '\361', '\362', '\363', '\364', '\365', '\366', '\367',
	'\370', '\371', '\372', '\333', '\334', '\335', '\336', '\337',
	'\340', '\341', '\342', '\343', '\344', '\345', '\346', '\347',
	'\350', '\351', '\352', '\353', '\354', '\355', '\356', '\357',
	'\360', '\361', '\362', '\363', '\364', '\365', '\366', '\367',
	'\370', '\371', '\372', '\373', '\374', '\375', '\376', '\377',
};

inline int genStrcasecmp ( const char *s1, const char *s2, char c )
{
    register unsigned char u1, u2;

    for (;;) {
		u1 = (unsigned char) *s1++;
		u2 = (unsigned char) *s2++;
		if (charmap[u1] != charmap[u2]) {
			return charmap[u1] - charmap[u2];
		}
		if ( u1 == c ) {
			return 0;
		}
    }
}
inline int genStrcasecmp ( const char *s1, const char *s2 )
{
    register unsigned char u1, u2;

    for (;;) {
		u1 = (unsigned char) *s1++;
		u2 = (unsigned char) *s2++;
		if (charmap[u1] != charmap[u2]) {
			return charmap[u1] - charmap[u2];
		}
		if ( u1 == '\0' ) {
			return 0;
		}
    }
}
inline int genStrncasecmp ( const char *s1, const char *s2, register size_t n )
{
	register unsigned char u1, u2;

	for (; n != 0; --n) {
		u1 = (unsigned char) *s1++;
		u2 = (unsigned char) *s2++;
		if (charmap[u1] != charmap[u2]) {
			return charmap[u1] - charmap[u2];
		}
		if (u1 == '\0') {
			return 0;
		}
	}
	return 0;
}
inline int genStrcasecmp ( const std::string& s1, const std::string& s2 )
{
	return genStrcasecmp ( s1.c_str (), s2.c_str () );
}
inline int genStrncasecmp ( const std::string& s1, const std::string& s2, register size_t n )
{
	return genStrncasecmp ( s1.c_str (), s2.c_str (), n );
}

struct IndexedInt {
	int number;
	int index;
	bool operator== ( const IndexedInt& rhs ) const
	{
		return number == rhs.number && index == rhs.index;
	}
};
typedef std::vector <IndexedInt> VectorIndexedInt;
typedef VectorIndexedInt::const_iterator VectorIndexedIntConstIterator;
typedef VectorIndexedInt::size_type VectorIndexedIntSizeType;

struct IndexedDouble {
	double number;
	int index;
};

struct IndexedConstCString {
	const char* name;
	int index;
};
typedef std::vector <IndexedConstCString> VectorIndexedConstCString;
typedef VectorIndexedConstCString::const_iterator VectorIndexedConstCStringConstIterator;
typedef VectorIndexedConstCString::size_type VectorIndexedConstCStringSizeType;

struct IndexedCString {
	char* name;
	int index;
};

struct IndexedString {
	std::string name;
	int index;
};

struct DoubleString {
	char* array;
	double number;
};

class SortDoubleStringAscending {
	public:
		int operator () ( const DoubleString& a, const DoubleString& b ) const
		{
			return ( a.number < b.number );
		}
};
class SortIndexedIntAscending {
	public:
		int operator () ( const IndexedInt& a, const IndexedInt& b ) const
		{
			return ( a.number == b.number ) ? a.index < b.index : a.number < b.number;
		}
		int operator () ( const IndexedInt& a, const int b ) const
		{
			return a.number < b;
		}
};
class SortIndexedDoubleAscending {
	public:
		int operator () ( const IndexedDouble& a, const IndexedDouble& b ) const
		{
			return ( a.number == b.number ) ? a.index < b.index : a.number < b.number;
		}
};
class SortIndexedStrcasecmpAscending {
	public:
		bool operator () ( const IndexedConstCString& a, const IndexedConstCString& b ) const
		{
			int ret = genStrcasecmp ( a.name, b.name );

			if ( ret == 0 ) {
				return ( a.index < b.index );
			}
			return ( ret < 0 );
		}	
		bool operator () ( const IndexedConstCString& a, const char* b ) const
		{
			return genStrcasecmp ( a.name, b ) < 0;
		}	
};
class SortIndexedStrcmpAscending {
	public:
		bool operator () ( const IndexedConstCString& a, const IndexedConstCString& b ) const
		{
			int ret = strcmp ( a.name, b.name );

			if ( ret == 0 ) {
				return ( a.index < b.index );
			}
			return ( ret < 0 );
		}	
};
class SortStrcmpAscending {
	public:
		bool operator () ( const char* a, const char* b ) const
		{
			return ( strcmp ( a, b ) < 0 );
		}	
};

class genStrcasecmpAscending {
public:
	bool operator () ( const char* a, const char* b ) const
	{
		return genStrcasecmp ( a, b ) < 0;
	}
	bool operator () ( std::string const& a, std::string const& b ) const
	{
		return genStrcasecmp ( a.c_str (), b.c_str () ) < 0;
	}
};

class genStrcasecmpAscendingDelim {
	char delim;
public:
	genStrcasecmpAscendingDelim ( char delim ) :
		delim ( delim ) {}
	bool operator () ( const char* a, const char* b ) const
	{
		return genStrcasecmp ( a, b, delim ) < 0;
	}
};

class genNumberStringAscending {
public:
	bool operator () ( const std::string& a, const std::string& b ) const
	{
		return atof ( a.c_str () ) < atof ( b.c_str () );
	}
};

class CheckStrcasecmpEquality {
	std::string ss;
public:
	CheckStrcasecmpEquality ( const std::string& ss ) :
		ss ( ss ) {}
	bool operator () ( const std::string& s )
	{
		return genStrcasecmp ( s.c_str (), ss.c_str () ) == 0;
	}
};

#endif /* ! __lgen_sort_h */
