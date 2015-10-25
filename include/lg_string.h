/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lg_string.h                                                   *
*                                                                             *
*  Created    : June 23rd 1996                                                *
*                                                                             *
*  Purpose    : Machine independent interface to string.h.                    *
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

#ifndef __lg_string_h
#define __lg_string_h

#include <cstring>

#include <string>
#include <lgen_define.h>

bool genStringIsFloat ( const std::string& str );
bool genStringIsPositiveFloat ( const std::string& str );
bool genStringIsPositiveNonZeroFloat ( const std::string& str );
bool genStringIsInteger ( const std::string& str );
bool genIsLower ( const std::string& str );
bool genIsUpper ( const std::string& str );
bool genAnyUpper ( const std::string& str );
std::string genToUpper ( const std::string& str );
std::string genToLower ( const std::string& str );
char* gen_strupr ( char* str );
std::string gen_strstrip ( const std::string& str );
std::string gen_strstripnum ( const std::string& str );
std::string gen_strstripchar ( const std::string& str, char c );
std::string gen_strstripchars ( const std::string& str, const std::string& s2 );
std::string gen_strtrim ( const std::string& str );
std::string gen_strtrim2 ( const std::string& str );
std::string gen_strstriptags ( const std::string& str, char ltag = '<', char rtag = '>' );
std::string gen_strstriptags2 ( const std::string& str, char ltag = '<', char rtag = '>' );
std::string genStrtrimSemiColon ( const std::string& str );
int gen_strerrcmp ( char* str1, char* str2, int maxerrs );
bool genEmptyString ( const std::string& str );
bool gen_strasciiprint ( char* str );
bool gen_strascii ( char* str );
int gen_strcharcount ( const std::string& str, char character );
bool gen_strcontains ( const std::string& str, char character );
int gen_longest_matching_prefix ( const char* str, char** prefix_list, int num_prefixes );
std::string gen_strrev ( const std::string& str );
std::string gen_itoa ( int number, const char* format = "%d" );
std::string gen_ftoa ( double number, const char* format );
std::string gen_ftoa ( double number, int precision );
bool isPrefix ( const std::string& str, const std::string& prefix );
bool isNoCasePrefix ( const std::string& str, const std::string& suffix );
bool isSuffix ( const std::string& str, const std::string& suffix );
bool isNoCaseSuffix ( const std::string& str, const std::string& suffix );
int genNumSubstrings ( const std::string& str, const std::string& substr );
StringVector genGetSubstrings ( const std::string& str, char delim );
int genNumChars ( const std::string& str, const std::string& chars );
std::string genStripSubstrings ( const std::string& str, const StringVector& subStr );
std::string genReplaceSubstrings ( const std::string& str, const std::string& subStr, const std::string& newStr );
std::string genEscapePath ( const std::string& path );
std::string genEscapeQuote ( const std::string& path );
std::string genTranslateCharacter ( const std::string& s, char c1, char c2 );
std::string genNextString ( const std::string& s, const std::string& delim, std::string::size_type& start, std::string::size_type& end );
int genNextInt ( const std::string& s, const std::string& delim, std::string::size_type& start, std::string::size_type& end );
double genNextDouble ( const std::string& s, const std::string& delim, std::string::size_type& start, std::string::size_type& end );
StringVector genGetSeparatedValues ( const std::string& s, const std::string& sep );
std::string genRandomString ( int len );

#endif /* ! __lg_string_h */
