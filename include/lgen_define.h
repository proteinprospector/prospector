/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lgen_define.h                                                 *
*                                                                             *
*  Created    : June 20th 1996                                                *
*                                                                             *
*  Purpose    : Some commonly used constants and macros.                      *
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

#ifndef __lgen_define_h
#define __lgen_define_h

template <class T>
inline const T& genMin(const T& a, const T& b) {
    return b < a ? b : a;
}

template <class T>
inline const T& genMax(const T& a, const T& b) {
    return  a < b ? b : a;
}

template <class T>
inline const T genAbsDiff(const T& a, const T& b) {
    return  (( a > b ) ? a - b : b - a );
}

template <class T>
inline const bool genIsOdd(const T& a) {
    return  ( ((int)a) & 1 );
}

template <class T>
inline const bool genIsEven(const T& a) {
    return  (!genIsOdd(a));
}

#include <deque>
#include <vector>
#include <list>
#include <string>
#include <stdlib.h>
#include <lgen_sort.h>

#define genArraySize(array) sizeof(array) / sizeof(array[0])

#ifdef VIS_C
typedef unsigned __int64 GENUINT64;
typedef __int64 GENINT64;
#else
#include <inttypes.h>
typedef uint64_t GENUINT64;
typedef int64_t GENINT64;
#endif

typedef std::deque <bool> BoolDeque;	// vector <bool> not recommended - see item 18,
										// Effective STL by Scott Meyers
typedef BoolDeque::size_type BoolDequeSizeType;

typedef std::string::size_type StringSizeType;

typedef std::vector <char> CharVector;
typedef CharVector::size_type CharVectorSizeType;

typedef std::vector <CharVector> CharVectorVector;
typedef CharVectorVector::size_type CharVectorVectorSizeType;

typedef std::vector <char*> CharPtrVector;
typedef CharPtrVector::size_type CharPtrVectorSizeType;
typedef CharPtrVector::iterator CharPtrVectorIterator;
typedef CharPtrVector::const_iterator CharPtrVectorConstIterator;
typedef std::vector <CharPtrVector> CharPtrVectorVector;

typedef std::vector <const char*> ConstCharPtrVector;
typedef ConstCharPtrVector::size_type ConstCharPtrVectorSizeType;
typedef ConstCharPtrVector::iterator ConstCharPtrVectorIterator;
typedef ConstCharPtrVector::const_iterator ConstCharPtrVectorConstIterator;
typedef std::vector <ConstCharPtrVector> ConstCharPtrVectorVector;

typedef std::vector <int> IntVector;
typedef IntVector::size_type IntVectorSizeType;
typedef IntVector::iterator IntVectorIterator;
typedef IntVector::const_iterator IntVectorConstIterator;
typedef std::vector <IntVector> IntVectorVector;
typedef IntVectorVector::size_type IntVectorVectorSizeType;
typedef IntVectorVector::iterator IntVectorVectorIterator;
typedef IntVectorVector::const_iterator IntVectorVectorConstIterator;

typedef std::vector <GENINT64> Int64Vector;
typedef std::vector <Int64Vector> Int64VectorVector;

typedef std::vector <unsigned int> UIntVector;
typedef UIntVector::size_type UIntVectorSizeType;
typedef UIntVector::iterator UIntVectorIterator;
typedef UIntVector::const_iterator UIntVectorConstIterator;
typedef std::vector <UIntVector> UIntVectorVector;

typedef std::vector <float> FloatVector;
typedef std::vector <FloatVector> FloatVectorVector;

typedef std::vector <double> DoubleVector;
typedef DoubleVector::size_type DoubleVectorSizeType;
typedef DoubleVector::iterator DoubleVectorIterator;

typedef std::vector <DoubleVector> DoubleVectorVector;
typedef DoubleVectorVector::size_type DoubleVectorVectorSizeType;

typedef std::vector <DoubleVectorVector> DoubleVectorVectorVector;
typedef DoubleVectorVectorVector::size_type DoubleVectorVectorVectorSizeType;

typedef std::vector <time_t> TimeVector;

typedef std::vector <std::string> StringVector;
typedef StringVector::size_type StringVectorSizeType;
typedef std::vector <std::string>::iterator StringVectorIterator;
typedef std::vector <StringVector> StringVectorVector;
typedef StringVectorVector::size_type StringVectorVectorSizeType;

typedef std::list <std::string> ListString;

#ifdef VIS_C
const char SLASH = '\\';
#else
const char SLASH = '/';
#endif

template <class T>
std::vector <T> reverseVector ( const std::vector <T>& vec )
{
	std::vector <T> v;
	std::copy ( vec.rbegin (), vec.rend (), std::back_inserter ( v ) );
	return v;
}
template <class T>
void deleteEntries ( const BoolDeque& flags, T& container )
{
	if ( !container.empty () ) {
		typename T::iterator it = container.begin ();
		for ( typename T::size_type i = 0 ; i < flags.size () ; i++ ) {
			if ( !flags [i] ) container.erase ( it );
			else it++;
		}
	}
}

typedef std::string::iterator StringIterator;
typedef std::string::const_iterator StringConstIterator;

typedef StringVector::iterator StringVectorIterator;
typedef StringVector::const_iterator StringVectorConstIterator;

typedef BoolDeque::iterator BoolDequeIterator;
typedef BoolDeque::const_iterator BoolDequeConstIterator;

#include <set>

typedef std::set <char> SetChar;
typedef SetChar::iterator SetCharIterator;
typedef SetChar::const_iterator SetCharConstIterator;

typedef std::set <int> SetInt;
typedef SetInt::iterator SetIntIterator;
typedef SetInt::const_iterator SetIntConstIterator;
typedef std::pair <SetIntIterator, bool> PairSetIntIteratorBool;
typedef std::vector <SetInt> VectorSetInt;
typedef VectorSetInt::iterator VectorSetIntIterator;
typedef VectorSetInt::const_iterator VectorSetIntConstIterator;

typedef std::set <unsigned int> SetUInt;
typedef SetUInt::iterator SetUIntIterator;
typedef SetUInt::const_iterator SetUIntConstIterator;

typedef std::set <double> SetDouble;
typedef SetDouble::iterator SetDoubleIterator;
typedef SetDouble::const_iterator SetDoubleConstIterator;
typedef std::pair <SetDoubleIterator, bool> PairSetDoubleIteratorBool;

typedef std::set <std::string> SetString;
typedef SetString::iterator SetStringIterator;
typedef SetString::const_iterator SetStringConstIterator;
typedef std::pair <SetStringIterator, bool> PairSetStringIteratorBool;

typedef std::set <std::string, genStrcasecmpAscending> NoCaseSetString;
typedef NoCaseSetString::iterator NoCaseSetStringIterator;
typedef NoCaseSetString::const_iterator NoCaseSetStringConstIterator;
typedef std::pair <NoCaseSetStringIterator, bool> PairNoCaseSetStringIteratorBool;

typedef std::set <StringVector> SetStringVector;
typedef SetStringVector::iterator SetStringVectorIterator;
typedef SetStringVector::const_iterator SetStringVectorConstIterator;

typedef std::set <int*> SetIntPtr;
typedef SetIntPtr::iterator SetIntPtrIterator;
typedef SetIntPtr::const_iterator SetIntPtrConstIterator;

typedef std::set <SetInt> SetSetInt;
typedef SetSetInt::iterator SetSetIntIterator;
typedef SetSetInt::const_iterator SetSetIntConstIterator;

#include <map>

typedef std::map <char, char> MapCharToChar;
typedef MapCharToChar::iterator MapCharToCharIterator;
typedef MapCharToChar::const_iterator MapCharToCharConstIterator;

typedef std::map <char, unsigned char> MapCharToUChar;
typedef MapCharToUChar::iterator MapCharToUCharIterator;
typedef MapCharToUChar::const_iterator MapCharToUCharConstIterator;

typedef std::map <char, std::string> MapCharToString;
typedef MapCharToString::iterator MapCharToStringIterator;
typedef MapCharToString::const_iterator MapCharToStringConstIterator;

typedef std::map <char, int> MapCharToInt;
typedef MapCharToInt::iterator MapCharToIntIterator;
typedef MapCharToInt::const_iterator MapCharToIntConstIterator;

typedef std::map <char, double> MapCharToDouble;
typedef MapCharToDouble::iterator MapCharToDoubleIterator;
typedef MapCharToDouble::const_iterator MapCharToDoubleConstIterator;

typedef std::map <char, StringVector> MapCharToStringVector;
typedef MapCharToStringVector::iterator MapCharToStringVectorIterator;
typedef MapCharToStringVector::const_iterator MapCharToStringVectorConstIterator;

typedef std::map <unsigned int, unsigned int> MapUIntToUInt;
typedef MapUIntToUInt::iterator MapUIntToUIntIterator;
typedef MapUIntToUInt::const_iterator MapUIntToUIntConstIterator;

typedef std::map <int, int> MapIntToInt;
typedef MapIntToInt::iterator MapIntToIntIterator;
typedef MapIntToInt::const_iterator MapIntToIntConstIterator;

typedef std::map <int, double> MapIntToDouble;
typedef MapIntToDouble::iterator MapIntToDoubleIterator;
typedef MapIntToDouble::const_iterator MapIntToDoubleConstIterator;

typedef std::map <std::string, char> MapStringToChar;
typedef MapStringToChar::iterator MapStringToCharIterator;
typedef MapStringToChar::const_iterator MapStringToCharConstIterator;

typedef std::map <std::string, char*> MapStringToCharPtr;
typedef MapStringToCharPtr::iterator MapStringToCharPtrIterator;
typedef MapStringToCharPtr::const_iterator MapStringToCharPtrConstIterator;

typedef std::map <std::string, unsigned int> MapStringToUInt;
typedef MapStringToUInt::iterator MapStringToUIntIterator;
typedef MapStringToUInt::const_iterator MapStringToUIntConstIterator;

typedef std::map <std::string, int> MapStringToInt;
typedef MapStringToInt::iterator MapStringToIntIterator;
typedef MapStringToInt::const_iterator MapStringToIntConstIterator;

typedef std::map <std::string, double> MapStringToDouble;
typedef MapStringToDouble::iterator MapStringToDoubleIterator;
typedef MapStringToDouble::const_iterator MapStringToDoubleConstIterator;

typedef std::map <std::string, double, genStrcasecmpAscending> MapNoCaseStringToDouble;
typedef MapNoCaseStringToDouble::iterator MapNoCaseStringToDoubleIterator;
typedef MapNoCaseStringToDouble::const_iterator MapNoCaseStringToDoubleConstIterator;

typedef std::map <std::string, std::string> MapStringToString;
typedef MapStringToString::iterator MapStringToStringIterator;
typedef MapStringToString::const_iterator MapStringToStringConstIterator;
typedef std::vector <MapStringToString> VectorMapStringToString;

typedef std::map <std::string, StringVector> MapStringToStringVector;
typedef MapStringToStringVector::iterator MapStringToStringVectorIterator;
typedef MapStringToStringVector::const_iterator MapStringToStringVectorConstIterator;

typedef std::map <int, IntVector> MapIntToIntVector;
typedef MapIntToIntVector::iterator MapIntToIntVectorIterator;
typedef MapIntToIntVector::const_iterator MapIntToIntVectorConstIterator;

typedef std::map <int, std::string> MapIntToString;
typedef MapIntToString::iterator MapIntToStringIterator;
typedef MapIntToString::const_iterator MapIntToStringConstIterator;

typedef std::map <double, int> MapDoubleToInt;
typedef MapDoubleToInt::iterator MapDoubleToIntIterator;
typedef MapDoubleToInt::const_iterator MapDoubleToIntConstIterator;

typedef std::map <std::string, std::streampos> MapStringToStreampos;
typedef MapStringToStreampos::iterator MapStringToStreamposIterator;
typedef MapStringToStreampos::const_iterator MapStringToStreamposConstIterator;

typedef std::pair <int, char> PairIntChar;
typedef std::vector <PairIntChar> VectorPairIntChar;
typedef VectorPairIntChar::size_type VectorPairIntCharSizeType;

typedef std::pair <char, std::string> PairCharString;
typedef std::vector <PairCharString> VectorPairCharString;
typedef VectorPairCharString::size_type VectorPairCharStringSizeType;

typedef std::pair <int, int> PairIntInt;
typedef std::vector <PairIntInt> VectorPairIntInt;
typedef VectorPairIntInt::size_type VectorPairIntIntSizeType;

typedef std::pair <double, double> PairDoubleDouble;

typedef std::pair <int, double> PairIntDouble;
typedef std::vector <PairIntDouble> VectorPairIntDouble;

typedef std::pair <std::string, std::string> PairStringString;
typedef std::vector <PairStringString> VectorPairStringString;
typedef VectorPairStringString::size_type VectorPairStringStringSizeType;

typedef std::pair <std::string, time_t> PairStringTime;
typedef std::vector <PairStringTime> VectorPairStringTime;

typedef std::pair <int, std::string> PairIntString;
typedef std::vector <PairIntString> VectorPairIntString;
typedef VectorPairIntString::iterator VectorPairIntStringIterator;
typedef VectorPairIntString::const_iterator VectorPairIntStringConstIterator;
typedef VectorPairIntString::size_type VectorPairIntStringSizeType;

typedef std::pair <std::string, bool> PairStringBool;
typedef std::pair <std::string, int> PairStringInt;
typedef std::pair <std::string, double> PairStringDouble;

typedef std::pair <int, IntVector> PairIntIntVector;

typedef std::pair <StringVector, std::string> PairStringVectorString;

typedef std::pair <StringVector, StringVector> PairStringVectorStringVector;

typedef std::pair <int, PairStringDouble> PairIntPairStringDouble;
typedef std::vector <PairIntPairStringDouble> VectorPairIntPairStringDouble;

typedef std::map <int, PairStringBool> MapIntToPairStringBool;
typedef MapIntToPairStringBool::iterator MapIntToPairStringBoolIterator;
typedef MapIntToPairStringBool::const_iterator MapIntToPairStringBoolConstIterator;

typedef std::map <int, PairIntInt> MapIntToPairIntInt;
typedef MapIntToPairIntInt::iterator MapIntToPairIntIntIterator;
typedef MapIntToPairIntInt::const_iterator MapIntToPairIntIntConstIterator;
typedef std::vector <MapIntToPairIntInt> VectorMapIntToPairIntInt;
typedef VectorMapIntToPairIntInt::const_iterator VectorMapIntToPairIntIntConstIterator;
typedef VectorMapIntToPairIntInt::size_type VectorMapIntToPairIntIntSizeType;

typedef std::map <int, VectorPairIntInt> MapIntToVectorPairIntInt;
typedef MapIntToVectorPairIntInt::iterator MapIntToVectorPairIntIntIterator;
typedef MapIntToVectorPairIntInt::const_iterator MapIntToVectorPairIntIntConstIterator;

typedef std::map <int, VectorPairIntString> MapIntToVectorPairIntString;
typedef MapIntToVectorPairIntString::iterator MapIntToVectorPairIntStringIterator;
typedef MapIntToVectorPairIntString::const_iterator MapIntToVectorPairIntStringConstIterator;

typedef std::map <std::string, PairIntDouble> MapStringToPairIntDouble;
typedef MapStringToPairIntDouble::iterator MapStringToPairIntDoubleIterator;
typedef MapStringToPairIntDouble::const_iterator MapStringToPairIntDoubleConstIterator;

typedef std::map <std::string, PairStringString> MapStringToPairStringString;
typedef MapStringToPairStringString::iterator MapStringToPairStringStringIterator;
typedef MapStringToPairStringString::const_iterator MapStringToPairStringStringConstIterator;

typedef std::map <PairIntInt, double> MapPairIntIntToDouble;
typedef MapPairIntIntToDouble::iterator MapPairIntIntToDoubleIterator;
typedef MapPairIntIntToDouble::const_iterator MapPairIntIntToDoubleConstIterator;

typedef std::map <PairIntInt, std::string> MapPairIntIntToString;
typedef MapPairIntIntToString::iterator MapPairIntIntToStringIterator;
typedef MapPairIntIntToString::const_iterator MapPairIntIntToStringConstIterator;

typedef std::map <PairStringInt, double> MapPairStringIntToDouble;
typedef MapPairStringIntToDouble::iterator MapPairStringIntToDoubleIterator;
typedef MapPairStringIntToDouble::const_iterator MapPairStringIntToDoubleConstIterator;

typedef std::map <PairStringString, int> MapPairStringStringToInt;
typedef MapPairStringStringToInt::iterator MapPairStringStringToIntIterator;
typedef MapPairStringStringToInt::const_iterator MapPairStringStringToIntConstIterator;

typedef std::multimap <std::string, int, genNumberStringAscending> MultiMapNumberStringToInt;
typedef MultiMapNumberStringToInt::iterator MultiMapNumberStringToIntIterator;
typedef MultiMapNumberStringToInt::const_iterator MultiMapNumberStringToIntConstIterator;

typedef std::set <PairIntChar> SetPairIntChar;
typedef SetPairIntChar::iterator SetPairIntCharIterator;
typedef std::pair <SetPairIntCharIterator, bool> PairSetPairIntCharIteratorBool;

typedef std::set <PairCharString> SetPairCharString;
typedef SetPairCharString::iterator SetPairCharStringIterator;
typedef std::pair <SetPairCharStringIterator, bool> PairSetPairCharStringIteratorBool;

typedef std::set <PairStringInt> SetPairStringInt;
typedef SetPairStringInt::iterator SetPairStringIntIterator;
typedef std::pair <SetPairStringIntIterator, bool> PairSetPairStringIntIteratorBool;

typedef std::set <PairIntString> SetPairIntString;
typedef SetPairIntString::iterator SetPairIntStringIterator;
typedef std::pair <SetPairIntStringIterator, bool> PairSetPairIntStringIteratorBool;

typedef std::set <PairStringString> SetPairStringString;
typedef SetPairStringString::iterator SetPairStringStringIterator;
typedef std::pair <SetPairStringStringIterator, bool> PairSetPairStringStringIteratorBool;

typedef std::vector <PairIntIntVector> VectorPairIntIntVector;
typedef VectorPairIntIntVector::size_type VectorPairIntIntVectorSizeType;
typedef VectorPairIntIntVector::const_iterator VectorPairIntIntVectorConstIterator;

typedef std::pair <bool, StringVector> PairBoolStringVector;

typedef std::pair <PairStringInt, PairStringInt> PairPairStringIntPairStringInt;

typedef std::map <PairStringString, PairStringString> MapPairStringStringPairStringString;
typedef MapPairStringStringPairStringString::const_iterator MapPairStringStringPairStringStringConstIterator;

typedef std::vector <PairStringInt> VectorPairStringInt;
typedef VectorPairStringInt::const_iterator VectorPairStringIntConstIterator;
typedef VectorPairStringInt::size_type VectorPairStringIntSizeType;

typedef std::vector <VectorPairStringInt> VectorVectorPairStringInt;
typedef VectorVectorPairStringInt::const_iterator VectorVectorPairStringIntConstIterator;
typedef VectorVectorPairStringInt::size_type VectorVectorPairStringIntSizeType;

typedef std::vector <VectorVectorPairStringInt> VectorVectorVectorPairStringInt;
typedef VectorVectorVectorPairStringInt::const_iterator VectorVectorVectorPairStringIntConstIterator;
typedef VectorVectorVectorPairStringInt::size_type VectorVectorVectorPairStringIntSizeType;

typedef std::map <std::string, MapCharToInt> MapStringToMapCharToInt;
typedef MapStringToMapCharToInt::const_iterator MapStringToMapCharToIntConstIterator;

class GenUnique {
	int i;
public:
	GenUnique ( int start = 0 ) { i = start; }
	int operator()() { return i++; }
};

inline PairStringString makePairStringString ( const std::string& s1, const std::string& s2 )
{
	return std::make_pair ( s1, s2 );
}

#endif /* ! __lgen_define */
