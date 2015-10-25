/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_formula.h                                                  *
*                                                                             *
*  Created    : December 21st 2000                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2000-2008) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_formula_h
#define __lu_formula_h

#include <vector>
#include <string>
#include <sstream>

#include <lgen_define.h>

class ParameterList;

template <class M>
class Formula {
	StringVector elementList;
	typename std::vector <M> multiplierList;
	int numElements;
	void add ( const Formula& rhs, bool subtract );

	StringVectorSizeType cur;
public:
	Formula ();
	Formula ( const std::string& formula );
	Formula ( const char* formula );
	void init ( const char* formula );
	Formula& operator= ( const Formula& rhs );
	Formula& operator+= ( const Formula& rhs );
	Formula& operator-= ( const Formula& rhs );
	void multiply ( M value );
	std::string getFormula ();

	void first () { cur = 0; }
	void next () { cur++; }
	bool isDone () const { return cur < elementList.size (); }
	const std::string& element () const { return elementList [cur];}
	const M multiplier () const { return multiplierList [cur];}
};
template <class M>
std::ostream& operator<< ( std::ostream& os, Formula <M>& f )
{
	bool begin = true;
 	for ( f.first () ; f.isDone () ; f.next () ) {
		if ( begin == false ) {
			os << " ";
		}
		os << f.element ();
		os << f.multiplier ();
		begin = false;
	}
	return os;
}
template <class M>
Formula <M>::Formula ()
{
	init ( "" );
}
template <class M>
Formula <M>::Formula ( const std::string& formula )
{
	init ( formula.c_str () );
}
template <class M>
Formula <M>::Formula ( const char* formula )
{
	init ( formula );
}
template <class M>
void Formula <M>::init ( const char* formula )
{
	while ( *formula ) {
		std::string element;
		while ( isdigit ( *formula ) ) {	// Allows isotope specification
			element += *formula++;
		}
		while ( isalpha ( *formula ) ) {
			element += *formula++;
		}
		elementList.push_back ( element );

		if ( isdigit ( *formula ) || *formula == '-' ) {
			char* next;
			multiplierList.push_back ( (M)strtod ( formula, &next ) );
			formula += next - formula;			// Wierd construction to allow formula to be const.
		}
		else multiplierList.push_back ( 1 );

		if ( *formula++ != ' ' ) break;
	}
	numElements = elementList.size ();
}
template <class M>
Formula <M>& Formula <M>::operator= ( const Formula <M>& rhs )
{
	if ( this != &rhs ) {
		elementList = rhs.elementList;
		multiplierList = rhs.multiplierList;
		numElements = rhs.numElements;
		cur = rhs.cur;
	}
	return *this;
}
template <class M>
Formula <M>& Formula <M>::operator+= ( const Formula& rhs )
{
	add ( rhs, false );
	return *this;
}

template <class M>
Formula <M>& Formula <M>::operator-= ( const Formula& rhs )
{
	add ( rhs, true );
	return *this;
}
template <class M> 
void Formula <M>::multiply ( M value )
{
	if ( value == 0 ) {
		elementList.clear ();
		multiplierList.clear ();
		numElements = 0;
	}
	else {
		for ( int i = 0 ; i < numElements ; i++ ) {
			multiplierList [i] *= value;
		}
	}
}
template <class M> 
void Formula <M>::add ( const Formula <M>& rhs, bool subtract )
{
	for ( int i = 0 ; i < rhs.numElements ; i++ ) {
		bool found = false;
		for ( int j = 0 ; j < numElements ; j++ ) {
			if ( rhs.elementList [i] == elementList [j] ) {
				multiplierList [j] += subtract ? - (rhs.multiplierList [i]) : rhs.multiplierList [i];
				found = true;
				break;
			}
		}
		if ( !found ) {
			elementList.push_back ( rhs.elementList [i] );
			if ( subtract )	multiplierList.push_back ( -rhs.multiplierList [i] );
			else			multiplierList.push_back ( rhs.multiplierList [i] );
			numElements++;
		}
	}
	int m = 0;
	for ( int k = 0 ; k < numElements ; k++ ) {
		if ( multiplierList [k] != 0 ) {
			elementList [m] = elementList [k];
			multiplierList [m] = multiplierList [k];
			m++;
		}
	}
	elementList.resize ( m );
	multiplierList.resize ( m );
	numElements = m;
}
template <class M>
std::string Formula <M>::getFormula ()
{
	std::stringstream formula;
	formula << *this;
	return ( formula.str () );
}
template <class M>
const Formula <M> operator+ ( const Formula <M>& lhs, const Formula <M>& rhs )
{
	return Formula <M> ( lhs ) += rhs;
}
template <class M>
const Formula <M> operator- ( const Formula <M>& lhs, const Formula <M>& rhs )
{
	return Formula <M> ( lhs ) -= rhs;
}


typedef Formula <int> ElementalFormula;
typedef std::vector <ElementalFormula> ElementalFormulaVector;
typedef ElementalFormulaVector::size_type ElementalFormulaVectorSizeType;

typedef Formula <int> AAFormula;

typedef std::map <char, ElementalFormula> MapCharToElementalFormula;
typedef MapCharToElementalFormula::const_iterator MapCharToElementalFormulaConstIterator;

typedef std::map <std::string, ElementalFormula> MapStringToElementalFormula;
typedef MapStringToElementalFormula::const_iterator MapStringToElementalFormulaConstIterator;
typedef MapStringToElementalFormula::iterator MapStringToElementalFormulaIterator;

typedef std::map <std::string, ElementalFormula, genStrcasecmpAscending> MapNoCaseStringToElementalFormula;
typedef MapNoCaseStringToElementalFormula::const_iterator MapNoCaseStringToElementalFormulaConstIterator;
typedef MapNoCaseStringToElementalFormula::iterator MapNoCaseStringToElementalFormulaIterator;


class ElementalFormulaInfo {
	ElementalFormula elementalFormula;
	std::string elementalFormulaString;
	std::string prefix;

	static std::string ELEMENTAL_FORMULA_INFO_DEFAULT;

	static std::string C;
	static std::string H;
	static std::string N;
	static std::string O;
	static std::string S;
	static std::string P;
	static std::string composition;
public:
	ElementalFormulaInfo ();
	ElementalFormulaInfo ( const std::string& prefix, const ParameterList* params );
	void putCGI ( std::ostream& os ) const;
	void putHiddenFormEntry ( std::ostream& os ) const;
	void putHiddenFormEntryJavascript ( std::ostream& os ) const;
	void printHTML ( std::ostream& os, const std::string& label ) const;
	ElementalFormula getElementalFormula () const { return elementalFormula; }
};

typedef std::vector <ElementalFormulaInfo> ElementalFormulaInfoVector;
typedef ElementalFormulaInfoVector::size_type ElementalFormulaInfoVectorSizeType;

#endif /* ! __lu_formula_h */
