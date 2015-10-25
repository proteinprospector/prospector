/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_formula.cpp                                                *
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
#include <lu_formula.h>
#include <lu_html_form.h>
#include <lu_cgi_val.h>
#include <lu_param_list.h>
using std::string;
using std::ostream;
using std::endl;

string ElementalFormulaInfo::ELEMENTAL_FORMULA_INFO_DEFAULT = "C2 H3 N1 O1";

string ElementalFormulaInfo::C = "_c";
string ElementalFormulaInfo::H = "_h";
string ElementalFormulaInfo::N = "_n";
string ElementalFormulaInfo::O = "_o";
string ElementalFormulaInfo::S = "_s";
string ElementalFormulaInfo::P = "_p";
string ElementalFormulaInfo::composition = "_composition";

ElementalFormulaInfo::ElementalFormulaInfo () :
	elementalFormula ( ELEMENTAL_FORMULA_INFO_DEFAULT ),
	elementalFormulaString ( elementalFormula.getFormula () ),
	prefix ( "" )
{
}
ElementalFormulaInfo::ElementalFormulaInfo ( const string& prefix, const ParameterList* params ) :
	prefix ( prefix )
{
	string value;
	if ( params->getValue ( prefix + composition, value ) )
		elementalFormula = value;
	else {
		string str;
		if ( params->getValue ( prefix + C, str ) && str != "0" ) elementalFormula += "C" + str;
		if ( params->getValue ( prefix + H, str ) && str != "0" ) elementalFormula += "H" + str;
		if ( params->getValue ( prefix + N, str ) && str != "0" ) elementalFormula += "N" + str;
		if ( params->getValue ( prefix + O, str ) && str != "0" ) elementalFormula += "O" + str;
		if ( params->getValue ( prefix + S, str ) && str != "0" ) elementalFormula += "S" + str;
		if ( params->getValue ( prefix + P, str ) && str != "0" ) elementalFormula += "P" + str;
	}
	elementalFormulaString = elementalFormula.getFormula ();
}
void ElementalFormulaInfo::putCGI ( ostream& os ) const
{
	if ( elementalFormulaString != "" ) printCGIString ( os, prefix + composition, elementalFormulaString );
}
void ElementalFormulaInfo::putHiddenFormEntry ( ostream& os ) const
{
	if ( elementalFormulaString != "" ) printHTMLFORMHidden ( os, prefix + composition, elementalFormulaString );
}
void ElementalFormulaInfo::putHiddenFormEntryJavascript ( ostream& os ) const
{
	if ( elementalFormulaString != "" ) printHTMLFORMJavascriptHidden ( os, prefix + composition, elementalFormulaString );
}
void ElementalFormulaInfo::printHTML ( ostream& os, const string& label ) const
{
	os << label << ": <b>" << elementalFormulaString << "</b><br />" << endl;
}
#ifdef FINISH
Formula::Formula ( const char* formula )
{
	while ( *formula ) {
		string element;
		while ( isdigit ( *formula ) ) {	// Allows isotope specification
			element += *formula++;
		}
		while ( isalpha ( *formula ) ) {
			element += *formula++;
		}
		elementList.push_back ( element );

		if ( isdigit ( *formula ) || *formula == '-' ) {
			char* next;
			multiplierList.push_back ( strtol ( formula, &next, 10 ) );
			formula += next - formula;			// Wierd construction to allow formula to be const.
		}
		else multiplierList.push_back ( 1 );

		if ( *formula++ != ' ' ) break;
	}
	numElements = elementList.size ();
}

Formula& Formula::operator+= ( const Formula& rhs )
{
	add ( rhs, false );
	return *this;
}

Formula& Formula::operator-= ( const Formula& rhs )
{
	add ( rhs, true );
	return *this;
}
void Formula::add ( const Formula& rhs, bool subtract )
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
	for ( int k = 0, m = 0 ; k < numElements ; k++ ) {
		if ( multiplierList [k] != 0 ) {
			elementList [m] = elementList [k];
			multiplierList [m] = multiplierList [k];
			m++;
		}
	}
	numElements = m;
}
string Formula::getFormula () const
{
	stringstream formula;
		 
 	for ( int i = 0 ; i < numElements ; i++ ) {							
		if ( i != 0 ) {
			formula << " ";
		}
		formula << elementList [i];
		formula << multiplierList [i];
	}
	return ( formula.str () );
}
const Formula operator+ ( const Formula& lhs, const Formula& rhs )
{
	return Formula ( lhs ) += rhs;
}

const Formula operator- ( const Formula& lhs, const Formula& rhs )
{
	return Formula ( lhs ) -= rhs;
}
#endif
