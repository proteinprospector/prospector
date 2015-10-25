/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_mass_seq.cpp                                               *
*                                                                             *
*  Created    : July 20th 1996                                                *
*                                                                             *
*  Purpose    : Functions concerned with peptide/protein formulae             *
*               conversions.                                                  *
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
#include <lg_io.h>
#include <lg_string.h>
#include <lgen_error.h>
#include <lu_aa_info.h>
#include <lu_mass_elem.h>
#include <lu_mass_seq.h>
#include <lu_mass.h>
using std::ostream;
using std::string;

ElementalFormula ProteinMW::proteinNTerminus;
ElementalFormula ProteinMW::proteinCTerminus;
DoubleVector ProteinMW::mw_amino_acid_wt;
double ProteinMW::mwTerminalWt;
int ProteinMW::mwPrecision = 0;


void ProteinMW::initialise ( const char* mw_aa, const char* mw_aa_file )
{
	proteinNTerminus = "H";
	proteinCTerminus = "O H";
	mwTerminalWt = formula_to_average_mass ( proteinNTerminus ) + formula_to_average_mass ( proteinCTerminus );
	mw_amino_acid_wt.resize ( AA_ARRAY_SIZE );
	for ( size_t i = 0 ; mw_aa [i] != '\0' ; i++ ) {
		char aa = mw_aa_file[i];
		ElementalFormula ef = AAInfo::getInfo ().getElemental ( aa );
		mw_amino_acid_wt[mw_aa[i]] = formula_to_average_mass ( ef );
	}
}
ProteinMW::ProteinMW ( const char* protein_string )
{
	int multiplier_list [AA_ARRAY_SIZE] = {0};

	while ( *protein_string ) multiplier_list [*protein_string++]++;

	mass = mwTerminalWt;

	for ( int i = 'A' ; i <= 'z' ; i++ ) {
		if ( multiplier_list [i] ) mass += multiplier_list [i] * mw_amino_acid_wt [i];
	}
}
ostream& operator<< ( ostream& os, const Mass& m )
{
	return m.print ( os );
}
ostream& Mass::printMass ( ostream& os, int precision ) const
{
	genPrint ( os, mass, precision, 8 );
	return os;
}
ostream& ProteinMW::print ( ostream& os ) const
{
	printMass ( os, mwPrecision );
	return os;
}
double peptide_formula_to_molecular_weight ( const char* peptideString )
{
	double molWt = terminal_wt;

	while ( *peptideString ) molWt += amino_acid_wt [*peptideString++];

	return molWt;
}
namespace {
void skipTag ( StringSizeType& i, const string& s, char ltag, char rtag )
{
	int ntag = 1;
	i++;
	for ( ; i < s.length () ; i++ ) {
		char a = s [i];
		if ( a == ltag ) ntag++;
		if ( a == rtag ) ntag--;
		if ( ntag == 0 ) return;
	}
}
string stripSequence ( const string& str, char ltag, char rtag )
{
	string newStr;
	bool tag = false;
	for ( StringSizeType i = 0 ; i < str.length () ; i++ ) {
		char s = str [i];
		if ( isdigit ( s ) || isspace ( s ) || s == '.' || s == '+' || s == '\'' || s == '\\' || s == '&' || s == ',' ) {
			// Ignore these
		}
		else if ( s == ltag ) {
			int start = i;
			skipTag ( i, str, ltag, rtag );
			newStr += str.substr ( start, i+1-start );
		}
		else newStr += s;
	}
	return newStr;
}
string checkSequence ( const string& pep, char ltag, char rtag )
{
	string s = gen_strstriptags2 ( pep, ltag, rtag );			// Check the sequence for invalid amino acids
	for ( StringSizeType i = 0 ; i < s.length () ; i++ ) {
		char c = s [i];
		if ( !aa_composition_mask [c] && c != 'B' && c != 'X' && c != 'Z' ) {
			return string ( 1, c );
		}
	}
	return "";
}
string convertToUpperIfLower ( const string& seq, char ltag, char rtag )
{
	if ( seq.empty () ) return "";
	for ( StringSizeType i = 0 ; i < seq.length () ; i++ ) {
		char c = seq [i];
		if ( c == ltag ) {	// Skip to rtag
			skipTag ( i, seq, ltag, rtag );
		}
		else if ( !islower ( c ) ) return seq;		// string is not all lower case return original string
	}
	// String is all lower case
	string newStr;
	for ( StringSizeType j = 0 ; j < seq.length () ; j++ ) {
		char c = seq [j];
		if ( c == ltag ) {
			int start = j;
			skipTag ( j, seq, ltag, rtag );
			newStr += seq.substr ( start, j+1-start );
		}
		else newStr += toupper ( c );
	}
	return newStr;
}
string convertSequence ( const string& s, char ltag, char rtag )
{
	string ret;
	for ( StringSizeType i = 0 ; i < s.length () ;  ) {
		char c = s [i];
		if ( c == ltag ) {
			int start = i;
			skipTag ( i, s, ltag, rtag );
			ret += s.substr ( start, i+1-start );
			i++;
		}
		else {
			char ch = AAInfo::getInfo ().getAA ( s.substr ( i, 3 ) );
			if ( ch == 0 || !aa_composition_mask [ch] ) {
				return "";
			}
			ret += ch;
			i += 3;
		}
	}
	return ret;
}
}
bool replaceBandZ ( string& s )
{
	bool changed = false;
	char ltag = '(';
	char rtag = ')';
	for ( StringSizeType i = 0 ; i < s.length () ; i++ ) {
		char c = s [i];
		if ( c == ltag ) {
			int start = i;
			skipTag ( i, s, ltag, rtag );
		}
		else if ( c == 'B' ) {
			s [i] = 'E';
			changed = true;
		}
		else if ( c == 'Z' ) {
			s [i] = 'Q';
			changed = true;
		}
	}
	return changed;
}
StringVector initSequence ( const string& seq )
{
	static char ltag = '(';
	static char rtag = ')';
	string pep = stripSequence ( seq, ltag, rtag );	// Delete illegal characters
	string invalidAA = checkSequence ( pep, ltag, rtag );		// This is the invalid AA to report
	bool flag = invalidAA.empty ();
	if ( !flag ) {
		pep = convertToUpperIfLower ( pep, ltag, rtag );
		string invAA2 = checkSequence ( pep, ltag, rtag );
		flag = invAA2.empty ();
		if ( !flag ) {
			pep = convertSequence ( pep, ltag, rtag );	// Check if user entered 3 letter code
			flag = !pep.empty ();
		}
	}
	if ( !flag ) {
		ErrorHandler::genError ()->error (
			"The input sequence contains the invalid amino acid: '" + invalidAA + "'" );
	}
	StringVector sv;
	for ( StringSizeType j = 0 ; j < pep.length () ; j++ ) {
		if ( j < pep.length () - 1 && pep [j+1] == '(' ) {
			string cur;
			cur += pep [j++];
			int bracket = 0;
			for ( ; j < pep.length () ; j++ ) {
				char a = pep [j];
				cur += a;
				if ( a == '(' ) bracket++;
				if ( a == ')' ) bracket--;
				if ( bracket == 0 ) break;
			}
			char aa;		// If the mod is a mutation replace it with the amino acid
			if ( ModCodeMap::instance ().getMutationFromPSIMod ( cur.substr ( 2, cur.length () - 3 ), aa ) )
				sv.push_back ( string ( 1, aa ) );
			else
				sv.push_back ( cur );
		}
		else
			sv.push_back ( string ( 1, pep [j] ) );
	}
	return sv;
}
