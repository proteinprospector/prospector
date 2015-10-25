/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_mass_aa.cpp                                                *
*                                                                             *
*  Created    : July 20th 1996                                                *
*                                                                             *
*  Purpose    : Gets the amino acid information from the params/aa.txt file.  *
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
#ifndef VIS_C
#include <stdexcept>
#endif
#include <lg_string.h>
#include <lu_pi.h>
#include <lu_getfil.h>
#include <lu_mass_elem.h>
#include <lu_aa_info.h>
using std::runtime_error;
using std::string;

AAInfo::AAInfo ()
{
	int numAA;
	char* info = getFileInfo ( MsparamsDir::instance ().getParamPath ( "aa.txt" ), '\n', AMINO_ACID_LINES_PER_ENTRY, true, &numAA );

	for ( int i = 0 ; i < numAA ; i++ ) {
		const char* name = ( i == 0 ) ? strtok ( info, "\n" ) : strtok ( NULL, "\n" );
		const string aaStr = strtok ( NULL, "\n" );
		const char aa = aaStr [0];
		aminoAcids.push_back ( aa );
		string lName = strtok ( NULL, "\n" );
		longNames [aa] = lName;
		longNameToAA [lName] = aa;
		ElementalFormula ef ( strtok ( NULL, "\n" ) );
		elementalFormula [aa] = ef;
		substituentA [aa] = strtok ( NULL, "\n" );
		substituentB [aa] = strtok ( NULL, "\n" );
		pKCTerm [aa] = convert_pk ( strtok ( NULL, "\n" ) );
		pKNTerm [aa] = convert_pk ( strtok ( NULL, "\n" ) );
		pKAcidicSC [aa] = convert_pk ( strtok ( NULL, "\n" ) );
		pKBasicSC [aa] = convert_pk ( strtok ( NULL, "\n" ) );
	}
	for ( int j = 0 ; j < numAA ; j++ ) {	// Must be done in a separate loop because of strtok
		const char aa = aminoAcids [j];
		monoisotopicMass [aa] = formula_to_monoisotopic_mass ( elementalFormula [aa] );
		averageMass [aa] = formula_to_average_mass ( elementalFormula [aa] );
	}
	delete [] info;
}

const int AAInfo::AMINO_ACID_LINES_PER_ENTRY = 10;

AAInfo& AAInfo::getInfo ()
{
	static AAInfo aaInfo;
	return aaInfo;
}
void AAInfo::addAminoAcid ( char aa, const ElementalFormula& formula )
{
	aminoAcids.push_back ( aa );
	longNames [aa] = string ( 1, aa );
	elementalFormula [aa] = formula;
	substituentA [aa] = "H1";
	substituentB [aa] = "0";
	pKCTerm [aa] = NO_PK_VALUE;
	pKNTerm [aa] = NO_PK_VALUE;
	pKAcidicSC [aa] = NO_PK_VALUE;
	pKBasicSC [aa] = NO_PK_VALUE;
	monoisotopicMass [aa] = formula_to_monoisotopic_mass ( elementalFormula [aa] );
	averageMass [aa] = formula_to_average_mass ( elementalFormula [aa] );
}
void AAInfo::addModAminoAcid ( char aa, const string& mod, const ElementalFormula& modFormula )
{
	string modAA ( 1, aa );
	modAA += '(' + mod + ')';
	ElementalFormula ef;
	if ( aa == 'n' ) {
		ef = "H";
		modAA = "term_" + modAA;
	}
	else if ( aa == 'c' ) {
		ef = "O H";
		modAA = "term_" + modAA;
	}
	else if ( aa == '.' || aa == 'B' || aa == 'Z' || aa == 'X' ) {
		ef = "";
	}
	else {
		ef = elementalFormula.find ( aa )->second;
	}
	ef += modFormula;
	modElementalFormula [modAA] = ef;
	modMonoisotopicMass [modAA] = formula_to_monoisotopic_mass ( ef );
	modAverageMass [modAA] = formula_to_average_mass ( ef );
}
double AAInfo::getModMonoisotopicMass ( const string& aa ) const
{
	MapNoCaseStringToDoubleConstIterator i = modMonoisotopicMass.find ( aa );
	if ( i != modMonoisotopicMass.end () ) return i->second;
	else {
		if ( genIsNumberStart ( aa [2] ) ) {
			return monoisotopicMass.find ( aa [0] )->second + atof ( aa.substr ( 2 ).c_str () );
		}
		else if ( gen_strcontains ( aa, '+' ) ) {
			double mod = 0.0;
			if ( getPlusSeparatedMod ( mod, aa, monoisotopicMass, modMonoisotopicMass ) ) {
				return mod;
			}
		}
		string err ( "Undefined amino acid " );
		err += aa;
		err += ".";
		throw runtime_error ( err );
	}
}
double AAInfo::getModAverageMass ( const string& aa ) const
{
	MapNoCaseStringToDoubleConstIterator i = modAverageMass.find ( aa );
	if ( i != modAverageMass.end () ) return i->second;
	else {
		if ( genIsNumberStart ( aa [2] ) ) {
			return averageMass.find ( aa [0] )->second + atof ( aa.substr ( 2 ).c_str () );
		}
		else if ( gen_strcontains ( aa, '+' ) ) {
			double mod = 0.0;
			if ( getPlusSeparatedMod ( mod, aa, averageMass, modAverageMass ) ) {
				return mod;
			}
		}
		std::string err ( "Undefined amino acid " );
		err += aa;
		err += ".";
		throw std::runtime_error ( err );
	}
}
bool AAInfo::getPlusSeparatedMod ( double& mod, const string& aa, const MapCharToDouble& mMap, const MapNoCaseStringToDouble& modMMap ) const
{
	double aaMass = mMap.find ( aa [0] )->second;
	StringVector sv = genGetSubstrings ( aa.substr ( 2, aa.length () - 3 ), '+' );		// Get a list of + separated mods
	string startMapping;				// startMapping is the new amino acid plus ( - eg MappingT would be "T("
	for ( StringVectorSizeType j = 0 ; j < sv.size () ; j++ ) {
		if ( isNoCasePrefix ( sv [j], "Mapping" ) ) {
			char newAA = sv [j][sv[j].length ()-1];
			startMapping = newAA + string ( "(" );
			aaMass = mMap.find ( newAA )->second;
		}
	}
	bool flag = false;
	int nAA = -1;
	for ( StringVectorSizeType k = 0 ; k < sv.size () ; k++ ) {
		string s;
		if ( startMapping.empty () ) {
			s = aa.substr ( 0, 2 ) + sv [k] + ")";
		}
		else if ( isNoCasePrefix ( sv [k], "Mapping" ) ) {								// Don't add the mod
			continue;
		}
		else {
			s = startMapping + sv [k] + ")";
		}
		nAA++;
		MapNoCaseStringToDoubleConstIterator i = modMMap.find ( s );
		if ( i == modMMap.end () ) {
			flag = false;
			break;
		}
		mod += i->second;
		flag = true;
	}
	if ( flag ) {
		mod -= nAA * aaMass;													// Only count the AA mass once
	}
	return flag;
}
double AAInfo::convert_pk ( const string& pkString )
{
	if ( pkString == "n/a" ) return NO_PK_VALUE;
	else return atof ( pkString.c_str () );
}
string AAInfo::convert_pk ( double pk )
{
	if ( pk == NO_PK_VALUE ) return ( "n/a" );
	else return ( gen_ftoa ( pk, "%.2f" ) );
}
StringVector AAInfo::aaListToElementalList ( const string& aaList )
{
	StringVector elementalList;

	for ( StringSizeType i = 0 ; i < aaList.length () ; i++ ) {
		elementalList.push_back ( getElementalString ( aaList [i] ) );
	}
	return ( elementalList );
}
int AAInfo::getUserAAElementNumber ( const string& element, char code )
{
	ElementalFormula ef = getElemental ( code );
	for ( ef.first () ; ef.isDone () ; ef.next () ) {
		if ( ef.element () == element ) return ef.multiplier ();
	}
	return 0;
}
