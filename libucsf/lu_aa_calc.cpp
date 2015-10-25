/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_aa_calc.cpp                                                *
*                                                                             *
*  Created    : September 17th 2001                                           *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifndef VIS_C
#include <stdexcept>
#endif
#include <lg_string.h>
#include <lgen_error.h>
#include <lu_aa_calc.h>
#include <lu_aa_info.h>
#include <lu_mass_elem.h>
#include <lu_mass_seq.h>
#include <lu_mass.h>
using std::string;
using std::runtime_error;

AACalculator::AACalculator ( bool monoisotopicFlag, const MapStringConstModPtr& constMods ) :
	localMonoisotopicFlag ( monoisotopicFlag ),
	nTermFormula ( "H" ),
	cTermFormula ( "O H" ),
	validTerminalFormula ( true )
{
	static string aa	= "ACDEFGHIKLMNPQRSTUVWXYhmstuvwxy";	/* NB. Rewrite the composition mask code after 32 aa's in this list */
	static string aaFile= "ACDEFGHIKLMNPQRSTUVWLYhmstuvwxy";
//-------------------------1234567890123456789012345678901
	if ( monoisotopicFlag ) massConvert = formula_to_monoisotopic_mass;
	else massConvert = formula_to_average_mass;
	nTerminusWt = massConvert ( nTermFormula );
	cTerminusWt = massConvert ( cTermFormula );

	for ( StringSizeType i = 0 ; i < aa.length () ; i++ ) {
		aminoAcidWt [aa[i]] = monoisotopicFlag	? AAInfo::getInfo ().getMonoisotopicMass ( aaFile [i] )
												: AAInfo::getInfo ().getAverageMass ( aaFile [i] );
		elemFormula [aa[i]] = AAInfo::getInfo ().getElementalString ( aaFile [i] );
	}
	bool validNTermFormula = true;
	bool validCTermFormula = true;
	for ( MapStringConstModPtrConstIterator j = constMods.begin () ; j != constMods.end () ; j++ ) {
		const ConstMod* cMod = (*j).second;
		string aaList = cMod->getAAList ();
		if ( aaList == "n" ) {
			nTermFormula = "H";
			if ( cMod->getValidFormula () ) {
				nTermFormula += cMod->getElementalFormula ();
				nTerminusWt = massConvert ( nTermFormula );
			}
			else {
				validNTermFormula = false;
				nTerminusWt = cMod->getMass () + massConvert ( nTermFormula );
			}
			nTerminusString = cMod->getLongName ();
		}
		else if ( aaList == "c" ) {
			cTermFormula = "O H";
			if ( cMod->getValidFormula () ) {
				cTermFormula += cMod->getElementalFormula ();
				cTerminusWt = massConvert ( cTermFormula );
			}
			else {
				validCTermFormula = false;
				cTerminusWt = cMod->getMass () + massConvert ( cTermFormula );
			}
			cTerminusString = cMod->getLongName ();
		}
		else {
			for ( StringVectorSizeType k = 0 ; k < aaList.length () ; k++ ) {
				char aaCode = aaList [k];
				elemFormula [aaCode] += cMod->getElementalFormula ();
				aminoAcidWt [aaCode] = massConvert ( elemFormula [aaCode] );
			}
		}
	}
	cationFormula = "H";
	cationWt = massConvert ( cationFormula ) - ELECTRON_REST_MASS;

	if ( validNTermFormula && validCTermFormula )
		terminalFormula = nTermFormula + cTermFormula + cationFormula;
	else
		validTerminalFormula = false;
	terminalWt = nTerminusWt + cTerminusWt + cationWt;
}
double AACalculator::calculatePeptideMW ( const string& peptide ) const
{
	double molWt = terminalWt;

	for ( StringConstIterator i = peptide.begin () ; i != peptide.end () ; i++ ) {
		molWt += aminoAcidWt [*i];
	}
	return ( molWt );
}
double AACalculator::calculatePeptideMW ( const StringVector& peptide ) const
{
	double molWt = terminalWt;

	for ( StringVectorConstIterator i = peptide.begin () ; i != peptide.end () ; i++ ) {
		if ( i->length () == 1 )
			molWt += aminoAcidWt [(*i)[0]];
		else
			molWt += getAminoAcidWt ( *i );
	}
	return molWt;
}
bool AACalculator::calculateStrippedElementalCompositionWithTerminii ( const string& peptide, const string& nterm, const string& cterm, const string& nloss, ElementalFormula& ef, bool falseIfMM ) const
{	// This function should only be called by MS-Display and Search Compare
	bool flag = true;
	string copyPeptide = peptide;
	if ( replaceBandZ ( copyPeptide ) ) {
		flag = false;
	}
	string sf;
	StringVector sv = initSequence ( copyPeptide );
	for ( StringVectorSizeType i = 0 ; i < sv.size () ; i++ ) {
		string s = sv [i];
		if ( s.length () == 1 ) {
			sf += s;
		}
		else {
			if ( genIsNumberStart ( s [2] ) ) {	// mass mod
				sf += s [0];
				if ( falseIfMM ) flag = false;
			}
			else 
				sf += s; 
		}
	}
	string nt = nterm;
	if ( !nt.empty () ) {
		if ( genIsNumberStart ( nt [0] ) ) {	// mass mod
			nt = "";
			if ( falseIfMM ) flag = false;
		}
	}
	string ct = cterm;
	if ( !ct.empty () ) {
		if ( genIsNumberStart ( ct [0] ) ) {	// mass mod
			ct = "";
			if ( falseIfMM ) flag = false;
		}
	}
	string nl = nloss;
	if ( !nl.empty () ) {
		if ( genIsNumberStart ( nl [0] ) ) {	// mass mod
			nl = "";
			if ( falseIfMM ) flag = false;
		}
	}
	ef = calculateFragElementalComposition ( initSequence ( sf ) );
	ef += "H";
	if ( !nt.empty () ) {
		ConstMod n ( nt + " (N-term)" );
		if ( !n.getValidFormula () ) throw AACalculatorNoElementalComposition ();
		ef += n.getElementalFormula ();
	}
	ef += "O H";
	if ( !ct.empty () ) {
		ConstMod c ( ct + " (C-term)" );
		if ( !c.getValidFormula () ) throw AACalculatorNoElementalComposition ();
		ef += c.getElementalFormula ();
	}
	ef += cationFormula;
	return flag;
}
ElementalFormula AACalculator::calculateElementalCompositionWithTerminii ( const string& peptide, const string& nterm, const string& cterm, const string& nloss ) const
{
	if ( !validTerminalFormula ) throw AACalculatorNoElementalComposition ();

	ElementalFormula elemFormula = calculateFragElementalComposition ( initSequence ( peptide ) );

	if ( nterm.empty () ) elemFormula += nTermFormula;
	else {
		ConstMod nt ( nterm + " (N-term)" );
		if ( !nt.getValidFormula () ) throw AACalculatorNoElementalComposition ();
		elemFormula += "H";
		elemFormula += nt.getElementalFormula ();
	}
	if ( cterm.empty () ) elemFormula += cTermFormula;
	else {
		ConstMod ct ( cterm + " (C-term)" );
		if ( !ct.getValidFormula () ) throw AACalculatorNoElementalComposition ();
		elemFormula += "O H";
		elemFormula += ct.getElementalFormula ();
	}
	if ( !nloss.empty () ) throw AACalculatorNoElementalComposition ();
	elemFormula += cationFormula;
	return ( elemFormula );
}
ElementalFormula AACalculator::calculateElementalComposition ( const string& peptide ) const
{
	if ( !validTerminalFormula ) throw AACalculatorNoElementalComposition ();

	ElementalFormula elemFormula = calculateFragElementalComposition ( peptide );

	elemFormula += terminalFormula;

	return ( elemFormula );
}
ElementalFormula AACalculator::calculateElementalComposition ( const StringVector& peptide ) const
{
	if ( !validTerminalFormula ) throw AACalculatorNoElementalComposition ();

	ElementalFormula elemFormula = calculateFragElementalComposition ( peptide );

	elemFormula += terminalFormula;

	return ( elemFormula );
}
ElementalFormula AACalculator::calculateFragElementalComposition ( const string& peptide ) const
{
	ElementalFormula ef;

	for ( StringConstIterator i = peptide.begin () ; i != peptide.end () ; i++ ) {
		ef += elemFormula [*i];
	}
	return ( ef );
}
ElementalFormula AACalculator::calculateFragElementalComposition ( const StringVector& peptide ) const
{
	ElementalFormula ef;

	for ( StringVectorConstIterator i = peptide.begin () ; i != peptide.end () ; i++ ) {
		string currentAA = *i;
		if ( currentAA.length () == 1 ) {
			ef += elemFormula [currentAA[0]];
		}
		else {
			try {
				ef += AAInfo::getInfo ().getModElementalString ( currentAA );
			}
			catch ( aaInfoNoElementalFormula ) {
				throw AACalculatorNoElementalComposition ();
			}
		}
	}
	return ( ef );
}
AAFormula AACalculator::calculateAAComposition ( const string& peptide )
{
	int multiplierList [AA_ARRAY_SIZE] = {0};

	for ( StringConstIterator a = peptide.begin () ; a != peptide.end () ; a++ ) {
		multiplierList [*a]++;
	}
	AAFormula aaFormula;
	for ( int i = 'A' ; i <= 'z' ; i++ ) {
		if ( multiplierList [i] ) {
			aaFormula += string ( 1, i ) + gen_itoa ( multiplierList [i] );
		}
	}
	return aaFormula;
}
double AACalculator::getAminoAcidWt ( const string& modAA ) const
{
	if ( localMonoisotopicFlag ) {
		double mass;
		try {
			mass = AAInfo::getInfo ().getModMonoisotopicMass ( modAA );
		}
		catch ( runtime_error e ) {
			ErrorHandler::genError ()->error ( e );
		}
		return mass;
	}
	else {
		return AAInfo::getInfo ().getModAverageMass ( modAA );
	}
}
