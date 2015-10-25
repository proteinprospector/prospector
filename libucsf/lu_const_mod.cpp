/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_const_mod.cpp                                              *
*                                                                             *
*  Created    : October 21st 2002                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2002-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lg_string.h>
#include <lgen_error.h>
#include <lgen_file.h>
#include <lgen_math.h>
#include <lu_const_mod.h>
#include <lu_mass.h>
#include <lu_getfil.h>
#include <lu_usermod.h>
#include <algorithm>
using std::string;
using std::sort;
using std::pair;

struct SingleMod {
	string name;
	string aaList;
	string longName;
	string formulaStr;
	ElementalFormula formula;
	SingleMod ( const string& name, const string& aaList, const string& longName, const string& formulaStr, const ElementalFormula& formula ) :
		name ( name ),
		aaList ( aaList ),
		longName ( longName ),
		formulaStr ( formulaStr ),
		formula ( formula ) {}
};
class sortSingleModByName {
public:
	int operator () ( const SingleMod* a, const SingleMod* b ) const
	{
		return genStrcasecmp ( a->name.c_str (), b->name.c_str () ) < 0;
	}
};

bool ConstMod::initialised = false;
MapStringToInt ConstMod::indexMap;
SingleModPtrVector ConstMod::singleMods;
SingleModPtrVector ConstMod::extraSCompSingleMods;
bool ConstMod::productFlag = false;

ConstMod::ConstMod ( const string& n )
{
	if ( initialised == false ) initialiseMods ();

	MapStringToIntConstIterator cur = indexMap.find ( n );
	if ( cur != indexMap.end () ) {
		int index = cur->second;
		aaList = singleMods [index]->aaList;
		longName = singleMods [index]->longName;
		formulaStr = singleMods [index]->formulaStr;
		formula = singleMods [index]->formula;
	}
	else {
		bool flag = false;
		if ( !n.empty () ) {
			if ( genIsNumberStart ( n [0] ) ) {
				int end = n.find ( '(' );
				longName = n.substr ( 0, end-1 );
				mass = atof ( longName.c_str () );
				if ( n.find ( "N-term" ) != string::npos ) aaList = "n";
				if ( n.find ( "C-term" ) != string::npos ) aaList = "c";
				formulaStr = "*";		// formula is unknown
				flag = true;
			}
		}
		if ( !flag )
			ErrorHandler::genError ()->error ( "The modification " + n + " does not occur in any of the usermod parameter files.\n" );
	}
}
void ConstMod::setMass ()
{
	if ( formulaStr != "*" ) {
		mass = massConvert ( formulaStr.c_str () );
	}
}
void ConstMod::initialiseMods () 
{
	SetString nameSet;
	SetString nameSet2;
	readUsermodFile ( "usermod.txt", nameSet, nameSet2 );			// This covers the situation where some mods are still defined in usermod.txt
	readUsermodFile ( "usermod_frequent.txt", nameSet, nameSet2 );
	readUsermodFile ( "usermod_glyco.txt", nameSet, nameSet2 );
	readUsermodFile ( "usermod_xlink.txt", nameSet, nameSet2 );
	readUsermodFile ( "usermod_unusual.txt", nameSet, nameSet2 );
	readUsermodFile ( "usermod_quant.txt", nameSet, nameSet2 );
	readUsermodFile ( "usermod_silac.txt", nameSet, nameSet2 );
	addExtraUsermods ( nameSet, nameSet2 );
	if ( productFlag ) readUsermodFile ( "usermod_msproduct.txt", nameSet, nameSet2 );
	sort ( singleMods.begin (), singleMods.end (), sortSingleModByName () );
	for ( SingleModPtrVectorSizeType j = 0 ; j < singleMods.size () ; j++ ) {
		indexMap [singleMods [j]->name] = j;
	}
	initialised = true;
}
void ConstMod::addExtraUsermods ( SetString& nameSet, SetString& nameSet2 )
{
	MapPairStringStringExtraUserModInfo mseumi = ExtraUserMods::instance ().getExtraUserMods ();
	for ( MapPairStringStringExtraUserModInfoConstIterator i = mseumi.begin () ; i != mseumi.end () ; i++ ) {
		processSingleUsermod ( (*i).first.first, (*i).second.getFormula (), (*i).first.second, nameSet, nameSet2 );
	}
}
void ConstMod::readUsermodFile ( const string& fName, SetString& nameSet, SetString& nameSet2 )
{
	string fullPath = MsparamsDir::instance ().getParamPath ( fName );
	if ( !genFileExists ( fullPath ) ) {							// It is not an error for the file not to exist
		return;
	}
	int numEntries;
	char* info = getFileInfo ( fullPath, '\n', 3, true, &numEntries );
	for ( int i = 0 ; i < numEntries ; i++ ) {
		string longName = ( i == 0 ) ? strtok ( info, "\n" ) : strtok ( NULL, "\n" );
		string formulaStr = strtok ( NULL, "\n" );
		string line3 = strtok ( NULL, "\n" );
		processSingleUsermod ( longName, formulaStr, line3, nameSet, nameSet2 );
	}
	delete [] info;
}
void ConstMod::processSingleUsermod ( const string& longName, const string& formulaStr, const string& line3, SetString& nameSet, SetString& nameSet2 )
{
	string aaList;
	bool sCompMod = false;
	string spec = line3;
	if ( isNoCasePrefix ( spec, "Protein N-term" ) || isNoCasePrefix ( spec, "N-term" ) ) {
		if ( !genStrcasecmp ( spec, "Protein N-term" ) || !genStrcasecmp ( spec, "N-term" ) ) {	// Independent of N-terminal amino acid
			spec = "N-term";
			aaList = "n";
		}
		else {		// This mod only required for Search Compare	// Dependent on N-terminal amino acid
			aaList = spec.substr ( spec.find ( '-' ) + 6 );
			spec = aaList;
			sCompMod = true;
		}
	}
	else if ( isNoCasePrefix ( spec, "Protein C-term" ) || isNoCasePrefix ( spec, "C-term" ) ) {
		if ( !genStrcasecmp ( spec, "Protein C-term" ) || !genStrcasecmp ( spec, "C-term" ) ) {	// Independent of C-terminal amino acid
			spec = "C-term";
			aaList = "c";
		}
		else {		// This mod only required for Search Compare	// Dependent on C-terminal amino acid
			aaList = spec.substr ( spec.find ( '-' ) + 6 );
			spec = aaList;
			sCompMod = true;
		}
	}
	else if ( isNoCasePrefix ( spec, "Uncleaved" ) ) {
		aaList = spec.substr ( 10 );
		sCompMod = true;
	}
	else if ( !genStrcasecmp ( spec, "Neutral loss" ) ) {
		aaList = ".";
		sCompMod = true;
	}
	else {
		aaList = spec;
	}
	string name = longName + " (" + spec + ")";
	if ( sCompMod ) {
		pair <SetStringIterator, bool> flag = nameSet2.insert ( name );
		if ( flag.second ) {
			extraSCompSingleMods.push_back ( new SingleMod ( name, aaList, longName, formulaStr, formulaStr ) );
		}
	}
	else {
		pair <SetStringIterator, bool> flag = nameSet.insert ( name );
		if ( flag.second ) {
			singleMods.push_back ( new SingleMod ( name, aaList, longName, formulaStr, formulaStr ) );
		}
	}
}
StringVector ConstMod::getNames ()
{
	if ( initialised == false ) initialiseMods ();
	StringVector sv;
	for ( SingleModPtrVectorSizeType i = 0 ; i < singleMods.size () ; i++ ) {
		sv.push_back ( singleMods [i]->name );
	}
	return sv;
}
StringVector ConstMod::getNLossNames ()
{
	if ( initialised == false ) initialiseMods ();
	StringVector sv;
	sv.push_back ( "" );
	for ( SingleModPtrVectorSizeType i = 0 ; i < extraSCompSingleMods.size () ; i++ ) {
		if ( extraSCompSingleMods [i]->aaList == "." ) {
			sv.push_back ( extraSCompSingleMods [i]->longName );
		}
	}
	return sv;
}
StringVector ConstMod::getNames ( const string& type )
{
	if ( initialised == false ) initialiseMods ();
	StringVector sv;
	sv.push_back ( "" );
	for ( SingleModPtrVectorSizeType i = 0 ; i < singleMods.size () ; i++ ) {
		if ( singleMods [i]->aaList == type ) {
			sv.push_back ( singleMods [i]->longName );
		}
	}
	return sv;
}
StringVector ConstMod::getNonTerminalNames ()
{
	if ( initialised == false ) initialiseMods ();
	StringVector sv;
	for ( SingleModPtrVectorSizeType i = 0 ; i < singleMods.size () ; i++ ) {
		string aas = singleMods [i]->aaList;
		if ( aas != "n" && aas != "c" ) {
			sv.push_back ( singleMods [i]->name );
		}
	}
	return sv;
}
