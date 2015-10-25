/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_usermod.cpp                                                *
*                                                                             *
*  Created    : December 23rd 1998                                            *
*                                                                             *
*  Purpose    : Gets the information on user defined modifications from the   *
*               params/usermod_*.txt files.                                   *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1998-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lg_string.h>
#include <lgen_file.h>
#include <lgen_error.h>
#include <lu_aa_info.h>
#include <lu_getfil.h>
#include <lu_mass.h>
#include <lu_html_form.h>
#include <lu_mass_elem.h>
#include <lu_usermod.h>
#include <lu_charge.h>
#include <lu_delim.h>
#include <lu_param_list.h>
#include <lu_cgi_val.h>
#include <lu_aa_calc.h>
#include <lu_fas_enz.h>
#include <lu_table.h>
#include <lu_srch_form.h>
using std::vector;
using std::string;
using std::ostream;
using std::endl;
using std::fill;
using std::find;
using std::sort;
using std::map;
using std::pair;
using std::set;
using std::make_pair;
using std::lower_bound;
using std::ostringstream;

struct SingleUserMod {
	string name;
	string aaList;
	string longName;
	string formulaStr;
	ElementalFormula formula;
	unsigned int mask;
	char terminalSpecificity;
	string menuString;
	static const char* typeOptions [];
	static const char* filenames [];
	static SetString actualTypeOptions;
	static string phosphoSTYMenuString;
	SingleUserMod ( const string& name, const string& aaList, const string& longName, const string& formulaStr, const ElementalFormula& formula, char terminalSpecificity, const string& typeString ) :
		name ( name ),
		aaList ( aaList ),
		longName ( longName ),
		formulaStr ( formulaStr ),
		formula ( formula ),
		terminalSpecificity ( terminalSpecificity ),
		menuString ( typeToMenuString ( typeString ) )
	{
		mask = string_to_mask ( aaList );
		actualTypeOptions.insert ( menuString );
		if ( name == "Phospho (STY)" ) {
			phosphoSTYMenuString = menuString;
		}
	}
	string getMenuString ()	const { return menuString; }
	string getSpecificity ()
	{
		string::size_type idx1 = name.rfind ( '(' );
		string::size_type idx2 = name.rfind ( ')' );
		return name.substr ( idx1+1, idx2-idx1-1 );
	}
	static string typeToMenuString ( const string& t )
	{
		for ( int i = 0 ; filenames [i] != 0 ; i++ ) {
			if ( t == filenames [i] ) return typeOptions [i];
		}
		return "Unknown";
	}
};
const char* SingleUserMod::typeOptions [] = { "Frequent",
												"Unusual",
												"Glycosylation",
												"Quant: SILAC",
												"Quant: Others",
												"Crosslinking",
												"MS-Product",
												0 };

const char* SingleUserMod::filenames [] = { "frequent",
											"unusual",
											"glyco",
											"silac",
											"quant",
											"xlink",
											"msproduct",
											0 };
SetString SingleUserMod::actualTypeOptions;
string SingleUserMod::phosphoSTYMenuString;


class sortSingleUserModByName {
public:
	int operator () ( const SingleUserMod* a, const SingleUserMod* b ) const
	{
		return genStrcasecmp ( a->name.c_str (), b->name.c_str () ) < 0;
	}
};

bool Usermod::initialised = false;
MapStringToInt Usermod::indexMap;
SingleUserModPtrVector Usermod::singleMods;
StringVector Usermod::names;
bool Usermod::productFlag = false;

Usermod::Usermod ( const string& outputString, const string& formulaStr, char terminalSpecificity, const string& aaList ) :
	outputString ( outputString ),
	formulaStr ( formulaStr ),
	formula ( formulaStr ),
	mask ( string_to_mask ( aaList ) ),
	terminalSpecificity ( terminalSpecificity ),
	aaList ( aaList )
{
}
Usermod::Usermod ( const string& n, bool dummy )		// Use short form of modification
{
	if ( initialised == false ) initialiseUsermod ();
	StringVectorConstIterator iter = lower_bound ( names.begin (), names.end (), n, genStrcasecmpAscending () );
	if ( iter != names.end () ) {
		int index = iter - names.begin ();
		aaList = singleMods [index]->aaList;
		outputString = singleMods [index]->longName;
		formula = singleMods [index]->formula;
		formulaStr = singleMods [index]->formulaStr;
		mask = singleMods [index]->mask;
		terminalSpecificity = singleMods [index]->terminalSpecificity;
		return;
	}
	else {
		std::cout << names [0] << "|" << n << "<br />" << std::endl;
		ErrorHandler::genError ()->error (
			"There is no entry in any of the usermod parameter files for the user defined modification: " + n + ".\n" +
			"Please notify the server administrator about this error message so the user defined modification: " +
			n + " can be used.\n" );
	}
}
Usermod::Usermod ( const string& n )
{
	if ( initialised == false ) initialiseUsermod ();

	MapStringToIntConstIterator cur = indexMap.find ( n );
	if ( cur != indexMap.end () ) {
		int index = cur->second;
		aaList = singleMods [index]->aaList;
		outputString = singleMods [index]->longName;
		formula = singleMods [index]->formula;
		formulaStr = singleMods [index]->formulaStr;
		mask = singleMods [index]->mask;
		terminalSpecificity = singleMods [index]->terminalSpecificity;
		return;
	}
	else {
		ErrorHandler::genError ()->error (
				"There is no entry in any of the usermod parameter files for the user defined modification: " + n + ".\n" +
				"Please notify the server administrator about this error message so the user defined modification: " +
				n + " can be used.\n" );
	}
}
/*
Entry in usermod parameter file:

Acetyl+Oxidation
C2 H2 O2
Protein N-term M

longName: Acetyl+Oxidation
formulaStr: C2 H2 O2
name: Acetyl+Oxidation (Protein N-term M)

*/
void Usermod::initialiseUsermod ()
{
	readUsermodFile ( "usermod.txt" );									// This covers the situation where some mods are still defined in usermod.txt
	for ( int i = 0 ; SingleUserMod::filenames [i] != "msproduct" ; i++ ) {
		readUsermodFile ( "usermod_" + string ( SingleUserMod::filenames [i] ) + ".txt" );
	}
	if ( productFlag ) readUsermodFile ( "usermod_msproduct.txt" );
	sort ( singleMods.begin (), singleMods.end (), sortSingleUserModByName () );
	for ( SingleUserModPtrVectorSizeType j = 0 ; j < singleMods.size () ; j++ ) {
		indexMap [singleMods [j]->name] = j;
	}
	for ( SingleUserModPtrVectorSizeType k = 0 ; k < singleMods.size () ; k++ ) {
		names.push_back ( singleMods [k]->name );
	}
	initialised = true;
}
void Usermod::readUsermodFile ( const string& fName )
{
	string fullPath = MsparamsDir::instance ().getParamPath ( fName );
	if ( !genFileExists ( fullPath ) ) {							// It is not an error for the file not to exist
		return;
	}
	GenCommentedIFStream ifs ( fullPath );
	string line;
	int phase = 1;
	string longName;
	string formulaStr;
	string line3;
	string type;					// If the file is usermod.txt then the type is undefined
	int idx1 = fName.find ( "_" );
	if ( idx1 != string::npos ) {
		int idx2 = fName.find ( "." );
		type = fName.substr ( idx1+1, idx2-idx1-1 );
	}
	else
		type = "unknown";
	while ( ifs.getUncommentedLine ( line ) ) {
		if ( phase == 1 )		longName = line;
		else if ( phase == 2 )	formulaStr = line;
		else if ( phase == 3 )	line3 = line;
		phase++;
		if ( phase == 4 ) {
			string name = longName + " (" + line3 + ")";
			string aaList;
			char terminalSpecificity;
			parseUsermodSpecificityLine ( line3, aaList, terminalSpecificity );
			singleMods.push_back ( new SingleUserMod ( name, aaList, longName, formulaStr, ElementalFormula ( formulaStr ), terminalSpecificity, type ) );
			phase = 1;
		}
	}
}
StringVector Usermod::getNames ( const string& type )
{
	if ( initialised == false ) initialiseUsermod ();
	StringVector sv;
	for ( int i = 0 ; i < singleMods.size () ; i++ ) {
		if ( type == singleMods [i]->menuString ) sv.push_back ( singleMods [i]->name );
	}
	return sv;
}
StringVector Usermod::getMSProdNames ()
{
	bool oldProductFlag = productFlag;
	productFlag = true;
	if ( initialised == false ) initialiseUsermod ();
	StringVector sv;
	SetString nameSet;
	for ( SingleModPtrVectorSizeType i = 0 ; i < singleMods.size () ; i++ ) {
		const SingleUserMod* sm = singleMods [i];
		string aas = sm->aaList;
		string name = sm->longName;
		for ( int j = 0 ; j < aas.size () ; j++ ) {
			char aa = aas [j];
			string n;
			if ( aa != 'n' && aa != 'c' && aa != '.' )	n = string ( 1, aa ) + "(" + name + ")";
			pair <SetStringIterator, bool> flag = nameSet.insert ( n );
			if ( flag.second ) {
				sv.push_back ( n );
			}
		}
	}
	sort ( sv.begin (), sv.end (), genStrcasecmpAscending () );
	productFlag = oldProductFlag;
	return sv;
}
StringVector Usermod::getSCompNames ()
{
	if ( initialised == false ) initialiseUsermod ();
	StringVector sv;
	SetString nameSet;
	for ( SingleModPtrVectorSizeType i = 0 ; i < singleMods.size () ; i++ ) {
		const SingleUserMod* sm = singleMods [i];
		string aas = sm->aaList;
		string name = sm->longName;
		for ( StringSizeType j = 0 ; j < aas.size () ; j++ ) {
			char aa = aas [j];
			string n;
			if ( aa == 'n' )		n = name + " (N-term)";
			else if ( aa == 'c' )	n = name + " (C-term)";
			else if ( aa == '.' )	n = name + " (Neutral loss)";
			else					n = name + " (" + string ( 1, aa ) + ")";
			pair <SetStringIterator, bool> flag = nameSet.insert ( n );
			if ( flag.second ) {
				sv.push_back ( n );
			}
		}
	}
	sort ( sv.begin (), sv.end (), genStrcasecmpAscending () );
	return sv;
}
void Usermod::parseUsermodSpecificityLine ( const string& line, string& aaList, char& terminalSpecificity )
{
	if ( !genStrcasecmp ( line, "Protein N-term" ) ) {
		aaList = "n";
		terminalSpecificity = 'n';
	}
	else if ( !genStrcasecmp ( line, "Protein C-term" ) ) {
		aaList = "c";
		terminalSpecificity = 'c';
	}
	else if ( !genStrcasecmp ( line, "N-term" ) ) {
		aaList = "n";
		terminalSpecificity = 'N';
	}
	else if ( !genStrcasecmp ( line, "C-term" ) ) {
		aaList = "c";
		terminalSpecificity = 'C';
	}
	else if ( !genStrcasecmp ( line, "Neutral loss" ) ) {
		aaList = ".";
		terminalSpecificity = '0';
	}
	else {
		if ( isPrefix ( line, "Protein N-term" ) ) {
			aaList = line.substr ( 15 );
			terminalSpecificity = 'n';
		}
		else if ( isPrefix ( line, "Protein C-term" ) ) {
			aaList = line.substr ( 15 );
			terminalSpecificity = 'c';
		}
		else if ( isPrefix ( line, "N-term" ) ) {
			aaList = line.substr ( 7 );
			terminalSpecificity = 'N';
		}
		else if ( isPrefix ( line, "C-term" ) ) {
			aaList = line.substr ( 7 );
			terminalSpecificity = 'C';
		}
		else if ( isPrefix ( line, "Uncleaved" ) ) {
			aaList = line.substr ( 10 );
			terminalSpecificity = 'e';
		}
		else {
			aaList = line;
			terminalSpecificity = '0';
		}
	}
}
void Usermod::initialiseUsermodAAInfo ( const StringVector& nStr )
{
	if ( initialised == false ) initialiseUsermod ();
	for ( StringVectorSizeType i = 0 ; i < nStr.size () ; i++ ) {
		string n = nStr [i];
		if ( !n.empty () && n [n.length ()-1] != ')' ) {	// MS-Tag, etc can have things like " - Rare" appended
			n = n.substr ( 0, n.find_last_of ( ")" ) + 1 );
		}
		MapStringToIntConstIterator cur = indexMap.find ( n );
		if ( cur != indexMap.end () ) {
			int index = cur->second;
			string aaList = singleMods [index]->aaList;
			string longName = singleMods [index]->longName;
			ElementalFormula& formula = singleMods [index]->formula;
			for ( StringSizeType j = 0 ; j < aaList.length () ; j++ ) {
				AAInfo::getInfo ().addModAminoAcid ( aaList [j], longName, formula );
			}
		}
		else {
			ErrorHandler::genError ()->error (
				"There is no entry in any of the usermod parameter files for the user defined modification: " + n + ".\n" +
				"<p>\n" +
				"Please notify the server administrator about this error message so the user defined modification: " +
				n + " can be used.\n" +
				"</p>\n" );
		}
	}
}
void Usermod::initialiseAllUsermodAAInfo ()
{
	if ( initialised == false ) initialiseUsermod ();
	for ( MapStringToIntConstIterator i = indexMap.begin () ; i != indexMap.end () ; i++ ) {
		int index = i->second;
		string aaList = singleMods [index]->aaList;
		string longName = singleMods [index]->longName;
		const ElementalFormula& formula = singleMods [index]->formula;
		for ( StringSizeType j = 0 ; j < aaList.length () ; j++ ) {
			AAInfo::getInfo ().addModAminoAcid ( aaList [j], longName, formula );
		}
	}
}
void Usermod::getAllUsermodInfo ( StringVector& namesList, StringVector& specificityList, StringVector& elementalFormulaStrList, StringVector& typeList )
{
	bool oldProductFlag = productFlag;
	productFlag = true;
	if ( initialised == false ) initialiseUsermod ();
	for ( MapStringToIntConstIterator i = indexMap.begin () ; i != indexMap.end () ; i++ ) {
		int index = i->second;
		namesList.push_back ( singleMods [index]->longName );
		specificityList.push_back ( singleMods [index]->getSpecificity () );
		elementalFormulaStrList.push_back ( singleMods [index]->formulaStr );
		typeList.push_back ( singleMods [index]->getMenuString () );
	}
	productFlag = oldProductFlag;
}
MapStringToInt Usermod::getN15Balance ()
{
	if ( initialised == false ) initialiseUsermod ();
	MapStringToInt n15Names;
	for ( StringVectorSizeType i = 0 ; i < singleMods.size () ; i++ ) {
		const string& n = singleMods [i]->longName;
		if ( isPrefix ( n, "Label:15N" ) ) {
			int pos = n.find ( "+" );
			if ( pos != string::npos ) {
				ElementalFormula& f = singleMods [i]->formula;
				int multiN = 0;
				int multi15N = 0;
				for ( f.first () ; f.isDone () ; f.next () ) {
					if ( f.element () == "N" ) {
						multiN = f.multiplier ();
					}
					if ( f.element () == "15N" ) {
						multi15N = f.multiplier ();
					}
				}
				int balance = multiN + multi15N;
				if ( balance > 0 ) {
					n15Names.insert ( make_pair ( n.substr ( pos+1 ), balance ) );
				}
			}
		}
	}
	return n15Names;
}
StringVector Usermod::getTypeOptions ()
{
	if ( initialised == false ) initialiseUsermod ();
	StringVector sv;
	for ( SetStringConstIterator i = SingleUserMod::actualTypeOptions.begin () ; i != SingleUserMod::actualTypeOptions.end () ; i++ ) {
		sv.push_back ( *i );
	}
	return sv;
}
string Usermod::getPhosphoSTYMenuString ()
{
	return SingleUserMod::phosphoSTYMenuString;
}
MultipleModification2::MultipleModification2 ( const ParameterList* pList )
{
	pList->getValue ( "mod_AA", modAANames );
	for ( StringVectorSizeType i = 0 ; i < modAANames.size () ; i++ ) {
		userMod.push_back ( new Usermod ( modAANames [i] ) );
	}
}
MultipleModification2::MultipleModification2 ( const ParameterList* pList, bool flag )
{
	StringVector names;
	string user1Name;
	string user2Name;
	string user3Name;
	string user4Name;
	pList->getValue ( "mod_AA", names );
	pList->getValue ( "user1_name", user1Name );
	pList->getValue ( "user2_name", user2Name );
	pList->getValue ( "user3_name", user3Name );
	pList->getValue ( "user4_name", user4Name );

	bool ox = false;
	for ( StringVectorSizeType ii = 0 ; ii < names.size () ; ii++ ) {
		if ( names [ii] == "Oxidation (M)" ) ox = true;
	}

	for ( StringVectorSizeType i = 0 ; i < names.size () ; i++ ) {
		string n = names [i];
		if ( n == "Peptide N-terminal Gln to pyroGlu" )	{
			modAANames.push_back ( "Gln->pyro-Glu (N-term Q)" );
			userMod.push_back ( new Usermod ( modAANames.back () ) );
		}
		else if ( n == "Oxidation of M" ) {
			modAANames.push_back ( "Oxidation (M)" );
			userMod.push_back ( new Usermod ( modAANames.back () ) );
		}
		else if ( n == "Protein N-terminus Acetylated" ) {
			modAANames.push_back ( "Acetyl (Protein N-term)" );
			userMod.push_back ( new Usermod ( modAANames.back () ) );
			modAANames.push_back ( "Met-loss+Acetyl (Protein N-term M)" );
			userMod.push_back ( new Usermod ( modAANames.back () ) );
			modAANames.push_back ( "Met-loss (Protein N-term M)" );
			userMod.push_back ( new Usermod ( modAANames.back () ) );
			if ( ox ) {
				modAANames.push_back ( "Acetyl+Oxidation (Protein N-term M)" );
				userMod.push_back ( new Usermod ( modAANames.back () ) );
			}
		}
		else if ( n == "Acrylamide Modified Cys" ) {
			modAANames.push_back ( "Propionamide (C)" );
			userMod.push_back ( new Usermod ( modAANames.back () ) );
		}
		else if ( n == "User Defined 1" ) {
			modAANames.push_back ( user1Name );
			userMod.push_back ( new Usermod ( user1Name ) );
		}
		else if ( n == "User Defined 2" ) {
			modAANames.push_back ( user2Name );
			userMod.push_back ( new Usermod ( user2Name ) );
		}
		else if ( n == "User Defined 3" ) {
			modAANames.push_back ( user3Name );
			userMod.push_back ( new Usermod ( user3Name ) );
		}
		else if ( n == "User Defined 4" ) {
			modAANames.push_back ( user4Name );
			userMod.push_back ( new Usermod ( user4Name ) );
		}
	}
}
void MultipleModification2::setUserMods ( const AAInitInfo& aaInitInfo ) const
{
	FragModContainer::setUserMods ( userMod, aaInitInfo );
}
void MultipleModification2::printHTML ( ostream& os ) const
{
	if ( modAANames.size () ) {
		os << "Considered modifications: | ";
		for ( StringVectorSizeType i = 0 ; i < modAANames.size () ; i++ ) {
			os << "<b>" << modAANames [i] << "</b> | ";
		}
		os << "<br />" << endl;
	}
}
void MultipleModification2::putCGI ( ostream& os ) const
{
	for ( StringVectorSizeType i = 0 ; i < modAANames.size () ; i++ ) {
		printCGIString ( os, "mod_AA", modAANames [i] );
	}
}
void MultipleModification2::putHiddenFormEntry ( ostream& os ) const
{
	for ( StringVectorSizeType i = 0 ; i < modAANames.size () ; i++ ) {
		printHTMLFORMHidden ( os, "mod_AA", modAANames [i] );
	}
}
void MultipleModification2::putHiddenFormJavascriptEntry ( ostream& os ) const
{
	for ( StringVectorSizeType i = 0 ; i < modAANames.size () ; i++ ) {
		printHTMLFORMJavascriptHidden ( os, "mod_AA", modAANames [i] );
	}
}
void MultipleModification2::copyToCGI ( ostream& os, const ParameterList* params )
{
	StringVector names;
	string user1Name;
	string user2Name;
	string user3Name;
	string user4Name;
	params->getValue ( "mod_AA", names );
	params->getValue ( "user1_name", user1Name );
	params->getValue ( "user2_name", user2Name );
	params->getValue ( "user3_name", user3Name );
	params->getValue ( "user4_name", user4Name );

	bool ox = false;
	for ( StringVectorSizeType ii = 0 ; ii < names.size () ; ii++ ) {
		if ( names [ii] == "Oxidation (M)" ) ox = true;
	}
	for ( StringVectorSizeType i = 0 ; i < names.size () ; i++ ) {
		string n = names [i];
		if ( n == "Peptide N-terminal Gln to pyroGlu" )	{
			printCGIString ( os, "mod_AA", "Gln->pyro-Glu (N-term Q)" );
		}
		else if ( n == "Oxidation of M" ) {
			printCGIString ( os, "mod_AA", "Oxidation (M)" );
		}
		else if ( n == "Protein N-terminus Acetylated" ) {
			printCGIString ( os, "mod_AA", "Acetyl (Protein N-term)" );
			printCGIString ( os, "mod_AA", "Met-loss+Acetyl (Protein N-term M)" );
			printCGIString ( os, "mod_AA", "Met-loss (Protein N-term M)" );
			if ( ox ) {
				printCGIString ( os, "mod_AA", "Acetyl+Oxidation (Protein N-term M)" );
			}
		}
		else if ( n == "Acrylamide Modified Cys" ) {
			printCGIString ( os, "mod_AA", "Propionamide (C)" );
		}
		else if ( n == "User Defined 1" ) {
			printCGIString ( os, "mod_AA", user1Name );
		}
		else if ( n == "User Defined 2" ) {
			printCGIString ( os, "mod_AA", user2Name );
		}
		else if ( n == "User Defined 3" ) {
			printCGIString ( os, "mod_AA", user3Name );
		}
		else if ( n == "User Defined 4" ) {
			printCGIString ( os, "mod_AA", user4Name );
		}
	}
}
MultipleModification::MultipleModification ( const ParameterList* pList )
{
	pList->getValue ( "mod_AA", modAANames );
	pList->getValue ( "user1_name", user1Name );
	pList->getValue ( "user2_name", user2Name );
	pList->getValue ( "user3_name", user3Name );
	pList->getValue ( "user4_name", user4Name );

	pyroglutamicAcidFlag = false;
	oxidationFlag = false;
	acetylationFlag = false;
	incompleteCysFlag = false;
	for ( StringVectorSizeType i = 0 ; i < modAANames.size () ; i++ ) {
		string name = modAANames [i];
		if ( name == "Peptide N-terminal Gln to pyroGlu" )	pyroglutamicAcidFlag = true;
		else if ( name == "Oxidation of M" )				oxidationFlag = true;
		else if ( name == "Protein N-terminus Acetylated" )	acetylationFlag = true;
		else if ( name == "Acrylamide Modified Cys" )		incompleteCysFlag = true;
		else if ( name == "User Defined 1" )				userMod.push_back ( new Usermod ( user1Name ) );
		else if ( name == "User Defined 2" )				userMod.push_back ( new Usermod ( user2Name ) );
		else if ( name == "User Defined 3" )				userMod.push_back ( new Usermod ( user3Name ) );
		else if ( name == "User Defined 4" )				userMod.push_back ( new Usermod ( user4Name ) );
		else												userMod.push_back ( new Usermod ( name ) );
	}
}
void MultipleModification::printHTML ( ostream& os ) const
{
	if ( modAANames.size () ) {
		os << "Considered modifications: | ";
		for ( StringVectorSizeType i = 0 ; i < modAANames.size () ; i++ ) {
			if ( modAANames [i] == "User Defined 1" ) os << "<b>" << user1Name << "</b> | ";
			if ( modAANames [i] == "User Defined 2" ) os << "<b>" << user2Name << "</b> | ";
			if ( modAANames [i] == "User Defined 3" ) os << "<b>" << user3Name << "</b> | ";
			if ( modAANames [i] == "User Defined 4" ) os << "<b>" << user4Name << "</b> | ";
			else os << "<b>" << modAANames [i] << "</b> | ";
		}
		os << "<br />" << endl;
	}
}
void MultipleModification::putCGI ( ostream& os ) const
{
	for ( StringVectorSizeType i = 0 ; i < modAANames.size () ; i++ ) {
		if ( modAANames [i] == "User Defined 1" ) printCGIString ( os, "mod_AA", user1Name );
		if ( modAANames [i] == "User Defined 2" ) printCGIString ( os, "mod_AA", user2Name );
		if ( modAANames [i] == "User Defined 3" ) printCGIString ( os, "mod_AA", user3Name );
		if ( modAANames [i] == "User Defined 4" ) printCGIString ( os, "mod_AA", user4Name );
		else printCGIString ( os, "mod_AA", modAANames [i] );
	}
}
void MultipleModification::putHiddenFormEntry ( ostream& os ) const
{
	for ( StringVectorSizeType i = 0 ; i < modAANames.size () ; i++ ) {
		if ( modAANames [i] == "User Defined 1" ) printHTMLFORMHidden ( os, "mod_AA", user1Name );
		if ( modAANames [i] == "User Defined 2" ) printHTMLFORMHidden ( os, "mod_AA", user2Name );
		if ( modAANames [i] == "User Defined 3" ) printHTMLFORMHidden ( os, "mod_AA", user3Name );
		if ( modAANames [i] == "User Defined 4" ) printHTMLFORMHidden ( os, "mod_AA", user4Name );
		else printHTMLFORMHidden ( os, "mod_AA", modAANames [i] );
	}
}
void MultipleModification::copyToCGI ( ostream& os, const ParameterList* params )
{
	params->copyToCGI ( os, "mod_AA" );
	params->copyToCGI ( os, "user1_name" );
	params->copyToCGI ( os, "user2_name" );
	params->copyToCGI ( os, "user3_name" );
	params->copyToCGI ( os, "user4_name" );
}
void MultipleModification::copyToHiddenFormEntry ( ostream& os, const ParameterList* params )
{
	params->copyToHiddenFormEntry ( os, "mod_AA" );
	params->copyToHiddenFormEntry ( os, "user1_name" );
	params->copyToHiddenFormEntry ( os, "user2_name" );
	params->copyToHiddenFormEntry ( os, "user3_name" );
	params->copyToHiddenFormEntry ( os, "user4_name" );
}
void MultipleModification::copyToHiddenFormJavascriptEntry ( ostream& os, const ParameterList* params )
{
	params->copyToHiddenFormJavascriptEntry ( os, "mod_AA" );
	params->copyToHiddenFormJavascriptEntry ( os, "user1_name" );
	params->copyToHiddenFormJavascriptEntry ( os, "user2_name" );
	params->copyToHiddenFormJavascriptEntry ( os, "user3_name" );
	params->copyToHiddenFormJavascriptEntry ( os, "user4_name" );
}
PossFragMods::PossFragMods ( const string& fragment, bool firstFragment )
{
	pyroglutamicAcid = 0;
	oxidation = 0;
	acetylation = 0;
	incompleteCys = 0;
	user.resize ( numUserMods );
	fill ( user.begin (), user.end (), 0 );
	unsigned int mask = aa_composition_mask [fragment [0]];
	if ( acetylationFlag && firstFragment && ( mask & acetylationMask ) ) acetylation = 1;
	if ( pyroglutamicAcidFlag ) if ( mask & pyroglutamicAcidMask ) pyroglutamicAcid = 1;
	for ( unsigned int i = 0 ; i < fragment.length () ; i++ ) {
		mask = aa_composition_mask [fragment [i]];
		if ( oxidationFlag )	if ( mask & oxidationMask )	oxidation++;
		if ( incompleteCysFlag )if ( mask & cysMask )		incompleteCys++;
		for ( int j = 0 ; j < numUserMods ; j++ ) {
			if ( mask & userMask [j] ) (user [j])++;
		}
	}
}

bool PossFragMods::pyroglutamicAcidFlag = false;
bool PossFragMods::oxidationFlag = false;
bool PossFragMods::acetylationFlag = false;
bool PossFragMods::incompleteCysFlag = false;
unsigned int PossFragMods::acetylationMask;
unsigned int PossFragMods::pyroglutamicAcidMask;
unsigned int PossFragMods::oxidationMask;
unsigned int PossFragMods::cysMask;
UIntVector PossFragMods::userMask;
int PossFragMods::numUserMods;
bool PossFragMods::initialised = false;

void PossFragMods::initialiseMasks ( const MultipleModification& multiMod )
{
	pyroglutamicAcidFlag= multiMod.getPyroglutamicAcidFlag ();
	oxidationFlag		= multiMod.getOxidationFlag ();
	acetylationFlag		= multiMod.getAcetylationFlag ();
	incompleteCysFlag	= multiMod.getIncompleteCysFlag ();

	acetylationMask		= string_to_mask ( "M" );
	pyroglutamicAcidMask= string_to_mask ( "Q" );
	oxidationMask		= string_to_mask ( "M" );
	cysMask				= string_to_mask ( "C" );
	numUserMods = multiMod.getNumUserMods ();
	for ( int i = 0 ; i < numUserMods ; i++ ) {
		userMask.push_back ( multiMod.getUserMask ( i ) );
	}
	initialised = true;
}

unsigned char FragModContainer::protNTermMask = 0x01;
unsigned char FragModContainer::protCTermMask = 0x02;
unsigned char FragModContainer::pepNTermMask = 0x04;
unsigned char FragModContainer::pepCTermMask = 0x08;
unsigned char FragModContainer::nLossMask = 0x10;
unsigned char FragModContainer::anyMask = 0x20;
IntVector FragModContainer::xLinkIdx;
MapCharToUChar FragModContainer::modAAs;
char FragModContainer::enzTermSpec;
double FragModContainer::minMod = 0.0;
double FragModContainer::minModMi = 0.0;
double FragModContainer::minModAv = 0.0;
double FragModContainer::maxMod = 0.0;
double FragModContainer::maxModMi = 0.0;
double FragModContainer::maxModAv = 0.0;
map <char, vector <pair <int, char> > > FragModContainer::modMap;
DoubleVector FragModContainer::userMod;
DoubleVector FragModContainer::userModMi;
DoubleVector FragModContainer::userModAv;
vector <ElementalFormula> FragModContainer::userModElemForm;
StringVector FragModContainer::userModOutputString;

FragModContainer::FragModContainer ( const string& fragment, bool firstFragment, bool lastFragment ) :
	firstFragment ( firstFragment ),
	lastFragment ( lastFragment ),
	nTermIdx ( -1 ),
	cTermIdx ( -1 )
{
	StringSizeType len = fragment.length ();
	for ( StringSizeType i = 0 ; i < len ; i++ ) {
		char aa = fragment [i];
		bool first = ( i == 0 );
		bool last = ( i == len-1 );
		MapCharToUCharConstIterator mIter = modAAs.find ( aa );
		if ( mIter != modAAs.end () ) {
			unsigned char mask = (*mIter).second;
			int flag = mask & anyMask;
			flag |= first && ( mask & pepNTermMask );
			flag |= first && firstFragment && ( mask & protNTermMask );
			flag |= last && ( mask & pepCTermMask );
			flag |= last && lastFragment && ( mask & protCTermMask );

			if ( flag ) {
				aasModified += aa;
				if ( first )	nTermIdx = 0;
				if ( last )		cTermIdx = aasModified.length () - 1;
			}
		}
	}
	MapCharToUCharConstIterator mIter2 = modAAs.find ( 'n' );
	if ( mIter2 != modAAs.end () ) {
		unsigned char mask = (*mIter2).second;
		int flag = mask & pepNTermMask;
		flag |= firstFragment && ( mask & protNTermMask );
		if ( flag ) aasModified += 'n';	// N-term
	}
	MapCharToUCharConstIterator mIter3 = modAAs.find ( 'c' );
	if ( mIter3 != modAAs.end () ) {
		unsigned char mask = (*mIter3).second;
		int flag = mask & pepCTermMask;
		flag |= lastFragment && ( mask & protCTermMask );
		if ( flag ) aasModified += 'c';	// C-term
	}
	if ( modAAs.find ( '.' ) != modAAs.end () ) aasModified += '.';	// Neutral loss
	curNTerm = false;
	curCTerm = false;
	curNLoss = false;
	modList.insert ( curMod );	// The empty mod
	modListSize = 1;
	consideredListSize = 1;
	getNextHit ( 0, 0 );
}
void FragModContainer::getNextHit ( int level, int start )
{
	VectorPairIntChar iv = modMap [aasModified[level]];						// int - mod index, char - term spec
	for ( VectorPairIntCharSizeType ii = start ; ii <= iv.size () ; ii++ ) {
		int i = ii - 1;
		bool termSpecCheck = true;
		bool multiTermCheck = true;
		bool prevNTerm = curNTerm;
		bool prevCTerm = curCTerm;
		bool prevNLoss = curNLoss;
		int index;
		if ( i != -1 ) {
			index = iv [i].first;
			char ts = iv [i].second;
			termSpecCheck = checkTermSpec ( ts, level );
			if ( termSpecCheck ) {
				multiTermCheck = checkMultiTerm ( ts );
				if ( multiTermCheck ) {
					if ( curMod.find ( index ) != curMod.end () )
						curMod [index] += 1;
					else
						curMod [index] = 1;
				}
			}
		}
		if ( termSpecCheck && multiTermCheck ) {
			if ( i != -1 ) {
				if ( !xLinkIdx.empty () ) {		// Consider crosslinks
					int numXL = 0;
					IntVector xLinkSums ( xLinkIdx.size (), 0 );
					MapIntToInt tempMod = curMod;
					for ( IntVectorSizeType x = 0 ; x < xLinkIdx.size () ; x++ ) {
						int idx = xLinkIdx [x];
						MapIntToIntConstIterator cur = tempMod.find ( idx );
						if ( cur != tempMod.end () ) {
							xLinkSums [x] = (*cur).second;
							numXL += (*cur).second;		// Add up crosslinks
							tempMod.erase ( idx );		// Erase this mod
						}
					}
					std::map <MapIntToInt, pair <int, IntVector> >::const_iterator cur = numXLMap.find ( tempMod );
					if ( cur != numXLMap.end () ) {
						if ( (*cur).second.first < numXL ) {
							numXLMap [tempMod].first = numXL;	// Set number of possible crosslinks for this set of mods
							numXLMap [tempMod].second = xLinkSums;	// Set number of possible crosslinks for this set of mods
						}
					}
					else {
						numXLMap [tempMod].first = numXL;
						numXLMap [tempMod].second = xLinkSums;
					}
					pair <std::set <MapIntToInt>::iterator, bool> flag = modList.insert ( tempMod );
					if ( flag.second ) {
						modListSize++;
						if ( modListSize > 50000 ) {
							ErrorHandler::genError ()->error ( "Too many modification combinations.\n" );
						}
					}
				}
				else {
					pair <std::set <MapIntToInt>::iterator, bool> flag = modList.insert ( curMod );
					if ( flag.second ) {
						modListSize++;
						if ( modListSize > 50000 ) {
							ErrorHandler::genError ()->error ( "Too many modification combinations.\n" );
						}
					}
				}
				consideredListSize++;
				if ( consideredListSize > 100000000 ) {		// 100,000,000
					ErrorHandler::genError ()->error ( "Too many modification combinations.\n" );
				}
			}
		}
		if ( level+1 < aasModified.length () ) {
			getNextHit ( level+1, ii == 0 ? 0 : 1 );	// This line causes a bug reported by Yongchang Qui  (see emails starting 2011-01-14. It can be fixed with the line getNextHit ( level+1, 0 ); However this causes issues for peptides with a lot of say Phospho sites
		}
		if ( termSpecCheck && multiTermCheck ) {
			curNTerm = prevNTerm;
			curCTerm = prevCTerm;
			curNLoss = prevNLoss;
			if ( i != -1 ) {
				curMod [index] -= 1;
				if ( curMod [index] == 0 ) {
					curMod.erase ( index );
				}
			}
		}
	}
}
bool FragModContainer::checkMultiTerm ( char ts )		// Checks that there is only a single N or C term modification
{
	if ( ts == 'N' || ts == 'n' ) {
		if ( curNTerm ) return false;
		else {
			curNTerm = true;
			return true;
		}
	}
	if ( ts == 'C' || ts == 'c' ) {
		if ( curCTerm ) return false;
		else {
			curCTerm = true;
			return true;
		}
	}
	if ( ts == '.' ) {
		if ( curNLoss ) return false;
		else {
			curNLoss = true;
			return true;
		}
	}
	return true;
}
bool FragModContainer::checkTermSpec ( char ts, int level )
{
	if ( ts == 'N' ) {									// Peptide N-term
		if ( aasModified [level] == 'n' ) return true;
		if ( level == nTermIdx )	return true;
		else 						return false;
	}
	else if ( ts == 'n' ) {								// Protein N-term
		if ( firstFragment ) {
			if ( aasModified [level] == 'n' ) return true;
			if ( level == nTermIdx )return true;
			else					return false;
		}
		else return false;
	}
	else if ( ts == 'C' ) {								// Peptide C-term
		if ( aasModified [level] == 'c' ) return true;
		if ( level == cTermIdx )	return true;
		else						return false;
	}
	else if ( ts == 'c' ) {								// Protein C-term
		if ( lastFragment ) {
			if ( aasModified [level] == 'c' ) return true;
			if ( level == cTermIdx )return true; 
			else					return false;
		}
		else return false;
	}
	else if ( ts == 'e' ) {
		if ( enzTermSpec == 'C' ) {
			if ( lastFragment )
				return true;
			else
				return ( level != cTermIdx );
		}
		else if ( enzTermSpec == 'N' ) {
			if ( firstFragment )
				return true;
			else
				return ( level != nTermIdx );
		}
	}
	return true;
}
double FragModContainer::getMass ( const MapIntToInt& mii )
{
	double mass = 0.0;
	for ( MapIntToInt::const_iterator i = mii.begin () ; i != mii.end () ; i++ ) {
		mass += (*i).second * userMod [(*i).first];
	}
	return mass;
}
double FragModContainer::getMonoisotopicMass ( const MapIntToInt& mii )
{
	double mass = 0.0;
	for ( MapIntToInt::const_iterator i = mii.begin () ; i != mii.end () ; i++ ) {
		mass += (*i).second * userModMi [(*i).first];
	}
	return mass;
}
double FragModContainer::getAverageMass ( const MapIntToInt& mii )
{
	double mass = 0.0;
	for ( MapIntToInt::const_iterator i = mii.begin () ; i != mii.end () ; i++ ) {
		mass += (*i).second * userModAv [(*i).first];
	}
	return mass;
}
ElementalFormula FragModContainer::getElementalFormula ( const MapIntToInt& mii )
{
	ElementalFormula ef;
	ElementalFormula efTemp;

	for ( MapIntToInt::const_iterator i = mii.begin () ; i != mii.end () ; i++ ) {
		int num = (*i).second;
		if ( num ) {
			efTemp = userModElemForm [(*i).first];
			efTemp.multiply ( num );
			ef += efTemp;
		}
	}
	return ef;
}
bool FragModContainer::getOxidation ( const MapIntToInt& mii )
{
	for ( MapIntToInt::const_iterator i = mii.begin () ; i != mii.end () ; i++ ) {
		if ( userModOutputString [(*i).first].find ( "Oxidation" ) != string::npos ) {
			return true;
		}
	}
	return false;
}
int FragModContainer::getNumOxidation ( const MapIntToInt& mii )
{
	for ( MapIntToInt::const_iterator i = mii.begin () ; i != mii.end () ; i++ ) {
		if ( userModOutputString [(*i).first].find ( "Oxidation" ) != string::npos ) {
			return (*i).second;
		}
	}
	return 0;
}
pair <int, IntVector> FragModContainer::countBridgeSites () const
{
	std::map <MapIntToInt, pair <int, IntVector> >::const_iterator it = numXLMap.find ( *cur );
	if ( it != numXLMap.end () ) {
		return (*it).second;
	}
	return make_pair ( 0, IntVector () );
}
bool FragModContainer::containsMods ( const MapIntToInt& mii, const MapStringToInt& mods )
{
	for ( MapStringToIntConstIterator i = mods.begin () ; i != mods.end () ; i++ ) {
		string m = (*i).first;
		int num = (*i).second;
		StringVectorIterator iter = find ( userModOutputString.begin (), userModOutputString.end (), m );
		if ( iter != userModOutputString.end () ) {
			int index = iter - userModOutputString.begin ();
			MapIntToIntConstIterator iter2 = mii.find ( index );
			if ( iter2 != mii.end () ) {
				int numMods = (*iter2).second;
				if ( numMods < num ) return false;// Not enough mods
			}
			else
				return false;// Mod not in current list
		}
		else
			return false;	// Mod not considered in search
	}
	return true;
}
bool FragModContainer::isExactMod ( const MapIntToInt& mii, const MapStringToInt& mods )
{
	int count = 0;
	for ( MapStringToIntConstIterator i = mods.begin () ; i != mods.end () ; i++ ) {
		string m = (*i).first;
		int num = (*i).second;
		StringVectorIterator iter = find ( userModOutputString.begin (), userModOutputString.end (), m );
		if ( iter != userModOutputString.end () ) {
			int index = iter - userModOutputString.begin ();
			MapIntToIntConstIterator iter2 = mii.find ( index );
			if ( iter2 != mii.end () ) {
				int numMods = (*iter2).second;
				if ( numMods != num ) return false;	// Num mods wrong
				else count++;						// Match
			}
			else
				return false;// Mod not in current list
		}
		else
			return false;	// Mod not considered in search
	}
	if ( count != mii.size () ) return false;		// All mods found but there are more
	return true;
}
bool FragModContainer::getModified ( const MapIntToInt& mii )
{
	return !mii.empty ();
}
void FragModContainer::printHTML ( ostream& os, const MapIntToInt& mii )
{
	if ( mii.empty () ) tableEmptyCell ( os );
	else {
		tableDataStart ( os );
			for ( MapIntToInt::const_iterator i = mii.begin () ; i != mii.end () ; i++ ) {
				os << (*i).second << userModOutputString [(*i).first] << " ";
			}
			os << endl;
		tableDataEnd ( os );
	}
}
void FragModContainer::printDelimited ( ostream& os, const MapIntToInt& mii )
{
	ostringstream ost;
	for ( MapIntToInt::const_iterator i = mii.begin () ; i != mii.end () ; i++ ) {
		ost << (*i).second;
		ost << userModOutputString [(*i).first];
		ost << "; ";
	}
	string s = ost.str ();
	if ( !s.empty () ) s = s.substr ( 0, s.length () - 2 ); // Delete trailing "; "
	delimitedCell ( os, s );
}
void FragModContainer::printXML ( ostream& os, const MapIntToInt& mii )
{
	for ( MapIntToInt::const_iterator i = mii.begin () ; i != mii.end () ; i++ ) {
		ParameterList::printXML ( os, "modification", (*i).second, userModOutputString [(*i).first] );
	}
}
void FragModContainer::setMassType ( bool flag )
{
	if ( flag ) {
		userMod = userModMi;
		minMod = minModMi;
		maxMod = maxModMi;
	}
	else {
		userMod = userModAv;
		minMod = minModAv;
		maxMod = maxModAv;
	}
}
unsigned char FragModContainer::termSpecToMask ( char ts )
{
	if ( ts == '0' ) return anyMask;
	if ( ts == 'e' ) return anyMask;
	if ( ts == 'N' ) return pepNTermMask;
	if ( ts == 'n' ) return protNTermMask;
	if ( ts == 'C' ) return pepCTermMask;
	if ( ts == 'c' ) return protCTermMask;
	if ( ts == '.' ) return nLossMask;
	return 0;
}
void FragModContainer::setUserMods ( const vector <Usermod*>& userMod, const AAInitInfo& aaInitInfo )	// This function can be called multiple times
{
	enzTermSpec = get_enzyme_terminal_specificity ();
	static SetString s;								// This is a set of distinct modification strings
	static map <char, set <pair <int, char> > > mm;	// This is a map from amino acids to a set of integer indicies to the modifications relevant to that amino acid
	static SetInt xx;
	for ( IntVectorSizeType i = 0 ; i < userMod.size () ; i++ ) {	// Get a list of unique mods
		const Usermod* u = userMod [i];
		string outStr = u->getOutputString ();
		bool xlinkFlag = isPrefix ( outStr, "Cross Link " );
		string aaList = u->getAAList ();
		char ts = u->getTerminalSpecificity ();
		if ( aaList == "." ) ts = '.';
		for ( StringSizeType j = 0 ; j < aaList.length () ; j++ ) {
			char aa = aaList [j];
			ConstMod* cm = aaInitInfo.getConstMod ( aa );
			ElementalFormula cmef;
			string fullOutStr = outStr;
			if ( cm ) {
				if ( xlinkFlag ) {
					ErrorHandler::genError ()->error ( "Constant modifications cannot be set for amino acids used in cross links.\n" );
				}
				cmef = cm->getElementalFormula ();
				fullOutStr = cm->getLongName () + "->" + outStr;
			}
			modAAs [aa] |= termSpecToMask ( ts ); 
			pair <SetStringIterator, bool> flag = s.insert ( fullOutStr );
			int idx;
			if ( flag.second ) {		// This is a new mod
				ElementalFormula ef = u->getElementalFormula ();
				ef -= cmef;
				userModElemForm.push_back ( ef );
				userModMi.push_back ( formula_to_monoisotopic_mass ( ef ) );
				userModAv.push_back ( formula_to_average_mass ( ef ) );
				minModMi = genMin ( minModMi, ( ( ts == '0' || ts == 'e' ) ? userModMi.back () * 10 : userModMi.back () ) );	// Allow up to 10 normal mods
				minModAv = genMin ( minModAv, ( ( ts == '0' || ts == 'e' ) ? userModAv.back () * 10 : userModAv.back () ) );
				maxModMi = genMax ( maxModMi, ( ( ts == '0' || ts == 'e' ) ? userModMi.back () * 10 : userModMi.back () ) );	// Allow up to 10 normal mods
				maxModAv = genMax ( maxModAv, ( ( ts == '0' || ts == 'e' ) ? userModAv.back () * 10 : userModAv.back () ) );
				userModOutputString.push_back ( fullOutStr );
				idx = userModOutputString.size () - 1;
			}
			else {
				idx = find ( userModOutputString.begin (), userModOutputString.end (), fullOutStr ) - userModOutputString.begin ();
			}
			mm [aa].insert ( make_pair ( idx, ts ) );
			if ( xlinkFlag ) xx.insert ( idx );
		}
	}
	modMap.clear ();
	for ( std::map <char, set <pair <int, char> > >::const_iterator k = mm.begin () ; k != mm.end () ; k++ ) {
		char aa = (*k).first;
		set <pair <int, char> > si = (*k).second;
		for ( std::set <pair <int, char> >::const_iterator m = si.begin () ; m != si.end () ; m++ ) {
			modMap [aa].push_back ( *m );
		}
	}
	xLinkIdx.clear ();
	for ( SetIntConstIterator p = xx.begin () ; p != xx.end () ; p++ ) {
		xLinkIdx.push_back ( *p );
	}
}
ExtraUserModInfo::ExtraUserModInfo ()
{
}
ExtraUserModInfo::ExtraUserModInfo ( const string& formula, const string& userAMass, const string& limit ) :
	formula ( formula ),
	userAMass ( userAMass ),
	limit ( limit )
{
}
ExtraUserMods::ExtraUserMods ()
{
}
ExtraUserMods& ExtraUserMods::instance ()
{
	static ExtraUserMods e;
	return e;
}
void ExtraUserMods::addUserMods ( const ParameterList* params )		// MS-Product
{
	SetString uniqLabels;
	for ( int i = 1 ; ; i++ ) {
		string sNum = gen_itoa ( i );
		string label = params->getStringValue ( "mod_" + sNum + "_label", "" );
		if ( !label.empty () ) {
			string comp = params->getStringValue ( "mod_" + sNum + "_composition", "" );
			if ( comp.empty () ) {
				ErrorHandler::genError ()->error ( "No elemental composition specified for Mod " + sNum + " Label.\n" );
			}
			else {
				string aaMod = params->getStringValue ( "mod_" + sNum + "_specificity", "" );
				if ( aaMod.empty () ) {
					ErrorHandler::genError ()->error ( "No specificity specified for Mod " + sNum + " Label.\n" );
				}
				pair <SetStringIterator, bool> flag = uniqLabels.insert ( label+aaMod );
				if ( !flag.second ) {
					ErrorHandler::genError ()->error ( "User defined variable modification definition repeated.\n" );
				}
				ExtraUserModInfo eumi ( comp, "", "" );
				extraUserMods [make_pair ( label, aaMod )] = eumi;
				string aaList;
				char terminalSpecificity;
				Usermod::parseUsermodSpecificityLine ( aaMod, aaList, terminalSpecificity );
				for ( int j = 0 ; j < aaList.length () ; j++ ) {
					char aa = aaList [j];
					AAInfo::getInfo ().addModAminoAcid ( aa, label, comp );
				}
			}
		}
		else {
			break;
		}
	}
}
void ExtraUserMods::addUserMods2 ( const ParameterList* params, int maxUserMods )		// MS-Tag. Search-Compare, MS-Viewer
{
	SetString uniqLabels;
	for ( int i = 1 ; i <= maxUserMods ; i++ ) {
		string sNum = gen_itoa ( i );
		string label = params->getStringValue ( "mod_" + sNum + "_label", "" );
		string aMass = params->getStringValue ( "mod_" + sNum + "_accurate_mass", "" );
		if ( !label.empty () || !aMass.empty () ) {
			string comp = params->getStringValue ( "mod_" + sNum + "_composition", "" );
			if ( comp.empty () && aMass.empty () ) {
				ErrorHandler::genError ()->error ( "No elemental composition or accurate mass specified for Mod " + sNum + " Label.\n" );
			}
			else {
				string aaMod = params->getStringValue ( "mod_" + sNum + "_specificity", "" );
				if ( aaMod.empty () ) {
					ErrorHandler::genError ()->error ( "No specificity specified for Mod " + sNum + " Label.\n" );
				}
				pair <SetStringIterator, bool> flag = uniqLabels.insert ( label+aaMod );
				if ( !flag.second ) {
					ErrorHandler::genError ()->error ( "User defined variable modification definition repeated.\n" );
				}
				string limit = params->getStringValue ( "mod_" + sNum + "_limit", "Common" );
				ExtraUserModInfo eumi ( comp, aMass, limit );
				extraUserMods [make_pair ( label, aaMod )] = eumi;
			}
		}
	}
}
void ExtraUserMods::copyToCGI ( ostream& os )
{
	int num = 1;
	for ( MapPairStringStringExtraUserModInfoConstIterator i = extraUserMods.begin () ; i != extraUserMods.end () ; i++ ) {
		string sNum = gen_itoa ( num );
		string form = (*i).second.getFormula ();
		if ( !form.empty () ) {		// Only copy if there is a formula
			printCGIString ( os, FormItemLabel::getName("mod_"+sNum, 1),		(*i).first.first );
			printCGIString ( os, FormItemSpecificity::getName("mod_"+sNum, 1),	(*i).first.second );
			printCGIString ( os, FormItemComposition::getName("mod_"+sNum, 1),	form );
		}
		num++;
	}
}
