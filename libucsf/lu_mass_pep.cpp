/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_mass_pep.cpp                                               *
*                                                                             *
*  Created    : July 20th 1996                                                *
*                                                                             *
*  Purpose    : Initialises some global arrays concerned with mass            *
*               calculations.                                                 *
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
#include <algorithm>
#include <lg_new.h>
#include <lg_string.h>
#define LUCSF_MASS_PEP_MAIN
#include <lu_aa_info.h>
#include <lu_const_mod.h>
#include <lu_mass_elem.h>
#include <lu_mass_seq.h>
#include <lu_mass.h>
#include <lu_html.h>
#include <lu_html_form.h>
#include <lu_inst.h>
#include <lu_mut_mtrx.h>
#include <lu_cgi_val.h>
#include <lu_param_list.h>
using std::string;
using std::vector;
using std::pair;
using std::ostream;
using std::runtime_error;
using std::stringstream;
using std::find;

string ModCodeMap::oneLetterCodes = "ABCDEFGHIKLMNPQRSTUVWXYZ";
string ModCodeMap::lowerCaseCodes = "hmsty";
const char* ModCodeMap::threeLetterCodes [] = {
	"Ala", "B", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys",
	"Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Sec", "Val",
	"Trp", "X", "Tyr", "Z", 0 };
const char* ModCodeMap::lowerCasePSIMods [] = {
	"M(Met->Hsl)", "M(Oxidation)", "S(Phospho)", "T(Phospho)", "Y(Phospho)", 0 };

ModCodeMap::ModCodeMap ()
{
	for ( string::size_type i = 0 ; i < oneLetterCodes.length () ; i++ ) {
		codeMap [oneLetterCodes [i]] = threeLetterCodes [i];
	}
	for ( string::size_type j = 0 ; j < lowerCaseCodes.length () ; j++ ) {
		codeMap2 [lowerCaseCodes [j]] = lowerCasePSIMods [j];
	}
	for ( string::size_type k = 0 ; k < oneLetterCodes.length () ; k++ ) {
		codeMap3 [threeLetterCodes [k]] = oneLetterCodes [k];
	}
}
ModCodeMap& ModCodeMap::instance ()
{
	static ModCodeMap d;
	return d;
}
string ModCodeMap::getModString ( char db, char mut )
{
	string s;
	if ( islower ( mut ) ) {
		s += codeMap2 [mut];
	}
	else {
		s += db;
		s += '(';
		s += codeMap [db];
		s += "->";
		s += codeMap [mut];
		s += ')';
	}
	return s;
}
bool ModCodeMap::getMutationFromPSIMod ( const string& psiMod, char& mutation )
{
	int end = psiMod.find ( "->" );
	if ( end != string::npos ) {
		string first = psiMod.substr ( 0, end );
		MapStringToCharConstIterator cur = codeMap3.find ( first );
		if ( cur != codeMap3.end () ) {
			string::size_type start = end+2;
			string::size_type len = psiMod.length () - start;
			if ( len > 0 ) {
				string second = psiMod.substr ( end+2, len );
				MapStringToCharConstIterator cur = codeMap3.find ( second );
				if ( cur != codeMap3.end () ) {
					mutation = (*cur).second;
					return true;
				}
			}
		}
	}
	return false;
}
void checkConstantAndVariableModCompatibility ( const ParameterList* params )
{
	StringVector constModNames = params->getStringVectorValue ( "const_mod" );
	StringVector userMods = params->getStringVectorValue ( "mod_AA" );
	StringVector msmsUserMods = params->getStringVectorValue ( "msms_mod_AA" );
	for ( int i = 0 ; i < constModNames.size () ; i++ ) {
		string cMod = constModNames [i];
		StringVectorIterator iter = find ( userMods.begin (), userMods.end (), cMod );
		StringVectorIterator iter2 = find ( msmsUserMods.begin (), msmsUserMods.end (), cMod );
		if ( iter != userMods.end () || iter2 != msmsUserMods.end () ) {
			throw runtime_error ( cMod + " is both a constant and variable modification." );
		}
	}
}
static int nom_amino_acid_wt [AA_ARRAY_SIZE];

AAInitInfo::AAInitInfo ( const ParameterList* params ) :
	cation			( "H" ),
	instrumentName	( params->getStringValue ( "instrument_name", "" ) )
{
	params->getValue ( "const_mod", constModNames );
	try {
		checkConstantAndVariableModCompatibility ( params );
	}
	catch ( runtime_error e ) {
		ErrorHandler::genError ()->error ( e );
	}
	SetChar aasUsed;
	for ( StringVectorSizeType i = 0 ; i < constModNames.size () ; i++ ) {
		ConstMod* cMod = new ConstMod ( constModNames [i] );
		string aaList = cMod->getAAList ();
		for ( StringVectorSizeType j = 0 ; j < aaList.length () ; j++ ) {
			char aa = aaList [j];
			pair <SetCharIterator, bool> flag = aasUsed.insert ( aa );
			if ( flag.second == false ) {
				ErrorHandler::genError ()->error (
					"Only a single constant modification can be chosen for each amino acid or terminus.\n" );
			}
		}
		constMods [aaList] = cMod;
		for ( int k = 0 ; k < aaList.length () ; k++ ) {
			cm [aaList[k]] = cMod;
		}
	}
	nTermName = params->getStringValue ( "nterm", "" );
	if ( !nTermName.empty () ) {
		ConstMod* c = new ConstMod ( nTermName + " (N-term)" );
		constMods ["n"] = c;
		cm ['n'] = c;
	}

	cTermName = params->getStringValue ( "cterm", "" );
	if ( !cTermName.empty () ) {
		ConstMod* c = new ConstMod ( cTermName + " (C-term)" );
		constMods ["c"] = c;
		cm ['c'] = c;
	}

	nLossName = params->getStringValue ( "nloss", "" );

	userAAElemForms.push_back ( ElementalFormulaInfo ( "user_aa", params ) );
	userAAElemForms.push_back ( ElementalFormulaInfo ( "user_aa_2", params ) );
	userAAElemForms.push_back ( ElementalFormulaInfo ( "user_aa_3", params ) );
	userAAElemForms.push_back ( ElementalFormulaInfo ( "user_aa_4", params ) );

	initialise_amino_acid_weights ( constMods, getUserAAElemForms (), true );
	initialiseInstrumentName ( instrumentName );

	for ( MapStringConstModPtrConstIterator k = constMods.begin () ; k != constMods.end () ; k++ ) {
		const ConstMod* cMod = (*k).second;
		string aaList = cMod->getAAList ();
		for ( StringVectorSizeType m = 0 ; m < aaList.length () ; m++ ) {
			PeptideSequence::addConstMod ( aaList [m], cMod->getLongName () );
		}
	}
}
MapCharToString AAInitInfo::getConstModMap () const
{
	MapCharToString mcs;
	for ( MapStringConstModPtrConstIterator i = constMods.begin () ; i != constMods.end () ; i++ ) {
		const ConstMod* cm = (*i).second;
		string aaList = cm->getAAList ();
		for ( int j = 0 ; j < aaList.size () ; j++ ) {
			mcs [aaList [j]] = cm->getLongName ();
		}
	}
	return mcs;
}
void AAInitInfo::putCGI ( ostream& os ) const
{
	for ( StringVectorSizeType j = 0 ; j < constModNames.size () ; j++ ) {
		printCGIString ( os, "const_mod", constModNames [j] );
	}
	printCGIString ( os, "instrument_name", instrumentName );
	for ( ElementalFormulaInfoVectorSizeType i = 0 ; i < userAAElemForms.size () ; i++ ) {
		userAAElemForms [i].putCGI ( os );
	}
}
void AAInitInfo::putHiddenFormEntryJavascript ( ostream& os ) const
{
	for ( StringVectorSizeType j = 0 ; j < constModNames.size () ; j++ ) {
		printHTMLFORMJavascriptHidden ( os, "const_mod", constModNames [j] );
	}
	printHTMLFORMJavascriptHidden ( os, "instrument_name", instrumentName );
	for ( ElementalFormulaInfoVectorSizeType i = 0 ; i < userAAElemForms.size () ; i++ ) {
		userAAElemForms [i].putHiddenFormEntryJavascript ( os );
	}
}
void AAInitInfo::copyToCGI ( ostream& os, const ParameterList* params, bool cysUnmodified )
{
	StringVector cModNames;
	params->getValue ( "const_mod", cModNames );
	if ( cysUnmodified ) {
		StringVectorIterator i = cModNames.begin ();
		for ( ; i != cModNames.end () ; i++ ) {
			ConstMod cm ( *i );
			if ( cm.getAAList ().find ( 'C' ) != string::npos ) break;
		}
		if ( i != cModNames.end () ) cModNames.erase ( i );
	}
	for ( StringVectorSizeType j = 0 ; j < cModNames.size () ; j++ ) {
		printCGIString ( os, "const_mod", cModNames [j] );
	}
	params->copyToCGI ( os, "instrument_name" );
	params->copyToCGI ( os, "user_aa_composition" );
	params->copyToCGI ( os, "user_aa_2_composition" );
	params->copyToCGI ( os, "user_aa_3_composition" );
	params->copyToCGI ( os, "user_aa_4_composition" );
}
void AAInitInfo::copyToHiddenFormEntry ( ostream& os, const ParameterList* params, bool cysUnmodified )
{
	StringVector cModNames;
	params->getValue ( "const_mod", cModNames );
	if ( cysUnmodified ) {
		StringVectorIterator i = cModNames.begin ();
		for ( ; i != cModNames.end () ; i++ ) {
			ConstMod cm ( *i );
			if ( cm.getAAList ().find ( 'C' ) != string::npos ) break;
		}
		if ( i != cModNames.end () ) cModNames.erase ( i );
	}
	for ( StringVectorSizeType j = 0 ; j < cModNames.size () ; j++ ) {
		printHTMLFORMHidden ( os, "const_mod", cModNames [j] );
	}
	params->copyToHiddenFormEntry ( os, "instrument_name" );
	params->copyToHiddenFormEntry ( os, "user_aa_composition" );
	params->copyToHiddenFormEntry ( os, "user_aa_2_composition" );
	params->copyToHiddenFormEntry ( os, "user_aa_3_composition" );
	params->copyToHiddenFormEntry ( os, "user_aa_4_composition" );
}
void AAInitInfo::copyToHiddenFormJavascriptEntry ( ostream& os, const ParameterList* params, bool cysUnmodified )
{
	StringVector cModNames;
	params->getValue ( "const_mod", cModNames );
	if ( cysUnmodified ) {
		StringVectorIterator i = cModNames.begin ();
		for ( ; i != cModNames.end () ; i++ ) {
			ConstMod cm ( *i );
			if ( cm.getAAList ().find ( 'C' ) != string::npos ) break;
		}
		if ( i != cModNames.end () ) cModNames.erase ( i );
	}
	for ( StringVectorSizeType j = 0 ; j < cModNames.size () ; j++ ) {
		printHTMLFORMJavascriptHidden ( os, "const_mod", cModNames [j] );
	}
	params->copyToHiddenFormJavascriptEntry ( os, "instrument_name" );
	params->copyToHiddenFormJavascriptEntry ( os, "user_aa_composition" );
	params->copyToHiddenFormJavascriptEntry ( os, "user_aa_2_composition" );
	params->copyToHiddenFormJavascriptEntry ( os, "user_aa_3_composition" );
	params->copyToHiddenFormJavascriptEntry ( os, "user_aa_4_composition" );
}
void AAInitInfo::printHTML ( ostream& os ) const
{
	for ( ElementalFormulaInfoVectorSizeType i = 0 ; i < userAAElemForms.size () ; i++ ) {
		if ( userAAElemForms [i].getElementalFormula ().getFormula () != "" ) {
			stringstream label;
			label << "User AA Formula " << i+1;
			userAAElemForms [i].printHTML ( os, label.str () );
		}
	}
	for ( StringVectorSizeType j = 0 ; j < constModNames.size () ; j++ ) {
		ParameterList::printHTML ( os, "Constant Modification", constModNames [j] );
	}
}
ElementalFormula AAInitInfo::getNTermModFormula () const
{
	MapStringConstModPtrConstIterator cur = constMods.find ( "n" );
	if ( cur != constMods.end () ) return (*cur).second->getElementalFormula ();
	else return "";
}
ElementalFormula AAInitInfo::getCTermModFormula () const
{
	MapStringConstModPtrConstIterator cur = constMods.find ( "c" );
	if ( cur != constMods.end () ) return (*cur).second->getElementalFormula ();
	else return "";
}

string MassInfo::MASS_TYPE = "parent_mass_convert";

string MassInfo::MASS_TYPE_DEFAULT = "monoisotopic";
bool MassInfo::MONOISOTOPIC_FLAG_DEFAULT = true;

MassInfo::MassInfo ( const ParameterList* params ) :
	massType		( params->getStringValue	( MASS_TYPE, MASS_TYPE_DEFAULT ) ),
	monoisotopicFlag( massType != "average" )
{
	modified_mass_convert ( monoisotopicFlag );
}
void MassInfo::putCGI ( ostream& os ) const
{
	printCGIString ( os, MASS_TYPE, massType );
}
void MassInfo::putCGIFragment ( ostream& os ) const
{
	if ( massType == "Par(mi)Frag(av)" )		printCGIString ( os, MASS_TYPE, "average" );
	else if ( massType == "Par(av)Frag(mi)" )	printCGIString ( os, MASS_TYPE, "monoisotopic" );
	else										printCGIString ( os, MASS_TYPE, massType );
}
void MassInfo::putHiddenFormEntry ( ostream& os ) const
{
	printHTMLFORMHidden ( os, MASS_TYPE, massType );
}
void MassInfo::copyToCGI ( ostream& os, const ParameterList* params )
{
	params->copyToCGI ( os, MASS_TYPE );
}
void MassInfo::copyToHiddenFormEntry ( ostream& os, const ParameterList* params )
{
	params->copyToHiddenFormEntry ( os, MASS_TYPE );
}
void MassInfo::copyToHiddenFormJavascriptEntry ( ostream& os, const ParameterList* params )
{
	params->copyToHiddenFormJavascriptEntry ( os, MASS_TYPE );
}
void MassInfo::printHTML ( ostream& os ) const
{
	os << "Masses are <b>" << massType << "</b><br />\n";
}

static void initialise ( const vector <ElementalFormula>& userAAElemFormula );
static void init_aa_arrays ();

static bool localMonoisotopicFlag;
static MapStringConstModPtr localConstMods;
static vector <ElementalFormula> localUserAAElemFormulae;
static int initialised = 0;
static double mi_amino_acid_wt [AA_ARRAY_SIZE];
static double av_amino_acid_wt [AA_ARRAY_SIZE];
static double mi_substituent_da [AA_ARRAY_SIZE];
static double mi_substituent_wa [AA_ARRAY_SIZE];
static double av_substituent_da [AA_ARRAY_SIZE];
static double av_substituent_wa [AA_ARRAY_SIZE];
static double mi_substituent_db [AA_ARRAY_SIZE];
static double mi_substituent_wb [AA_ARRAY_SIZE];
static double av_substituent_db [AA_ARRAY_SIZE];
static double av_substituent_wb [AA_ARRAY_SIZE];

void initialise_amino_acid_weights ( const MapStringConstModPtr& constMods, const vector <ElementalFormula>& userAAElemFormulae, bool monoisotopicFlag )
{
	if ( monoisotopicFlag ) massConvert = formula_to_monoisotopic_mass;
	else massConvert = formula_to_average_mass;
	if ( initialised == 0 ) initialise ( userAAElemFormulae );
	localConstMods = constMods;
	localUserAAElemFormulae = userAAElemFormulae;
	localMonoisotopicFlag = monoisotopicFlag;
	h1 = massConvert ( "H" );
	bool validNTerminusFormula = true;
	ElementalFormula nTermElemForm ( "H" );
	n_terminus_wt = h1;
	MapStringConstModPtrConstIterator cur = constMods.find ( "n" );
	if ( cur != constMods.end () ) {
		const_cast <ConstMod*> ((*cur).second)->setMass ();
		n_terminus_wt += (*cur).second->getMass ();
		validNTerminusFormula = (*cur).second->getValidFormula ();
		if ( validNTerminusFormula ) {
			nTermElemForm += (*cur).second->getElementalFormula ();
		}
	}
	o_h =  massConvert ( "O H" );
	bool validCTerminusFormula = true;
	ElementalFormula cTermElemForm ( "O H" );
	c_terminus_wt = o_h;
	MapStringConstModPtrConstIterator cur2 = constMods.find ( "c" );
	if ( cur2 != constMods.end () ) {
		const_cast <ConstMod*> ((*cur2).second)->setMass ();
		c_terminus_wt += (*cur2).second->getMass ();
		validCTerminusFormula = (*cur2).second->getValidFormula ();
		if ( validCTerminusFormula ) {
			cTermElemForm += (*cur2).second->getElementalFormula ();
		}
	}
	string cation = "H";
	cation_wt = massConvert ( cation.c_str () ) - ELECTRON_REST_MASS;

	terminal_wt = n_terminus_wt + c_terminus_wt + cation_wt;
	if ( validNTerminusFormula && validCTerminusFormula ) {
		ElementalFormula tempFormula = nTermElemForm;
		tempFormula += cTermElemForm;
		tempFormula += cation;
		terminal_formula = gen_new_string ( tempFormula.getFormula ().c_str () );
	}
	internal_n_terminus_wt = massConvert ( "H" );
	if ( monoisotopicFlag ) {
		amino_acid_wt = mi_amino_acid_wt;
		substituent_da = mi_substituent_da;
		substituent_wa = mi_substituent_wa;
		substituent_db = mi_substituent_db;
		substituent_wb = mi_substituent_wb;
		nominal_mass_offset = MONOISOTOPIC_NOMINAL_MASS_OFFSET;
	}
	else {
		amino_acid_wt = av_amino_acid_wt;
		substituent_da = av_substituent_da;
		substituent_wa = av_substituent_wa;
		substituent_db = av_substituent_db;
		substituent_wb = av_substituent_wb;
		nominal_mass_offset = AVERAGE_NOMINAL_MASS_OFFSET;
	}
	for ( int j = 0 ; j < AA_ARRAY_SIZE ; j++ ) {
		aaArrayMT.push_back ( static_cast <MassType> (amino_acid_wt [j] * MASS_TYPE_MULTIPLIER) );
		substituentDaMT.push_back ( static_cast <MassType> (substituent_da [j] * MASS_TYPE_MULTIPLIER) );
		substituentDbMT.push_back ( static_cast <MassType> (substituent_db [j] * MASS_TYPE_MULTIPLIER) );
		substituentWaMT.push_back ( static_cast <MassType> (substituent_wa [j] * MASS_TYPE_MULTIPLIER) );
		substituentWbMT.push_back ( static_cast <MassType> (substituent_wb [j] * MASS_TYPE_MULTIPLIER) );
	}
	cnbr_homoserine_lactone_mod = amino_acid_wt ['h'] - amino_acid_wt ['M'];

	h2_o = massConvert ( "H2 O" );
	n_h3 = massConvert ( "N H3" );
	h2 = massConvert ( "H2" );
	c_o = massConvert ( "C O" );

	double c13 = 13.003354838 - 12.0;

	a_tag_offset			= c_o;
	a_h2o_tag_offset		= massConvert ( "C O2 H2" );
	a_nh3_tag_offset		= massConvert ( "C O N H3" );
	a_h3po4_tag_offset		= massConvert ( "H3 P O4" ) + massConvert ( "C O" );
	b_tag_offset			= 0.0;
	b_h2o_tag_offset		= h2_o;
	b_plus_h2o_tag_offset	= - h2_o;
	b_nh3_tag_offset		= n_h3;
	b_soch4_tag_offset		= massConvert ( "S O C H4" );
	b_h3po4_tag_offset		= massConvert ( "H3 P O4" );
	c_tag_offset			= - n_h3;
	cPlus1DaTagOffset		= c_tag_offset - c13;
	cPlus2DaTagOffset		= cPlus1DaTagOffset - c13;
	cMinus1DaTagOffset		= c_tag_offset + h1;
	d_tag_offset			= 0.0;

	c2_h2_n_o = massConvert ( "C2 H2 N O" );
	s_o_c_h4 = massConvert ( "S O C H4" );
	h3_p_o4 = massConvert ( "H3 P O4" );
	double c2_h3_n_o = massConvert ( "C2 H3 N O" );
	v_tag_offset		= - c2_h3_n_o;
	w_tag_offset		= - h2;
	x_tag_offset		= - c_o;
	y_tag_offset		= - h2;
	y_nh3_tag_offset	= n_h3 - h2;
	y_h2o_tag_offset	= h2_o - h2;
	y_soch4_tag_offset	= massConvert ( "S O C H4" ) - h2;
	y_h3po4_tag_offset	= massConvert ( "H3 P O4" ) - h2;
	Y_tag_offset		= 0.0;
	z_tag_offset		= n_h3 - h1 - h2;
	zPlus1DaTagOffset	= z_tag_offset - h1;
	zPlus2DaTagOffset	= zPlus1DaTagOffset - c13;
	zPlus3DaTagOffset	= zPlus2DaTagOffset - c13;

	for ( MapStringConstModPtrConstIterator i = constMods.begin () ; i != constMods.end () ; i++ ) {
		const ConstMod* cMod = (*i).second;
		string aaList = cMod->getAAList ();
		if ( aaList == "n" || aaList == "c" ) continue;
		for ( StringVectorSizeType j = 0 ; j < aaList.length () ; j++ ) {
			char aaCode = aaList [j];
			amino_acid_wt [aaCode] = massConvert ( AAInfo::getInfo ().getElementalString ( aaCode ).c_str () ) + massConvert ( cMod->getElementalFormula ().getFormula ().c_str () );
			nom_amino_acid_wt [aaCode] = (int) amino_acid_wt [aaCode];

			double substituent = massConvert ( "H1" );
			double n_c2_h4 = massConvert ( "N C2 H4" );
			substituent_da [aaCode] = substituent ? substituent + n_c2_h4 - amino_acid_wt [aaCode] : 0.0;
			double c3_h_o = massConvert ( "C3 H O" );
			substituent_wa [aaCode] = substituent ? substituent + c3_h_o - amino_acid_wt [aaCode] : 0.0;
			substituent = massConvert ( "0" );
			substituent_db [aaCode] = substituent ? substituent + n_c2_h4 - amino_acid_wt [aaCode] : 0.0;
			substituent_wb [aaCode] = substituent ? substituent + c3_h_o - amino_acid_wt [aaCode] : 0.0;

			aaArrayMT [aaCode] = static_cast <MassType> (amino_acid_wt [aaCode] * MASS_TYPE_MULTIPLIER);
			substituentDaMT [aaCode] = static_cast <MassType> (substituent_da [aaCode] * MASS_TYPE_MULTIPLIER);
			substituentDbMT [aaCode] = static_cast <MassType> (substituent_db [aaCode] * MASS_TYPE_MULTIPLIER);
			substituentWaMT [aaCode] = static_cast <MassType> (substituent_wa [aaCode] * MASS_TYPE_MULTIPLIER);
			substituentWbMT [aaCode] = static_cast <MassType> (substituent_wb [aaCode] * MASS_TYPE_MULTIPLIER);
		}
	}
	phosphorylation_mask	= string_to_mask ( "st" );
	oxidized_m_mask			= string_to_mask ( "m" );
}
void modified_mass_convert ( bool monoisotopicFlag )
{
	initialise_amino_acid_weights ( localConstMods, localUserAAElemFormulae, monoisotopicFlag );
}
unsigned int string_to_mask ( const string& str )
{
	unsigned int mask = 0;

	for ( StringSizeType i = 0 ; i < str.length () ; i++ ) {
		mask |= aa_composition_mask [str [i]];
	}
	return mask;
}
unsigned int string_to_mask ( const char* str, int len )
{
	unsigned int mask = 0;

	for ( int i = 0 ; i < len ; i++ ) {
		mask |= aa_composition_mask [str [i]];
	}
	return mask;
}
int sumBasicResidues ( const string& fragment )
{
	int num = 0;

	for ( StringSizeType i = 0 ; i < fragment.length () ; i++ ) {
		if ( aa_composition_mask [fragment [i]] & instInf->getPosChargeBearingMask () ) num++; 
	}
	return num;
}
int sumBasicResidues ( const StringVector& fragment )
{
	int num = 0;

	for ( StringVectorSizeType i = 0 ; i < fragment.size () ; i++ ) {
		const string& f = fragment [i];
		if ( aa_composition_mask [f[0]] & instInf->getPosChargeBearingMask () ) num++;
		else {
			if ( f.size () > 8 && isPrefix ( f.substr ( 2 ), "Cation" ) ) num++;
		}
	}
	return num;
}
int sumResidues ( const string& fragment, unsigned int mask )
{
	int num = 0;

	for ( StringSizeType i = 0 ; i < fragment.length () ; i++ ) {
		if ( aa_composition_mask [fragment [i]] & mask ) num++; 
	}
	return num;
}
static void initialise ( const vector <ElementalFormula>& userAAElemFormula )
{
	static char mw_aa []		= "ABCDEFGHIKLMNPQRSTUVWXYZstuvwxy";
	static char mw_aa_file []	= "ADCDEFGHIKLMNPQRSTUVWLYEstuvwxy";
	unsigned int i;

	AAInfo::getInfo ().addAminoAcid ( 'u', userAAElemFormula [0] );
	AAInfo::getInfo ().addAminoAcid ( 'v', userAAElemFormula [1] );
	AAInfo::getInfo ().addAminoAcid ( 'w', userAAElemFormula [2] );
	AAInfo::getInfo ().addAminoAcid ( 'x', userAAElemFormula [3] );

	init_aa_arrays ();

	ProteinMW::initialise ( mw_aa, mw_aa_file );

	for ( i = 0 ; i < mw_aa [i] != '\0' ; i++ ) {
		char aa = mw_aa_file[i];
		pk_c_term		[mw_aa[i]] = AAInfo::getInfo ().getCTermPK ( aa );
		pk_n_term		[mw_aa[i]] = AAInfo::getInfo ().getNTermPK ( aa );
		pk_acidic_sc	[mw_aa[i]] = AAInfo::getInfo ().getAcidicSCPK ( aa );
		pk_basic_sc		[mw_aa[i]] = AAInfo::getInfo ().getBasicSCPK ( aa );
	}
	initialised = 1;
}
static char aa []			= "ACDEFGHIKLMNPQRSTUVWXYhmstuvwxy";// Also defined in lu_tag_frag.cpp /* NB. Rewrite the composition mask code after 32 aa's in this list */
static char aa_file []		= "ACDEFGHIKLMNPQRSTUVWLYhmstuvwxy";
//-----------------------------1234567890123456789012345678901

static void init_aa_arrays ()
{
	for ( unsigned int i = 0 ; aa [i] != '\0' ; i++ ) {
		double substituent;
		mi_amino_acid_wt	[aa[i]] = AAInfo::getInfo ().getMonoisotopicMass ( aa_file [i] );
		av_amino_acid_wt	[aa[i]] = AAInfo::getInfo ().getAverageMass ( aa_file [i] );
		nom_amino_acid_wt	[aa[i]] = (int)mi_amino_acid_wt [aa_file [i]];

		substituent = formula_to_monoisotopic_mass ( AAInfo::getInfo ().getSubstituentA ( aa_file[i] ).c_str () );
		mi_substituent_da [aa[i]] = substituent ? substituent + formula_to_monoisotopic_mass ( "N C2 H4" ) - mi_amino_acid_wt [aa[i]] : 0.0;
		mi_substituent_wa [aa[i]] = substituent ? substituent + formula_to_monoisotopic_mass ( "C3 H O" ) - mi_amino_acid_wt [aa[i]] : 0.0;
		substituent = formula_to_average_mass ( AAInfo::getInfo ().getSubstituentA ( aa_file[i] ).c_str () );
		av_substituent_da [aa[i]] = substituent ? substituent + formula_to_average_mass ( "N C2 H4" ) - av_amino_acid_wt [aa[i]] : 0.0;
		av_substituent_wa [aa[i]] = substituent ? substituent + formula_to_average_mass ( "C3 H O" ) - av_amino_acid_wt [aa[i]] : 0.0;

		substituent = formula_to_monoisotopic_mass ( AAInfo::getInfo ().getSubstituentB ( aa_file[i] ).c_str () );
		mi_substituent_db [aa[i]] = substituent ? substituent + formula_to_monoisotopic_mass ( "N C2 H4" ) - mi_amino_acid_wt [aa[i]] : 0.0;
		mi_substituent_wb [aa[i]] = substituent ? substituent + formula_to_monoisotopic_mass ( "C3 H O" ) - mi_amino_acid_wt [aa[i]] : 0.0;
		substituent = formula_to_average_mass ( AAInfo::getInfo ().getSubstituentB ( aa_file[i] ).c_str () );
		av_substituent_db [aa[i]] = substituent ? substituent + formula_to_average_mass ( "N C2 H4" ) - av_amino_acid_wt [aa[i]] : 0.0;
		av_substituent_wb [aa[i]] = substituent ? substituent + formula_to_average_mass ( "C3 H O" ) - av_amino_acid_wt [aa[i]] : 0.0;

		aa_composition_mask [aa[i]] = ( i == 0 ) ? 1 : aa_composition_mask [aa[i-1]] * 2;
	}
}
double getNTermMassShift ( const string& modAA )
{
	if ( localMonoisotopicFlag )
		return AAInfo::getInfo ().getModMonoisotopicMass ( "term_" + modAA ) - n_terminus_wt;
	else
		return AAInfo::getInfo ().getModAverageMass ( "term_" + modAA ) - n_terminus_wt;
}
double getCTermMassShift ( const string& modAA )
{
	if ( localMonoisotopicFlag )
		return AAInfo::getInfo ().getModMonoisotopicMass ( "term_" + modAA ) - c_terminus_wt;
	else
		return AAInfo::getInfo ().getModAverageMass ( "term_" + modAA ) - c_terminus_wt;
}
double getNeutralLossMassShift ( const string& modAA )
{
	if ( localMonoisotopicFlag ) {
		return AAInfo::getInfo ().getModMonoisotopicMass ( modAA );
	}
	else {
		return AAInfo::getInfo ().getModAverageMass ( modAA );
	}
}
void modifyAAInfo ( char aaCode, const string& modAA )
{
	if ( localMonoisotopicFlag ) {
		mi_amino_acid_wt [aaCode] = AAInfo::getInfo ().getModMonoisotopicMass ( modAA );
	}
	else {
		av_amino_acid_wt [aaCode] = AAInfo::getInfo ().getModAverageMass ( modAA );
	}
}
double getAminoAcidWt ( const string& modAA )
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
double getMaxDSubstituent ()
{
	double ret = 0.0;
	for ( int i = 0 ; i < AA_ARRAY_SIZE ; i++ ) {
		ret = genMax ( ret, -substituent_da [i] );
		ret = genMax ( ret, -substituent_db [i] );
	}
	return ( ret );
}
double getMaxWSubstituent ()
{
	double ret = 0.0;
	for ( int i = 0 ; i < AA_ARRAY_SIZE ; i++ ) {
		ret = genMax ( ret, -substituent_wa [i] );
		ret = genMax ( ret, -substituent_wb [i] );
	}
	return ( ret );
}
int* get_protein_int_mass_array ( const char* protein )
{
	static IntVector proteinMassArray;

	proteinMassArray.clear ();

	while ( *protein ) proteinMassArray.push_back ( nom_amino_acid_wt [*protein++] );
	proteinMassArray.push_back ( 0 );

	return ( &proteinMassArray [0] );
}
double* get_protein_double_mass_array ( const char* protein )
{
	static DoubleVector proteinMassArray;

	proteinMassArray.clear ();

	while ( *protein ) proteinMassArray.push_back ( amino_acid_wt [*protein++] );
	proteinMassArray.push_back ( 0 );

	return ( &proteinMassArray [0] );
}
