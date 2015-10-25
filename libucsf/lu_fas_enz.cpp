/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_fas_enz.cpp                                                *
*                                                                             *
*  Created    : July 20th 1996                                                *
*                                                                             *
*  Purpose    : Enzymatic and chemical digest cleavage functions.             *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1996-2009) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lg_new.h>
#include <lgen_error.h>
#define LUCSF_FAS_ENZ_MAIN
#include <lu_param_list.h>
#include <lu_fas_enz.h>
#include <lu_mass.h>
#include <lu_composit.h>
#include <lu_getfil.h>
#include <lu_html_form.h>
#include <lu_cgi_val.h>
using std::string;
using std::ostream;
using std::endl;

string EnzymeParameters::ENZYME = "enzyme";
string EnzymeParameters::ALLOW_NON_SPECIFIC = "allow_non_specific";
string EnzymeParameters::MISSED_CLEAVAGES = "missed_cleavages";
string EnzymeParameters::END_TERMINUS = "end_terminus";
string EnzymeParameters::STRIPPING_TERMINAL = "stripping_terminal";
string EnzymeParameters::START_STRIP = "start_strip";
string EnzymeParameters::END_STRIP = "end_strip";
string EnzymeParameters::COMP_ION = "comp_ion";
string EnzymeParameters::COMP_MASK_TYPE = "comp_mask_type";

string const EnzymeParameters::NO_ENZYME = "No enzyme";
string EnzymeParameters::ENZYME_DEFAULT = "Trypsin";
string EnzymeParameters::ALLOW_NON_SPECIFIC_DEFAULT = "at 0 termini";
int EnzymeParameters::MISSED_CLEAVAGES_DEFAULT = 1;
bool EnzymeParameters::END_TERMINUS_DEFAULT = false;
char EnzymeParameters::STRIPPING_TERMINAL_DEFAULT = 'N';
int EnzymeParameters::START_STRIP_DEFAULT = 2;
int EnzymeParameters::END_STRIP_DEFAULT = 4;
string EnzymeParameters::COMP_MASK_TYPE_DEFAULT = "AND";

EnzymeParameters::EnzymeParameters ( const ParameterList* params, bool noNoEnzyme ) :
	enzyme				( initEnzyme ( params, noNoEnzyme ) ),
	allowNonSpecific	( params->getStringValue	( ALLOW_NON_SPECIFIC, ALLOW_NON_SPECIFIC_DEFAULT ) ),
	missedCleavages		( params->getIntValue		( MISSED_CLEAVAGES, MISSED_CLEAVAGES_DEFAULT ) ),
	endTerminus			( params->getBoolValue		( END_TERMINUS, END_TERMINUS_DEFAULT ) ),
	strippingTerminal	( params->getCharValue		( STRIPPING_TERMINAL, STRIPPING_TERMINAL_DEFAULT ) ),
	startStrip			( params->getIntValue		( START_STRIP, START_STRIP_DEFAULT ) ),
	endStrip			( params->getIntValue		( END_STRIP, END_STRIP_DEFAULT ) ),
	compIons			( params->getStringVectorValue ( COMP_ION ) ),
	compMaskType		( params->getStringValue		( COMP_MASK_TYPE, COMP_MASK_TYPE_DEFAULT ) ),
	compMask			( CompositionSearch::calculateIncludeMask ( compIons ) )
{
}
string EnzymeParameters::initEnzyme ( const ParameterList* params, bool noNoEnzyme )
{
	string value = params->getStringValue ( ENZYME, ENZYME_DEFAULT );
	if ( noNoEnzyme && value == NO_ENZYME ) return ENZYME_DEFAULT;
	else return value;
}
void EnzymeParameters::printHTML ( ostream& os ) const
{
	ParameterList::printHTML ( os, "Digest Used", enzyme );
	if ( allowNonSpecific != ALLOW_NON_SPECIFIC_DEFAULT ) {
		ParameterList::printHTML ( os, "Allow Non-Specific", allowNonSpecific );
	}
	ParameterList::printHTML ( os, "Max. # Missed Cleavages", missedCleavages );
	if ( endTerminus ) {
		os << "End Terminus Parameters<br />" << endl;
		ParameterList::printHTML ( os, "Stripping Terminal", strippingTerminal );
		ParameterList::printHTML ( os, "Start Strip", startStrip );
		ParameterList::printHTML ( os, "End Strip", endStrip );
	}
	if ( !compIons.empty () ) {
		ParameterList::printHTMLContainer ( os, "Composition Include", compIons );
		ParameterList::printHTML ( os, "Composition Mask Type", compMaskType );
	}
}
void EnzymeParameters::putCGI ( ostream& os ) const
{
	printCGIString ( os, ENZYME, enzyme );
	printCGI ( os, MISSED_CLEAVAGES, missedCleavages );

	if ( endTerminus ) {
		printCGI ( os, END_TERMINUS, endTerminus );
		printCGI ( os, STRIPPING_TERMINAL, strippingTerminal );
		printCGI ( os, START_STRIP, startStrip );
		printCGI ( os, END_STRIP, endStrip );
	}
	if ( !compIons.empty () ) {
		for ( StringVectorSizeType i = 0 ; i < compIons.size () ; i++ ) {
			printCGI ( os, COMP_ION, compIons [i] );
		}
		printCGI ( os, COMP_MASK_TYPE, compMaskType );
	}
}
void EnzymeParameters::putNoEnzymeCGI ( ostream& os ) const
{
	printCGIString ( os, "enzyme", "No enzyme" );
	if ( !compIons.empty () ) {
		for ( StringVectorSizeType i = 0 ; i < compIons.size () ; i++ ) {
			printCGI ( os, COMP_ION, compIons [i] );
		}
		printCGI ( os, COMP_MASK_TYPE, compMaskType );
	}
}
void EnzymeParameters::putHiddenFormEntryJavascript ( ostream& os ) const
{
	printHTMLFORMJavascriptHidden ( os, ENZYME, enzyme );
	printHTMLFORMJavascriptHidden ( os, MISSED_CLEAVAGES, missedCleavages );

	if ( endTerminus ) {
		printHTMLFORMJavascriptHidden ( os, END_TERMINUS, endTerminus );
		printHTMLFORMJavascriptHidden ( os, STRIPPING_TERMINAL, strippingTerminal );
		printHTMLFORMJavascriptHidden ( os, START_STRIP, startStrip );
		printHTMLFORMJavascriptHidden ( os, END_STRIP, endStrip );
	}
	if ( !compIons.empty () ) {
		for ( StringVectorSizeType i = 0 ; i < compIons.size () ; i++ ) {
			printHTMLFORMJavascriptHidden ( os, COMP_ION, compIons [i] );
		}
		printHTMLFORMJavascriptHidden ( os, COMP_MASK_TYPE, compMaskType );
	}
}
void EnzymeParameters::putNoEnzymeHiddenFormEntryJavascript ( ostream& os ) const
{
	printHTMLFORMJavascriptHidden ( os, "enzyme", "No enzyme" );
	if ( !compIons.empty () ) {
		for ( StringVectorSizeType i = 0 ; i < compIons.size () ; i++ ) {
			printHTMLFORMJavascriptHidden ( os, COMP_ION, compIons [i] );
		}
		printHTMLFORMJavascriptHidden ( os, COMP_MASK_TYPE, compMaskType );
	}
}
void EnzymeParameters::copyToCGI ( ostream& os, const ParameterList* params, bool noNoEnzyme )
{
	printCGIString ( os, ENZYME, initEnzyme ( params, noNoEnzyme ) );
	params->copyToCGI ( os, MISSED_CLEAVAGES );
	if ( params->copyToCGI ( os, END_TERMINUS ) ) {
		params->copyToCGI ( os, STRIPPING_TERMINAL );
		params->copyToCGI ( os, START_STRIP );
		params->copyToCGI ( os, END_STRIP );
	}
	params->copyToCGI ( os, COMP_ION );
	params->copyToCGI ( os, COMP_MASK_TYPE );
}
void EnzymeParameters::copyToHiddenFormEntry ( ostream& os, const ParameterList* params, bool noNoEnzyme )
{
	printHTMLFORMHidden ( os, ENZYME, initEnzyme ( params, noNoEnzyme ) );
	params->copyToHiddenFormEntry ( os, MISSED_CLEAVAGES );
	if ( params->copyToHiddenFormEntry ( os, END_TERMINUS ) ) {
		params->copyToHiddenFormEntry ( os, STRIPPING_TERMINAL );
		params->copyToHiddenFormEntry ( os, START_STRIP );
		params->copyToHiddenFormEntry ( os, END_STRIP );
	}
	params->copyToHiddenFormEntry ( os, COMP_ION );
	params->copyToHiddenFormEntry ( os, COMP_MASK_TYPE );
}
void EnzymeParameters::copyToHiddenFormJavascriptEntry ( ostream& os, const ParameterList* params, bool noNoEnzyme )
{
	printHTMLFORMJavascriptHidden ( os, ENZYME, initEnzyme ( params, noNoEnzyme ) );
	params->copyToHiddenFormJavascriptEntry ( os, MISSED_CLEAVAGES );
	if ( params->copyToHiddenFormJavascriptEntry ( os, END_TERMINUS ) ) {
		params->copyToHiddenFormJavascriptEntry ( os, STRIPPING_TERMINAL );
		params->copyToHiddenFormJavascriptEntry ( os, START_STRIP );
		params->copyToHiddenFormJavascriptEntry ( os, END_STRIP );
	}
	params->copyToHiddenFormJavascriptEntry ( os, COMP_ION );
	params->copyToHiddenFormJavascriptEntry ( os, COMP_MASK_TYPE );
}
string DigestTable::getBreakMask ( const string& enzymeName ) const
{
	MapStringToStringConstIterator cur = breakMaskTable.find ( enzymeName );
	if ( cur == breakMaskTable.end () )
		DigestTable::errorHandler ( enzymeName );
	return (*cur).second;
}
string DigestTable::getExcludeMask ( const string& enzymeName ) const
{
	MapStringToStringConstIterator cur = excludeMaskTable.find ( enzymeName );
	if ( cur == excludeMaskTable.end () )
		DigestTable::errorHandler ( enzymeName );
	return (*cur).second;
}
char DigestTable::getSpecificity ( const string& enzymeName ) const
{
	MapStringToCharConstIterator cur = digestSpecificityTable.find ( enzymeName );
	if ( cur == digestSpecificityTable.end () )
		DigestTable::errorHandler ( enzymeName );
	return (*cur).second;
}
void DigestTable::errorHandler ( const string& enzymeName )
{
	string errMessage = "Invalid enzyme name: " + enzymeName;
	ErrorHandler::genError ()->error ( errMessage );
}
DigestTable::DigestTable ()
{
	GenCommentedIFStream ifs ( MsparamsDir::instance ().getParamPath ( "enzyme.txt" ) );
	string line;
	int phase = 1;
	while ( ifs.getUncommentedLine ( line ) ) {
		if ( phase == 1 )		names.push_back ( line );
		else if ( phase == 2 )	breakMaskTable [names.back ()] = line;
		else if ( phase == 3 )	excludeMaskTable [names.back ()] = line;
		else if ( phase == 4 )	digestSpecificityTable [names.back ()] = line [0];
		if ( phase == 4 )	phase = 1;
		else				phase++;
	}
	GenCommentedIFStream ifs2 ( MsparamsDir::instance ().getParamPath ( "enzyme_comb.txt" ) );
	while ( ifs2.getUncommentedLine ( line ) ) {
		names.push_back ( line );
	}
}
DigestTable& DigestTable::instance ()
{
	static DigestTable digestTable;
	return digestTable;
}

static UIntVector break_mask_array;
static StringVector break_aas_array;
static UIntVector exclude_mask_array;
static StringVector exclude_aas_array;
static CharVector digestSpecificityArray;
static int numDigests = 0;

static char digestSpecificity = '0';
static string break_aas;
static unsigned int break_mask;
static unsigned int exclude_mask;

static IntVector& calc_fasta_c_term_fragments ( const string& peptideFormula );
static IntVector& calc_fasta_n_term_fragments ( const string& peptideFormula );
static IntVector& calc_fasta_multi_digest_fragments ( const string& peptideFormula );

void init_fasta_enzyme_function ( const string& enzyme )
{
	static bool called = false;
	if ( called ) return;
	called = true;
	static DigestTable& digestTable = DigestTable::instance ();

	if ( enzyme == "No enzyme" ) {
		ErrorHandler::genError ()->error ( "No enzyme not valid option for this program." );
		return;
	}
	StringVector enzymeNameTable;
	cnbr_digest = false;
	StringSizeType start = 0;
	StringSizeType end = 0;
	for ( ; ; ) {
		end = enzyme.find_first_of ( "/", start );
		string enzymeName = enzyme.substr ( start, end-start );
		if ( enzymeName == "CNBr" ) cnbr_digest = true;
		enzymeNameTable.push_back ( enzymeName );
		if ( end == string::npos ) break;
		start = end + 1;
	}
	numDigests = enzymeNameTable.size ();

	for ( int i = 0 ; i < numDigests ; i++ ) {
		string br_aas = digestTable.getBreakMask (enzymeNameTable[i]);
		break_aas_array.push_back ( br_aas );
		break_mask_array.push_back ( string_to_mask ( br_aas ) );

		string exclude_aas = digestTable.getExcludeMask (enzymeNameTable[i]);
		exclude_aas_array.push_back ( exclude_aas );
		if ( exclude_aas != "-" ) exclude_mask_array.push_back ( string_to_mask ( exclude_aas ) );
		else exclude_mask_array.push_back ( 0 );

		digestSpecificityArray.push_back ( digestTable.getSpecificity (enzymeNameTable[i]) );
	}	
	if ( numDigests == 1 ) {
		break_mask = break_mask_array [0];
		break_aas = break_aas_array [0];
		exclude_mask = exclude_mask_array [0];
		digestSpecificity = digestSpecificityArray [0];
		if ( digestSpecificity == 'C' ) enzyme_fragmenter = calc_fasta_c_term_fragments;
		else enzyme_fragmenter = calc_fasta_n_term_fragments;
	}
	else {
		int i;
		for ( i = 1 ; i < numDigests ; i++ ) {
			if ( digestSpecificityArray [i] != digestSpecificityArray [0] ) break;
		}
		if ( i == numDigests ) digestSpecificity = digestSpecificityArray [0];	/* All digest specificities must be the same for digestSpecificity to be set */

		enzyme_fragmenter = calc_fasta_multi_digest_fragments;
	}
}
char get_enzyme_terminal_specificity ()
{
	return digestSpecificity;
}
static IntVector& calc_fasta_c_term_fragments ( const string& peptideFormula )
{
	static IntVector cleavageIndex;
	int numAA = peptideFormula.length ();
	cleavageIndex.reserve ( numAA );
	cleavageIndex.clear ();

	int numFragments = 0;
	if ( numAA ) {
		int penultimateAA = numAA - 1;

		for ( int i = 0 ; i < penultimateAA ; i++ ) {
			if ( aa_composition_mask [peptideFormula [i]] & break_mask ) {
				if ( (aa_composition_mask [peptideFormula [i+1]] & exclude_mask) == 0 ) {
					cleavageIndex.push_back ( i );
				}
			}
		}
		cleavageIndex.push_back ( penultimateAA );
	}
	return ( cleavageIndex );
}
static IntVector& calc_fasta_n_term_fragments ( const string& peptideFormula )
{
	static IntVector cleavageIndex;
	int numAA = peptideFormula.length ();
	cleavageIndex.reserve ( numAA );
	cleavageIndex.clear ();

	if ( numAA ) {
		int penultimateAA = numAA - 1;

		for ( int i = 0 ; i < penultimateAA ; i++ ) {
			if ( aa_composition_mask [peptideFormula [i+1]] & break_mask ) {
				if ( (aa_composition_mask [peptideFormula [i]] & exclude_mask ) == 0 ) {
					cleavageIndex.push_back ( i );
				}
			}
		}
		cleavageIndex.push_back ( penultimateAA );
	}
	return ( cleavageIndex );
}
static IntVector& calc_fasta_multi_digest_fragments ( const string& peptideFormula )
{
	static IntVector cleavageIndex;
	int numAA = peptideFormula.length ();
	cleavageIndex.reserve ( numAA );
	cleavageIndex.clear ();

	if ( numAA ) {
		int penultimateAA = numAA - 1;
		for ( int i = 0 ; i < penultimateAA ; i++ ) {
			for ( int j = 0 ; j < numDigests ; j++ ) {
				char break_test_aa = ( digestSpecificityArray [j] == 'C' ) ? peptideFormula [i] : peptideFormula [i+1];
				if ( aa_composition_mask [break_test_aa] & break_mask_array [j] ) {
					char exclude_test_aa = ( digestSpecificityArray [j] == 'C' ) ? peptideFormula [i+1] : peptideFormula [i];
					if ( (aa_composition_mask [exclude_test_aa] & exclude_mask_array [j]) == 0 ) {
						cleavageIndex.push_back ( i );
						break;
					}
				}
			}
		}
		cleavageIndex.push_back ( penultimateAA );
	}
	return ( cleavageIndex );
}
DoubleVector& get_cleaved_masses ( const string& protein, const IntVector& cleavageIndex )
{
	static DoubleVector cleavedMassArray ( 36000 );
	StringSizeType numAA = protein.length ();
	if ( numAA > cleavedMassArray.size () ) cleavedMassArray.resize ( numAA );

	for ( StringSizeType i = 0, j = 0 ; i < numAA ; ) {
		double mass = 0.0;
		StringSizeType index = cleavageIndex [j];
		while ( i <= index ) {
			mass += amino_acid_wt [protein[i++]];
		}
		cleavedMassArray [j++] = mass;
	}
	return cleavedMassArray;
}
/*
This function is similar to the one above except that it stops adding up the
mass of a fragment if it is beyond limit for efficiency reasons. It is useful
where you are looking for a fragment whose mass you know to be less than limit.
*/
DoubleVector& get_cleaved_masses_to_limit ( const string& protein, const IntVector& cleavageIndex, double limit )
{
	static DoubleVector cleavedMassArray ( 36000 );
	StringSizeType numAA = protein.length ();
	if ( numAA > cleavedMassArray.size () ) cleavedMassArray.resize ( numAA );

	for ( StringSizeType i = 0, j = 0 ; i < numAA ; ) {
		double mass = 0.0;
		StringSizeType index = cleavageIndex [j];
		while ( i <= index ) {
			if ( mass <= limit ) mass += amino_acid_wt [protein[i]];
			i++;
		}
		cleavedMassArray [j++] = mass;
	}
	return cleavedMassArray;
}
