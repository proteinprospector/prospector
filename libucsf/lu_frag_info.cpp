/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_frag_info.cpp                                              *
*                                                                             *
*  Created    : July 21st 1996                                                *
*                                                                             *
*  Purpose    : Function to calculate information on enzyme fragments.        *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1996-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lu_aa_calc.h>
#include <lu_fas_enz.h>
#include <lu_frag_info.h>
#include <lu_prod_par.h>
#include <lu_param_list.h>
#include <lu_mass_seq.h>
#include <lu_mut_mtrx.h>
#include <lu_table.h>
#include <lu_delim.h>
using std::string;
using std::ostream;
using std::endl;

EnzymeFragmentContainer::EnzymeFragmentContainer ( const string& protein, const EnzymeParameters& enzymeParameters ) :
	protLen ( protein.length () )
{
	init_fasta_enzyme_function ( enzymeParameters.getEnzyme () );

	IntVector& cleavageIndicies = enzyme_fragmenter ( protein.c_str () );
	int numEnzymeFragments = cleavageIndicies.size ();
	unsigned int mask = enzymeParameters.getCompMask ();
	string compMaskType = enzymeParameters.getCompMaskType ();

	int missedCleavages = enzymeParameters.getMissedCleavages ();
	bool endTerminus = enzymeParameters.getEndTerminus ();
	int startStrip = endTerminus ? enzymeParameters.getStartStrip () : 0;
	int endStrip = endTerminus ? enzymeParameters.getEndStrip () : 0;
	char strippingTerminal = endTerminus ? enzymeParameters.getStrippingTerminal () : 'N';

	for ( int i = 0 ; i < numEnzymeFragments ; i++ ) {
		int startIndex = ( i == 0 ) ? 0 : cleavageIndicies [i-1] + 1;
		for ( int j = 0 ; i + j < numEnzymeFragments && j <= missedCleavages ; j++ ) {
			int endIndex = cleavageIndicies [i+j];
			for ( int k = startStrip ; k <= endStrip ; k++ ) {
				int start = startIndex + ( ( strippingTerminal == 'N' ) ? k : 0 );
				int end = endIndex - ( ( strippingTerminal == 'C' ) ? k : 0 );
				if ( end > start ) {
					EnzymeFragment ef ( protein, start, end, j );
					if ( mask ) {
						if ( compMaskType == "AND" ) {
							if ( ef.containsAAMaskAnd ( mask ) ) {
								enzymeFragmentList.push_back ( ef );
							}
						}
						else if ( compMaskType == "OR" ) {
							if ( ef.containsAAMaskOr ( mask ) ) {
								enzymeFragmentList.push_back ( ef );
							}
						}
					}
					else
						enzymeFragmentList.push_back ( ef );
				}
			}
		}
	}
}
EnzymeFragment::EnzymeFragment ( const string& protein, int start, int end, int missedCleavages ) :
	fragment ( protein, start, end - start + 1 ),
	startAA ( start + 1 ),
	endAA ( end + 1 ),
	missedCleavages ( missedCleavages )
{
	int len = end - start + 1;
	if ( cnbr_digest && fragment [len-1] == 'M' ) fragment [len-1] = 'h';

	previousAA = ( start == 0 ) ? '-' : ( cnbr_digest && protein [start-1] == 'M' ) ? 'h' : protein [start-1];
	nextAA = ( end == protein.length () - 1 ) ? '-' : protein [end+1];
	maxCharge = sumBasicResidues ( fragment ) + 1;
}
int EnzymeFragment::getNumPossMod ( char aa, bool nTermLinkPossible ) const
{
	int numPossMod = gen_strcharcount ( fragment, aa );
	if ( nextAA != '-' ) {									// If not the C-terminal fragment
		if ( fragment [fragment.length () - 1] == aa ) {	// This is a digest cleavage point
			numPossMod -= 1;
		}
	}
	if ( nTermLinkPossible && previousAA == '-' ) numPossMod += 1;	// N-terminal link possible
	return numPossMod;
}
double EnzymeFragment::getMass () const
{
	return peptide_formula_to_molecular_weight ( fragment.c_str () );
}
ElementalFormula EnzymeFragment::getElementalFormula ( const AACalculator& aaCalc ) const
{
	return aaCalc.calculateElementalComposition ( fragment );
}
bool EnzymeFragment::containsAAMaskAnd ( unsigned int mask ) const
{
	return ( mask & string_to_mask ( fragment.c_str () ) ) == mask;
}
bool EnzymeFragment::containsAAMaskOr ( unsigned int mask ) const
{
	return ( mask & string_to_mask ( fragment.c_str () ) ) != 0;
}
void EnzymeFragment::printHeaderHTML ( ostream& os )
{
	tableHeader ( os, "Start" );
	tableHeader ( os, "End" );
	tableHeader ( os, "Missed<br />Cleavages" );
	tableHeader ( os, "Sequence", "", "left" );
}
void EnzymeFragment::printHeaderDelimited ( ostream& os )
{
	delimitedHeader ( os, "Start AA" );
	delimitedHeader ( os, "End AA" );
	delimitedHeader ( os, "Missed Cleavages" );
	delimitedHeader ( os, "Previous AA" );
	delimitedHeader ( os, "Sequence" );
	delimitedHeader ( os, "Next AA" );
}
void EnzymeFragment::printHTML ( ostream& os, bool hideLinks, const MSProductLink& productLink, int maxCharge, char aa ) const
{
	PeptideSequence ps ( fragment );
	tableDataStart ( os, "", "right" );
		os << startAA << endl;
	tableDataEnd ( os );
	tableDataStart ( os, "", "right" );
		os << endAA << endl;
	tableDataEnd ( os );
	tableDataStart ( os, "", "center" );
		os << missedCleavages << endl;
	tableDataEnd ( os );
	tableDataStart ( os, "", "left", true );
		os << "(" << previousAA << ")";
		productLink.write1 ( os, ps, hideLinks, maxCharge, aa );
		os << "(" << nextAA << ")" << endl;
	tableDataEnd ( os );
}
void EnzymeFragment::printDelimited ( ostream& os ) const
{
	PeptideSequence ps ( fragment );
	delimitedCell ( os, startAA );
	delimitedCell ( os, endAA );
	delimitedCell ( os, missedCleavages );
	delimitedCell ( os, previousAA );
	delimitedCell ( os, ps.getPeptideString () );
	delimitedCell ( os, nextAA );
}
void EnzymeFragment::printXML ( ostream& os ) const
{
	PeptideSequence ps ( fragment );
	ParameterList::printXML ( os, "start_aa", startAA );
	ParameterList::printXML ( os, "end_aa", endAA );
	ParameterList::printXML ( os, "missed_cleavages", missedCleavages );
	ParameterList::printXML ( os, "previous_aa", previousAA );
	ParameterList::printXML ( os, "database_sequence", ps.getPeptideString () );
	ParameterList::printXML ( os, "next_aa", nextAA );
}
