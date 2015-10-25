/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_comp_srch.cpp                                              *
*                                                                             *
*  Created    : August 3rd 1996                                               *
*                                                                             *
*  Purpose    : A tool that can be used to complete the possible amino acid   *
*               composition of a peptide given a parent mass and partial      *
*               composition from the immonium ion data contained in a user's  *
*               tandem mass spectrum.                                         *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1996-2010) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifndef VIS_C
#include <stdexcept>
#endif
#include <lu_mass_conv.h>
#include <lu_comp_srch.h>
#include <lu_param_list.h>
#include <lu_charge.h>
#include <lu_comp_par.h>
#include <lu_aa_info.h>
#include <lu_aa_calc.h>
#include <lu_html.h>
#include <lu_elem_comp.h>
#include <lu_iso_par.h>
#include <lu_table.h>
using std::string;
using std::ostream;
using std::endl;
using std::sort;
using std::runtime_error;

MSCompSearch::MSCompSearch ( const MSCompParameters& params ) :
	MSProgram ( params ),
	compParams ( params ),
	cTermFormula ( "" )
{
}
void MSCompSearch::printBodyHTML ( ostream& os )
{
	bool monoisotopicFlag = compParams.getMonoisotopicFlag ();
	string combinationType = compParams.getCombinationType ();
	Combination_hit* combinations;
	char** elementalCombinations;
	int numCombinations;
	double effectiveParentMass;
	int maxHitsExceeded = 0;
	string fullAAList = compParams.getAAList ();

	AAInitInfo aaInitInfo = compParams.getAAInitInfo ();
	aaCalc = new AACalculator ( monoisotopicFlag, aaInitInfo.getConstMods () );
	cTermFormula = aaInitInfo.getCTermModFormula ();

	Peak parentPeak ( compParams.getParentMass (), compParams.getParentMassTolerance ()->getTolerance ( compParams.getParentMass (), compParams.getParentCharge () ), compParams.getParentCharge (), 100.0, getAdductMass ( monoisotopicFlag ) );
	effectiveParentMass = subtract_known_composition ( parentPeak.getMass (), compParams.getCompIons () );

	PeakContainer peakContainer;
	compositionSearch = new CompositionSearch ( compParams.getCompSearchParams (), peakContainer );

	startJavascript ( os );
	MSIsotopeLink isotopeLink;
	isotopeLink.printHTML ( os );
	endJavascript ( os );

	outputHTMLHeader ( os, &parentPeak, effectiveParentMass );

	StringVector rit = compParams.getRequestedIonTypes ();
	for ( StringVectorSizeType i = 0 ; i < rit.size () ; i++ ) {
		string ionType = rit [i];
		double effectiveIonMass = effectiveParentMass + ionTypeToOffset ( ionType );
		if ( effectiveIonMass > 1600 ) {
			ErrorHandler::genError ()->error ( "The effective parent mass (the ion mass minus the mass of the known amino acids) is too large.\n" );
		}
		if ( combinationType  == "Elemental" ) {
			maxHitsExceeded = calculate_elem_combinations ( AAInfo::getInfo ().aaListToElementalList ( fullAAList ), effectiveIonMass, parentPeak.getTolerance (), monoisotopicFlag, terminal_formula, compParams.getMaxReportedHits (), &elementalCombinations, &numCombinations );
		}
		else if ( combinationType == "Any Elemental" ) {
			StringVector s;
			maxHitsExceeded = calculate_elem_combinations ( s, effectiveIonMass, parentPeak.getTolerance (), monoisotopicFlag, "", compParams.getMaxReportedHits (), &elementalCombinations, &numCombinations );
		}
		else {
			try {
				maxHitsExceeded = calculate_combinations ( fullAAList, effectiveIonMass, parentPeak.getTolerance (), compParams.getMaxReportedHits (), combinationType != "Amino Acid", &combinations, &numCombinations );
			}
			catch ( runtime_error e ) {
				ErrorHandler::genError ()->error ( e );
			}
		}
		HitsContainer hits;
		if ( combinationType == "Amino Acid" )			hits = getHitsInformationAminoAcid			( combinations, numCombinations, ionType );
		if ( combinationType == "Peptide Elemental" )	hits = getHitsInformationPeptideElemental	( combinations, numCombinations, ionType );
		if ( combinationType == "Elemental" || combinationType == "Any Elemental" ) {
			hits = getHitsInformationElemental ( elementalCombinations, numCombinations, ionType );
		}
		outputHTMLTable ( os, compParams, hits, &parentPeak, maxHitsExceeded, ionType, isotopeLink );
	}
}
void MSCompSearch::printBodyXML ( ostream& os )
{
}
void MSCompSearch::outputHTMLHeader ( ostream& os, Peak* parentPeak, double effectiveParentMass )
{
	parentPeak->printHTML ( os, "Submitted Ion", instInf->getParentPeakPrecision () );

	compositionSearch->printHTML ( os );
	os << "MS-Comp's effective parent mass: <b>";
	genPrint ( os, effectiveParentMass, instInf->getParentPeakPrecision ().getMassDecimalPlaces () );
	os << "</b><br />";
//print_combination_search_aa_info ();

	os << "<hr />" << endl;
}
void MSCompSearch::outputHTMLTable ( ostream& os, const MSCompParameters& params, HitsContainer& hits, Peak* parentPeak, int maxHitsExceeded, const string& ionType, const MSIsotopeLink& isotopeLink )
{
	AAFormula knownAAComposition = AACalculator::calculateAAComposition ( compositionSearch->getKnownCompositionSequence () );
	string combinationType = params.getCombinationType ();
	double offset = ionTypeToOffset ( ionType );
	int i;
	int numHits = hits.size ();

	ParameterList::printHTML ( os, "Ion Type", ionType );
	if ( maxHitsExceeded )
		os << "<b>The maximum number of compositions has been exceeded.</b><br />" << endl;
	else {
		if ( numHits == 0 )
			os << "<b>No matching AA compositions.</b><br />" << endl;
		else {
			if ( combinationType == "Amino Acid" ) {
				ParameterList::printHTML ( os, "Unique AA Compositions (Combinations)", numHits );
				double totalPermutations = 0.0;
				for ( i = 0 ; i < numHits ; i++ ) {
					totalPermutations += hits [i].getNumPermutations ();
				}
				os << "Unique Peptide Sequences (Permutations): <b>";
				genPrint ( os, totalPermutations, 0 );
				os << "</b><br />";
			}
			ParameterList::printHTML ( os, "Unique Elemental Compositions", getNumUniqElemCompositions ( hits ) );

			os << "<table>" << endl;

			for ( i = 0 ; i < numHits ; i++ ) {
				if ( i == 0 ) CompHit::printHTMLHeader ( os, combinationType, params, ionType );
				hits [i].printHTML ( os, combinationType, isotopeLink, knownAAComposition, parentPeak, offset, params );
			}
			os << "</table>" << endl;
		}
	}
}
HitsContainer MSCompSearch::getHitsInformationAminoAcid ( Combination_hit* combinations, int numCombinations, const string& ionType )
{
	HitsContainer hits;
	for ( int i = 0 ; i < numCombinations ; i++ ) {
		string fullSequence = combinationPlusCombination ( combinations [i].array, compositionSearch->getKnownCompositionSequence () );
		if ( compositionSearch->doCompositionSearch ( fullSequence ) ) {
			hits.push_back ( CompHit ( aaCalc->calculatePeptideMW ( fullSequence ), fullSequence ) );
		}
	}
	sort ( hits.begin (), hits.end (), sort_comp_hits () );

	for ( HitsContainerSizeType j = 0 ; j < hits.size () ; j++ ) {
		if ( j == 0 || ( genAbsDiff ( hits [j].pepMW, hits [j-1].pepMW ) > FORMULA_DIFF_TOL ) ) 
			hits [j].elemComp = adjustedElementalFormula ( aaCalc->calculateElementalComposition ( hits [j].getSequence () ), ionType );
		else
			hits [j].elemComp = hits [j-1].elemComp;
	}
	return ( hits );
}
HitsContainer MSCompSearch::getHitsInformationPeptideElemental ( Combination_hit* combinations, int numCombinations, const string& ionType )
{
	HitsContainer hits;

	for ( int i = 0 ; i < numCombinations ; i++ ) {
		string fullSequence = combinationPlusCombination ( combinations [i].array, compositionSearch->getKnownCompositionSequence () );
		ElementalFormula ef = adjustedElementalFormula ( aaCalc->calculateElementalComposition ( fullSequence ), ionType );
		hits.push_back ( CompHit ( ef, aaCalc->calculatePeptideMW ( fullSequence ) ) );
	}
	sort ( hits.begin (), hits.end (), sort_comp_hits () );

	return ( hits );
}
HitsContainer MSCompSearch::getHitsInformationElemental ( char** elementalCombinations, int numCombinations, const string& ionType )
{
	double ( *mass_convert ) (const char*);
	bool monoisotopicFlag = compParams.getMonoisotopicFlag ();
	if ( monoisotopicFlag ) mass_convert = formula_to_monoisotopic_mass;
	else mass_convert = formula_to_average_mass;
	ElementalFormula knownElementalComposition = aaCalc->calculateFragElementalComposition ( compositionSearch->getKnownCompositionSequence () );
	HitsContainer hits;

	for ( int i = 0 ; i < numCombinations ; i++ ) {
		ElementalFormula elemComp = knownElementalComposition;
		elemComp += elementalCombinations [i];
		double mw = mass_convert ( elemComp.getFormula().c_str () );
		ElementalFormula ef = adjustedElementalFormula ( elemComp, ionType );
		hits.push_back ( CompHit ( mw, ef ) );
	}
	sort ( hits.begin (), hits.end (), sort_comp_hits () );

	return ( hits );
}
int MSCompSearch::getNumUniqElemCompositions ( const HitsContainer& hits )
{
	int numHits = hits.size ();
	if ( numHits == 0 ) return 0;
	int numUniqueElemComp = 1;
	for ( int i = 1; i < numHits ; i++ ) {
		if ( genAbsDiff ( hits [i].pepMW, hits [i-1].pepMW ) > FORMULA_DIFF_TOL ) numUniqueElemComp++;
	}
	return ( numUniqueElemComp );
}
double MSCompSearch::ionTypeToOffset ( const string& ionType )
{
	double offset;
	double corr = terminal_wt - n_terminus_wt;
	if ( ionType == "MH+" )		offset = 0.0;
	if ( ionType == "a" )		offset = a_tag_offset			+ corr + ELECTRON_REST_MASS;
	if ( ionType == "a-H2O" )	offset = a_h2o_tag_offset		+ corr + ELECTRON_REST_MASS;
	if ( ionType == "a-NH3" )	offset = a_nh3_tag_offset		+ corr + ELECTRON_REST_MASS;
	if ( ionType == "a-H3PO4" )	offset = a_h3po4_tag_offset		+ corr + ELECTRON_REST_MASS;
	if ( ionType == "b" )		offset = b_tag_offset			+ corr + ELECTRON_REST_MASS;
	if ( ionType == "b-H2O" )	offset = b_h2o_tag_offset		+ corr + ELECTRON_REST_MASS;
	if ( ionType == "b+H2O" )	offset = b_plus_h2o_tag_offset	+ corr + ELECTRON_REST_MASS;
	if ( ionType == "b-NH3" )	offset = b_nh3_tag_offset		+ corr + ELECTRON_REST_MASS;
	if ( ionType == "b-SOCH4" )	offset = b_soch4_tag_offset		+ corr + ELECTRON_REST_MASS;
	if ( ionType == "b-H3PO4" )	offset = b_h3po4_tag_offset		+ corr + ELECTRON_REST_MASS;
	if ( ionType == "c" )		offset = c_tag_offset			+ corr + ELECTRON_REST_MASS;

	if ( ionType == "y" )		offset = y_tag_offset + h2;
	if ( ionType == "y-H2O" )	offset = y_h2o_tag_offset + h2;
	if ( ionType == "y-NH3" )	offset = y_nh3_tag_offset + h2;
	if ( ionType == "y-SOCH4" )	offset = y_soch4_tag_offset + h2;
	if ( ionType == "y-H3PO4" )	offset = y_h3po4_tag_offset + h2;

	return offset;
}

CompHit::CompHit ( double mw, const string& sequence ) :
	sequence ( sequence ),
	pepMW ( mw ),
	aaComp ( AACalculator::calculateAAComposition ( sequence ) ),
	numPermutations ( calc_num_permutations_from_combination ( sequence ) )
{
}
CompHit::CompHit ( ElementalFormula& elemComp, double mw ) :
	elemComp ( elemComp ),
	pepMW ( mw )
{
}
CompHit::CompHit ( double mw, ElementalFormula& elemComp ) :
	pepMW ( mw ),
	elemComp ( elemComp ),
	doubleBondEquivalent ( formula_to_double_bond_equivalent ( elemComp ) )
{
}
void CompHit::printHTMLHeader ( ostream& os, const string& combinationType, const MSCompParameters& params, const string& ionType )
{
	tableRowStart ( os );
		if ( combinationType == "Amino Acid" ) {
			tableHeader ( os, "AA<br /> Composition" );
			tableHeader ( os, "Missing <br />AA Composition<br />(Effective Parent Mass)" );
		}
		tableHeader ( os, "Elemental<br /> Composition" );
		tableHeaderStart ( os );
			os << " Peptide Mass<br /> (" << ionType << ")" << endl;
		tableHeaderEnd ( os );
		tableHeaderStart ( os );
			os << "Mass<br />Error<br />(" << params.getParentMassTolerance ()->getUnitsString () << ")" << endl;
		tableHeaderEnd ( os );
		if ( combinationType == "Elemental" || combinationType == "Any Elemental" ) {
			tableHeader ( os, "Double<br />Bond<br />Equivalent<br />" );
		}
		if ( combinationType == "Amino Acid" ) {
			tableHeader ( os, "Number of<br />Sequence<br />Permutations" );
		}
	tableRowEnd ( os );
}
void CompHit::printHTML ( ostream& os, const string& combinationType, const MSIsotopeLink& isotopeLink, const AAFormula& knownAAComposition, const Peak* parentPeak, double offset, const MSCompParameters& params )
{
	tableRowStart ( os );
		if ( combinationType == "Amino Acid" ) {
			tableDataStart ( os );
				os << aaComp << endl;
			tableDataEnd ( os );
			AAFormula unknownAAComp = aaComp;
			unknownAAComp -= knownAAComposition;
			tableDataStart ( os );
				os << unknownAAComp << endl;
			tableDataEnd ( os );
		}
		tableDataStart ( os );
			isotopeLink.write ( os, elemComp, params.getParentCharge () );
		tableDataEnd ( os );
		tableDataStart ( os, "", "right" );
			genPrint ( os, pepMW - offset, instInf->getParentPeakPrecision ().getMassDecimalPlaces () );
			os << endl;
		tableDataEnd ( os );
		tableDataStart ( os, "", "right" );
			genPrint ( os, params.getParentMassTolerance ()->getError ( parentPeak->getMass (), pepMW - offset, parentPeak->getCharge () ), 2 );
			os << endl;
		tableDataEnd ( os );
		if ( combinationType == "Elemental" || combinationType == "Any Elemental" ) {
			tableDataStart ( os, "", "right" );
				genPrint ( os, doubleBondEquivalent, 1 );
				os << endl;
			tableDataEnd ( os );
		}
		if ( combinationType == "Amino Acid" ) {
			tableDataStart ( os, "", "right" );
				genPrint ( os, numPermutations, 0 );
				os << endl;
			tableDataEnd ( os );
		}
	tableRowEnd ( os );
}
ElementalFormula MSCompSearch::adjustedElementalFormula ( const ElementalFormula& elementalFormula, const string& ionType )
{
	ElementalFormula adjustedFormula ( elementalFormula );

	if ( ionType != "MH+" && ionType != "y" ) {

		if ( ionType == "y-H2O" )	adjustedFormula -= "H2 O";
		else if ( ionType == "y-NH3" )	adjustedFormula -= "N H3";
		else if ( ionType == "y-SOCH4" )adjustedFormula -= "S O C H4";
		else if ( ionType == "y-H3PO4" )adjustedFormula -= "H3 P O4";
		else {
			ElementalFormula terminusAdjustFormula = "H";
			terminusAdjustFormula += cTermFormula;
			terminusAdjustFormula += "O H";
			adjustedFormula -= terminusAdjustFormula;

			if ( ionType == "b+H2O" )	adjustedFormula += "H2 O";
			if ( ionType == "c" )		adjustedFormula += "N H3";
			if ( ionType == "a" )		adjustedFormula -= "C O";
			if ( ionType == "a-H2O" )	adjustedFormula -= "C O2 H2";
			if ( ionType == "a-NH3" )	adjustedFormula -= "C O N H3";
			if ( ionType == "a-H3PO4" )	adjustedFormula -= "H3 C P O5";
			//if ( ionType == "b" )	do nothing
			if ( ionType == "b-H2O" )	adjustedFormula -= "H2 O";
			if ( ionType == "b-NH3" )	adjustedFormula -= "N H3";
			if ( ionType == "b-SOCH4" )	adjustedFormula -= "S O C H4";
			if ( ionType == "b-H3PO4" )	adjustedFormula -= "H3 P O4";
		}
	}
	return ( adjustedFormula );
}
void MSCompSearch::printParamsBodyHTML ( ostream& os ) const
{
	compParams.printHTML ( os );
}
