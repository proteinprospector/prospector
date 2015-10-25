/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_mass_elem.cpp                                              *
*                                                                             *
*  Created    : July 20th 1996                                                *
*                                                                             *
*  Purpose    : Gets the elemental information from the params/elements.txt   *
*               file. Functions for manipulating elemental formulae.          *
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
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <lu_mass_elem.h>
#include <lu_getfil.h>
#include <lu_param_list.h>
using std::string;
using std::fill;
using std::find;
using std::find_if;
using std::ostream;
using std::runtime_error;

static void replaceElement ( const string& hElement, const string& element, int nominalLight, int nominalHeavy, double percent );

IsotopePurity::IsotopePurity ( const ParameterList* params ) :
	percentC13 ( params->getDoubleValue ( "percent_C13", 100.0 ) ),
	percentN15 ( params->getDoubleValue ( "percent_N15", 100.0 ) ),
	percentO18 ( params->getDoubleValue ( "percent_O18", 100.0 ) )
{
	if ( percentC13 != 100.0 ) replaceElement ( "13C", "C", 12, 13, percentC13 );
	if ( percentN15 != 100.0 ) replaceElement ( "15N", "N", 14, 15, percentN15 );
	if ( percentO18 != 100.0 ) replaceElement ( "18O", "O", 16, 18, percentO18 );
}
void IsotopePurity::copyToCGI ( ostream& os, const ParameterList* params )
{
	if ( params->getDoubleValue ( "percent_C13", 100.0 ) != 100.0 ) params->copyToCGI ( os, "percent_C13" );
	if ( params->getDoubleValue ( "percent_N15", 100.0 ) != 100.0 ) params->copyToCGI ( os, "percent_N15" );
	if ( params->getDoubleValue ( "percent_O18", 100.0 ) != 100.0 ) params->copyToCGI ( os, "percent_O18" );
}
static bool initialised = false;
static StringVector	elementTable;
static IntVector	valencyTable;
static DoubleVectorVector isotopeMassArray;
static DoubleVectorVector isotopeAbundanceArray;
static DoubleVector	monoisotopicMassArray;
static DoubleVector	idealMonoisotopicMassArray;
static DoubleVector	nominalMassArray;
static DoubleVector	averageMassArray;

static double formula_to_numeric ( const string& formula, const DoubleVector& numericalArray );
static double formula_to_numeric ( ElementalFormula& elemForm, const DoubleVector& numericalArray );
static double formula_to_percent_element ( ElementalFormula& elemForm, const string& element, const DoubleVector& numericalArray );
static void initialise_mass_calculations ();
static inline void initialise ();

static inline void initialise ()
{
	if ( initialised == false ) initialise_mass_calculations ();
}
static double formula_to_numeric ( const string& formula, const DoubleVector& numericalArray )
{
	ElementalFormula elemForm ( formula );
	return formula_to_numeric ( elemForm, numericalArray );
}
static double formula_to_numeric ( ElementalFormula& elemForm, const DoubleVector& numericalArray )
{
	int* multiplier_list = formula_to_multiplier_list ( elemForm );
	double numeric_value = 0.0;

	for ( StringVectorSizeType i = 0 ; i < elementTable.size () ; i++ ) {
		numeric_value += numericalArray [i] * multiplier_list [i];
	}
	delete [] multiplier_list;

	return numeric_value;
}
double formula_to_monoisotopic_mass ( const char* formula )
{
	initialise ();
	return formula_to_numeric ( formula, monoisotopicMassArray );
}
double formula_to_monoisotopic_mass ( ElementalFormula& elemForm )
{
	initialise ();
	return formula_to_numeric ( elemForm, monoisotopicMassArray );
}
double formula_to_ideal_monoisotopic_mass ( const char* formula )
{
	initialise ();
	return formula_to_numeric ( formula, idealMonoisotopicMassArray );
}
double formula_to_nominal_mass ( const char* formula )
{
	initialise ();
	return formula_to_numeric ( formula, nominalMassArray );
}
double formula_to_nominal_mass ( ElementalFormula& elemForm )
{
	initialise ();
	return formula_to_numeric ( elemForm, nominalMassArray );
}
double formula_to_average_mass ( const char* formula )
{
	initialise ();
	return formula_to_numeric ( formula, averageMassArray );
}
double formula_to_average_mass ( ElementalFormula& elemForm )
{
	initialise ();
	return formula_to_numeric ( elemForm, averageMassArray );
}
static double formula_to_percent_element ( ElementalFormula& elemForm, const string& element, const DoubleVector& numericalArray )
{
	int* multiplier_list = formula_to_multiplier_list ( elemForm );
	double single_value = 0.0;
	double total_value = 0.0;

	for ( StringVectorSizeType i = 0 ; i < elementTable.size () ; i++ ) {
		if ( elementTable [i] == element ) single_value = numericalArray [i] * multiplier_list [i];
		total_value += numericalArray [i] * multiplier_list [i];
	}
	delete [] multiplier_list;

	return single_value / total_value;
}
double formula_to_percent_element_monoisotopic ( ElementalFormula& formula, const char* element )
{
	initialise ();

	return formula_to_percent_element ( formula, element, monoisotopicMassArray );
}
double formula_to_percent_element_average ( ElementalFormula& formula, const char* element )
{
	initialise ();

	return formula_to_percent_element ( formula, element, averageMassArray );
}
double formula_to_double_bond_equivalent ( ElementalFormula& f )
{
	double doubleBondEquivalent = 1.0;

	for ( f.first () ; f.isDone () ; f.next () ) {
		for ( StringVectorSizeType i = 0 ; i < elementTable.size () ; i++ ) {
			if ( elementTable [i] == f.element () ) {
				int numAtoms = f.multiplier ();
				switch ( valencyTable [i] ) {
					case 1:
						doubleBondEquivalent -= numAtoms / 2.0;
						break;
					case 3:
						doubleBondEquivalent += numAtoms / 2.0;
						break;
					case 4:
						doubleBondEquivalent += numAtoms;
						break;
				}
			}
		}
	}
	return doubleBondEquivalent;
}
double num_atoms_to_double_bond_equivalent ( int* num_atoms )
{
	double double_bond_equivalent = 1.0;

	for ( StringVectorSizeType i = 0 ; i < elementTable.size () ; i++ ) {
		if ( num_atoms [i] ) {
			switch ( valencyTable [i] ) {
				case 1:
					double_bond_equivalent -= num_atoms [i] / 2.0;
					break;
				case 3:
					double_bond_equivalent += num_atoms [i] / 2.0;
					break;
				case 4:
					double_bond_equivalent += num_atoms [i];
					break;
			}
		}
	}
	return double_bond_equivalent;
}
void get_elemental_info ( StringVector* pass_element_table, DoubleVectorVector* pass_isotope_mass, DoubleVectorVector* pass_isotope_abundance )
{
	initialise ();

	*pass_element_table = elementTable;
	*pass_isotope_mass = isotopeMassArray;
	*pass_isotope_abundance = isotopeAbundanceArray;
}
void get_element_table ( StringVector* pass_element_table )
{
	initialise ();

	*pass_element_table = elementTable;
}
int* formula_to_multiplier_list ( ElementalFormula& ef )
{
	int* element_multiplier_list = new int [elementTable.size ()];
	fill ( element_multiplier_list, element_multiplier_list + elementTable.size (), 0 );

	for ( ef.first () ; ef.isDone () ; ef.next () ) {
		for ( StringVectorSizeType i = 0 ; i < elementTable.size () ; i++ ) {
			if ( elementTable [i] == ef.element () ) {
				element_multiplier_list [i] += ef.multiplier ();
				break;
			}
		}
	}
	return element_multiplier_list;
}
static void initialise_mass_calculations ()
{
	int numValidElements;
	char* elementalInfo = getParamsFileInfo ( "elements.txt", &numValidElements );

	for ( int i = 0 ; i < numValidElements ; i++ ) {
		elementTable.push_back ( ( i == 0 ) ? strtok ( elementalInfo, " " ) : strtok ( NULL, " " ) );
		valencyTable.push_back ( atoi ( strtok ( NULL, " " ) ) );
		int numIsotopes = atoi ( strtok ( NULL, " " ) );
		DoubleVector iMass;
		DoubleVector iAbundance;

		double averageMass = 0.0;
		for ( int j = 0 ; j < numIsotopes ; j++ ) {
			double mass = atof ( strtok ( NULL, " " ) );
			double abundance = atof ( ( j == numIsotopes - 1 ) ?
				strtok ( NULL, "\n" ) : strtok ( NULL, " " ) );
			iMass.push_back ( mass );
			iAbundance.push_back ( abundance );
			averageMass += mass * abundance;
		}
		isotopeMassArray.push_back ( iMass );
		isotopeAbundanceArray.push_back ( iAbundance );
		averageMassArray.push_back ( averageMass );
		monoisotopicMassArray.push_back ( iMass [0] );
		idealMonoisotopicMassArray.push_back ( iMass [0] );
		nominalMassArray.push_back ( floor ( iMass [0] + 0.5 ) );
	}
	initialised = true;
}
static void replaceElementMasses ( const string& element, const DoubleVector& masses, const DoubleVector& abundances )
{
	StringVectorIterator iter = find ( elementTable.begin (), elementTable.end (), element );
	if ( iter == elementTable.end () ) {
		throw runtime_error ( "Element not available in element.txt file." );
	}
	else {
		int index = iter - elementTable.begin ();
		double averageMass = 0.0;
		for ( DoubleVectorSizeType i = 0 ; i < masses.size () ; i++ ) {
			averageMass += masses [i] * abundances [i];
		}
		isotopeMassArray [index] = masses;
		isotopeAbundanceArray [index] = abundances;
		averageMassArray [index] = averageMass;
		monoisotopicMassArray [index] = masses [0];
		nominalMassArray [index] = floor ( masses [0] + 0.5 );
	}
}
static void addElementToTable ( const string& element, int valency, const DoubleVector& masses, const DoubleVector& abundances )
{
	elementTable.push_back ( element );
	valencyTable.push_back ( valency );

	double averageMass = 0.0;
	for ( DoubleVectorSizeType i = 0 ; i < masses.size () ; i++ ) {
		averageMass += masses [i] * abundances [i];
	}
	isotopeMassArray.push_back ( masses );
	isotopeAbundanceArray.push_back ( abundances );
	averageMassArray.push_back ( averageMass );
	monoisotopicMassArray.push_back ( masses [0] );
	idealMonoisotopicMassArray.push_back ( masses [0] );
	nominalMassArray.push_back ( floor ( masses [0] + 0.5 ) );
}
class CheckMass {
	int iMass;
public:
	CheckMass ( int iMass ) :
		iMass ( iMass ) {}
	bool operator () ( double dMass )
	{
		return iMass == floor ( dMass + 0.5 );
	}
};
static double getAccurateElementMass ( const string& element, int intMass )
{
	StringVectorIterator iter = find ( elementTable.begin (), elementTable.end (), element );
	if ( iter != elementTable.end () ) {
		DoubleVector isoMasses = isotopeMassArray [iter - elementTable.begin ()];
		DoubleVectorIterator iter2 = find_if ( isoMasses.begin (), isoMasses.end (), CheckMass ( intMass ) );
		if ( iter2 != isoMasses.end () ) {
			return isoMasses [iter2-isoMasses.begin ()];
		}
		else throw runtime_error ( "Requested mass not available in element.txt file." );
	}
	else throw runtime_error ( "Element not available in element.txt file." );
}
static void replaceElement ( const string& hElement, const string& element, int nominalLight, int nominalHeavy, double percent )
{
	DoubleVector masses;
	DoubleVector abundances;
	masses.push_back ( getAccurateElementMass ( element, nominalHeavy ) );
	masses.push_back ( getAccurateElementMass ( element, nominalLight ) );
	double ratio = percent / 100.0;
	abundances.push_back ( ratio );
	abundances.push_back ( 1.0 - ratio );
	replaceElementMasses ( hElement, masses, abundances );
}
