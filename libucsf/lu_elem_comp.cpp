/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_elem_comp.cpp                                              *
*                                                                             *
*  Created    : September 30th 1997                                           *
*                                                                             *
*  Purpose    : Elemental composition searches.                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1997-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <cmath>
#include <cstdio>
#include <lg_new.h>
#include <lu_mass_elem.h>

static bool initialised = false;
static char** hit;
static double continue_mass;
static double start_mass;
static double end_mass;
static StringVector elementTable;
static DoubleVector max_percent;
static DoubleVector min_percent;
static IntVector max_num_atoms;
static IntVector min_num_atoms;
static int* num_atoms;
static int num_valid_elements;
static int num;
static int max_hits;
static int max_hits_exceeded;
static int nitrogen_index;
static DoubleVector wt;
static double max_dbe;
static double min_dbe;
static int dbe_checks;

static void get_combinations ( int i, double sum );
static void add_hit ();
static void initialise_elem_combs ();
static inline void initialise ();

static inline void initialise ()
{
	if ( initialised == false ) initialise_elem_combs ();
}
int calculate_elem_combinations ( const StringVector& superElementList, double mass, double tolerance, bool monoisotopicFlag, const char* start_formula, int maxReportedHits, char*** combinations, int* num_combinations )
{
	double ( *mass_convert ) (const char*);
	if ( monoisotopicFlag ) mass_convert = formula_to_monoisotopic_mass;
	else mass_convert = formula_to_average_mass;
	int i, j;
	double initial_mass;
	double max_dbe_per_mass;
	double min_dbe_per_mass;
	double ( *calc_percent_element ) ( ElementalFormula& formula, const char* );
	int numSuperElements = superElementList.size ();

	if ( monoisotopicFlag ) calc_percent_element = formula_to_percent_element_monoisotopic;
	else calc_percent_element = formula_to_percent_element_average;

	initialise ();

	for ( i = 0 ; i < num_valid_elements ; i++ ) {
		wt [i] = mass_convert ( elementTable [i].c_str () );
	}

	initial_mass = mass_convert ( start_formula );

	for ( i = 0 ; i < num_valid_elements ; i++ ) {
		if ( numSuperElements == 0 ) {
			min_percent [i] = 0.0;
			max_percent [i] = 100.0;
		}
		else {
			for ( j = 0, max_percent [i] = 0.0, min_percent [i] = 100.0 ; j < numSuperElements ; j++ ) {
				ElementalFormula superElemForm ( superElementList [j] );
				double percent = calc_percent_element ( superElemForm, elementTable [i].c_str () );
				max_percent [i] = genMax ( max_percent [i], percent );
				min_percent [i] = genMin ( min_percent [i], percent );
			}
		}
		max_num_atoms [i] = (int) ( ( mass - initial_mass + tolerance ) * ( max_percent [i] ) / wt [i] );
		min_num_atoms [i] = (int) ( ( mass - initial_mass - tolerance ) * ( min_percent [i] ) / wt [i] );
	}

	dbe_checks = ( numSuperElements != 0 );
	for ( i = 0 ; i < numSuperElements ; i++ ) {
		ElementalFormula superElemForm ( superElementList [i] );
		double dbe_per_mass = ( formula_to_double_bond_equivalent ( superElemForm ) - 1.0 ) / mass_convert ( superElementList [i].c_str () );
		max_dbe_per_mass = ( i == 0 ) ? dbe_per_mass : genMax ( max_dbe_per_mass, dbe_per_mass );
		min_dbe_per_mass = ( i == 0 ) ? dbe_per_mass : genMin ( min_dbe_per_mass, dbe_per_mass );
	}
	ElementalFormula startElemForm ( start_formula );
	if ( numSuperElements ) {
		max_dbe = ( mass + tolerance - initial_mass ) * max_dbe_per_mass;
		min_dbe = ( mass - tolerance - initial_mass ) * min_dbe_per_mass;
		max_dbe += formula_to_double_bond_equivalent ( startElemForm );
		min_dbe += formula_to_double_bond_equivalent ( startElemForm );
	}

	num_atoms = formula_to_multiplier_list ( startElemForm );
	for ( i = 0 ; i < num_valid_elements ; i++ ) {
		initial_mass += wt [i] * min_num_atoms [i];
		max_num_atoms [i] += num_atoms [i];
		min_num_atoms [i] += num_atoms [i];
		num_atoms [i] = min_num_atoms [i];
	}

	hit = new char* [maxReportedHits + 2];
	max_hits = maxReportedHits;
	max_hits_exceeded = 0;

	num = 0;
	start_mass = mass - tolerance;
	end_mass = mass + tolerance;

	continue_mass = mass + tolerance - floor ( mass_convert ( "H" ) );

	if ( initial_mass >= start_mass && initial_mass <= end_mass )
		add_hit ();
	else
		get_combinations ( 0, initial_mass );

	*combinations = hit;
	*num_combinations = num;

	delete [] num_atoms;

	return ( max_hits_exceeded );
}
static void get_combinations ( int i, double sum )
{
	for ( ; i < num_valid_elements ; i++ ) {
		if ( num_atoms [i] >= max_num_atoms [i] ) continue;
		(num_atoms [i])++;
		sum += wt [i];
		if ( sum < continue_mass ) get_combinations ( i, sum );
		else if ( sum < end_mass && sum > start_mass ) add_hit ();
		sum -= wt [i];
		(num_atoms [i])--;
	}
}
static void add_hit ()
{
	char temp_string [100];
	int nominal_mass;
	int num_nitrogens;
	int i;

	if ( num < max_hits ) {
		for ( i = 0, temp_string [0] = 0 ; i < num_valid_elements ; i++ ) {
			if ( num_atoms [i] ) {
				if ( *temp_string ) strcat ( temp_string, " " );
				sprintf ( temp_string, "%s%s%d", temp_string, elementTable [i].c_str (), num_atoms [i] );
			}
		}
		nominal_mass = (int) formula_to_nominal_mass ( temp_string );
		num_nitrogens = num_atoms [nitrogen_index];
		if ( ( genIsEven ( nominal_mass ) && genIsOdd ( num_nitrogens ) ) || ( genIsOdd ( nominal_mass ) && genIsEven ( num_nitrogens ) ) ) {
			if ( dbe_checks ) {
				double dbe = num_atoms_to_double_bond_equivalent ( num_atoms );
				if ( dbe > min_dbe && dbe < max_dbe )
					hit [num++] = gen_new_string ( temp_string );
			}
			else hit [num++] = gen_new_string ( temp_string );
		}
	}
	else max_hits_exceeded = 1;	
}
static void initialise_elem_combs ()
{
	get_element_table ( &elementTable );
	num_valid_elements = elementTable.size ();

	for ( int i = 0 ; i < num_valid_elements ; i++ ) {
		if ( elementTable [i] == "N" ) nitrogen_index = i;
	}

	wt.resize ( num_valid_elements );
	max_percent.resize ( num_valid_elements );
	min_percent.resize ( num_valid_elements );
	max_num_atoms.resize ( num_valid_elements );
	min_num_atoms.resize ( num_valid_elements );

	initialised = true;
}
