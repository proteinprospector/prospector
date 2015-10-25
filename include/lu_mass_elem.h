/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_mass_elem.h                                                *
*                                                                             *
*  Created    : May 30th 2001                                                 *
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
*  Copyright (2001-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_mass_elem_h
#define __lu_mass_elem_h

#include <lgen_define.h>
#include <lu_formula.h>

class IsotopePurity {
	double percentC13;
	double percentN15;
	double percentO18;
public:
	IsotopePurity ( const ParameterList* params );
	static void copyToCGI ( std::ostream& os, const ParameterList* params );
	double getPercentC13 () const { return percentC13; }
	double getPercentN15 () const { return percentN15; }
	double getPercentO18 () const { return percentO18; }
};

double formula_to_monoisotopic_mass ( const char* formula );
double formula_to_monoisotopic_mass ( ElementalFormula& elemForm );
double formula_to_ideal_monoisotopic_mass ( const char* formula );
double formula_to_nominal_mass ( const char* formula );
double formula_to_nominal_mass ( ElementalFormula& elemForm );
double formula_to_average_mass ( const char* formula );
double formula_to_average_mass ( ElementalFormula& elemForm );
double formula_to_percent_element_monoisotopic ( ElementalFormula& formula, const char* element );
double formula_to_percent_element_average ( ElementalFormula& formula, const char* element );
double formula_to_double_bond_equivalent ( ElementalFormula& formula );
double num_atoms_to_double_bond_equivalent ( int* num_atoms );
void get_elemental_info ( StringVector* pass_element_table, DoubleVectorVector* pass_isotope_mass, DoubleVectorVector* pass_isotope_abundance );
void get_element_table ( StringVector* pass_element_table );
int* formula_to_multiplier_list ( ElementalFormula& ef );

#endif /* ! __lu_mass_elem_h */
