 /*****************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_aa_calc.h                                                  *
*                                                                             *
*  Created    : September 18th 2001                                           *
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

#ifndef __lu_aa_calc_h
#define __lu_aa_calc_h

#include <string>
#include <lu_formula.h>
#include <lu_const_mod.h>

struct AACalculatorNoElementalComposition {};

class AACalculator {
	enum { AA_ARRAY_SIZE = 128 };
	bool localMonoisotopicFlag;
	double aminoAcidWt [AA_ARRAY_SIZE];
	typedef double ( *MassConvert ) (ElementalFormula&);
	MassConvert massConvert;
	ElementalFormula elemFormula [AA_ARRAY_SIZE];
	bool validTerminalFormula;
	ElementalFormula nTermFormula;
	ElementalFormula cTermFormula;
	ElementalFormula cationFormula;
	ElementalFormula terminalFormula;
	double nTerminusWt;
	double cTerminusWt;
	double cationWt;
	double terminalWt;
	std::string nTerminusString;
	std::string cTerminusString;
	double getAminoAcidWt ( const std::string& modAA ) const;
public:
	AACalculator ( bool monoisotopicFlag, const MapStringConstModPtr& constMods );
	double calculatePeptideMW ( const std::string& peptide ) const;
	double calculatePeptideMW ( const StringVector& peptide ) const;
	bool calculateStrippedElementalCompositionWithTerminii ( const std::string& peptide, const std::string& nterm, const std::string& cterm, const std::string& nloss, ElementalFormula& ef, bool falseIfMM = true ) const;
	ElementalFormula calculateElementalCompositionWithTerminii ( const std::string& peptide, const std::string& nterm, const std::string& cterm, const std::string& nloss ) const;
	ElementalFormula calculateElementalComposition ( const std::string& peptide ) const;
	ElementalFormula calculateElementalComposition ( const StringVector& peptide ) const;
	ElementalFormula calculateFragElementalComposition ( const std::string& peptide ) const;
	ElementalFormula calculateFragElementalComposition ( const StringVector& peptide ) const;
	static AAFormula calculateAAComposition ( const std::string& peptide );
	double getNTerminusWt () const { return nTerminusWt; }
	double getCTerminusWt () const { return cTerminusWt; }
	std::string getNTerminusString () const { return nTerminusString; }
	std::string getCTerminusString () const { return cTerminusString; }
};

#endif /* ! __lu_aa_calc_h */
