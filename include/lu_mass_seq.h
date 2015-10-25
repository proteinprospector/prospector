/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_mass_seq.h                                                 *
*                                                                             *
*  Created    : May 30th 2001                                                 *
*                                                                             *
*  Purpose    : Functions concerned with peptide/protein formulae             *
*               conversions.                                                  *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2010) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_mass_seq_h
#define __lu_mass_seq_h

#include <lgen_define.h>
#include <lu_formula.h>

class Mass {
protected:
	double mass;
public:
	Mass () {};
	Mass ( double m ) : mass ( m )  {}
	Mass& operator= ( double rhs ) { mass = rhs; return *this; }
	double getMass () const {return mass; }
	virtual std::ostream& print ( std::ostream& os ) const = 0;
	std::ostream& printMass ( std::ostream& os, int precision ) const;
};
inline bool operator< ( const Mass& lhs, const Mass& rhs )
	{ return lhs.getMass () < rhs.getMass (); }
std::ostream& operator<< ( std::ostream& os, const Mass& m );

class ProteinMW : public Mass {
	static ElementalFormula proteinNTerminus;
	static ElementalFormula proteinCTerminus;
	static double mwTerminalWt;
	static DoubleVector mw_amino_acid_wt;
public:
	ProteinMW () {};
	ProteinMW ( const char* protein_string );
	std::ostream& print ( std::ostream& os ) const;
	static void initialise ( const char* mw_aa, const char* mw_aa_file );
	static int mwPrecision;
};

double peptide_formula_to_molecular_weight ( const char* peptideString );
bool replaceBandZ ( std::string& s );
StringVector initSequence ( const std::string& seq );

#endif /* ! __lu_mass_seq_h */
