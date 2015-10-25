/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_patt_par.h                                                 *
*                                                                             *
*  Created    : June 14th 2001                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2008) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_patt_par_h
#define __lu_patt_par_h

#include <lu_srch_par.h>

class MSPatternParameters : public MSSearchParameters {

	std::string regularExpression;
	int maxAASubstitutions;
	int maxPeptideHits;

	std::string massType;
	bool monoisotopicFlag;
public:
	MSPatternParameters ( const ParameterList* params );

	std::string getRegularExpression () const { return regularExpression; }
	int getMaxAASubstitutions () const { return maxAASubstitutions; }
	int getMaxPeptideHits () const { return maxPeptideHits; }

	std::string getMassType () const { return massType; }
	bool getMonoisotopicFlag () const { return monoisotopicFlag; }

	void printHTML ( std::ostream& os ) const;
};

#endif /* ! __lu_patt_par_h */
