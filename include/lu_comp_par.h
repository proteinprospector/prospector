/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_comp_par.h                                                 *
*                                                                             *
*  Created    : June 9th 2001                                                 *
*                                                                             *
*  Purpose    :                                                               *
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

#ifndef __lu_comp_par_h
#define __lu_comp_par_h

#include <lu_prog_par.h>
#include <lu_tol.h>
#include <lu_comb_perm.h>
#include <lu_composit.h>
#include <lu_mass.h>

class MSCompParameters : public MSProgramParameters {
	AAInitInfo aaInitInfo;
	double parentMass;
	int parentCharge;
	std::string combinationType;
	ConsideredAA consideredAA; 
	int maxReportedHits;
	ToleranceInfo parentMassTolerance;
	MassInfo massInfo;
	CompositionSearchParameters compSearchParams;
	StringVector requestedIonTypes;
public:
	MSCompParameters ( const ParameterList* params );

	AAInitInfo getAAInitInfo () const { return aaInitInfo; }
	double getParentMass () const { return parentMass; }
	int getParentCharge () const { return parentCharge; }
	std::string getCombinationType () const { return combinationType; }
	std::string getAAList () const { return consideredAA.getAAList (); }
	int getMaxReportedHits () const { return maxReportedHits; }
	Tolerance* getParentMassTolerance () const { return parentMassTolerance.getTolerance (); }
	bool getMonoisotopicFlag () const { return massInfo.getMonoisotopicFlag (); }
	CompositionSearchParameters getCompSearchParams () const { return compSearchParams; }
	StringVector getCompIons () const { return compSearchParams.getCompIons (); }
	StringVector getRequestedIonTypes () const { return requestedIonTypes; }

	void printHTML ( std::ostream& os ) const;
};

#endif /* ! __lu_comp_par_h */
