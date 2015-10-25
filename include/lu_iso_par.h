/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_iso_par.h                                                  *
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
*  Copyright (2001-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_iso_par_h
#define __lu_iso_par_h

#include <lu_mass_elem.h>
#include <lu_prog_par.h>
#include <lu_mass.h>

class IsotopePeakStats;
class IsotopeProfile;

class MSIsotopeParameters : public MSProgramParameters {
	AAInitInfo aaInitInfo;
	StringVector distributionType;
	StringVector outputString;
	IntVector parentCharge;
	StringVector intensityString;
	std::vector <IsotopePeakStats*> isotopePeakStats;
	bool displayGraph;
	bool detailedReport;
	std::string profileType;
	double resolution;
	IsotopePurity isotopePurity;
	IsotopeProfile* isotopeProfile;
	DoubleVector intensityVector;
public:
	MSIsotopeParameters ( const ParameterList* params );
	void initIsotopePeakStats ( const ParameterList* params );
	IsotopeProfile* initIsotopeProfile ();

	AAInitInfo getAAInitInfo () const { return aaInitInfo; }
	double getResolution () const { return resolution; }
	std::vector <IsotopePeakStats*> getIsotopePeakStats () const { return isotopePeakStats; }
	bool getDisplayGraph () const { return displayGraph; }
	bool getDetailedReport () const { return detailedReport; }
	IsotopeProfile* getIsotopeProfile () const { return isotopeProfile; }

	void printHTML ( std::ostream& os ) const;
};

class MSIsotopeLink {
	void putCGI ( std::ostream& os ) const;
public:
	MSIsotopeLink () {}
	void write ( std::ostream& os, ElementalFormula& elementalComposition, const int charge ) const;
	void printHTML ( std::ostream& os ) const;
};


#endif /* ! __lu_iso_par_h */
