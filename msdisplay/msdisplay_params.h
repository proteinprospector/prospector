/******************************************************************************
*                                                                             *
*  Program    : msdisplay_params.h                                            *
*                                                                             *
*  Filename   : msdisplay_params.h                                            *
*                                                                             *
*  Created    : December 20th 2002                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2002-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __msdisplay_params_h
#define __msdisplay_params_h

#include <lu_prog_par.h>
#include <lu_mass.h>
#include <lu_mass_elem.h>
#include <lu_parent.h>
#include <lu_quan_multi.h>
#include <lu_quan_ratio.h>
#include <lu_spec_id.h>

class LinkInfo;

class MSDisplayParameters : public MSProgramParameters {
	AAInitInfo aaInitInfo;
	std::string quanType;
	std::string rawType;
	IsotopePurity isotopePurity;
	PurityCorrection* purityCorrection;
	double resolution;
	int numPeaks;
	std::string project;
	SpecID specID;
	double mOverZ;
	std::string rtIntervalStart;
	std::string rtIntervalEnd;
	double snrThreshold;
	double areaThreshold;
	double intensityThreshold;
	int charge;
	std::string formula;
	std::string nterm;
	std::string cterm;
	std::string nloss;
	std::string formula2;
	std::string nterm2;
	std::string cterm2;
	LinkInfo* linkInfo;
	QuanPeptide qPeptide;

	double displayStartMass;
	double displayEndMass;
	bool displayAllMasses;
	std::string msType;
	std::string msmsInfo;
	double reporterIonWindow;
public:
	MSDisplayParameters ( const ParameterList* params );
	const AAInitInfo& getAAInitInfo () const { return aaInitInfo; }
	std::string getQuanType () const { return quanType; }
	std::string getRawType () const { return rawType; }
	IsotopePurity getIsotopePurity () const { return isotopePurity; }
	double getResolution () const { return resolution; }
	int getNumPeaks () const { return numPeaks; }
	std::string getProject () const { return project; }
	SpecID getSpecID () const { return specID; }
	double getMOverZ () const { return mOverZ; }
	std::string getRTIntervalStart () const { return rtIntervalStart; }
	std::string getRTIntervalEnd () const { return rtIntervalEnd; }
	double getSNRThreshold () const { return snrThreshold; }
	double getAreaThreshold () const { return areaThreshold; }
	double getIntensityThreshold () const { return intensityThreshold; }
	int getCharge () const { return charge; }
	std::string getFormula () const { return formula; }
	std::string getNTerm () const { return nterm; }
	std::string getCTerm () const { return cterm; }
	std::string getNLoss () const { return nloss; }
	std::string getFormula2 () const { return formula2; }
	std::string getNTerm2 () const { return nterm2; }
	std::string getCTerm2 () const { return cterm2; }
	bool isXlink () const { return !formula2.empty (); }
	std::string getBridgeFormula () const;
	const QuanPeptide& getQPeptide () const { return qPeptide; }

	double getDisplayStartMass () const { return displayStartMass; }
	double getDisplayEndMass () const { return displayEndMass; }
	bool getDisplayAllMasses () const { return displayAllMasses; }
	int getMSType () const { return msType == "MS" ? 0 : 1; }
	std::string getMSMSInfo () const { return msmsInfo; }
};

#endif /* ! __msdisplay_params_h */
