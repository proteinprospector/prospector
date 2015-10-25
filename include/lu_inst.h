/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_inst.h                                                     *
*                                                                             *
*  Created    : November 3rd 1997                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1997-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_inst_h
#define __lu_inst_h

#include <lgen_define.h>
#include <lu_fragmentation.h>

#ifdef LUCSF_INST_MAIN
#define LUCSF_INST_EXTERN
#else
#define LUCSF_INST_EXTERN extern
#endif

class PeakPrecision {
	int massDecimalPlaces;
	int errorSigFig;
	int intensitySigFig;
public:
	void setMassDecimalPlaces	( int n ) { massDecimalPlaces = n; }
	void setErrorSigFig			( int n ) { errorSigFig = n; }
	void setIntensitySigFig		( int n ) { intensitySigFig = n; }
	int getMassDecimalPlaces () const { return massDecimalPlaces; }
	int getErrorSigFig () const { return errorSigFig; }
	int getIntensitySigFig () const { return intensitySigFig; }
};

class InstrumentInfo {
	std::string name;
	int parentPrecision;
	int parentErrorSigFig;
	int parentIntensitySigFig;
	int fragmentPrecision;
	int fragmentErrorSigFig;
	int fragmentIntensitySigFig;
	PeakPrecision parentPeakPrecision;
	PeakPrecision fragmentPeakPrecision;

	double quanTolerance;

	MSMSFragmentation* frag;

	bool allowIncorrectCharge;
	
	void setPrecisions ();
public:
	InstrumentInfo ( const std::string& instrumentName );
	~InstrumentInfo ();

	std::string getName () const { return name; }
	std::string getFragName () const { return frag->getFragName (); }
	unsigned int getAmmoniaLossMask () const;
	unsigned int getWaterLossMask () const;
	unsigned int getPosChargeBearingMask () const;
	unsigned int getDIonExcludeMask () const;
	unsigned int getWIonExcludeMask () const;
	unsigned int getVIonExcludeMask () const;
	StringVector getIonTypes () const;
	DoubleVector getLossMasses () const;

	PeakPrecision getParentPeakPrecision () const { return parentPeakPrecision; }
	PeakPrecision getFragmentPeakPrecision () const { return fragmentPeakPrecision; }

	double getMaximumInternalIonMass () const;

	ScoreType getUnmatchedScore () const;

	ScoreType getImmoniumScore () const;
	ScoreType getRelatedIonScore () const;
	ScoreType getMScore () const;

	ScoreType getAScore () const;
	ScoreType getALossScore () const;
	ScoreType getAPhosLossScore () const;
	ScoreType getBScore () const;
	ScoreType getBPlusH2OScore () const;
	ScoreType getBLossScore () const;
	ScoreType getBPhosLossScore () const;
	ScoreType getCLadderScore () const;
	ScoreType getCPlus2DaScore () const;
	ScoreType getCPlus1DaScore () const;
	ScoreType getCScore () const;
	ScoreType getCMinus1DaScore () const;
	ScoreType getDScore () const;

	ScoreType getVScore () const;
	ScoreType getWScore () const;
	ScoreType getXScore () const;
	ScoreType getNLadderScore () const;
	ScoreType getYScore () const;
	ScoreType getYLossScore () const;
	ScoreType getYPhosLossScore () const;
	ScoreType getBigYScore () const;
	ScoreType getZScore () const;
	ScoreType getZPlus1DaScore () const;
	ScoreType getZPlus2DaScore () const;
	ScoreType getZPlus3DaScore () const;

	ScoreType getBP2Score () const;
	ScoreType getBP2LossScore () const;
	ScoreType getBP2PhosLossScore () const;

	ScoreType getYP2Score () const;
	ScoreType getYP2LossScore () const;
	ScoreType getYP2PhosLossScore () const;

	ScoreType getInternalAScore () const;
	ScoreType getInternalBScore () const;
	ScoreType getInternalLossScore () const;

	ScoreType getMH3PO4Score () const;
	ScoreType getMSOCH4Score () const;

	bool getScoring () const;

	bool getAllowIncorrectCharge ()	const { return allowIncorrectCharge; }

	bool getChargeReducedFragmentation () const;
	bool getETDTypeFragmentation () const;
};

LUCSF_INST_EXTERN InstrumentInfo* instInf;

void initialiseInstrumentName ( const std::string& instrumentName );
void resetInstrumentName ( const std::string& instrumentName );

class InstrumentList {
	StringVector names;
	void initialise ();
	InstrumentList ();
public:
	static InstrumentList& instance ();
	StringVector getNames () const { return names; }
};

#endif /* ! __lu_inst_h */
