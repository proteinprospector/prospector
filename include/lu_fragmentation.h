/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_fragmentation.h                                            *
*                                                                             *
*  Created    : May 14th 2007                                                 *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2007-2009) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_fragmentation_h
#define __lu_fragmentation_h

#include <string>
#include <lg_io.h>
#include <lgen_define.h>

class ParameterList;

typedef int ScoreType;
static const double SCORE_TYPE_MULTIPLIER = 1000.0;	//Should be 1 for floating point
static const ScoreType DEFAULT_HIT_SCORE = static_cast<ScoreType> ( SCORE_TYPE_MULTIPLIER );
static const ScoreType DEFAULT_MISS_SCORE = 0;
typedef std::vector <ScoreType> ScoreTypeVector;

class MSMSFragmentation {
	std::string fragName;
	std::string ammoniaLoss;
	std::string waterLoss;
	std::string posChargeBearing;
	std::string dIonExclude;
	std::string vIonExclude;
	std::string wIonExclude;
	unsigned int ammoniaLossMask;
	unsigned int waterLossMask;
	unsigned int posChargeBearingMask;
	unsigned int dIonExcludeMask;
	unsigned int vIonExcludeMask;
	unsigned int wIonExcludeMask;

	ScoreType unmatchedScore;

	ScoreType immoniumScore;
	ScoreType relatedIonScore;
	ScoreType mScore;

	ScoreType aScore;
	ScoreType aLossScore;
	ScoreType aPhosLossScore;
	ScoreType bScore;
	ScoreType bPlusH2OScore;
	ScoreType bLossScore;
	ScoreType bPhosLossScore;
	ScoreType cLadderScore;
	ScoreType cPlus2DaScore;
	ScoreType cPlus1DaScore;
	ScoreType cScore;
	ScoreType cMinus1DaScore;
	ScoreType dScore;

	ScoreType vScore;
	ScoreType wScore;
	ScoreType xScore;
	ScoreType nLadderScore;
	ScoreType yScore;
	ScoreType yLossScore;
	ScoreType yPhosLossScore;
	ScoreType bigYScore;
	ScoreType zScore;
	ScoreType zPlus1DaScore;
	ScoreType zPlus2DaScore;
	ScoreType zPlus3DaScore;

	ScoreType bp2Score;
	ScoreType bp2LossScore;
	ScoreType bp2PhosLossScore;

	ScoreType yp2Score;
	ScoreType yp2LossScore;
	ScoreType yp2PhosLossScore;

	ScoreType internalAScore;
	ScoreType internalBScore;
	ScoreType internalLossScore;

	ScoreType MH3PO4Score;
	ScoreType MSOCH4Score;

	bool scoring;

	StringVector ionTypes;

	double maximumInternalIonMass;

	bool chargeReducedFragmentation;

	void setMasks ();
	bool getScoreValue ( const ParameterList& params, const std::string& name, ScoreType& value );
public:
	MSMSFragmentation ( const std::string& fragName );

	std::string getFragName () const { return fragName; }

	unsigned int getAmmoniaLossMask ()		const { return ammoniaLossMask; }
	unsigned int getWaterLossMask ()		const { return waterLossMask; }
	unsigned int getPosChargeBearingMask ()	const { return posChargeBearingMask; }
	unsigned int getDIonExcludeMask ()		const { return dIonExcludeMask; }
	unsigned int getWIonExcludeMask ()		const { return wIonExcludeMask; }
	unsigned int getVIonExcludeMask ()		const { return vIonExcludeMask; }

	ScoreType getUnmatchedScore ()		const { return unmatchedScore; }

	ScoreType getImmoniumScore ()		const { return immoniumScore; }
	ScoreType getRelatedIonScore ()		const { return relatedIonScore; }
	ScoreType getMScore ()				const { return mScore; }

	ScoreType getAScore ()				const { return aScore; }
	ScoreType getALossScore ()			const { return aLossScore; }
	ScoreType getAPhosLossScore ()		const { return aPhosLossScore; }
	ScoreType getBScore ()				const { return bScore; }
	ScoreType getBPlusH2OScore ()		const { return bPlusH2OScore; }
	ScoreType getBLossScore ()			const { return bLossScore; }
	ScoreType getBPhosLossScore ()		const { return bPhosLossScore; }
	ScoreType getCLadderScore ()		const { return cLadderScore; }
	ScoreType getCPlus2DaScore ()		const { return cPlus2DaScore; }
	ScoreType getCPlus1DaScore ()		const { return cPlus1DaScore; }
	ScoreType getCScore ()				const { return cScore; }
	ScoreType getCMinus1DaScore ()		const { return cMinus1DaScore; }
	ScoreType getDScore ()				const { return dScore; }

	ScoreType getVScore ()				const { return vScore; }
	ScoreType getWScore ()				const { return wScore; }
	ScoreType getXScore ()				const { return xScore; }
	ScoreType getNLadderScore ()		const { return nLadderScore; }
	ScoreType getYScore ()				const { return yScore; }
	ScoreType getYLossScore ()			const { return yLossScore; }
	ScoreType getYPhosLossScore ()		const { return yPhosLossScore; }
	ScoreType getBigYScore ()			const { return bigYScore; }
	ScoreType getZScore ()				const { return zScore; }
	ScoreType getZPlus1DaScore ()		const { return zPlus1DaScore; }
	ScoreType getZPlus2DaScore ()		const { return zPlus2DaScore; }
	ScoreType getZPlus3DaScore ()		const { return zPlus3DaScore; }

	ScoreType getBP2Score ()			const { return bp2Score; }
	ScoreType getBP2LossScore ()		const { return bp2LossScore; }
	ScoreType getBP2PhosLossScore ()	const { return bp2PhosLossScore; }

	ScoreType getYP2Score ()			const { return yp2Score; }
	ScoreType getYP2LossScore ()		const { return yp2LossScore; }
	ScoreType getYP2PhosLossScore ()	const { return yp2PhosLossScore; }

	ScoreType getInternalAScore ()		const { return internalAScore; }
	ScoreType getInternalBScore ()		const { return internalBScore; }
	ScoreType getInternalLossScore ()	const { return internalLossScore; }

	ScoreType getMH3PO4Score ()	const { return MH3PO4Score; }
	ScoreType getMSOCH4Score ()	const { return MSOCH4Score; }

	bool getScoring ()				const { return scoring; }
	
	StringVector getIonTypes ()		const { return ionTypes; }

	double getMaximumInternalIonMass () const { return maximumInternalIonMass; }
	bool getChargeReducedFragmentation () const { return chargeReducedFragmentation; }
	DoubleVector getLossMasses ();
};

class FragmentationParams : public GenNameValueStream {
	FragmentationParams ();
public:
	static FragmentationParams& instance ( bool reset = false );
	static time_t getLastModifiedTime ();
};

#endif /* ! __lu_fragmentation_h */
