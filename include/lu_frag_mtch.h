/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_frag_mtch.h                                                *
*                                                                             *
*  Created    : June 20th 1996                                                *
*                                                                             *
*  Purpose    : msfit style search algorithms.                                *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1996-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_frag_mtch_h
#define __lu_frag_mtch_h

#include <vector>

#include <lu_prog.h>
#include <lu_charge.h>

struct HitStats;
class MSFitParameters;
class ModificationTable;

class MowseScore {
	DoubleVectorVector mowseScore;
	double mowsePFactor;
	DoubleVector* mowse;
	bool statsNeedUpdating;
	static const char SCORE_MATCH;
	static const double MAX_MOWSE_PROTEIN_MASS;
	static const double MOWSE_PROTEIN_AMU_BIN;
	static const int MAX_MOWSE_PROTEIN_BIN_INDEX;
	static const int NUM_MOWSE_PROTEIN_BINS;

	static const double MAX_MOWSE_PEPTIDE_MASS;
	static const double MOWSE_PEPTIDE_AMU_BIN;
	static const int NUM_MOWSE_PEPTIDE_BINS;
	void assembleMowseStats ();
public:
	MowseScore ( double mowsePFactor );
	void setMowseArray ( double proteinMW );
	void accumulateMowseScore ( double fragmentMass, bool singleCleavage );
	void calculateMowseScores ( HitStats* hs, double proteinMW, const DoubleVector& peakMass );
};

class MSFitSearch {
	static int maxFitPeaks;
protected:
	static const char SCORE_MATCH;
	static const char NO_SCORE_MATCH;

	PeakContainer peaks;
	int numPeaks;
	DoubleVector peakMassLowerBound;
	DoubleVector tolerance;
	DoubleVector peakMass;
	CharVector massMatched;
	double lowMass;
	double highMass;
	int numScoreMatches;
	MowseScore* mowseScore;
	CharVector mowseMissedCleavages;
	unsigned int compMask;
	bool compMaskTypeAnd;
	bool compMaskTypeOr;
	int missedCleavages;
	bool checkComposition ( const char* fragment, int len );
	void init ( const MSFitParameters& params );
public:
	MSFitSearch ( MSDataPoint* dataSet, const MSFitParameters& params );
	virtual ~MSFitSearch ();
	PeakContainer getPeaks () const { return ( peaks );};
	HitStats getProteinFragmentStats ();
	void calculateMowseScores ( HitStats* hs, double protein_mw );
	virtual int matchFragments ( char* protein, const IntVector& cleavageIndex );
	virtual ModificationTable* getModificationTable () const { return 0; }
};
typedef std::vector <MSFitSearch*> MSFitSearchPtrVector;
typedef MSFitSearchPtrVector::size_type MSFitSearchPtrVectorSizeType;

class MSFitModifiedSearch : public MSFitSearch {
	bool pyroglutamicAcidFlag;
	bool user1Flag;
	bool user2Flag;
	bool user3Flag;
	bool user4Flag;
	bool oxidationFlag;
	bool acetylationFlag;
	bool incompleteCysFlag;
	double highModifiedMass;
	unsigned int acetylationMask;
	unsigned int pyroglutamicAcidMask;
	unsigned int oxidationMask;
	unsigned int cysMask;
	unsigned int user1Mask;
	unsigned int user2Mask;
	unsigned int user3Mask;
	unsigned int user4Mask;
	double pyroglutamicAcidMod;
	double acetylationMod;
	double oxidationMod;
	double incompleteCysMod;
	double user1Mod;
	double user2Mod;
	double user3Mod;
	double user4Mod;
	void init ( const MSFitParameters& params );
public:
	MSFitModifiedSearch ( MSDataPoint* dataSet, const MSFitParameters& params );
	~MSFitModifiedSearch ();
	int matchFragments ( char* protein, const IntVector& cleavageIndex );
};

class MSFitAllowErrorsSearch : public MSFitSearch {
	static ModificationTable* modificationTable;
	static double maxParentError;
	int minParentIonMatches;
	char enzymeTerminalSpecificity;
	void init ( const MSFitParameters& params );
public:
	MSFitAllowErrorsSearch ( MSDataPoint* dataSet, const MSFitParameters& params );
	~MSFitAllowErrorsSearch ();
	int matchFragments ( char* protein, const IntVector& cleavageIndex );
	ModificationTable* getModificationTable () const { return modificationTable; }
};

MSFitSearch* getMSFitSearch ( MSDataPoint* dataSet, const MSFitParameters& params );

#endif /* ! __lu_frag_mtch_h */
