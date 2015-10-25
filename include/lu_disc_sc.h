/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_disc_sc.h                                                  *
*                                                                             *
*  Created    : December 4th 2003                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2003-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_disc_sc_h
#define __lu_disc_sc_h

#include <string>

class DiscriminantScore {
	double offsetCoeff;
	double bestScoreCoeff;
	double scoreDiffCoeff;
	double maxBestScore;

	double offsetCoeff2;
	double bestScoreCoeff2;
	double expectCoeff2;
	double maxBestScore2;
public:
	DiscriminantScore ( const std::string& discScoreName );
	double calculateDiscriminantScore ( double bestScore, double scoreDiff, double expect ) const;
	static const double MIN_DISC_SCORE;
	double getBestScoreCoeff () const { return bestScoreCoeff; }
	double getOffsetCoeff () const { return offsetCoeff; }
	double getMaxBestScore () const { return maxBestScore; }
	double getScoreDiffCoeff () const { return scoreDiffCoeff; }

	double getBestScoreCoeff2 () const { return bestScoreCoeff2; }
	double getOffsetCoeff2 () const { return offsetCoeff2; }
	double getMaxBestScore2 () const { return maxBestScore2; }
	double getExpectCoeff2 () const { return expectCoeff2; }
};

class DiscriminantScoreInstrumentList {
	StringVector names;
	void initialise ();
	DiscriminantScoreInstrumentList ();
public:
	static DiscriminantScoreInstrumentList& instance ();
	StringVector getNames () const { return names; }
};


#endif /* ! __lu_disc_sc_h */
