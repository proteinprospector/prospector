/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_disc_sc.cpp                                                *
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
*  Copyright (2003-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <cmath>
#include <limits>
#include <lgen_error.h>
#include <lu_getfil.h>
#include <lu_disc_sc.h>
#include <lu_param_list.h>
using std::string;
using std::getline;

const double DiscriminantScore::MIN_DISC_SCORE = -std::numeric_limits<double>::max();

DiscriminantScore::DiscriminantScore ( const string& discScoreName ) :
	offsetCoeff		( MIN_DISC_SCORE ),
	bestScoreCoeff	( MIN_DISC_SCORE ),
	scoreDiffCoeff	( MIN_DISC_SCORE ),
	maxBestScore	( MIN_DISC_SCORE ),

	offsetCoeff2	( MIN_DISC_SCORE ),
	bestScoreCoeff2	( MIN_DISC_SCORE ),
	expectCoeff2	( MIN_DISC_SCORE ),
	maxBestScore2	( MIN_DISC_SCORE )
{
	GenIFStream fromFile ( MsparamsDir::instance ().getParamPath ( "disc_score.txt" ) );
	string line;
	bool flag = false;
	while ( getline ( fromFile, line ) ) {
		if ( line == discScoreName ) {
			ParameterList params ( fromFile );
			params.getValue ( "offset",		offsetCoeff );
			params.getValue ( "best_score", bestScoreCoeff );
			params.getValue ( "score_diff", scoreDiffCoeff );
			params.getValue ( "max_best_score", maxBestScore );
			flag = true;
			break;
		}
	}
	if ( !flag ) ErrorHandler::genError ()->error ( "Invalid or unspecified discriminant score name (file: disc_score.txt).\nFunction: DiscriminantScore.\n" );
	flag = false;
	GenIFStream fromFile2 ( MsparamsDir::instance ().getParamPath ( "disc_score2.txt" ) );
	while ( getline ( fromFile2, line ) ) {
		if ( line == discScoreName ) {
			ParameterList params ( fromFile2 );
			params.getValue ( "offset",		offsetCoeff2 );
			params.getValue ( "best_score", bestScoreCoeff2 );
			params.getValue ( "expectation", expectCoeff2 );
			params.getValue ( "max_best_score", maxBestScore2 );
			flag = true;
		}
	}
	if ( !flag ) ErrorHandler::genError ()->error ( "Invalid or unspecified discriminant score name (file: disc_score2.txt).\nFunction: DiscriminantScore.\n" );
}
double DiscriminantScore::calculateDiscriminantScore ( double bestScore, double scoreDiff, double expect ) const
{
	double discScore;
	if ( expect == -1.0 || expect == 0.0 ) {
		if ( maxBestScore != MIN_DISC_SCORE ) {
			bestScore = genMin ( bestScore, maxBestScore );
		}
		discScore = offsetCoeff;
		if ( bestScoreCoeff != MIN_DISC_SCORE ) discScore += bestScore * bestScoreCoeff;
		if ( scoreDiffCoeff != MIN_DISC_SCORE ) discScore += scoreDiff * scoreDiffCoeff;
	}
	else {
		if ( maxBestScore2 != MIN_DISC_SCORE ) {
			bestScore = genMin ( bestScore, maxBestScore2 );
		}
		discScore = offsetCoeff2;
		if ( bestScoreCoeff2 != MIN_DISC_SCORE ) discScore += bestScore * bestScoreCoeff2;
		if ( expectCoeff2 != MIN_DISC_SCORE )	discScore += -log10 ( expect ) * expectCoeff2;
	}
	return discScore;
}
DiscriminantScoreInstrumentList::DiscriminantScoreInstrumentList ()
{
	initialise ();
}
void DiscriminantScoreInstrumentList::initialise ()
{
	GenCommentedIFStream ifs ( MsparamsDir::instance ().getParamPath ( "disc_score.txt" ) );
	string line;
	while ( ifs.getUncommentedLine ( line ) ) {
		names.push_back ( line );
		ParameterList params ( ifs );
	}
}
DiscriminantScoreInstrumentList& DiscriminantScoreInstrumentList::instance ()
{
	static DiscriminantScoreInstrumentList discriminantScoreInstrumentList;
	return discriminantScoreInstrumentList;
}
