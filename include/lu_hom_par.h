/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_hom_par.h                                                  *
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
*  Copyright (2001-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_hom_par_h
#define __lu_hom_par_h

#include <lu_srch_par.h>
#include <lu_tol.h>

class MSHomologyParameters : public MSSearchParameters {

	int minMatches;
	std::string scoreMatrix;
	int maxPeptideHits;

	ToleranceInfo productMassTolerance;

	int previousAA;
	int nextAA;

	StringVector possibleSequences;
	IntVector sequenceSet;
	IntVector maxSeqErrors;
public:
	MSHomologyParameters ( const ParameterList* params );

	int getMinMatches () const { return minMatches; }
	std::string getScoreMatrix () const { return scoreMatrix; }
	int getMaxPeptideHits () const { return maxPeptideHits; }

	Tolerance* getProductMassTolerance () const { return productMassTolerance.getTolerance (); }

	int getPreviousAA () const { return previousAA; }
	int getNextAA () const { return nextAA; }

	StringVector getPossibleSequences () const { return possibleSequences; }
	IntVector getSequenceSet () const { return sequenceSet; }
	IntVector getMaxSeqErrors () const { return maxSeqErrors; }

	void printHTML ( std::ostream& os ) const;
};

#endif /* ! __lu_hom_par_h */
