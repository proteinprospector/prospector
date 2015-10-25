/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_seq_exp.h                                                  *
*                                                                             *
*  Created    : October 19th 1999                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1999-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_seq_exp_h
#define __lu_seq_exp_h

#include <string>
#include <vector>

class Tolerance;

struct TreeList {
	int num_permutations;
	std::string* permutation;
};

class SequenceExpression {
	const Tolerance* massTolerance;
	std::string originalExp;
	char* baseExp;
	char** expList;
	int expListLength;
	TreeList* treeList;
	char* sequence;
	std::vector <std::string> sequenceList;
	void compile ( const std::string& expression );
	char* nullTerminate ( char* str, int open, int close );
	char* getTreeExpression ( const char* expression );
	void get_next_sequence ( TreeList tL, std::string rest, int index );
	std::string getMassCompositions ( double mass ) const;
public:
	SequenceExpression ( const std::string& seqExpression, const Tolerance* massTolerance );
	~SequenceExpression ();
	std::vector <std::string> getSequenceList () {return sequenceList;};
};

#endif /* ! __lu_seq_exp_h */
