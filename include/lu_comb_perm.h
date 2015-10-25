/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_comb_perm.h                                                *
*                                                                             *
*  Created    : January 6th 1997                                              *
*                                                                             *
*  Purpose    : Functions concerned with combinations and permutations of     *
*               amino acids.                                                  *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1997-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_comb_perm_h
#define __lu_comb_perm_h

#include <lgen_sort.h>

class ParameterList;

struct Combination_hit {
	char* array;
	double sum;
};
const double FORMULA_DIFF_TOL = 0.000000001;

std::string init_comb_search_aa_list ( const std::string& possibleAminoAcids, const std::string& aaExclude, const std::string& aaAdd );
int calculate_combinations ( const std::string& aa, double mass, double tolerance, int maxReportedHits, bool elemental_flag, Combination_hit** combinations, int* num_combinations );
double calc_num_permutations_from_combination ( const std::string& combination );
std::string combinationPlusCombination ( const std::string& comb1, const std::string& comb2 );
double subtract_known_composition ( double parentMass, const StringVector& comp_ions );
void permutations_from_combination ( char** perm_sequences, int* num_perm_sequences, const std::string& combination, char n_term, char c_term );
double get_sequence_probability ( char* sequence, int err );
void free_permutations ();

class CombinationHit {
	double mass;
	char* sequence;
	int sequenceLength;
	friend class SortPossibleCombinationsAscending;
public:
	CombinationHit ( char* sequence, double mass );
	~CombinationHit ();
	char* getSequence () const { return sequence; }
	int getSequenceLength () const { return sequenceLength; }
	double getMass () const { return mass; }
};

class SortPossibleCombinationsAscending {
	public:
		int operator () ( const CombinationHit* a, const CombinationHit* b ) const
		{
			return a->mass < b->mass;
		}
};

typedef std::vector <CombinationHit*> CombinationHitContainer;
typedef CombinationHitContainer::size_type CombinationHitContainerSizeType;
typedef CombinationHitContainer::iterator CombinationHitIterator;

class PossibleCombinations {
	std::string aaList;
	int aaListLength;
	int maxLength;
	double maxMass;
	bool permutationFlag;
	char* sequence;
	int numLevels;
	CombinationHitContainer hits;
	void getNext ( int startAA, int level, double mass );
public:
	PossibleCombinations ( const std::string& aaList, int maxLength, double maxMass, bool permutationFlag = false );
	~PossibleCombinations ();
	CombinationHitContainer getHits () const { return hits; }
	double getMinMass () const { return hits.front ()->getMass (); }
	double getMaxMass () const { return hits.back ()->getMass (); }
};

class ConsideredAA {
	std::string possibleAminoAcids;
	std::string aaExclude;
	std::string aaAdd;
	std::string aaList;
public:
	ConsideredAA ( const ParameterList* params );
	std::string getAAList () const { return aaList; }
	void printHTML ( std::ostream& os ) const;
};

#endif /* ! __lu_comb_perm_h */
