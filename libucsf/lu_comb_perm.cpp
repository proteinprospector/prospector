/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_comb_perm.cpp                                              *
*                                                                             *
*  Created    : January 3rd 1997                                              *
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
*  Copyright (1997-2009) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifndef VIS_C
#include <stdexcept>
#endif
#include <nr.h>
#include <lg_new.h>
#include <lg_memory.h>
#include <lg_stdlib.h>
#include <lgen_error.h>
#include <lu_mass.h>
#include <lu_comb_perm.h>
#include <lu_param_list.h>
using std::ostream;
using std::string;
using std::sort;
using std::copy;
using std::set_difference;
using std::endl;
using std::runtime_error;

static const int MAX_AA_PERMUTATIONS = 500000000;
static const int MAX_NUM_UNIQUE_AA = 64;
static double wt [MAX_NUM_UNIQUE_AA];
static double wt_diff [MAX_NUM_UNIQUE_AA];
static char* listAA;
static char array [200];
static Combination_hit* hit;
static double start_mass;
static double end_mass;
static int len;
static int num;
static int max_hits;
static int max_hits_exceeded;
static int amount [MAX_NUM_UNIQUE_AA];
static int len_var_combination;
static int len_combination;
static int num_uniq_aa;
static char uniq_aa [MAX_NUM_UNIQUE_AA];
static char* permutation;
static char* permutation_ptr;
static void ( *add_hit ) ( int, double );

static string combinationMinusCombination ( const string& comb1, const string& comb2 );
static void get_combinations ( int level, int start_letter, double sum );
static void add_hit_super_element ( int level, double sum );
static void add_hit_elemental ( int level, double sum );
static void get_permutations ( int level );
static int sort_list_by_mass ( const void* i, const void* j );

string init_comb_search_aa_list ( const string& possibleAminoAcids, const string& aaExclude, const string& aaAdd )
{
	return ( combinationPlusCombination ( combinationMinusCombination ( possibleAminoAcids, aaExclude ), aaAdd ) );
}
int calculate_combinations ( const string& aa, double mass, double tolerance, int maxReportedHits, bool elemental_flag, Combination_hit** combinations, int* num_combinations )
{
	int i;

	listAA = gen_new_string ( aa.c_str () );
	len = strlen ( listAA );
	gen_qsort ( (char*)listAA, len, sizeof (char), sort_list_by_mass );

	hit = new Combination_hit [maxReportedHits + 2];
	max_hits = maxReportedHits;
	max_hits_exceeded = 0;

	for ( i = 0 ; i < len ; i++ ) {
		wt [i] = amino_acid_wt [listAA[i]];
		if ( wt [i] < amino_acid_wt ['G'] - 1.0 ) {
			throw runtime_error ( "The mass of amino acid " + string ( 1, listAA[i] ) + " is too low to be included in the search." );
		}
	}
	for ( i = 0 ; i < len ; i++ ) {
		wt_diff [i] = ( i == len - 1 ) ? 0.0 : wt[i+1] - wt[i];
	}
	num = 0;
	start_mass = mass - tolerance;
	end_mass = mass + tolerance;

	if ( elemental_flag ) add_hit = add_hit_elemental;
	else add_hit = add_hit_super_element;

	if ( terminal_wt >= start_mass && terminal_wt <= end_mass )
		add_hit ( 0, terminal_wt );
	else
		get_combinations ( 0, 0, terminal_wt );

	*combinations = hit;
	*num_combinations = num;

	return ( max_hits_exceeded );
}
static void get_combinations ( int level, int i, double sum )
{
	sum += wt [i];
	while ( i < len ) {
		if ( sum > end_mass ) break;
		array [level] = listAA [i];
		if ( sum < start_mass ) get_combinations ( level + 1, i, sum );
		else {
			if ( num < max_hits ) add_hit ( level, sum );
			else max_hits_exceeded = 1;
		}
		sum += wt_diff [i++];
	}
}
static void add_hit_super_element ( int level, double sum )
{
	array [level+1] = 0;
	hit [num].sum = sum;
	hit [num].array = gen_new_string ( array );
	num++;
}
static void add_hit_elemental ( int level, double sum )
{
	int i;

	for ( i = 0 ; i < num ; i++ ) {
		if ( genAbsDiff ( hit [i].sum, sum ) < FORMULA_DIFF_TOL ) return;
	}
	add_hit_super_element ( level, sum );	
}
double calc_num_permutations_from_combination ( const string& combination )
{
	int len = combination.length ();
	double num_perms;
	int i, j;

	num_perms = factrl ( len );

	for ( i = 1, j = 1 ; i < len ; i++ ) {
		if ( combination [i] != combination [i-1] ) {
			if ( j > 1 ) {
				num_perms /= factrl ( j );
				j = 1;
			}	
		}
		else j++;
	}
	num_perms /= factrl ( j );
	return ( num_perms );
}
string combinationPlusCombination ( const string& comb1, const string& comb2 )
{
	string newComb = comb1;
	newComb += comb2;
	sort ( newComb.begin (), newComb.end () );
	return ( newComb );
}
string combinationMinusCombination ( const string& comb1, const string& comb2 )
{
	string newComb;

	for ( StringSizeType i = 0 ; i < comb1.length () ; i++ ) {
		for ( StringSizeType j = 0 ; j < comb2.length () ; j++ ) {
			if ( comb2 [j] == comb1 [i] ) goto label;
		}
		newComb += comb1 [i];
		label:;
	}
	return ( newComb );
}
double subtract_known_composition ( double parentMass, const StringVector& comp_ions )
{
	for ( StringVectorSizeType i = 0 ; i < comp_ions.size () ; i++ ) {
		if ( comp_ions [i].length () == 1 ) {
			parentMass -= amino_acid_wt [comp_ions[i][0]];
		}
	}
	return ( parentMass );
}
void permutations_from_combination ( char** perm_sequences, int* num_perm_sequences, const string& combination, char n_term, char c_term )
{
	char n_term_seq [2];
	char c_term_seq [2];
	int len_n_term = n_term ? 1 : 0;
	int len_c_term = c_term ? 1 : 0;

	len_combination = combination.length ();

	string tempCombination = combination;

	if ( n_term ) {
		n_term_seq [0] = n_term;
		n_term_seq [1] = 0;
		tempCombination = combinationMinusCombination ( tempCombination, n_term_seq );
	}
	if ( c_term ) {
		c_term_seq [0] = c_term;
		c_term_seq [1] = 0;
		tempCombination = combinationMinusCombination ( tempCombination, c_term_seq );
	}
	len_var_combination = tempCombination.length ();

	double n_perms = calc_num_permutations_from_combination ( tempCombination );
	
	if ( n_perms > MAX_AA_PERMUTATIONS ) {
		sprintf ( gen_error_message, "Your data generates too many different amino acid permutations to proceed. The current maximum is %d.", MAX_AA_PERMUTATIONS );
		ErrorHandler::genError ()->error ( gen_error_message );
	}
	int num_permutations = (int)n_perms;

	gen_bzero ( (char*)amount, MAX_NUM_UNIQUE_AA * sizeof (int) );

	if ( c_term ) array [len_n_term + len_var_combination] = c_term;
	if ( n_term ) array [0] = n_term;

	array [len_combination] = 0;

	permutation = new char [num_permutations * ( len_combination + 1 )];
	permutation_ptr = permutation;

	num_uniq_aa = 0;
	for ( int i = 0 ; i < len_var_combination ; i++ ) {
		if ( i == 0 || tempCombination [i] != tempCombination [i-1] ) {
			uniq_aa [num_uniq_aa++] = tempCombination [i];
		}
		amount [num_uniq_aa-1]++;
	}
	get_permutations ( len_n_term );

	*perm_sequences = permutation;
	*num_perm_sequences = num_permutations;
}
double get_sequence_probability ( char* sequence, int err )
{
	int len = strlen ( sequence );
	static double aa_prob [] = {
		0.073686,		// A
		0.000000,		// B
		0.017806,       // C
		0.051768,		// D
		0.062030,		// E
		0.040020,		// F
		0.069358,		// G
		0.022315,		// H
		0.056327,		// I
		0.000000,		// J
		0.057993,		// K
		0.092001,		// L
		0.023345,		// M
		0.046234,		// N
		0.000000,		// O
		0.050189,		// P
		0.040621,		// Q
		0.051965,		// R
		0.073744,		// S
		0.059710,		// T
		0.000000,		// U
		0.064316,		// V
		0.013343,		// W
		0.000585,		// X
		0.032643,		// Y
		0.000000 };		// Z
	double seq_prob = 1.0;
	for ( int i = 0 ; i < len ; i++ ) {
		seq_prob *= aa_prob [sequence [i] - 'A'];
	}
	seq_prob = 1 / seq_prob;
	double effective_num_aa = exp ( log ( seq_prob ) / len );
	double bc = bico ( len, err );
	double err_prob = pow ( effective_num_aa - 1, err );
	double prob = seq_prob / ( bc * err_prob );

	return ( prob );
}
static void get_permutations ( int level )
{
	level++;
	for ( int i = 0 ; i < num_uniq_aa ; i++ ) {
		if ( amount [i] ) {
			amount [i] -= 1;
			array [level-1] = uniq_aa [i];
			if ( level < len_var_combination ) get_permutations ( level );
			else {
				strcpy ( permutation_ptr, array );
				permutation_ptr += len_combination + 1;
			}
			amount [i] += 1;
		}
	}
}
void free_permutations ()
{
	delete [] permutation;
}
static int sort_list_by_mass ( const void* i, const void* j )
{
	double ret = amino_acid_wt [*(char*)i] - amino_acid_wt [*(char*)j];

	if ( ret > 0.0 ) return ( 1 );
	if ( ret < 0.0 ) return ( -1 );
	return ( 0 );
}
/* Function not used - but might be useful in the future */
/*
static void calc_num_combinations ( int start_aa, int level )
{
	for ( int i = start_aa ; i < aa_list_length ; i++ ) {
		level += 1;
		if ( level < num_levels ) calc_num_combinations ( i, level );
		if ( level == num_levels ) {
			num_hits++;
		}
		level -= 1;
	}
}
*/
#ifdef FINISH
static void get_next_sequence ( char* array, int index )
{
	int i, j;

	for ( i = 0, j = 0 ; i < length ; i++ ) {
		switch ( array [i] ) {
			case '{':
				permutations_from_combination ( perm_sequences, num_perm_sequences, combination, 0, 0 );
				for ( j = 0 ; j < num_perm_sequences ; j++ ) {
					get_next_sequence ( perm_sequences [j], int level );
				}
				break;
			default:
				sequence [index++] = array [i];
				break;
		}
		if ( index == length ) {
			sequence [index] = 0;
			perms [num_sequences] = gen_new_string ( sequence );
			num_sequences++;
		}
	}
}
#endif
CombinationHit::CombinationHit ( char* sequence, double mass )
{
	this->sequence = gen_new_string ( sequence );
	this->mass = mass;
	this->sequenceLength = strlen ( sequence );
}
CombinationHit::~CombinationHit ()
{
	delete [] sequence;
}
PossibleCombinations::PossibleCombinations ( const string& aaList, int maxLength, double maxMass, bool permutationFlag ) :
	aaList ( aaList ),
	aaListLength ( aaList.length () ),
	maxLength ( maxLength ),
	maxMass ( maxMass ),
	permutationFlag ( permutationFlag )
{
	sequence = new char [maxLength + 1];

	for ( int i = 1 ; i <= maxLength ; i++ ) {
		numLevels = i;
		getNext ( 0, 0, 0.0 );
	}
	sort ( hits.begin (), hits.end (), SortPossibleCombinationsAscending () );
}
PossibleCombinations::~PossibleCombinations ()
{
	delete [] sequence;
	for ( CombinationHitIterator i = hits.begin () ; i != hits.end () ; i++ ) {
		delete *i;
	}
}
void PossibleCombinations::getNext ( int startAA, int level, double mass )
{
	int start = permutationFlag ? 0 : startAA;
	for ( int i = start ; i < aaListLength ; i++ ) {
		sequence [level] = aaList [i];
		sequence [level+1] = 0;
		mass += amino_acid_wt [aaList[i]];
		level += 1;
		if ( level < numLevels ) getNext ( i, level, mass );
		if ( level == numLevels && mass <= maxMass ) {
			hits.push_back ( new CombinationHit ( sequence, mass ) );
		}
		sequence [level] = 0;
		mass -= amino_acid_wt [aaList[i]];
		level -= 1;
	}
}
ConsideredAA::ConsideredAA ( const ParameterList* params )
{
	possibleAminoAcids = "ACDEFGHIKLMNPQRSTUVWY";

	CharVector cvvalue;
	if ( params->getValue ( "aa_exclude", cvvalue ) ) {
		aaExclude.resize ( cvvalue.size () );
		copy ( cvvalue.begin (), cvvalue.end (), aaExclude.begin () );
	}
	sort ( aaExclude.begin (), aaExclude.end () );

	if ( params->getValue ( "aa_add", cvvalue ) ) {
		aaAdd.resize ( cvvalue.size () );
		copy ( cvvalue.begin (), cvvalue.end (), aaAdd.begin () );
	}
	sort ( aaAdd.begin (), aaAdd.end () );

	aaList.resize ( possibleAminoAcids.length () ); 
	set_difference ( possibleAminoAcids.begin (), possibleAminoAcids.end (), aaExclude.begin (), aaExclude.end (), aaList.begin () );
	aaList.resize ( strlen(aaList.c_str ()) ); 
	aaList += aaAdd;
	sort ( aaList.begin (), aaList.end () );
}
void ConsideredAA::printHTML ( ostream& os ) const
{
	if ( aaExclude.size () ) {
		os << "Absent Amino Acids: <b>";
		for ( StringSizeType i = 0 ; i < aaExclude.size () ; i++ ) os << aaExclude [i];
		os << "</b><br />" << endl;
	}
	if ( aaAdd.size () ) {
		os << "Modified Amino Acids Considered: <b>";
		for ( StringSizeType i = 0 ; i < aaAdd.size () ; i++ ) os << aaAdd [i];
		os << "</b><br />";
	}
}
