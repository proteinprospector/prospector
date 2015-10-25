/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_seq_exp.cpp                                                *
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
*  Copyright (1997-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <cstdio>
#include <nr.h>
#include <lg_new.h>
#include <lg_string.h>
#include <lgen_error.h>
#include <lu_mass.h>
#include <lu_comb_perm.h>
#include <lu_seq_exp.h>
#include <lu_tol.h>
using std::string;
using std::sort;
using std::next_permutation;

static const int max_combinations = 50;
static const int max_bracket = 65;

SequenceExpression::SequenceExpression ( const string& seqExpression, const Tolerance* massTolerance ) :
	massTolerance ( massTolerance )
{
	compile ( seqExpression );
}
SequenceExpression::~SequenceExpression ()
{
}
void SequenceExpression::compile ( const string& expression )
{
	int i, j, k;
	originalExp = "";
	for ( const char* ptr = expression.c_str () ; *ptr ; ) {
		if ( isdigit ( *ptr ) ) {
			char* next;
			originalExp += getMassCompositions ( strtod ( ptr, &next ) );
			ptr += next - ptr;
		}
		else originalExp += *ptr++;
	}
	string::size_type ind = originalExp.find_first_of ( "." );	// This is illegal if not in a number	
	if ( ind != string::npos ) {
		ErrorHandler::genError ()->error (
			"Possible sequences field contains the illegal character: '" + string ( 1, originalExp[ind] ) + "'" );
	}
	char* baseExp = gen_new_string ( originalExp.c_str () );
	int baseExpLength = originalExp.length ();

	expList = new char* [gen_strcharcount ( baseExp, '{' ) + gen_strcharcount ( baseExp, '[' )];
	for ( i = 0, j = 0 ; i < baseExpLength ; i++ ) {
		if ( baseExp [i] == '{' ) {
			expList [j] = gen_new_string ( nullTerminate ( baseExp + i, '{', '}' ) );
		}
		else if ( baseExp [i] == '[' ) {
			expList [j] = gen_new_string ( nullTerminate ( baseExp + i, '[', ']' ) );
		}
		else continue;
		bool present = false;
		for ( k = 0 ; k < j ; k++ ) {
			if ( !strcmp ( expList [j], expList [k] ) ) present = true;
		}
		if ( present == false ) j++;
	}
	expListLength = j;

	treeList = new TreeList [expListLength+1];
	treeList [0].permutation = new string [1];
	treeList [0].num_permutations = 1;
	treeList [0].permutation [0] = getTreeExpression ( originalExp.c_str () );
	for ( i = 0 ; i < expListLength ; i++ ) {
		if ( expList [i][0] == '[' ) {
			char* temp_exprn = getTreeExpression ( expList [i] + 1 ); 
			treeList [i+1].permutation = new string [gen_strcharcount ( temp_exprn, '|' ) + 1];
			treeList [i+1].permutation [0] = strtok ( temp_exprn, "|]" );
			treeList [i+1].num_permutations = 1;
			for ( j = 1 ; ; j++ ) {
				char* exprn = strtok ( NULL, "|]" );
				if ( exprn == NULL ) break;
				treeList [i+1].permutation [j] = exprn;
				treeList [i+1].num_permutations++;
			}
		}
		if ( expList [i][0] == '{' ) {
			char* exprn = &expList [i][1];
			int len = strlen ( exprn );
			sort ( exprn, exprn + len );
			treeList [i+1].permutation = new string [(int)factrl (len)];
			j = 0;
			do {
				treeList [i+1].permutation [j++] = exprn;
				treeList [i+1].num_permutations = j;
			} while ( next_permutation ( exprn, exprn + len ) );
		}
	}
	for ( i = 0 ; i < expListLength ; i++ ) {
		delete [] expList [i];
	}
	delete [] expList;
	int index = 0;

	sequence = new char [baseExpLength+1];
	get_next_sequence ( treeList [0], "", index );
	delete [] sequence;

	for ( i = 0 ; i < expListLength + 1 ; i++ ) {
		delete [] treeList [i].permutation;
	}
	delete [] treeList;
}
char* SequenceExpression::nullTerminate ( char* str, int open, int close )
{
	for ( int i = 0, j = 0 ; ; i++ ) {
		if ( str [i] == open ) j++;
		if ( str [i] == close ) j--;
		if ( j == 0 ) {
			str [i] = 0;
			return ( str );
		}
	}
}
char* SequenceExpression::getTreeExpression ( const char* expression )
{
	int expressionLength = strlen ( expression );
	int i, j, k;
	char* treeExpression = new char [expressionLength+1];

	for ( i = 0, j = 0 ; i < expressionLength ; i++ ) {
		char code = expression [i];
		switch ( code ) {
			case '[':
			case '{':
				for ( k = 0 ; k < expListLength ; k++ ) {
					char* exprn = expList [k];
					int len = strlen ( exprn );
					if ( !strncmp ( expression + i, exprn, len ) ) {
						if ( k + 1 == max_bracket ) {
							ErrorHandler::genError ()->error ( "The sequence " + originalExp + " is too complex." );
						}
						treeExpression [j++] = k + 1;
						i += len;
						break;
					}
				}
				break;
			default:
				treeExpression [j++] = code;
				break;
		}
	}
	treeExpression [j] = 0;
 
	return ( treeExpression );
}
void SequenceExpression::get_next_sequence ( TreeList tL, string rest, int index )
{
	for ( int i = 0 ; i < tL.num_permutations ; i++ ) {
		string temp_sequence = tL.permutation [i] + rest;
		int len = temp_sequence.length ();
		int j, k;
		for ( j = 0, k = index ; j < len ; j++ ) {
			char code = temp_sequence [j];
			if ( code < max_bracket ) {
				get_next_sequence ( treeList [code], temp_sequence.substr ( j + 1, len - 1), k );
				goto lab;
			}
			else {
				sequence [k++] = code;
			}
		}
		sequence [k] = 0;
		sequenceList.push_back ( sequence );
		lab:;
	}
}
string SequenceExpression::getMassCompositions ( double mass ) const
{
	Combination_hit* combinations;
	int numCombinations;
	if ( calculate_combinations ( "ACDEFGHIKLMNPQRSTUVWY", mass + terminal_wt, massTolerance->getTolerance ( mass + terminal_wt), max_combinations, false, &combinations, &numCombinations ) ) {
		sprintf ( gen_error_message, "The mass %f leads to more than %d amino acid combinations.", mass, max_combinations );
		ErrorHandler::genError ()->error ( gen_error_message );
	}
	if ( numCombinations == 0 ) {
		sprintf ( gen_error_message, "The mass %.4f does not correspond to any amino acid combinations.", mass );
		ErrorHandler::genError ()->error ( gen_error_message );
	}
	string compositions = "";
	for ( int i = 0 ; i < numCombinations ; i++ ) {
		compositions = compositions + '{' + combinations [i].array + '}';
		if ( i != numCombinations - 1 ) compositions += '|';
	}
	return ( compositions );
}
