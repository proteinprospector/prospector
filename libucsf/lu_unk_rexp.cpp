/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_unk_rexp.cpp                                               *
*                                                                             *
*  Created    : April 21st 1997                                               *
*                                                                             *
*  Purpose    : Code to search the Unknome hit list in MS-Homology.           *
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
#include <lg_new.h>
#include <lg_string.h>
#include <lgen_error.h>
#include <lu_unk_rexp.h>
#include <lu_seq_exp.h>
using std::vector;
using std::copy;
using std::string;

static void initialise_unknome_search ();

static vector <char*> hits;
static vector <int> seqset;
static vector <int> maxerrs;
static int num_hits;
static int check_len;

static int* node;
static int maxnumerrs;
static int maxlen;
static int minlen;
static int* sum_errors;

static IntVector hit_number (1);
static int hit_number_array_size = 1;
static int num_errors;
static int nhits;

static const char* loc1;
static const char* loc2;
#include <iostream>
void initialise_regular_expression_hit_list ( const StringVector& reg_exp_hits, const IntVector& set, const IntVector& max_allow_errors, const Tolerance* massTolerance, int check_peptide_length )
{
	if ( reg_exp_hits.empty () ) {
		ErrorHandler::genError ()->error ( "No sequences entered.\n" );
	}
	for ( StringVectorSizeType i = 0 ; i < reg_exp_hits.size () ; i++ ) {
		SequenceExpression seq_exp ( reg_exp_hits [i], massTolerance );
		StringVector seq_list = seq_exp.getSequenceList ();
		for ( StringVectorSizeType j = 0 ; j < seq_list.size () ; j++ ) {
			hits.push_back ( gen_new_string ( seq_list [j].c_str () ) );
			seqset.push_back ( set [i] );
			maxerrs.push_back ( max_allow_errors [i] );
			string s = hits.back ();
			if ( s.length () < 4 ) {
				ErrorHandler::genError ()->error ( "One or more of the possible sequences is too short.\n" );
			}
			if ( s.length () <= maxerrs.back () + 3 ) {
				string err;
				if ( reg_exp_hits.size () == 1 )
					err = "The maximum number of errors is too large for the specified sequence.\n";
				else
					err = "The maximum number of errors is too large for the specified sequence for spectrum " + gen_itoa (i+1) + ".\n";
				ErrorHandler::genError ()->error ( err );
			}
		}
	}
	num_hits = hits.size ();
	check_len = check_peptide_length;
	initialise_unknome_search ();
}
static void initialise_unknome_search ()
{
	int i, j;

	node = new int [num_hits];
	node [0] = 0;
	minlen = maxlen = strlen ( hits [0] );
	for ( i = 1 ; i < num_hits ; i++ ) {
		for ( j = 0 ; hits [i] [j] == hits [i-1] [j] ; j++ );
		node [i] = j;
		minlen = genMin ( (int) strlen ( hits [i] ), minlen );
		maxlen = genMax ( (int) strlen ( hits [i] ), maxlen );
	}
	maxnumerrs = 0;
	for ( i = 0 ; i < num_hits ; i++ ) {
		maxnumerrs = genMax ( maxerrs [i], maxnumerrs );
	}
	sum_errors = new int [maxlen];
}
bool multi_unknome_hits_present ( char* offset )
{
	static int same_string = 0;
	const char* loc = check_len ? loc2 : loc1 + 1;

	if ( same_string ) {
		if ( unknome_hit_present ( loc ) ) {		/* Called after a hit has already been found in a reading frame. */
			return ( true );
		}
	}
	else {
		if ( unknome_hit_present ( offset ) ) {		/* Called at the start of a reading frame. */
			same_string = 1;
			return ( true );
		}
	}
	same_string = 0;								/* End of reading frame, no more hits. */
	return ( false );
}
bool unknome_hit_present ( const char* offset )
{
	int len = strlen ( offset );

	nhits = 0;
	if ( len < minlen ) return ( false );	/* Shorter than any hit. */
	if ( check_len ) len = 1;			/* Used for digest specific fragments */
	num_errors = maxnumerrs + 1;
	int bestLen = maxlen;
	for ( int i = 0 ; i < len ; i++ ) {
		loc1 = offset + i;				/* loc1 is the start of the section of protein under examination */
		sum_errors [0] = 0;
		for ( int j = 0, k = 0 ; j < num_hits ; j++ ) {
			int nod = node [j];
			if ( k >= nod ) {
				const char* ep = &(hits [j][nod]);	/* ep is the pointer in the hit string */
				const char* lp = loc1 + nod;		/* lp is the pointer in the protein string */
				int nerrs = sum_errors [nod];

				for ( k = nod ; ; k++ ) {
					if ( *ep == 0 ) {
						loc2 = lp;
						int curLen = loc2 - loc1;
						if ( nerrs < num_errors ) {	// If there are several hits for a particular protein alignment
							hit_number [0] = j;
							nhits = 1;				// only the one (or more) with the least errors is reported.
							num_errors = nerrs;
							bestLen = curLen;
						}
						else if ( nhits && ( nerrs == num_errors ) ) {
							if ( curLen < bestLen ) {	// only the one (or more) which is shortest is reported.
								hit_number [0] = j;
								nhits = 1;
								num_errors = nerrs;
								bestLen = curLen;
							}
							else if ( nhits && ( bestLen == curLen ) ) {
								if ( nhits == hit_number_array_size ) {
									hit_number_array_size *= 2;
									hit_number.resize ( hit_number_array_size );
								}
								hit_number [nhits++] = j;
							}
						}
						break;
					}
					if ( *lp == 0 ) break;
					sum_errors [k] = nerrs;
					if ( *ep++ != *lp++ ) {
						if ( ++nerrs > maxerrs [j] ) break;
					}
				}
			}
		}
		if ( num_errors <= maxnumerrs ) return ( true );
	}
	return ( false );
}
bool get_unknome_rexp_results ( char** measured_sequence, int* sequence_set, int* num_substitutions, bool reset )
{
	static int i;

	if ( reset ) i = 0;

	for ( ; i < nhits ; ) {
		*measured_sequence = hits [hit_number [i]];
		*num_substitutions = num_errors;
		*sequence_set = seqset [hit_number [i]];
		i++;
		return ( true );
	}
	return ( false );
}
char* get_regular_expression_match ()
{
	static vector <char> expression;
	int length = loc2 - loc1;
	expression.resize ( length + 1 );
	copy ( loc1, loc2, expression.begin () );
	expression [length] = 0;
	return ( &expression [0] );
}
const char* get_regular_expression_loc1 ()
{
	return ( loc1 );
}
const char* get_regular_expression_loc2 ()
{
	return ( loc2 );
}
