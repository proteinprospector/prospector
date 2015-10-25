/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_unk_rexp.h                                                 *
*                                                                             *
*  Created    : April 22nd 1997                                               *
*                                                                             *
*  Purpose    : Code to search the Unknome hit list in MS-Pattern.            *
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

#ifndef __lu_unk_rexp_h
#define __lu_unk_rexp_h

#include <lgen_define.h>

class Tolerance;

void initialise_regular_expression_hit_list ( const StringVector& reg_exp_hits, const IntVector& set, const IntVector& max_allow_errors, const Tolerance* massTolerance, int check_peptide_length );
bool multi_unknome_hits_present ( char* offset );
bool unknome_hit_present ( const char* offset );
bool get_unknome_rexp_results ( char** measured_sequence, int* sequence_set, int* num_substitutions, bool reset );
char* get_regular_expression_match ();
const char* get_regular_expression_loc1 ();
const char* get_regular_expression_loc2 ();

#endif /* ! __lu_unk_rexp_h */
