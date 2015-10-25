/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lgen_math.h                                                   *
*                                                                             *
*  Created    : January 3rd 1997                                              *
*                                                                             *
*  Purpose    : Maths utilities.                                              *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1997-2008) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lgen_math_h
#define __lgen_math_h

#include <lgen_define.h>

double gen_binomial_probability ( int num_independent_trials, double prob_of_success, int num_successes );
double gen_multinomial_coeff ( int num_independent_trials, int* num_successes, int num_classes );
double gen_multinomial_probability_product ( int num_independent_trials, double* prob_of_success, int* num_successes, int num_classes );
double gen_multinomial_probability ( int num_independent_trials, double* prob_of_success, int* num_successes, int num_classes );
double gen_poisson_probability ( int num_independent_trials, double probability_of_success, int num_successes );
double gen_array_sum ( double* array, int array_len );
int gen_integer_series_sum ( int number );
int gen_round_up_int_divide ( int top, int bottom );
int genNearestIndex ( const DoubleVector& dv, double num );
#ifdef VIS_C
#define isnan(x) ((x) != (x))
#endif
inline bool genIsNumberStart ( int c )
{
	return isdigit ( c ) || c == '-' || c == '+';
}

#endif /* ! __lgen_math_h */
