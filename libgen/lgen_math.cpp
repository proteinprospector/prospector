/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lgen_math.cpp                                                 *
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
#include <nr.h>
using std::lower_bound;

double gen_binomial_probability ( int num_independent_trials, double prob_of_success, int num_successes )
{
	return ( bico ( num_independent_trials, num_successes ) * pow ( prob_of_success, (double)num_successes ) * pow ( ( 1.0 - prob_of_success ), (double)( num_independent_trials - num_successes ) ) );
}
double gen_multinomial_coeff ( int num_independent_trials, int* num_successes, int num_classes )
{
	int sum_of_n;
	double sum_of_n_factorial;
	double product_of_n_factorial;
	int i;

	for ( i = 0, sum_of_n = 0 ; i < num_classes ; i++ ) {
		sum_of_n += num_successes [i];
	}
	sum_of_n_factorial = factrl ( sum_of_n );
	for ( i = 0, product_of_n_factorial = 1.0 ; i < num_classes ; i++ ) {
		if ( num_successes [i] > 1 ) product_of_n_factorial *= factrl ( num_successes [i] );
	}
	return ( bico ( num_independent_trials, sum_of_n ) * sum_of_n_factorial / product_of_n_factorial );
}
double gen_multinomial_probability_product ( int num_independent_trials, double* prob_of_success, int* num_successes, int num_classes )
{
	double sum_of_probabilities;
	double product_of_probabilities;
	int sum_of_n;
	int i;

	for ( i = 0, sum_of_n = 0 ; i < num_classes ; i++ ) {
		sum_of_n += num_successes [i];
	}
	if ( sum_of_n > num_independent_trials ) return ( 0.0 );
	for ( i = 0, sum_of_probabilities = 0.0 ; i < num_classes ; i++ ) {
		sum_of_probabilities += prob_of_success [i];
	}
	for ( i = 0, product_of_probabilities = 1.0 ; i < num_classes ; i++ ) {
		if ( num_successes [i] > 0 ) product_of_probabilities *= pow ( prob_of_success [i], (double) num_successes [i] );
	}
	product_of_probabilities *= pow ( ( 1.0 - sum_of_probabilities ), (double)( num_independent_trials - sum_of_n ) );

	return ( product_of_probabilities );
}
double gen_multinomial_probability ( int num_independent_trials, double* prob_of_success, int* num_successes, int num_classes )
{
	double multinomial_coefficient;
	double product_of_probabilities = gen_multinomial_probability_product ( num_independent_trials, prob_of_success, num_successes, num_classes );

	if ( product_of_probabilities == 0.0 ) return ( 0.0 );

	multinomial_coefficient = gen_multinomial_coeff ( num_independent_trials, num_successes, num_classes );

	return ( multinomial_coefficient * product_of_probabilities );
}
double gen_poisson_probability ( int num_independent_trials, double probability_of_success, int num_successes )
{
	double mean = num_independent_trials * probability_of_success;

	return ( pow ( mean, num_successes ) * exp ( - mean ) / factrl ( num_successes ) );
}
double gen_array_sum ( double* array, int array_len )
{
	double sum = 0.0;

	if ( array != NULL ) for ( int i = 0 ; i < array_len ; i++ ) sum += array [i];

	return ( sum );
}
int gen_integer_series_sum ( int number )
{
	int sum = 0;

	for ( int i = 1 ; i <= number ; i++ ) sum += i;

	return ( sum );
}
int gen_round_up_int_divide ( int top, int bottom )
{
	return ( ( top / bottom ) + ( ( top % bottom ) ? 1 : 0 ) );
}
int genNearestIndex ( const DoubleVector& dv, double num )	// Finds the index of the element in dv nearest to num
{
	int ind = lower_bound ( dv.begin (), dv.end (), num ) - dv.begin ();
	if ( ind == 0 || num == dv [ind] ) return ind;
	else {
		double diff1 = genAbsDiff ( num, dv [ind] );
		double diff2 = genAbsDiff ( num, dv [ind-1] );
		return diff2 < diff1 ? ind-1 : ind;
	}
}
