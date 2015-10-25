/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_pi.cpp                                                     *
*                                                                             *
*  Created    : April 30th 1998                                               *
*                                                                             *
*  Purpose    : Functions to calculate the pI of a protein.                   *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1998-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <nr.h>
#include <lg_io.h>
#include <lg_memory.h>
#include <lu_pi.h>
#include <lu_mass.h>
using std::string;
using std::ostream;

static double ph_to_charge ( double ph );

static int multiplier_list [AA_ARRAY_SIZE];
static double pk_n_terminus_aa;
static double pk_c_terminus_aa;
static double ln10;

int ProteinPI::piPrecision = 1;
ProteinPI::ProteinPI ( const string& proteinString )
{
	static int array_size = AA_ARRAY_SIZE * sizeof (int);
	int len = proteinString.length ();
	
	if ( len == 0 ) {
		pi = PI_NOT_CALCULATED;
	}
	else {
		pk_n_terminus_aa = pk_n_term [proteinString [0]];
		pk_c_terminus_aa = pk_c_term [proteinString [len-1]];
		if ( pk_n_terminus_aa == NO_PK_VALUE || pk_c_terminus_aa == NO_PK_VALUE ) {
			pi = PI_NOT_CALCULATED;
		}
		else {
			ln10 = log ( 10.0 );
			gen_bzero ( (char*)multiplier_list, array_size );

			for ( int i = 0 ; i < len ; i++ ) {
				multiplier_list [proteinString[i]]++;
			}
			pi = zbrent ( ph_to_charge, 1.0, 14.0, 0.01 );

			pi = nrerrno == 0 ? pi : PI_NOT_CALCULATED;
		}
	}
}
ostream& operator<< ( ostream& os, const ProteinPI& ppi )
{
	double pi = ppi.getProteinPI ();
	if ( pi == PI_NOT_CALCULATED ) os << "----";
	else {
		genPrint ( os, pi, ProteinPI::piPrecision );
	}
	return os;
}
static double ph_to_charge ( double ph )
{
	double charge = 0.0;
	double alpha;
	int i;

	for ( i = 'A' ; i <= 'z' ; i++ ) {
		if ( multiplier_list [i] ) {
			if ( pk_basic_sc [i] != NO_PK_VALUE ) {
				alpha = exp ( ln10 * ( ph - pk_basic_sc [i] ) );
				charge += multiplier_list [i] / ( alpha + 1.0 );
			}
			else if ( pk_acidic_sc [i] != NO_PK_VALUE ) {
				alpha = exp ( ln10 * ( ph - pk_acidic_sc [i] ) );
				charge -= multiplier_list [i] * alpha / ( 1.0 + alpha );
			}
		}
	}
	alpha = exp ( ln10 * ( ph - pk_n_terminus_aa ) );
	charge += 1.0 / ( alpha + 1.0 );								/* N terminus aa */

	alpha = exp ( ln10 * ( ph - pk_c_terminus_aa ) );
	charge -= alpha / ( 1.0 + alpha );								/* C terminus aa */

	return ( charge );
}
