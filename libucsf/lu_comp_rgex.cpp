/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_comp_rgex.cpp                                              *
*                                                                             *
*  Created    : July 21st 1996                                                *
*                                                                             *
*  Purpose    : Regular expression style composition searches.                *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1996-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lg_memory.h>
#include <lu_mass.h>

static int peptide_include [AA_ARRAY_SIZE];
static char* peptide_include_alternatives [AA_ARRAY_SIZE];
static int peptide_exclude [AA_ARRAY_SIZE];

int init_composition_search ( char* include, char* exclude )
{
	int i, j, k;
	int count;
	char* next;

	gen_bzero ( (char*)peptide_include, AA_ARRAY_SIZE * sizeof (int) );
	gen_bzero ( (char*)peptide_include_alternatives, AA_ARRAY_SIZE * sizeof (int) );
	while ( *include ) {
		while ( isspace ( *include ) ) include++;
		if ( *include == '[' ) {
			include++;
			for ( count = 0 ; include [count] != ']' ; count++ ) {
				if ( isupper ( include [count] ) == 0 ) return ( -1 );
			}
			for ( i = 0 ; i < count ; i++ ) {
				if ( isdigit ( include [count+1] ) ) {
					peptide_include [include[i]] += strtol ( include + count + 1, &next, 10 );
				}
				else {
					peptide_include [include[i]] += 1;
				}
				peptide_include_alternatives [include[i]] = new char [count + 1];
				for ( j = 0, k = 0 ; j < count ; j++ ) {
					if ( include [i] != include [j] ) {
						peptide_include_alternatives [include[i]][k++] = include [j];
					}
				}
				peptide_include_alternatives [include[i]][k] = 0;
			}
			include = ( isdigit ( include [count+1] ) ) ? next : include + count + 1;
		}
		else {
			if ( isupper ( *include ) ) {
				if ( isdigit ( include [1] ) ) {
					peptide_include [*include] += strtol ( ++include, &include, 10 );
				}
				else {
					peptide_include [*include++] += 1;
				}
			}
			else {
				return ( -1 );
			}
		}
		while ( isspace ( *include ) ) include++;
	}

	gen_bzero ( (char*)peptide_exclude, AA_ARRAY_SIZE * sizeof (int) );
	while ( *exclude ) {
		while ( isspace ( *exclude ) ) exclude++;
		if ( isupper ( *exclude ) ) {
			peptide_exclude [*exclude++] = 1;
		}
		while ( isspace ( *exclude ) ) exclude++;
	}

	return ( 0 );
}
int composition_search ( char* peptide )
{
	int temp_peptide_include [AA_ARRAY_SIZE];
	int i;

	gen_bcopy ( (char*)peptide_include, (char*)temp_peptide_include, AA_ARRAY_SIZE * sizeof (int) );
	while ( *peptide ) {
		temp_peptide_include [*peptide]--;
		if ( peptide_include_alternatives [*peptide] ) {
			/*decrements the matched character */
			temp_peptide_include [*peptide]--;
			for ( i = 0 ; peptide_include_alternatives [*peptide][i] != 0 ; i++ ) {
				/*decrements the alternative */
				temp_peptide_include [peptide_include_alternatives [*peptide][i]]--;
			}
		}
		if ( peptide_exclude [*peptide] ) return ( 0 );
		peptide++;
	}
	for ( i = 'A' ; i <= 'Z' ; i++ ) {
		if ( temp_peptide_include [i] > 0 ) return ( 0 );
	}
	return ( 1 );
}
