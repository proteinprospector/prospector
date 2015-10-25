/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_mat_score.cpp                                              *
*                                                                             *
*  Created    : June 1st 1999                                                 *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1999-2008) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lu_getfil.h>
using std::string;
using std::fill;

static const int MATRIX_SIZE = 26;
static int*** matrix;
static int** chosenMatrix = NULL;
static StringVector names;
static bool initialised = false;

static void initialise_mat_score ();
static inline void initialise ();

static inline void initialise ()
{
	if ( initialised == false ) initialise_mat_score ();
}
void set_score_matrix ( const string& score_matrix )
{
	initialise ();
	for ( StringVectorSizeType i = 0 ; i < names.size () ; i++ ) {
		if ( score_matrix == names [i] ) {
			chosenMatrix = matrix [i];
			return;
		}
	}
}
int matrix_score ( const char* qseq, const char* dbseq )
{
	int score = 0;
	while ( *qseq ) {
		int i = *qseq++ - 'A';
		int j = *dbseq++ - 'A';
		if ( i < 0 || j < 0 ) continue;
		if ( i < MATRIX_SIZE && j < MATRIX_SIZE ) score += chosenMatrix [i][j];
	}
	return ( score );
}
bool matrixMatch ( const string& seq1, const string& seq2 )
{
	int len1 = seq1.length ();
	int len2 = seq2.length ();
	if ( len1 == len2 ) {
		for ( int i = 0 ; i < len1 ; i++ ) {
			const char s1 = seq1 [i];
			const char s2 = seq2 [i];
			if ( s1 != s2 ) {
				if		( s1 == 'I' && s2 == 'L' ) continue;
				else if	( s1 == 'L' && s2 == 'I' ) continue;
				else if	( s1 == 'K' && s2 == 'Q' ) continue;
				else if	( s1 == 'Q' && s2 == 'K' ) continue;
				return false;
			}
		}
		return true;
	}
	return false;
}
StringVector getScoreMatrixNameList ()
{
	initialise ();
	return names;
}
static void initialise_mat_score ()
{
	int numEntries;
	int j, k;
	int score [MATRIX_SIZE];

	char* info = getFileInfo ( MsparamsDir::instance ().getParamPath ( "mat_score.txt" ), '>', 1, true, &numEntries );
	matrix = new int** [numEntries];
	for ( int i = 0 ; i < numEntries ; i++ ) {
		names.push_back ( ( i == 0 ) ? strtok ( info, "\n" ) : strtok ( NULL, "\n" ) );
		StringVector lines;
		for ( ; ; ) {
			string line = strtok ( NULL, "\n" );
			if ( line == ">" ) break;
			lines.push_back ( line );
		}
		int numLines = lines.size ();
		matrix [i] = new int* [MATRIX_SIZE];
		for ( j = 0 ; j < MATRIX_SIZE ; j++ ) {
			matrix [i][j] = new int [MATRIX_SIZE];
			fill ( matrix [i][j], matrix [i][j] + MATRIX_SIZE, 0 );
		}
		char* order = new char [numLines];
		for ( j = 0 ; j < numLines ; j++ ) {
			sscanf ( lines [j].c_str (), "%c", order + j );
		}
		for ( j = 0 ; j < numLines ; j++ ) {
			char dummy;
			sscanf ( lines [j].c_str (), "%c %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",
				&dummy,
				score, score+1, score+2, score+3, score+4, score+5, score+6, score+7, score+8, score+9, score+10, score+11, score+12,
				score+13, score+14, score+15, score+16, score+17, score+18, score+19, score+20, score+21, score+22, score+23, score+24, score+25 );
			for ( k = 0 ; k <= j ; k++ ) {
				int row = order[j]-'A';
				int col = order[k]-'A';
				matrix [i][row][col] = matrix [i][col][row] = score [k];
			}
		}
		delete [] order;
	}
	initialised = true;
}
