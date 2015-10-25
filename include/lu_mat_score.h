/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_mat_score.h                                                *
*                                                                             *
*  Created    : September 23rd 1999                                           *
*                                                                             *
*  Purpose    :                                                               *
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

#ifndef __lu_mat_score_h
#define __lu_mat_score_h

#include <string>
#include <lgen_define.h>

void set_score_matrix ( const std::string& score_matrix );
int matrix_score ( const char* qseq, const char* dbseq );
bool matrixMatch( const std::string& seq1, const std::string& seq2 );
StringVector getScoreMatrixNameList ();

#endif /* ! __lu_mat_score_h */
