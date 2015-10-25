/******************************************************************************
*                                                                             *
*  Library    : liboracle                                                     *
*                                                                             *
*  Filename   : lo_define.h                                                   *
*                                                                             *
*  Created    :                                                               *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2002-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifndef __lo_define_h
#define __lo_define_h

// Macro simplified repetative call to OCIDefineByPos
// Note the preincrement on i, which will force us to start at index 1because ocidefinebypos skips position 0, 
// but doesn't return an error. It fails later when you try to fetch with some silly seemingly unrelated error
//
#define OCIDEFINEBYPOS(VALUE, OCITYPE) OCIDefineByPos(	cursor, &pDefn[i], phErr, ++i, &VALUE, sizeof(VALUE), \
														OCITYPE, (dvoid *) 0, (ub2 *) 0, (ub2 *) 0, OCI_DEFAULT)
#endif /* ! __lo_define_h */
