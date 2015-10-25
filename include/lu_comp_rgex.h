/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_comp_rgex.h                                                *
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

#ifndef __lu_comp_rgex_h
#define __lu_comp_rgex_h

int init_composition_search ( char* include, char* exclude );
int composition_search ( char* peptide );

#endif /* ! __lu_comp_rgex_h */
