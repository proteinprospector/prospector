/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_elem_comp.h                                                *
*                                                                             *
*  Created    : 30th September 1997                                           *
*                                                                             *
*  Purpose    : Elemental composition searches.                               *
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

#ifndef __lu_elem_comp_h
#define __lu_elem_comp_h

int calculate_elem_combinations ( const StringVector& superElementList, double mass, double tolerance, bool monoisotopicFlag, const char* start_formula, int maxReportedHits, char*** combinations, int* num_combinations );

#endif /* ! __lu_elem_comp_h */
