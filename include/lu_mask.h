/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_mask.h                                                     *
*                                                                             *
*  Created    : September 13th 2001                                           *
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
*  Copyright (2001-2009) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_mask_h
#define __lu_mask_h

#include <string>
#include <vector>

class CompositionMask {
	unsigned int mask;
public:
	CompositionMask ( const std::string& str );
	bool contains ( const CompositionMask& m ) const { return ( mask & m.mask ) != 0; }
};

typedef std::vector <CompositionMask> CompositionMaskVector;

#endif /* ! __lu_mask_h */
