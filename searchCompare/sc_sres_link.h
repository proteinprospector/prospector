/******************************************************************************
*                                                                             *
*  Program    : searchCompare                                                 *
*                                                                             *
*  Filename   : sc_sres_link.h                                                *
*                                                                             *
*  Created    : July 15th 2003                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2003-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __sc_sres_link_h
#define __sc_sres_link_h

#include <lu_prog.h>

class SResLink {
	std::string urlProg;
	mutable StringVector idFilterList;
	bool xlink;
	void putCGI ( std::ostream& os ) const;
public:
	SResLink ( bool rawForwarding, bool xlink = false );
	void write ( std::ostream& os, const std::string& accessionNumber, const std::string& id, const std::string& rank ) const;
	void printHTML ( std::ostream& os ) const;
};

#endif /* ! __sc_sres_link_h */
