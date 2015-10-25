/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_acc_link.h                                                 *
*                                                                             *
*  Created    : September 5th 2001                                            *
*                                                                             *
*  Purpose    : Functions for dealing with accession number links.            *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_acc_link_h
#define __lu_acc_link_h

#include <vector>
#include <string>
#include <ostream>
#include <lgen_define.h>

class AccessionNumberLinkInfo {
	int ind;
	StringVector filenamePrefix;
	StringVector accessionNumberLink;
	static StringVector aLink;
	int longestMatchingPrefix ( const std::string& str ) const;
	std::string getAccessionNumberLink ( const std::string& database ) const;
	static int num;
	void init ( const std::string& filename );
public:
	AccessionNumberLinkInfo ();
	AccessionNumberLinkInfo ( const std::string& linksName );
	static void write ( std::ostream& os, const std::string& aNum, int num, bool link );
	void printHTML ( std::ostream& os, const std::string& database ) const;
	void write2 ( std::ostream& os, const std::string& aNum, bool link ) const;
};

#endif /* ! __lu_acc_link_h */
