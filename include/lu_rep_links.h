/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_rep_links.h                                                *
*                                                                             *
*  Created    : April 29th 2014                                               *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2014-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_rep_links_h
#define __lu_rep_links_h

#include <string>
#include <ostream>
#include <lgen_define.h>

class ParameterList;

class ReportLinkProgram {
	std::string linkProgram;
	std::string linkName;
public:
	ReportLinkProgram ( const std::string& linkProgram, const std::string& linkName );
	void printHTML ( std::ostream& os, const ParameterList* pList ) const;
	void write ( std::ostream& os, const std::string& p, const std::string& l ) const;
};

class ReportLinks {
	std::string url;
	int rowsPerPage;
	int pageNumber;
	VectorPairStringString vpss;
	PairStringString getPair ( const std::string& form, const std::string& name ) const;
public:
	ReportLinks ( const std::string& url, int rowsPerPage, int pageNumber, int numDataLines );
	void printHTML ( std::ostream& os, const ReportLinkProgram& link ) const;
};



#endif /* ! __lu_rep_links_h */
