/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_pp_param.cpp                                               *
*                                                                             *
*  Created    : July 17th 2003                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2003-2009) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <cstdlib>
#include <lu_pp_param.h>
using std::string;

namespace PPParamParser {

static string getParameter ( const string& line, const string& param )
{
	int start = line.find ( param );
	if ( start == string::npos ) {
		return "";
	}
	else {
		start += param.length ();
		int end = line.find ( "$", start );
		return line.substr ( start, end - start );
	}
}

int getIntValue ( const string& str, const string& tag )
{
	return atoi ( getParameter ( str, tag ).c_str () );
}
double getDoubleValue ( const string& str, const string& tag )
{
	return atof ( getParameter ( str, tag ).c_str () );
}
string getStringValue ( const string& str, const string& tag )
{
	return getParameter ( str, tag );
}

}
