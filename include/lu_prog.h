/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_prog.h                                                     *
*                                                                             *
*  Created    : June 12th 2001                                                *
*                                                                             *
*  Purpose    : Functions to define argruments for msdigest, msproduct and    *
*               msisotope when called from a results page.                    *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_prog_h
#define __lu_prog_h

#include <ostream>
#include <string>

class ParameterList;

class ProgramLink {
protected:
	static const ParameterList* params;
public:
	static void openLink ( std::ostream& os, const std::string& linkName, int num, bool newWindow = true );
	static void closeLink ( std::ostream& os );
	static std::string getURLStart ( const std::string& program, const std::string& type = ".cgi" );
	static const ParameterList* getParams () { return params; };
	static void setParams ( const ParameterList* p ) { params = p; };
	static void putForm ( std::ostream& os, const std::string& variableName, const std::string& script );
	static void putFormLink ( std::ostream& os, const std::string& variableName, const std::string& text );
};

#endif /* ! __lu_prog_h */
