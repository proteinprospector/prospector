/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_tofil_par.h                                                *
*                                                                             *
*  Created    : October 18th 2001                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_tofil_par_h
#define __lu_tofil_par_h

#include <string>
#include <ostream>

class ParameterList;

class ToFileParameters {
	bool toFile;
	std::string outputPath;
	std::string outputFilename;
	std::string scriptName;
	std::string scriptType;
	bool script;
public:
	ToFileParameters ( const ParameterList* params, const std::string& programName );
	bool getToFile () const { return toFile; }
	std::string getOutputFilename () const { return outputFilename; }
	bool isRandomSearch () const;
	std::string getOutputPath ( const std::string& suffix ) const;
	void printOutputFileLink ( std::ostream& os, const std::string& linkText, const std::string& fileSuffix ) const;
	bool getScript () const { return script; }
	std::string getScriptType () const { return scriptType; }
	std::string getScriptName () const { return scriptName; }
};

#endif /* ! __lu_tofil_par_h */
