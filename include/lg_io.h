/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lg_io.h                                                       *
*                                                                             *
*  Created    : July 25rd 2000                                                *
*                                                                             *
*  Purpose    : Error checking interface to io library.                       *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2000-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lg_io_h
#define __lg_io_h

#include <fstream>
#include <iostream>
#include <string>
#include <lgen_define.h>
#include <sstream>

class GenIFStream : public std::ifstream {
public:
	GenIFStream ( const std::string& fullPath, std::ios_base::openmode mode = std::ios_base::in );
	GenIFStream ();
	void open ( const std::string& fullPath, std::ios_base::openmode mode = std::ios_base::in );
};
class GenCommentedIFStream : public GenIFStream {
public:
	GenCommentedIFStream ( const std::string& fullPath );
	bool getUncommentedLine ( std::string& line );
};

class GenNameValueStream : public GenCommentedIFStream {
	MapStringToStringVector params;
public:
	GenNameValueStream ( const std::string& fullPath );
	bool getValue ( const std::string& name, std::string& value );
	bool getValue ( const std::string& name, int& value );
	bool getValue ( const std::string& name, unsigned int& value );
	bool getValue ( const std::string& name, double& value );
	bool getValue ( const std::string& name, StringVector& value ) const;
	bool getBoolValue ( const std::string& name, bool defaultVal = false ) const;
	int getIntValue ( const std::string& name, int defaultVal = 0 ) const;
	unsigned int getUIntValue ( const std::string& name, unsigned int defaultVal = 0 ) const;
	double getDoubleValue ( const std::string& name, double defaultVal = 0.0 ) const;
	std::string getStringValue ( const std::string& name, const std::string& defaultVal = "" ) const;
	StringVector getStringVectorValue ( const std::string& name, const StringVector& defaultVal = StringVector () ) const;
	StringVector getNameList () const;
};

class GenOFStream : public std::ofstream {
public:
	GenOFStream ( const std::string& fullPath, std::ios_base::openmode mode = std::ios_base::out );
	void open ( const std::string& fullPath, std::ios_base::openmode mode = std::ios_base::out );
};

void genPrint ( std::ostream& os, double number, int precision = -1, int width = -1 );
void genPrintSigFig ( std::ostream& os, double number, int precision );

std::istream& genUniversalGetLine ( std::istream& is, std::string& s );

#endif /* ! __lg_io_h */
