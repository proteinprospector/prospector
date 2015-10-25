/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_repos_info.h                                               *
*                                                                             *
*  Created    : November 30th 2009                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2009-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_repos_info_h
#define __lu_repos_info_h

#include <lgen_define.h>
#include <lgen_xml.h>

class InstrumentDir : public GenNameValueStream {	// Old style inst_dir.txt file - just associates a directory with an instrument type
	InstrumentDir ();
public:
	~InstrumentDir ();
	static InstrumentDir& instance ();
	std::string getInstrumentType ( const StringVector& fullFileList );
};

class InstrumentType {
	std::string fileSuffix;
	MapStringToStringVector params;
public:
	InstrumentType () {}
	InstrumentType ( const std::string& fileSuffix, const MapStringToStringVector& params ) :
		fileSuffix ( fileSuffix ),
		params ( params ) {}
	std::string getFileSuffix () const { return fileSuffix; }
	MapStringToStringVector getParams () const { return params; }
};

typedef std::map <PairStringInt, InstrumentType> MapPairStringIntInstrumentType;
typedef MapPairStringIntInstrumentType::const_iterator MapPairStringIntInstrumentTypeConstIterator;

typedef std::map <std::string, MapPairStringIntInstrumentType> MapStringMapPairStringIntInstrumentType;
typedef MapStringMapPairStringIntInstrumentType::const_iterator MapStringMapPairStringIntInstrumentTypeConstIterator;

class RepositoryInfo : public PPExpat {
	MapStringMapPairStringIntInstrumentType instruments;
	MapStringToString rawDir;
	MapStringToString rawType;
	static MapPairStringIntInstrumentType instTypes;
	static int index;
	static std::string directoryName;
	static std::string rDir;
	static std::string rType;
	static std::string typeName;
	static std::string suffix;
	static MapStringToStringVector params;
	static bool parametersFlag;
	static std::string paramName;
	static std::string value;
	static std::string paramFilePath;
	void startElement ( const char* name, const char** attributes );
	void characterDataHandler ( const char* str, int len );
	void endElement ( const char* name );
	RepositoryInfo ();
	std::pair <std::string, const InstrumentType*> getInstrumentTypeFromPath ( const std::string& filepath ) const;
	std::string getInstrumentType ( const std::string& filepath ) const;
	MapStringToStringVector getInstrumentParams ( const std::string& filepath ) const;
public:
	~RepositoryInfo ();
	static RepositoryInfo& instance ();
	static bool exists ();
	std::string getInstrumentType ( const StringVector& fullFileList ) const;
	MapStringToStringVector getInstrumentParams ( const StringVector& fullFileList ) const;
	std::string getCentroidSuffix ( const std::string& filepath ) const;
	std::string getAdjustedRawPath ( const std::string& filepath ) const;
	std::string getRawType ( const std::string& filepath ) const;
	StringVector getNames () const;
	MapPairStringIntInstrumentType getInstrumentInfo ( const std::string& name ) const;
};

#endif /* ! __lu_repos_info_h */
