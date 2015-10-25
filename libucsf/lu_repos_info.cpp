/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_repos_info.cpp                                             *
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
#ifdef VIS_C
#pragma warning( disable : 4503 )	// Disable warnings about decorated name length being too long
#endif
#ifndef VIS_C
#include <stdexcept>
#endif
#include <lg_string.h>
#include <lgen_file.h>
#include <lgen_error.h>
#include <lu_getfil.h>
#include <lu_repos_info.h>

using std::string;
using std::runtime_error;
using std::pair;
using std::make_pair;

InstrumentDir::InstrumentDir () :
	GenNameValueStream ( MsparamsDir::instance ().getParamPath ( "inst_dir.txt" ) )
{
}
InstrumentDir::~InstrumentDir () {}
InstrumentDir& InstrumentDir::instance ()
{
	static InstrumentDir id;
	return id;
}
string InstrumentDir::getInstrumentType ( const StringVector& fullFileList )
{
	string inst;
	string oldInst;
	for ( StringVectorSizeType i = 0 ; i < fullFileList.size () ; i++ ) {
		int end = fullFileList [i].find ( "/" );
		string instDir = fullFileList [i].substr ( 0, end );
		if ( !getValue ( instDir, inst ) ) {
			throw runtime_error ( "Instrument directory " + instDir + " not in inst_dir.txt parameter file." );
		}
		if ( i != 0 ) {
			if ( inst != oldInst ) {
				throw runtime_error ( "Files from different instrument types can't be used in the same project." );
			}
		}
		else oldInst = inst;
	}
	return inst;
}
MapPairStringIntInstrumentType RepositoryInfo::instTypes;
int RepositoryInfo::index = 1;
string RepositoryInfo::directoryName;
string RepositoryInfo::rDir;
string RepositoryInfo::rType;
string RepositoryInfo::typeName;
string RepositoryInfo::suffix;
MapStringToStringVector RepositoryInfo::params;
bool RepositoryInfo::parametersFlag = false;
string RepositoryInfo::paramName;
string RepositoryInfo::value;
string RepositoryInfo::paramFilePath = MsparamsDir::instance ().getParamPath ( "repository.xml" );

RepositoryInfo::RepositoryInfo ()
{
	parametersFlag = false;
	try {
		parseXMLFromFile ( paramFilePath );
	}
	catch ( runtime_error e ) {
		ErrorHandler::genError ()->error ( e );
	}
}
RepositoryInfo::~RepositoryInfo () {}
RepositoryInfo& RepositoryInfo::instance ()
{
	static RepositoryInfo ri;
	return ri;
}
bool RepositoryInfo::exists ()
{
	return genFileExists ( paramFilePath );
}
void RepositoryInfo::startElement ( const char* name, const char** attributes )
{
	if ( parametersFlag ) {
		paramName = name;
		value = "";
	}
	else if ( !strcmp ( name, "parameters" ) ) {
		parametersFlag = true;
	}
	else if ( !strcmp ( name, "directory" ) ) {
		getAttributeValue ( attributes, "name", directoryName );
		getAttributeValue ( attributes, "raw_dir", rDir );
		getAttributeValue ( attributes, "raw_type", rType );
	}
	else if ( !strcmp ( name, "type" ) ) {
		suffix = "";
		typeName = "";
		params.clear ();
		getAttributeValue ( attributes, "name", typeName );
		getAttributeValue ( attributes, "suffix", suffix );
	}
	else if ( !strcmp ( name, "instrument" ) ) {
		instTypes.clear ();
		directoryName = "";
		rDir = "";
		rType = "";
	}
}
void RepositoryInfo::characterDataHandler ( const char* str, int len )
{
	if ( parametersFlag ) value.append ( str, len );
}
void RepositoryInfo::endElement ( const char* name )
{
	if ( !strcmp ( name, "parameters" ) ) {
		parametersFlag = false;
	}
	else if ( parametersFlag ) {
		params [paramName].push_back ( value );
	}
	else if ( !strcmp ( name, "type" ) ) {
		instTypes [make_pair ( typeName, index )] = InstrumentType ( suffix, params );
		index++;
	}
	else if ( !strcmp ( name, "instrument" ) ) {
		instruments [directoryName] = instTypes;
		if ( !rDir.empty () ) rawDir [directoryName] = rDir;
		if ( !rType.empty () ) rawType [directoryName] = rType;
	}
}
pair <string, const InstrumentType*> RepositoryInfo::getInstrumentTypeFromPath ( const string& filepath ) const
{
	string shortFilename = genShortFilenameFromPath ( filepath );	// File name with no suffix
	pair <string, const InstrumentType*> ret ( "", 0 );
	string instDir = filepath.substr ( 0, filepath.find ( "/" ) );
	MapStringMapPairStringIntInstrumentTypeConstIterator cur = instruments.find ( instDir );
	if ( cur != instruments.end () ) {
		const MapPairStringIntInstrumentType& msti = (*cur).second;
		for ( MapPairStringIntInstrumentTypeConstIterator i = msti.begin () ; i != msti.end () ; i++ ) {
			string suffix = (*i).second.getFileSuffix ();
			if ( suffix.empty () ) ret = make_pair ( (*i).first.first, &(*i).second );	// If there is a no suffix condition this is the default
			if ( isNoCaseSuffix ( shortFilename, suffix ) ) return make_pair ( (*i).first.first, &(*i).second );
		}
	}
	else {
		throw runtime_error ( "Instrument directory " + instDir + " not in the repository.xml parameter file." );
	}
	if ( ret.second == 0 ) {
		throw runtime_error ( "Suffix match not found in the repository.xml parameter file." );
	}
	return ret;
}
string RepositoryInfo::getInstrumentType ( const string& filepath ) const
{
	pair <string, const InstrumentType*> p = getInstrumentTypeFromPath ( filepath );
	return p.first;
}
MapStringToStringVector RepositoryInfo::getInstrumentParams ( const string& filepath ) const
{
	pair <string, const InstrumentType*> p = getInstrumentTypeFromPath ( filepath );
	return p.second->getParams ();
}
string RepositoryInfo::getInstrumentType ( const StringVector& fullFileList ) const
{
	string inst;
	string oldInst;
	for ( StringVectorSizeType i = 0 ; i < fullFileList.size () ; i++ ) {
		inst = getInstrumentType ( fullFileList [i] );
		if ( i != 0 ) {
			if ( inst != oldInst ) {
				throw runtime_error ( "Files from different instrument types can't be used in the same project." );
			}
		}
		else oldInst = inst;
	}
	return inst;
}
MapStringToStringVector RepositoryInfo::getInstrumentParams ( const StringVector& fullFileList ) const
{
	MapStringToStringVector params;
	for ( StringVectorSizeType i = 0 ; i < fullFileList.size () ; i++ ) {
		params = getInstrumentParams ( fullFileList [i] );
		if ( !params.empty () ) return params;
	}
	return params;
}
string RepositoryInfo::getCentroidSuffix ( const string& filepath ) const
{
	pair <string, const InstrumentType*> p = getInstrumentTypeFromPath ( filepath );
	return p.second->getFileSuffix ();
}
string RepositoryInfo::getAdjustedRawPath ( const string& filepath ) const
{
	string instDir = filepath.substr ( 0, filepath.find ( "/" ) );
	MapStringToStringConstIterator cur = rawDir.find ( instDir );
	if ( cur != rawDir.end () )
		return genReplaceSubstrings ( filepath, instDir, (*cur).second );
	else
		return filepath;
}
string RepositoryInfo::getRawType ( const string& filepath ) const
{
	string instDir = filepath.substr ( 0, filepath.find ( "/" ) );
	MapStringToStringConstIterator cur = rawType.find ( instDir );
	if ( cur != rawType.end () )
		return (*cur).second;
	else
		return "";
}
StringVector RepositoryInfo::getNames () const
{
	StringVector sv;
	for ( MapStringMapPairStringIntInstrumentTypeConstIterator i = instruments.begin () ; i != instruments.end () ; i++ ) {
		sv.push_back ( (*i).first );
	}
	return sv;
}
MapPairStringIntInstrumentType RepositoryInfo::getInstrumentInfo ( const string& name ) const
{
	MapStringMapPairStringIntInstrumentTypeConstIterator cur = instruments.find ( name );
	return (*cur).second;
}
