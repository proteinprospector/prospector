/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_param_list.h                                               *
*                                                                             *
*  Created    : January 31st 2001                                             *
*                                                                             *
*  Purpose    : Functions to read in parameters.                              *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_param_list_h
#define __lu_param_list_h

#include <iostream>
#include <string>

#include <lg_io.h>
#include <lgen_define.h>
#include <lu_xml.h>

class CgiEnvironmentVariables;
class GenElapsedTime;

class ParameterList {
	std::string programBinaryName;
	MapStringToStringVector params;
	MapStringToString fileContents;
	std::string logFileName;
	GenElapsedTime* elapsedTime;
	static void parseCookie ( const std::string& c );
	void parseLine ( const std::string& line );
	void parseQueryString ( const std::string& queryString );
	void parseXMLStringFromFile ( const std::string& file );
	void parseXMLString ( const std::string& xmlString );
	void parseMultipartQueryEntry ( const std::string& contentType, const std::string& contentLength, const std::string& serverProtocol );
	bool writeUploadFile ( const std::string& filepath, const std::string& dataStart );
	static std::string getTempFileFullPath ( const std::string& filename, bool tempFlag );
	static void writeStdinToLog ( unsigned int cl );
	static std::string virtualDir;
	static std::string serverName;
	static std::string serverPort;
	static std::string server;
	static std::string remoteAddr;
	static std::string method;
	static MapStringToString cookie;
	static const size_t BLOCK_SIZE;
	static std::string batchFileEscapeURL ( const std::string& url );
public:
	ParameterList ( int argc, char** argv );
	ParameterList ( const std::string& filename );
	ParameterList ( const std::string& filename, bool paramFile );
	ParameterList ( std::istream& is );
	ParameterList ( const std::string& str, bool flag1, bool flag2 );
	ParameterList ( const std::string& str, bool flag1, bool flag2, bool flag3 );
	ParameterList ( const std::string& filename, bool flag1, bool flag2, bool flag3, bool flag4 );
	ParameterList ( const std::string& programBinaryName, bool flag1, bool flag2, bool flag3, bool flag4, bool flag5 );
	~ParameterList ();
	void writeLogStart ( const CgiEnvironmentVariables& cgiEnvVar );
	void writeLogError ( const std::string& err ) const;
	void writeLogEnd () const;
	static std::string getVirtualDir () { return virtualDir; }
	static std::string getServerName () { return serverName; }
	static std::string getServerPort () { return serverPort; }
	static void setVirtualDir ( const std::string& v ) { virtualDir = v; }
	static void setServerName ( const std::string& n ) { serverName = n; }
	static void setServerPort ( const std::string& p ) { serverPort = p; }
	static std::string getServer () { return server; }
	static std::string getMethod () { return method; }
	static std::string getRemoteAddr () { return remoteAddr; }
	std::string getProgramBinaryName () const { return programBinaryName; }
	bool getValue ( const std::string& name, bool& value ) const;
	bool getValue ( const std::string& name1, const std::string& name2, bool& value ) const;
	bool getValue ( const std::string& name, char& value ) const;
	bool getValue ( const std::string& name, int& value ) const;
	bool getValue ( const std::string& name, double& value ) const;
	bool getValue ( const std::string& name, std::string& value ) const;
	bool getValue ( const std::string& name, const char*& value ) const;
	bool getValue ( const std::string& name, CharVector& value ) const;
	bool getValue ( const std::string& name, StringVector& value ) const;
	bool getBoolValue ( const std::string& name, bool defaultVal = false ) const;
	int getIntValue ( const std::string& name, int defaultVal = 0 ) const;
	unsigned int getUIntValue ( const std::string& name, unsigned int defaultVal = 0 ) const;
	char getCharValue ( const std::string& name, char defaultVal = '\0' ) const;
	std::string getStringValue ( const std::string& name, const std::string& defaultVal = "" ) const;
	std::string getStrippedStringValue ( const std::string& name, const std::string& defaultVal = "" ) const;
	std::string getFileStringValue ( const std::string& name, const std::string& defaultVal ) const;
	double getDoubleValue ( const std::string& name, double defaultVal = 0.0 ) const;
	IntVector getIntVectorValue ( const std::string& name ) const;
	StringVector getUniqueStringVectorValue ( const std::string& name ) const;
	StringVector getStringVectorValue ( const std::string& name ) const;
	StringVector getPQStringVectorValue ( const std::string& name ) const;
	BoolDeque getBoolDequeValue ( const std::string& name, int size ) const;
	std::string getDirectoryValue ( const std::string& name, const std::string& defaultVal = "" ) const;
	static std::string getCookieValue ( const std::string& name, const std::string& defaultVal = "" );
	std::string getFileContents ( const std::string& name, const std::string& defaultVal = "" ) const;
	StringVector getNameList ();
	int size () { return params.size (); }
	bool empty () { return params.size () == 0; }

	void removeName ( const std::string& name );
	void addName ( const std::string& name, const std::string& value );
	void addOrReplaceName ( const std::string& name, int value );
	void addOrReplaceName ( const std::string& name, const std::string& value );
	void addOrReplaceName ( const std::string& name, const StringVector& value );
	void setValue ( const std::string& name, const std::string& value );
	void copyName ( const ParameterList* plist, const std::string& name );
	void copyFileContents ( const ParameterList* plist, const std::string& name );
	void appendParameters ( const ParameterList* newParams );
	void appendParameters ( const ParameterList& newParams );
	bool isEqual ( const ParameterList* plist, const std::string& name ) const;

	void HTMLParameters ( std::ostream& os ) const;
	void XMLParameterFile ( const std::string& filename ) const;
	void XMLParameters ( std::ostream& os ) const;
	void getParams ( StringVector& sv, StringVectorVector& svv ) const;
	void pepXMLParameters ( std::ostream& os, int ntab ) const;
	VectorXMLOutputItemPtr getMZIdentMLUserParameters () const;
	void perlFileParameters ( std::ostream& os ) const;
	void batchFileParameters ( std::ostream& os ) const;
	std::string getURL () const;
	bool copyToCookie ( std::ostream& os, const std::string& cookieName ) const;
	void copyToCGI ( std::ostream& os ) const;
	bool copyToCGI ( std::ostream& os, const std::string& name ) const;
	bool copyToCGI ( std::ostream& os, const std::string& name, int num ) const;
	void copyToHiddenFormEntry ( std::ostream& os ) const;
	bool copyToHiddenFormEntry ( std::ostream& os, const std::string& name ) const;
	bool copyToHiddenFormEntry ( std::ostream& os, const std::string& name, int num ) const;
	std::string getCommandLineNVPair ( const std::string& name ) const;
	bool copyToHiddenFormJavascriptEntry ( std::ostream& os, const std::string& name ) const;
	bool copyToHiddenFormJavascriptEntry ( std::ostream& os, const std::string& name, int num ) const;

	template <class T>
		static void printXML ( std::ostream& os, const std::string& name, const T& value )
		{ os << "<" << name << ">" << value << "</" << name << ">" << std::endl;}
	static void printXML ( std::ostream& os, const std::string& name, int number, const std::string& value )
		{ os << "<" << name << ">" << number << value << "</" << name << ">" << std::endl;}
	template <class T>
		static void printXMLNC ( std::ostream& os, const std::string& name, T& value )
		{ os << "<" << name << ">" << value << "</" << name << ">" << std::endl;}
	static void printDoubleXMLFixed ( std::ostream& os, const std::string& name, double value, int precision )
	{
		os << "<" << name << ">";
		genPrint ( os, value, precision );
		os << "</" << name << ">";
		os << std::endl;
	}
	static void printDoubleXMLSigFig ( std::ostream& os, const std::string& name, double value, int precision )
	{
		os << "<" << name << ">";
		genPrintSigFig ( os, value, precision );
		os << "</" << name << ">";
		os << std::endl;
	}
	static void printXMLContainer ( std::ostream& os, const std::string& name, const IntVector& value );
	static void printXMLContainer ( std::ostream& os, const std::string& name, const DoubleVector& value );
	static void printXMLContainer ( std::ostream& os, const std::string& name, const StringVector& value );

	template <class T>
		static void printHTML ( std::ostream& os, const std::string& name, const T& value )
	{ os << name << ": <b>" << value << "</b><br />" << std::endl; }
	template <class T>
		static void printHTMLRange ( std::ostream& os, const std::string& name, const T& value1, const T& value2 )
	{ os << name << ": <b>" << value1 << " - " << value2 << "</b><br />" << std::endl;}
	template <class T>
		static void printHTMLNC ( std::ostream& os, const std::string& name, T& value )
	{ os << name << ": <b>" << value << "</b><br />" << std::endl;}
	static void printDoubleHTMLFixed ( std::ostream& os, const std::string& name, double value, int precision )
	{
		os << name << ": <b>";
		genPrint ( os, value, precision );
		os << "</b><br />";
		os << std::endl;
	}
	static void printDoubleHTMLFixedRange ( std::ostream& os, const std::string& name, double value1, double value2, int precision )
	{
		os << name << ": <b>";
		genPrint ( os, value1, precision );
		os << " - ";
		genPrint ( os, value2, precision );
		os << "</b><br />";
		os << std::endl;
	}
	static void printDoubleHTMLSigFig ( std::ostream& os, const std::string& name, double value, int precision, const std::string& units = "" )
	{
		os << name << ": <b>";
		genPrintSigFig ( os, value, precision );
		if ( units != "" ) os << " " << units;
		os << "</b><br />";
		os << std::endl;
	}
	static void printDoubleHTMLSigFigRange ( std::ostream& os, const std::string& name, double value1, double value2, int precision, const std::string& units = "" )
	{
		os << name << ": <b>";
		genPrintSigFig ( os, value1, precision );
		os << " - ";
		genPrintSigFig ( os, value2, precision );
		if ( units != "" ) os << " " << units;
		os << "</b><br />";
		os << std::endl;
	}
	static void printHTMLContainer ( std::ostream& os, const std::string& name, const StringVector& value );
	static std::string spaceToPlus ( const std::string& str1 );
	static std::string convertValue ( const std::string& oldValue );
};

#ifndef VIS_C
void initParamSigtermHandler ( ParameterList* pList );
#endif

#endif /* ! __lu_param_list_h */
