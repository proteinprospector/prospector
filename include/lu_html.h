/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_html.h                                                     *
*                                                                             *
*  Created    : July 20th 1996                                                *
*                                                                             *
*  Purpose    : Functions for printing out commonly used elements of HTML     *
*               reports.                                                      *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1996-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_html_h
#define __lu_html_h

#include <string>
#include <vector>
#include <lgen_define.h>
#include <lgen_error.h>

void init_html ( std::ostream& os, const std::string& title, const std::string& manual = "", bool reset = false );
void init_html_premature_stop ( const std::string& programName, bool stopButton, const StringVector& filesToDelete = StringVector () );
void printAbortFunctions ( std::ostream& os );
void printAbortButton ( std::ostream& os, const std::string& searchJobKey );
void htmlPrintProteinSequence ( std::ostream& os, const std::string& protein, const IntVector& cleavageIndicies, const CharVector& coverageMap, bool compressedDisplay );
void printHTMLLink ( std::ostream& os, const std::string& outputPath, const std::string& filename, const std::string& label );
void print_charge_superscript ( std::ostream& os, int charge );
void printNumHits ( std::ostream& os, const std::string& programName, int numHits, int maxReportedHits );
void startJavascript ( std::ostream& os, const std::string& version = "" );
void endJavascript ( std::ostream& os );
void htmlPreformattedTextFile ( std::ostream& os, const std::string& fileName );

std::string get_msproduct_charge_color ( int charge );

static std::string mh_plus_html	= "MH<sup>+</sup>";

class HTMLErrorHandler : public ErrorHandler {
	void messageDisplay ( const std::string& messageString );
	void errorDisplay ( const std::string& errorString, bool endProgram );
public:
	HTMLErrorHandler ();
};

void basicJavascriptFunctions ( std::ostream& os );
void updateJavascriptList ( const std::string& jsFilename, const StringVector& nameList, const std::string& selectedName );

class ExpandableJavascriptBlock {
	std::string title;
	std::string idName;
	std::string linkName;
	std::string defaultToggle;
	static int blockNumber;
	static std::string idNamePrefix;
	static bool functionsPrinted;
	static void printFunctions ( std::ostream& os );
public:
	ExpandableJavascriptBlock ( const std::string& title, bool open = false );
	void printHeader ( std::ostream& os );
	void printFooter ( std::ostream& os );
};

class UpdatingJavascriptMessage {
	static std::string outstandingSpanID;
	int index;
	std::string id;
	std::string spanID;
public:
	UpdatingJavascriptMessage ();
	void writeMessage ( std::ostream& os, int num );
	void writeMessage ( std::ostream& os, const std::string& message );
	void startWriteMessage ( std::ostream& os );
	void addMessage ( std::ostream& os, const std::string& message ) const;
	void endWriteMessage ( std::ostream& os );
	void deletePreviousMessage ( std::ostream& os ) const;
	static void deleteOutstandingSpanID ( std::ostream& os );
};

class CheckboxArraySettingJavascript {
	static int num;
	static void printFunctions ( std::ostream& os, const BoolDeque& flagList );
public:
	CheckboxArraySettingJavascript ( std::ostream& os, const std::string& str, const std::string& buttonLabel, const BoolDeque& flagList );
};

class CheckboxSettingJavascript {
	StringVector names;
	static bool functionsPrinted;
	void printButton ( std::ostream& os, const std::string& buttonLabel, bool flag ) const;
public:
	CheckboxSettingJavascript ( const std::string& str );
	CheckboxSettingJavascript ( const StringVector& names );
	void print ( std::ostream& os ) const;
	static void printFunctions ( std::ostream& os );
};

void refreshJavascript ( std::ostream& os, int timeout, const std::string& url = "", bool manualRedirectMessage = false, bool automaticUpdateMessage = false );
void back2Javascript ( std::ostream& os );
bool setCookie ( std::ostream& os, const std::string& name, const std::string& value, bool persistant );
void deleteCookie ( std::ostream& os, const std::string& name );
void setValueFromCookie ( std::ostream& os, const std::string& name );
std::string getElementalFormulaHTML ( const std::string& formula );
void compressedHTMLStyle ( std::ostream& os );

#endif /* ! __lu_html_h */
