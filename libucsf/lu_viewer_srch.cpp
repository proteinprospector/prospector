/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_viewer_srch.cpp                                            *
*                                                                             *
*  Created    : March 8th 2011                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2011-2015) The Regents of the University of California.         *
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
#include <iostream>
#include <lg_io.h>
#include <lg_string.h>
#include <lg_time.h>
#include <lgen_file.h>
#include <lgen_uncompress.h>
#include <lu_ambiguity.h>
#include <lu_blib.h>
#include <lu_msp.h>
#include <lu_cgi_val.h>
#include <lu_const_mod.h>
#include <lu_file_type.h>
#include <lu_getfil.h>
#include <lu_html.h>
#include <lu_msf.h>
#include <lu_param_list.h>
#include <lu_prod_par.h>
#include <lu_table.h>
#include <lu_viewer_form.h>
#include <lu_viewer_srch.h>
#include <lu_repository.h>
#include <lu_rep_links.h>
#include <lu_xml.h>
#include <lu_delim.h>
#include <lu_xml_data.h>

using std::runtime_error;
using std::ostream;
using std::string;
using std::stable_sort;
using std::remove_if;
using std::vector;
using std::map;
using std::cout;
using std::endl;
using std::ios;
using std::find;
using std::find_if;
using std::fill;
using std::unique;
using std::ostringstream;
using namespace FileTypes;

class SortRowsByColumn {
	IntVector sortLevel;
	CharVector sortOrderType;
	CharVector sortOrderDirection;
public:
	SortRowsByColumn ( int numCols, const IntVector& sortLevel, const StringVector& sot, const StringVector& sod ) :
		sortLevel ( sortLevel )
	{
		for ( int i = 0 ; i < sortLevel.size () ; i++ ) {
			if ( sortLevel [i] > numCols ) {
				throw runtime_error ( "Sort column is greater than the number of columns in the report." );
			}
			sortOrderType.push_back ( ( sot [i] == "Alphabetic" ) ? 'A' : 'N' );
			sortOrderDirection.push_back ( ( sod [i] == "Ascending" ) ? 'A' : 'D' );
		}
	}
	bool operator () ( const StringVector& a, const StringVector& b ) const
	{
		return sortColumns ( a, b, 0 );
	}
	bool sortColumns ( const StringVector& a, const StringVector& b, int level ) const
	{
		if ( level < sortLevel.size () ) {
			if ( sortOrderType [level] == 'A' ) {
				string a1 = a [sortLevel [level]-1];
				string b1 = b [sortLevel [level]-1];
				if ( a1 == b1 ) {
					return sortColumns ( a, b, level+1 );
				}
				else {
					if ( sortOrderDirection [level] == 'A' )
						return a1 < b1;
					else
						return a1 > b1;
				}
			}
			else {
				double a1 = atof ( a [sortLevel [level]-1].c_str () );
				double b1 = atof ( b [sortLevel [level]-1].c_str () );
				if ( a1 == b1 ) {
					return sortColumns ( a, b, level+1 );
				}
				else {
					if ( sortOrderDirection [level] == 'A' )
						return a1 < b1;
					else
						return a1 > b1;
				}
			}
		}
		else
			return false;
	}
};
enum FilterTypes {
    Equals = 1,
    NotEqualTo = 2,
    GreaterThanAlphabetic = 3,
    GreaterThanNumeric = 4,
    LessThanAlphabetic = 5,
    LessThanNumeric = 6,
    Contains = 7,
    DoesNotContain = 8,
    Prefix = 9,
    Suffix = 10
};
class FilterRowsByColumn {
	IntVector filterColumn;
	StringVectorVector filterValue;
	IntVector filterType;
public:
	FilterRowsByColumn ( int numCols, const IntVector& filterColumn, const StringVector& ft, const StringVectorVector& filterValue ) :
		filterColumn ( filterColumn ),
		filterValue ( filterValue )
	{
		for ( int i = 0 ; i < filterColumn.size () ; i++ ) {
			if ( filterColumn [i] > numCols ) {
				throw runtime_error ( "Filter column is greater than the number of columns in the report." );
			}
			int f = 0;
			if		( ft [i] == "Equals" )					f = Equals;
			else if	( ft [i] == "Not Equal To" )			f = NotEqualTo;
			else if	( ft [i] == "Greater Than Alphabetic" )	f = GreaterThanAlphabetic;
			else if	( ft [i] == "Greater Than Numeric" )	f = GreaterThanNumeric;
			else if	( ft [i] == "Less Than Alphabetic" )	f = LessThanAlphabetic;
			else if	( ft [i] == "Less Than Numeric" )		f = LessThanNumeric;
			else if	( ft [i] == "Contains" )				f = Contains;
			else if	( ft [i] == "Does Not Contain" )		f = DoesNotContain;
			else if	( ft [i] == "Prefix" )					f = Prefix;
			else if	( ft [i] == "Suffix" )					f = Suffix;
			filterType.push_back ( f );
		}
	}
	bool operator () ( const StringVector& a )										// return true will remove
	{
		for ( IntVectorSizeType i = 0 ; i < filterColumn.size () ; i++ ) {
			string col = a [filterColumn [i]-1];
			StringVector v = filterValue [i];
			if ( filterType [i] == NotEqualTo ) {
				for ( StringVectorSizeType j = 0 ; j < v.size () ; j++ ) {
					string val = v [j];
					if ( col == val ) return true;
				}
			}
			else if ( filterType [i] == DoesNotContain ) {
				for ( StringVectorSizeType j = 0 ; j < v.size () ; j++ ) {
					string val = v [j];
					if ( col.find ( val ) != string::npos ) return true;
				}
			}
			else {
				bool retain = false;
				for ( StringVectorSizeType j = 0 ; j < v.size () ; j++ ) {
					string val = v [j];
					switch ( filterType [i] ) {
						case Equals:
							if ( col == val ) retain = true;
							break;
						case GreaterThanAlphabetic:
							if ( col > val ) retain = true;
							break;
						case GreaterThanNumeric:
							if ( atof ( col.c_str () ) > atof ( val.c_str () ) ) retain = true;
							break;
						case LessThanAlphabetic:
							if ( col < val ) retain = true;
							break;
						case LessThanNumeric:
							if ( atof ( col.c_str () ) < atof ( val.c_str () ) ) retain = true;
							break;
						case Contains:
							if ( col.find ( val ) != string::npos ) retain = true;
							break;
						case Prefix:
							if ( isNoCasePrefix ( col, val ) ) retain = true;
							break;
						case Suffix:
							if ( isNoCaseSuffix ( col, val ) ) retain = true;
							break;
						default:
							break;
					}
					if ( retain ) break;
				}
				if ( !retain ) return true;
			}
		}
		return false;
	}
};

class TestReplicates {
	IntVector testColumns;
public:
	TestReplicates ( const IntVector& testColumns ) :
		testColumns ( testColumns )
	{
	}
	bool operator () ( const StringVector& a, const StringVector& b )
	{
		for ( IntVectorSizeType i = 0 ; i < testColumns.size () ; i++ ) {
			int tc = testColumns[i]-1;
			if ( a [tc] != b [tc] ) return false;
		}
		return true;
	}
};
enum SpectrumIdentifierTypes {
    SpecIDPPRT = 1,
    SpecIDTitle = 2,
    SpecIDSpecNo = 3,
    SpecIDMZ = 4,
    SpecIDScanNo = 5,
    SpecIDScanNoMinusOne = 6
};

namespace {

string getSpectrumIdentifierTypeString ( int sit )
{
	if ( sit == SpecIDPPRT )				return "Protein Prospector RT";
	else if ( sit == SpecIDTitle )			return "Scan Title (Mascot/X!Tandem)";
	else if ( sit == SpecIDSpecNo )			return "Spectrum Number";
	else if ( sit == SpecIDMZ )				return "m/z";
	else if ( sit == SpecIDScanNo )			return "Scan Number";
	else if ( sit == SpecIDScanNoMinusOne )	return "Scan Number - 1";
	return "";
}

}

class ViewerMSTagLink {
	std::string dataFilename;
	const MSViewerForm* msfv;
	static int num;
	int index;
	std::string getLinkName () const
	{
		if ( num == 0 ) return "vTagLink";
		else			return "vTagLink" + gen_itoa ( index );
	}
	void putFormStart ( ostream& os, const string& variableName ) const;
	void putFormEnd ( ostream& os ) const;
public:
	ViewerMSTagLink ( const string& dataFilename, const MSViewerForm* msfv );
	void putDataFilename ( ostream& os ) const;
	void putForm1 ( ostream& os, const string& variableName, const SpecID& specID ) const;
	void putForm2 ( ostream& os, const string& variableName, const string& title ) const;
	void putForm3 ( ostream& os, const string& variableName, const string& index ) const;
	void putForm4 ( ostream& os, const string& variableName, const string& mOverZ ) const;
	void putForm5 ( ostream& os, const string& variableName, const string& scanNumber ) const;
	static void putFormLink ( ostream& os, const string& variableName, const string& str );
};

class ViewerMSViewerLink {
public:
	ViewerMSViewerLink () {}
	void write ( ostream& os, const string& searchKey ) const;
	void write2 ( ostream& os, const string& searchKey ) const;
	void printHTML ( ostream& os ) const;
};

void ViewerMSViewerLink::printHTML ( ostream& os ) const
{
	os << "viewerLink" << "=\"";
	os << ProgramLink::getURLStart ( "mssearch" );
	os << "?";
	printCGIString ( os, "search_name", "msviewer" );
	printCGIString ( os, "report_title", "MS-Viewer" );
	printCGIString ( os, "rows_per_page", "100" );
	os << "\";\n";
}
void ViewerMSViewerLink::write ( ostream& os, const string& searchKey ) const
{
	ProgramLink::openLink ( os, "viewerLink", -1 );
	printCGIString ( os, "search_key", searchKey );
	os << "\\\">";
	os << searchKey;
	ProgramLink::closeLink ( os );
}
void ViewerMSViewerLink::write2 ( ostream& os, const string& searchKey ) const
{
	ProgramLink::openLink ( os, "viewerLink", -1 );
	printCGIString ( os, "search_key", searchKey );
	printCGIString ( os, "delete", "1" );
	os << "\\\">";
	os << "Delete";
	ProgramLink::closeLink ( os );
}

int ViewerMSTagLink::num = 0;
ViewerMSTagLink::ViewerMSTagLink ( const string& dataFilename, const MSViewerForm* msfv ) :
	dataFilename ( dataFilename ),
	msfv ( msfv )
{
	num++;
	index = num;
}
void ViewerMSTagLink::putDataFilename ( ostream& os ) const
{
	startJavascript ( os );
	os << "function printMSTagLinkItems" << index << " () {" << endl;
		printHTMLFORMJavascriptHidden2 ( os, "data_filename", dataFilename );
	os << "}" << endl;
	endJavascript ( os );
}
void ViewerMSTagLink::putFormStart ( ostream& os, const string& variableName ) const
{
	static int idx = 1;
	if ( idx == 1 ) {		// These are items which are passed in but unchanged between calls
		startJavascript ( os );
		os << "function printMSTagLinkItems () {" << endl;

			msfv->printTagParamsHTMLJavascriptHidden ( os );
			printHTMLFORMJavascriptHidden ( os, "form", "mstagfromviewer" );
			printHTMLFORMJavascriptHidden ( os, "data_source", "Data From File" );
		os << "}" << endl;
		endJavascript ( os );
	}
	idx++;
	ProgramLink::putForm ( os, variableName, "msform.cgi" );
	startJavascript ( os );
	os << "\tprintMSTagLinkItems" << index << " ();" << endl;		// Data file name
	os << "\tprintMSTagLinkItems ();" << endl;						// Generic tag parameters
	endJavascript ( os );
}
void ViewerMSTagLink::putFormEnd ( ostream& os ) const
{
	os << "</form>" << endl;
}
void ViewerMSTagLink::putForm1 ( ostream& os, const string& variableName, const SpecID& specID ) const
{
	putFormStart ( os, variableName );
	specID.printHTMLHidden ( os );
	putFormEnd ( os );
}
void ViewerMSTagLink::putForm2 ( ostream& os, const string& variableName, const string& title ) const
{
	putFormStart ( os, variableName );
	printHTMLFORMHidden ( os, "title", title );
	putFormEnd ( os );
}
void ViewerMSTagLink::putForm3 ( ostream& os, const string& variableName, const string& index ) const
{
	putFormStart ( os, variableName );
	printHTMLFORMHidden ( os, "index", index );
	putFormEnd ( os );
}
void ViewerMSTagLink::putForm4 ( ostream& os, const string& variableName, const string& mOverZ ) const
{
	putFormStart ( os, variableName );
	printHTMLFORMHidden ( os, "m_over_z", mOverZ );
	putFormEnd ( os );
}
void ViewerMSTagLink::putForm5 ( ostream& os, const string& variableName, const string& scanNumber ) const
{
	putFormStart ( os, variableName );
	printHTMLFORMHidden ( os, "scan_number", scanNumber );
	putFormEnd ( os );
}
void ViewerMSTagLink::putFormLink ( ostream& os, const string& variableName, const string& str )
{
	ProgramLink::putFormLink ( os, variableName, str );
}

struct ViewerDataSet {
	string searchKey;
	string dataSet;
	string results;
	GenTime date;

	ViewerDataSet ( const string& searchKey, const string& dataSet, const string& results, const GenTime& date );
};

ViewerDataSet::ViewerDataSet ( const string& searchKey, const string& dataSet, const string& results, const GenTime& date ) :
	searchKey ( searchKey ),
	dataSet ( dataSet ),
	results ( results ),
	date ( date )
{
}

class sortViewerDataSet {
public:
	int operator () ( const ViewerDataSet& a, const ViewerDataSet& b ) const
	{
		return ( a.date > b.date );
	}
};

MSViewerSearch::MSViewerSearch ( const MSViewerParameters& params ) :
	MSProgram ( params ),
	viewerParams ( params ),
	ppSpectrumColumnNumber ( 0 )
{
	init_html ( cout, "MS-Viewer Report" );
	if ( viewerParams.getSaveSettings () ) {
		viewerParams.saveDataSet ();
	}
	else if ( viewerParams.getDeleteFlag () ) {
		viewerParams.deleteDataSet ();
	}
	else {
		processPeakListAndResults ();
		if ( viewerParams.getCommandLine () ) {
			viewerParams.saveDataSet ();
		}
	}
}
MSViewerSearch::~MSViewerSearch ()
{
}
void MSViewerSearch::processPeakListAndResults ()
{
	processPeakListFile ();
	string ff = viewerParams.getResultsFpath ();
	if ( genIsDirectory ( ff ) ) {
		int numPepXML = viewerParams.getNumPepXML ();
		int numMSF = viewerParams.getNumMSF ();
		if ( numMSF ) processMSF ( true );
		if ( numPepXML ) processPepXMLResultsFile ( true );
	}
	else if ( isPepXMLFile ( ff ) ) {
		processPepXMLResultsFile ( false );
	}
	else if ( isPrideXMLFile ( ff ) ) {
		processPrideXMLResultsFile ();
	}
	else if ( isSQLiteFile ( ff ) ) {
		bool blibFlag = false;
		bool msfFlag = false;
		{
			SQLiteIdentify sqlID ( ff );
			sqlID.testID ();
			blibFlag = sqlID.isBLib ();
			msfFlag = sqlID.isMSF ();
		}
		if ( blibFlag )
			processBLibLibrary ();
		else if ( msfFlag )
			processMSF ( false );
		else
			throw runtime_error ( "Unknown SQLite format." ); 
	}
	else if ( isMSPFile ( ff ) ) {
		processLibrary ( "msp" );
	}
	else if ( isSPTXTFile ( ff ) ) {
		processLibrary ( "sptxt" );
	}
	else {
		processResultsFile ();
	}
}
void MSViewerSearch::getMenuColumnNumbers ()
{
	string s = viewerParams.getSpectrumIdentifier ();
	if ( s == "Protein Prospector RT" )				spectrumIdentifier = SpecIDPPRT;
	else if ( s == "Scan Title (Mascot/X!Tandem)" )	spectrumIdentifier = SpecIDTitle;
	else if ( s == "Spectrum Number" )				spectrumIdentifier = SpecIDSpecNo;
	else if ( s == "m/z" )							spectrumIdentifier = SpecIDMZ;
	else if ( s == "Scan Number" )					spectrumIdentifier = SpecIDScanNo;
	else if ( s == "Scan Number - 1" )				spectrumIdentifier = SpecIDScanNoMinusOne;
	modificationReporting = viewerParams.getModificationReporting ();
	if ( modificationReporting == "Variable Mods In Peptide" || modificationReporting == "Variable Mods Column" ) {
		constMods = viewerParams.getConstMods ();
	}
	fractionColumnNumber = viewerParams.getFractionColumnNumber ();
	scanIDColumnNumber	= viewerParams.getScanIDColumnNumber ();
	peptideColumnNumber	= viewerParams.getPeptideColumnNumber ();
	peptide2ColumnNumber= 0;
	linkerColumnNumber	= 0;
	zColumnNumber		= viewerParams.getZColumnNumber ();
	constantModsColumnNumber = 0;
	variableModsColumnNumber = 0;
	allModsColumnNumber = 0;
	if ( modificationReporting == "All Mods (2 Columns)" ) {
		constantModsColumnNumber = viewerParams.getConstantModColumnNumber ();
	}
	else if ( modificationReporting == "Variable Mods Column" || modificationReporting == "All Mods (2 Columns)" ) {
		variableModsColumnNumber = viewerParams.getVariableModColumnNumber ();
	}
	else if ( modificationReporting == "All Mods (1 Column)" ) {
		allModsColumnNumber = viewerParams.getAllModColumnNumber ();
	}
}
int MSViewerSearch::getColumnNumber ( const string& id ) const
{
	const StringVector& h = headerLines [0];
	StringVectorConstIterator svi = find_if ( h.begin (), h.end (), CheckStrcasecmpEquality ( id ) );
	if ( svi != h.end () )	return svi - h.begin () + 1;
	else {
		ErrorHandler::genError ()->error ( "Can't find column header '" + id + "'." );
		return -1;
	}
}
int MSViewerSearch::getOptionalColumnNumber ( const string& id ) const
{
	const StringVector& h = headerLines [0];
	StringVectorConstIterator svi = find ( h.begin (), h.end (), id );
	if ( svi != h.end () )	return svi - h.begin () + 1;
	else					return 0;
}
void MSViewerSearch::getColumnNumbers ()
{
	string spectrumIdentifierHeader;
	string fractionHeader;
	string scanIDHeader;
	string peptideHeader;
	string chargeHeader;
	string modificationsHeader;
	viewerParams.getScriptParameters2 ( spectrumIdentifierHeader, fractionHeader, scanIDHeader, peptideHeader, chargeHeader, modificationsHeader );

	if ( spectrumIdentifierHeader == "Scan Title" )				spectrumIdentifier = SpecIDTitle;
	else if ( spectrumIdentifierHeader == "PP RT" )				spectrumIdentifier = SpecIDPPRT;
	else if ( spectrumIdentifierHeader == "Spectrum Number" )	spectrumIdentifier = SpecIDSpecNo;
	else if ( spectrumIdentifierHeader == "m/z" )				spectrumIdentifier = SpecIDMZ;
	else if ( spectrumIdentifierHeader == "Scan Number" )		spectrumIdentifier = SpecIDScanNo;
	else if ( spectrumIdentifierHeader == "Scan Number - 1" )	spectrumIdentifier = SpecIDScanNoMinusOne;
	fractionColumnNumber	= fractionHeader == "N/A" ? 0 : getColumnNumber ( fractionHeader );
	scanIDColumnNumber		= getColumnNumber ( scanIDHeader );
	peptideColumnNumber		= getColumnNumber ( peptideHeader );
	peptide2ColumnNumber	= 0;
	linkerColumnNumber		= 0;
	zColumnNumber			= getColumnNumber ( chargeHeader );
	modificationReporting	= "Variable Mods Column";
	constantModsColumnNumber= 0;
	variableModsColumnNumber= getColumnNumber ( modificationsHeader );
	allModsColumnNumber		= 0;
}
void MSViewerSearch::getProspectorColumnNumbers ()
{
	spectrumIdentifier		= SpecIDPPRT;
	fractionColumnNumber	= getOptionalColumnNumber ( "Fraction" );

	scanIDColumnNumber	= getColumnNumber ( "RT" );
	ppSpectrumColumnNumber = getOptionalColumnNumber ( "Spectrum" );
	peptideColumnNumber	= getOptionalColumnNumber ( "DB Peptide" );
	if ( !peptideColumnNumber )
		peptideColumnNumber	= getColumnNumber ( "Peptide" );
	peptide2ColumnNumber= 0;
	linkerColumnNumber	= 0;
	zColumnNumber		= getColumnNumber ( "z" );

	allModsColumnNumber	= getOptionalColumnNumber ( "Mods" );

	if ( allModsColumnNumber )
		modificationReporting = "All Mods (1 Column)";
	else {
		constantModsColumnNumber = getOptionalColumnNumber ( "Constant Mods" );
		variableModsColumnNumber = getOptionalColumnNumber ( "Variable Mods" );
		if ( constantModsColumnNumber && variableModsColumnNumber ) {
			modificationReporting = "All Mods (2 Columns)";
		}
		else {
			constMods = viewerParams.getConstMods ();
			if ( variableModsColumnNumber )
				modificationReporting = "Variable Mods Column";
			else {
				if ( !constMods.empty () )
					modificationReporting = "Variable Mods In Peptide";
				else
					modificationReporting = "All Mods In Peptide";
			}
		}
	}
}
void MSViewerSearch::getProspectorXLColumnNumbers ()
{
	spectrumIdentifier		= SpecIDPPRT;
	fractionColumnNumber	= getOptionalColumnNumber ( "Fraction" );

	scanIDColumnNumber	= getColumnNumber ( "RT" );
	ppSpectrumColumnNumber = getOptionalColumnNumber ( "Spectrum" );
	peptideColumnNumber	= getColumnNumber ( "Peptide 1" );
	peptide2ColumnNumber= getColumnNumber ( "Peptide 2" );
	linkerColumnNumber	= getOptionalColumnNumber ( "Linker" );
	zColumnNumber		= getColumnNumber ( "z" );
	constMods = viewerParams.getConstMods ();
	if ( !constMods.empty () )
		modificationReporting = "Variable Mods In Peptide";
	else
		modificationReporting = "All Mods In Peptide";
}
void MSViewerSearch::processResultsFile ()
{
	resultsFPath = viewerParams.getResultsFpath ();

	int numTitleLines;
	int numHeaderLines;
	string separator;
	bool script = false;
	if ( viewerParams.getProspector () ) {
		separator = "\t";
		numHeaderLines = 1;
	}
	else if ( viewerParams.getProspectorXL () ) {
		separator = "\t";
		numTitleLines = 0;
		numHeaderLines = 1;
	}
	else if ( viewerParams.getOther () ) {
		numTitleLines = viewerParams.getNumTitleLines ();
		numHeaderLines = viewerParams.getNumHeaderLines ();
		separator = viewerParams.getSeparator ();
	}
	else {
		viewerParams.getScriptParameters ( numTitleLines, numHeaderLines, separator );
		script = true;
	}

	string line;
	GenIFStream istr ( resultsFPath );
	if ( viewerParams.getProspector () ) {
		for ( ; ; ) {
			genUniversalGetLine ( istr, line );
			if ( line.compare ( 0, 11, "Search Name" ) ) break;
			titleLines.push_back ( line );
		}
	}
	else {
		for ( int a = 0 ; a < numTitleLines ; a++ ) {
			genUniversalGetLine ( istr, line );
			titleLines.push_back ( line );
		}
	}
	for ( int b = 0 ; b < numHeaderLines ; b++ ) {
		genUniversalGetLine ( istr, line );
		headerLines.push_back ( getColumns ( line, separator ) );
	}
	int maxCols = 0;
	while ( genUniversalGetLine ( istr, line ) ) {
		rows.push_back ( getColumns ( line, separator ) );
		maxCols = genMax ( maxCols, (int)rows.back ().size () );
	}
	istr.close ();
	extraTitleLines = 0;
	int i = 0;
	if ( numHeaderLines == 0 ) {
		for ( ; i < rows.size () ; i++ ) {
			if ( rows [i].size () != maxCols ) {
				extraTitleLines++;
			}
			else break;
		}
	}
	for ( ; i < rows.size () ; i++ ) {
		if ( rows [i].size () != maxCols ) {
			rows [i].resize ( maxCols );
		}
		if ( rows [i].size () < 3 ) {
			ErrorHandler::genError ()->error ( "Illegal results file format." );
		}
	}
	numCols = maxCols;

	setColumnNumbers ();
	filterRows ();
	sortRows ();
	removeReplicates ();
	setSpecialColumns ();
	if ( script ) setParams ( false, separator, "" );
}
void MSViewerSearch::setColumnNumbers ()
{
	if ( viewerParams.getProspector () )		getProspectorColumnNumbers ();
	else if ( viewerParams.getProspectorXL () )	getProspectorXLColumnNumbers ();
	else if ( viewerParams.getOther () )		getMenuColumnNumbers ();
	else										getColumnNumbers ();
}
void MSViewerSearch::sortRows ()
{
	IntVector sortLevel = viewerParams.getSortLevel ();
	StringVector sortOrderType = viewerParams.getSortOrderType ();
	StringVector sortOrderDirection = viewerParams.getSortOrderDirection ();
	if ( !sortLevel.empty () && !rows.empty () ) {
		try {
			stable_sort ( rows.begin (), rows.end (), SortRowsByColumn ( numCols, sortLevel, sortOrderType, sortOrderDirection ) );
		}
		catch ( runtime_error e ) {
			ErrorHandler::genError ()->error ( e );
		}
	}
}
void MSViewerSearch::filterRows ()
{
	IntVector filterColumn = viewerParams.getFilterColumn ();
	StringVector filterType = viewerParams.getFilterType ();
	StringVectorVector filterValue = viewerParams.getFilterValue ();
	if ( !filterColumn.empty () && !rows.empty () ) {
		try {
			int i = remove_if ( rows.begin (), rows.end (), FilterRowsByColumn ( numCols, filterColumn, filterType, filterValue ) ) - rows.begin ();
			rows.erase ( rows.begin () + i, rows.end () );
		}
		catch ( runtime_error e ) {
			ErrorHandler::genError ()->error ( e );
		}
	}
}
void MSViewerSearch::removeReplicates ()
{
	IntVector replicateTest = viewerParams.getReplicateTest ();
	if ( !replicateTest.empty () && !rows.empty () ) {
		rows.erase ( unique ( rows.begin (), rows.end (), TestReplicates ( replicateTest ) ), rows.end () );
	}
}
void MSViewerSearch::setSpecialColumns ()
{
	setFractionName ();
	setScanID ();
	setPeptide ();
	setPeptide2 ();
	setLinker ();
	setCharge ();
	setConstantMods ();
	setVariableMods ();
	setAllMods ();
}
void MSViewerSearch::setFractionName ()
{
	if ( fractionColumnNumber ) {
		if ( fractionColumnNumber > numCols ) {
			ErrorHandler::genError ()->error ( "Fraction column number is greater than the number of columns in the report.\n" );
		}
		for ( int i = extraTitleLines ; i < rows.size () ; i++ ) {
			fractionName.push_back ( rows [i][fractionColumnNumber-1] );
		}
	}
	else {
		if ( peakListFractionNames.empty () ) {
			for ( int i = extraTitleLines ; i < rows.size () ; i++ ) {
				fractionName.push_back ( "" );
			}
		}
		else {
			for ( int i = extraTitleLines ; i < rows.size () ; i++ ) {
				fractionName.push_back ( peakListFractionNames [0] );
			}
		}
	}
}
void MSViewerSearch::setScanID ()
{
	if ( scanIDColumnNumber ) {
		if ( scanIDColumnNumber > numCols ) {
			ErrorHandler::genError ()->error ( "scan id column number is greater than the number of columns in the report.\n" );
		}
		for ( int i = extraTitleLines ; i < rows.size () ; i++ ) {
			scanID.push_back ( rows [i][scanIDColumnNumber-1] );
			if ( ppSpectrumColumnNumber ) {
				string specNum = rows [i][ppSpectrumColumnNumber-1];
				if ( specNum != "1" ) {
					scanID.back () += "-1-" + specNum;
				}
			}
		}
	}
}
void MSViewerSearch::setPeptide ()
{
	if ( peptideColumnNumber ) {
		if ( peptideColumnNumber > numCols ) {
			ErrorHandler::genError ()->error ( "Peptide column number is greater than the number of columns in the report.\n" );
		}
		for ( int i = extraTitleLines ; i < rows.size () ; i++ ) {
			peptide.push_back ( rows [i][peptideColumnNumber-1] );
		}
	}
}
void MSViewerSearch::setPeptide2 ()
{
	if ( peptide2ColumnNumber ) {
		if ( peptide2ColumnNumber > numCols ) {
			ErrorHandler::genError ()->error ( "Peptide 2 column number is greater than the number of columns in the report.\n" );
		}
		for ( int i = extraTitleLines ; i < rows.size () ; i++ ) {
			peptide2.push_back ( rows [i][peptide2ColumnNumber-1] );
		}
	}
}
void MSViewerSearch::setLinker ()
{
	if ( linkerColumnNumber ) {
		if ( linkerColumnNumber > numCols ) {
			ErrorHandler::genError ()->error ( "Linker column number is greater than the number of columns in the report.\n" );
		}
		for ( int i = extraTitleLines ; i < rows.size () ; i++ ) {
			linker.push_back ( rows [i][linkerColumnNumber-1] );
		}
	}
}
void MSViewerSearch::setCharge ()
{
	if ( zColumnNumber ) {
		if ( zColumnNumber > numCols ) {
			ErrorHandler::genError ()->error ( "Charge column number is greater than the number of columns in the report.\n" );
		}
		for ( int i = extraTitleLines ; i < rows.size () ; i++ ) {
			charge.push_back ( atoi ( rows [i][zColumnNumber-1].c_str () ) );
		}
	}
}
void MSViewerSearch::setConstantMods ()
{
	if ( modificationReporting == "All Mods (2 Columns)" ) {
		if ( constantModsColumnNumber ) {
			if ( constantModsColumnNumber > numCols ) {
				ErrorHandler::genError ()->error ( "Constant Mods column number is greater than the number of columns in the report.\n" );
			}
			for ( int i = extraTitleLines ; i < rows.size () ; i++ ) {
				constantMods.push_back ( rows [i][constantModsColumnNumber-1] );
			}
		}
	}
}
void MSViewerSearch::setVariableMods ()
{
	if ( modificationReporting == "Variable Mods Column" || modificationReporting == "All Mods (2 Columns)" ) {
		if ( variableModsColumnNumber ) {
			if ( variableModsColumnNumber > numCols ) {
				ErrorHandler::genError ()->error ( "Variable Mods column number is greater than the number of columns in the report.\n" );
			}
			for ( int i = extraTitleLines ; i < rows.size () ; i++ ) {
				variableMods.push_back ( rows [i][variableModsColumnNumber-1] );
			}
		}
	}
}
void MSViewerSearch::setAllMods ()
{
	if ( modificationReporting == "All Mods (1 Column)" ) {
		if ( allModsColumnNumber ) {
			if ( allModsColumnNumber > numCols ) {
				ErrorHandler::genError ()->error ( "Mods column number is greater than the number of columns in the report.\n" );
			}
			for ( int i = extraTitleLines ; i < rows.size () ; i++ ) {
				allMods.push_back ( rows [i][allModsColumnNumber-1] );
			}
		}
	}
}
void MSViewerSearch::processPepXMLResultsFile ( bool multi )
{
	resultsFPath = viewerParams.getResultsFpath ();
	string dir;
	string shortFName;
	string fracName;
	if ( multi ) {
		dir = resultsFPath;									// A directory for the peak lists
		shortFName = viewerParams.getResName ();			// The name of the upload file. This is where the peaks will go.
		FileList fList ( resultsFPath, "", "", false );
		StringVector n1 = fList.getNameList ();
		for ( StringVectorSizeType i = 0, j = 0 ; i < n1.size () ; i++ ) {
			string n = n1 [i];
			if ( n == "." || n == ".." || n == "params.xml" ) continue;
			string f = resultsFPath + SLASH + n;		// A pepXML file
			if ( genIsDirectory ( f ) ) continue;
			{
				if ( !peakListFractionNames.empty () ) {
					if ( isFileType ( n, PEPXML ) ) fracName = n.substr ( 0, n.length () - 8 );
					else							fracName = genShortFilenameFromPath ( n );
				}
				StringVectorVector hL;
				PPExpatPepXML ppepx ( fracName, hL, rows, zColumnNumber, fractionColumnNumber, scanIDColumnNumber, peptideColumnNumber, variableModsColumnNumber );
				ppepx.parseXMLFromFile ( f );
				if ( headerLines.empty () ) {
					headerLines = hL;
				}
				else {
					if ( headerLines != hL ) {
						ErrorHandler::genError ()->message ( "The header lines from the different msf files don't match.\n" );
					}
				}
			}
			genUnlink ( f );
		}
	}
	else {
		dir = genDirectoryFromPath ( resultsFPath );
		shortFName = genShortFilenameFromPath ( resultsFPath );

		if ( !peakListFractionNames.empty () ) fracName = peakListFractionNames [0];
		PPExpatPepXML ppepx ( fracName, headerLines, rows, zColumnNumber, fractionColumnNumber, scanIDColumnNumber, peptideColumnNumber, variableModsColumnNumber );
		ppepx.parseXMLFromFile ( resultsFPath );
	}
	for ( int a = 0 ; a < rows.size () ; a++ ) {
		fractionName.push_back ( rows [a][fractionColumnNumber-1] );
	}
	string resFile = dir + SLASH + shortFName + ".txt";
	writeResultsFile ( resFile );
	if ( multi ) {
		const_cast <MSViewerParameters*> (&viewerParams)->setResultsFpath ( resFile );
	}
	else {
		genUnlink ( resultsFPath );
	}
	spectrumIdentifier = SpecIDScanNo;
	setParams ( false, "\t", ".txt" );

	if ( !rows.empty () ) {
		numCols = rows [0].size ();
	}
	extraTitleLines = 0;
	setColumnNumbers ();
	sortRows ();
	setSpecialColumns ();
}
void MSViewerSearch::processPrideXMLResultsFile ()
{
	resultsFPath = viewerParams.getResultsFpath ();
	fractionColumnNumber = 0;
	MSMSDataPointVector msmsDataPointList;
	SpecID specID;
	specID.setFraction ( -1 );

	StringVector constantHeader;
	SetString variableHeader;
	VectorMapStringToString rowsVariable;
	PPExpatPrideXML ppepx ( constantHeader, variableHeader, rows, rowsVariable, msmsDataPointList, -1, SpectrumRange (), specID, -9999, IntVector () );
	ppepx.parseXMLFromFile ( resultsFPath );
	makeHeaderLines ( constantHeader, variableHeader );
	makeRowLines ( variableHeader, rowsVariable );

	IntVector scans;
	int numRows = rows.size ();
	peptide2ColumnNumber = 0;
	linkerColumnNumber = 0;
	scanIDColumnNumber = 4;
	peptideColumnNumber = 5;
	variableModsColumnNumber = 6;
	zColumnNumber = 8;
	for ( int a = 0 ; a < numRows ; a++ ) {
		scans.push_back ( atoi ( rows [a][scanIDColumnNumber-1].c_str () ) );
		scanID.push_back ( rows [a][scanIDColumnNumber-1] );
		peptide.push_back ( rows [a][peptideColumnNumber-1] );
		variableMods.push_back ( rows [a][variableModsColumnNumber-1] );
		charge.push_back ( atoi ( rows [a][zColumnNumber-1].c_str () ) );
	}
	stable_sort ( scans.begin (), scans.end () );
	scans.erase ( unique ( scans.begin (), scans.end () ), scans.end () );

	string dir = genDirectoryFromPath ( resultsFPath );
	string shortFName = genShortFilenameFromPath ( resultsFPath );

	string resFile = dir + SLASH + shortFName + ".txt";
	writeResultsFile ( resFile );
	genUnlink ( resultsFPath );

	writeMGFFile ( dir, shortFName, msmsDataPointList, scans );
	peakListFractionNames.push_back ( shortFName );
	peakListCentroidFiles.push_back ( shortFName + ".mgf" );

	spectrumIdentifier = SpecIDPPRT;
	setParams ( true, "\t", ".txt" );

	peakListFPath = viewerParams.getPeakListFpath ();
	if ( !peakListFractionNames.empty () ) {
		fractionName.resize ( scanID.size () );
		fill ( fractionName.begin (), fractionName.end (), peakListFractionNames [0] );
	}
	if ( !rows.empty () ) {
		numCols = rows [0].size ();
	}
	extraTitleLines = 0;
}
void MSViewerSearch::setParams ( bool peakList, const string& separator, const string& suffix )
{
	resultsFPath = viewerParams.getResultsFpath ();
	string dir = genDirectoryFromPath ( resultsFPath );
	string shortFName = genShortFilenameFromPath ( resultsFPath );
	string newSuffix = suffix;
	if ( newSuffix.empty () )	newSuffix = "." + genSuffixFromPath ( resultsFPath );
	string resFile = dir + SLASH + shortFName + newSuffix;

	string newSeparator = separator;
	if ( separator == "\t" )	newSeparator = "Tab Delimited";
	else						newSeparator = "CSV";
	string parFile = dir + SLASH + "params.xml";
	ParameterList pList ( parFile, false, false, false, false );

	string spectrumIDString = getSpectrumIdentifierTypeString ( spectrumIdentifier );

	pList.addOrReplaceName ( "results_filepath", resFile );
	if ( peakList ) pList.addOrReplaceName ( "peak_list_filepath", dir + SLASH + shortFName );
	pList.addOrReplaceName ( "results_file_format", "Other" );
	pList.addOrReplaceName ( "spectrum_identifier", spectrumIDString );
	pList.addOrReplaceName ( "column_num_scan_id", scanIDColumnNumber );
	pList.addOrReplaceName ( "column_num_peptide", peptideColumnNumber );
	pList.addOrReplaceName ( "column_num_variable_mod", variableModsColumnNumber );
	pList.addOrReplaceName ( "column_num_z", zColumnNumber );
	pList.addOrReplaceName ( "column_num_fraction", fractionColumnNumber );
	pList.addOrReplaceName ( "column_separator", newSeparator );
	pList.addOrReplaceName ( "num_header_lines", 1 );
	pList.XMLParameterFile ( parFile );

	const_cast <MSViewerParameters*> (&viewerParams)->setResultsFpath ( resFile );
	if ( peakList ) const_cast <MSViewerParameters*> (&viewerParams)->setPeakListFpath ( dir + SLASH + shortFName );
	const_cast <MSViewerParameters*> (&viewerParams)->setResultsFileFormat ( "Other" );
	const_cast <MSViewerParameters*> (&viewerParams)->setSpectrumIdentifier ( spectrumIDString );
	const_cast <MSViewerParameters*> (&viewerParams)->setScanIDColumnNumber ( scanIDColumnNumber );
	const_cast <MSViewerParameters*> (&viewerParams)->setPeptideColumnNumber ( peptideColumnNumber );
	const_cast <MSViewerParameters*> (&viewerParams)->setVariableModColumnNumber ( variableModsColumnNumber );
	const_cast <MSViewerParameters*> (&viewerParams)->setZColumnNumber ( zColumnNumber );
	const_cast <MSViewerParameters*> (&viewerParams)->setFractionColumnNumber ( fractionColumnNumber );
	const_cast <MSViewerParameters*> (&viewerParams)->setSeparator ( newSeparator );
	const_cast <MSViewerParameters*> (&viewerParams)->setNumHeaderLines ( 1 );
}
void MSViewerSearch::processBLibLibrary ()
{
	resultsFPath = viewerParams.getResultsFpath ();
	string dir = genDirectoryFromPath ( resultsFPath );
	string shortFName = genShortFilenameFromPath ( resultsFPath );
	{
		BlibRead br ( resultsFPath );
		br.read ( dir + SLASH + shortFName, headerLines, rows );
		peakListFractionNames = br.getPeakListFractionNames ();
		peakListCentroidFiles = br.getPeakListCentroidFileNames ();
	}
	peptide2ColumnNumber = 0;
	linkerColumnNumber	= 0;
	peptideColumnNumber = 1;
	variableModsColumnNumber = 2;
	zColumnNumber = 4;
	fractionColumnNumber = 10;
	scanIDColumnNumber = 11;

	for ( int a = 0 ; a < rows.size () ; a++ ) {
		scanID.push_back ( rows [a][scanIDColumnNumber-1] );
		peptide.push_back ( rows [a][peptideColumnNumber-1] );
		variableMods.push_back ( rows [a][variableModsColumnNumber-1] );
		charge.push_back ( atoi ( rows [a][zColumnNumber-1].c_str () ) );
		fractionName.push_back ( rows [a][fractionColumnNumber-1] );
	}

	string resFile = dir + SLASH + shortFName + ".txt";
	writeResultsFile ( resFile );
	genUnlink ( resultsFPath );

	spectrumIdentifier = SpecIDScanNo;
	setParams ( true, "\t", ".txt" );

	peakListFPath = viewerParams.getPeakListFpath ();

	if ( !rows.empty () ) {
		numCols = rows [0].size ();
	}
	extraTitleLines = 0;
}
void MSViewerSearch::processLibrary ( const string& type )
{
	resultsFPath = viewerParams.getResultsFpath ();
	string dir = genDirectoryFromPath ( resultsFPath );
	string shortFName = genShortFilenameFromPath ( resultsFPath );
	if ( type == "msp" ) {
		MSPRead msp ( resultsFPath );
		msp.read ( dir + SLASH + shortFName, headerLines, rows );
		peakListFractionNames = msp.getPeakListFractionNames ();
		peakListCentroidFiles = msp.getPeakListCentroidFileNames ();
	}
	else if ( type == "sptxt" ) {
		SPTXTRead sptxt ( resultsFPath );
		sptxt.read ( dir + SLASH + shortFName, headerLines, rows );
		peakListFractionNames = sptxt.getPeakListFractionNames ();
		peakListCentroidFiles = sptxt.getPeakListCentroidFileNames ();
	}
	peptide2ColumnNumber = 0;
	linkerColumnNumber	= 0;
	peptideColumnNumber = 2;
	variableModsColumnNumber = 4;
	zColumnNumber = 6;
	scanIDColumnNumber = 7;
	fractionColumnNumber = 0;

	for ( int a = 0 ; a < rows.size () ; a++ ) {
		scanID.push_back ( rows [a][scanIDColumnNumber-1] );
		peptide.push_back ( rows [a][peptideColumnNumber-1] );
		variableMods.push_back ( rows [a][variableModsColumnNumber-1] );
		charge.push_back ( atoi ( rows [a][zColumnNumber-1].c_str () ) );
		//fractionName.push_back ( rows [a][fractionColumnNumber-1] );
	}

	string resFile = dir + SLASH + shortFName + ".txt";
	writeResultsFile ( resFile );
	genUnlink ( resultsFPath );

	spectrumIdentifier = SpecIDPPRT;
	setParams ( true, "\t", ".txt" );

	peakListFPath = viewerParams.getPeakListFpath ();
	if ( !peakListFractionNames.empty () ) {
		fractionName.resize ( scanID.size () );
		fill ( fractionName.begin (), fractionName.end (), peakListFractionNames [0] );
	}
	if ( !rows.empty () ) {
		numCols = rows [0].size ();
	}
	extraTitleLines = 0;
}
void MSViewerSearch::processMSF ( bool multi )
{
	resultsFPath = viewerParams.getResultsFpath ();		// The full path of the directory containing the msf files (multi) or the file (single)
	MapCharToString mcrs = viewerParams.getConstMods ();// What's on the menu
	MSFRead::setConstMods ( mcrs );						// Sets the const mods from what's on the menu
	string dir;
	string shortFName;
	if ( multi ) {
		dir = resultsFPath;									// A directory for the peak lists
		shortFName = viewerParams.getResName ();			// The name of the upload file. This is where the peaks will go.
		FileList fList ( resultsFPath, "", "", false );
		StringVector n1 = fList.getNameList ();
		for ( StringVectorSizeType i = 0 ; i < n1.size () ; i++ ) {
			string n = n1 [i];
			if ( n != "." && n != ".." && n != "params.xml" ) {
				string f = resultsFPath + SLASH + n;		// An msf file
				{
					MSFRead msf ( f );
					msf.read ( dir + SLASH + shortFName, rows );
					StringVectorVector hL;
					msf.readHeader ( hL );					// Header has to be read after as the score columns are extracted from the database.
					if ( headerLines.empty () ) {
						headerLines = hL;
					}
					else {
						if ( headerLines != hL ) {
							ErrorHandler::genError ()->message ( "The header lines from the different msf files don't match.\n" );
						}
					}
					StringVector plfn = msf.getPeakListFractionNames ();
					peakListFractionNames.insert ( peakListFractionNames.end (), plfn.begin (), plfn.end () );
					StringVector plcf = msf.getPeakListCentroidFileNames ();
					peakListCentroidFiles.insert ( peakListCentroidFiles.end (), plcf.begin (), plcf.end () );
				}
				genUnlink ( f );
			}
		}
	}
	else {
		dir = genDirectoryFromPath ( resultsFPath );					// The directory from the above path
		shortFName = genShortFilenameFromPath ( resultsFPath );			// The msf file name without the msf suffix. This is the directory where the peak lists will go.
		{
			MSFRead msf ( resultsFPath );
			msf.read ( dir + SLASH + shortFName, rows );
			msf.readHeader ( headerLines );					// Header has to be read after as the score columns are extracted from the database.
			peakListFractionNames = msf.getPeakListFractionNames ();
			peakListCentroidFiles = msf.getPeakListCentroidFileNames ();
		}
	}
	peptide2ColumnNumber = 0;
	linkerColumnNumber	= 0;
	zColumnNumber = 2;
	peptideColumnNumber = 3;
	variableModsColumnNumber = 4;
	fractionColumnNumber = 5;
	scanIDColumnNumber = 6;
	constMods = MSFRead::getConstModMap ();
	if ( mcrs.empty () ) const_cast <MSViewerParameters*> (&viewerParams)->setConstMod ( MSFRead::getMsfConstModsStr () );	// If the menu was empty set it to what's in the msf file
	for ( int a = 0 ; a < rows.size () ; a++ ) {
		scanID.push_back ( rows [a][scanIDColumnNumber-1] );
		peptide.push_back ( rows [a][peptideColumnNumber-1] );
		variableMods.push_back ( rows [a][variableModsColumnNumber-1] );
		charge.push_back ( atoi ( rows [a][zColumnNumber-1].c_str () ) );
		fractionName.push_back ( rows [a][fractionColumnNumber-1] );
	}
	string resFile = dir + SLASH + shortFName + ".txt";
	writeResultsFile ( resFile );
	if ( multi ) {
		const_cast <MSViewerParameters*> (&viewerParams)->setResultsFpath ( resFile );
	}
	else {
		genUnlink ( resultsFPath );
	}
	spectrumIdentifier = SpecIDPPRT;
	setParams ( true, "\t", ".txt" );
	peakListFPath = viewerParams.getPeakListFpath ();
	if ( !rows.empty () ) {
		numCols = rows [0].size ();
	}
	extraTitleLines = 0;
}
void MSViewerSearch::makeHeaderLines ( const StringVector& constantHeader, const SetString& variableHeader )
{
	headerLines.push_back ( constantHeader );
	for ( SetStringConstIterator i = variableHeader.begin () ; i != variableHeader.end () ; i++ ) {
		headerLines [0].push_back ( *(i) );
	}
}
void MSViewerSearch::makeRowLines ( const SetString& variableHeader, const VectorMapStringToString rowsVariable )
{
	for ( StringVectorVectorSizeType i = 0 ; i < rows.size () ; i++ ) {
		for ( SetStringConstIterator j = variableHeader.begin () ; j != variableHeader.end () ; j++ ) {
			MapStringToStringConstIterator cur = rowsVariable [i].find ( *(j) );
			if ( cur != rowsVariable [i].end () ) {
				rows [i].push_back ( (*cur).second );
			}
			else {
				rows [i].push_back ( "" );
			}
		}
	}
}
void MSViewerSearch::writeResultsFile ( const string& resFile ) const
{
	GenOFStream ofs ( resFile );
	for ( int i = 0 ; i < headerLines.size () ; i++ ) {
		int siz = headerLines [i].size ();
		for ( int j = 0 ; j < siz ; j++ ) {
			ofs << headerLines [i][j];
			if ( j != siz - 1 ) ofs << '\t';
		}
		ofs << '\n';
	}
	int numRows = rows.size ();
	for ( int x = 0 ; x < numRows ; x++ ) {
		int siz = rows [x].size ();
		for ( int y = 0 ; y < siz ; y++ ) {
			ofs << rows [x][y];
			if ( y != siz - 1 ) ofs << '\t';
		}
		ofs << '\n';
	}
	ofs.close ();
}
void MSViewerSearch::writeMGFFile ( const string& dir, const string& shortFName, MSMSDataPointVector& msmsDataPointList, const IntVector& scans ) const
{
	genCreateDirectory ( dir + SLASH + shortFName );
	string mgfFileBase = dir + SLASH + shortFName + SLASH + shortFName;
	string mgfFile = mgfFileBase + ".mgf";
	GenOFStream osf ( mgfFile );
	int m = 0;
	for ( MSMSDataPointVectorSizeType i = 0 ; i < msmsDataPointList.size () ; i++ ) {
		const MSMSDataPoint& msms = msmsDataPointList [i];
		if ( atoi ( msms.getMSMSInfo ().c_str () ) == scans [m] ) {
			m++;
			int z = msms.getPrecursorCharge ();
			PPProject::writeMGFScan ( osf, msmsDataPointList [i], z );
		}
	}
	osf.close ();
}
void MSViewerSearch::processPeakListFile ()
{
	peakListFPath = viewerParams.getPeakListFpath ();
	if ( peakListFPath.empty () ) return;
	PPProject ppp ( peakListFPath );
	if ( ppp.initialised () ) {
		peakListFractionNames = ppp.getFractionNames ();
		peakListCentroidFiles = ppp.getCentroidFiles ();
	}
	else {
		string err = ppp.getErrMessage ();
		genUnlinkDirectory ( peakListFPath );
		if ( err.empty () )
			throw runtime_error ( "Upload file has an invalid format." );
		else
			throw runtime_error ( err );
	}
}
StringVector MSViewerSearch::getColumns ( const string& line, const string& separator )
{
	StringVector cols;
	string::size_type start = 0;
	string::size_type end = 0;
	bool append = false;
	string s2;
	for ( ; ; ) {
		end = line.find_first_of ( separator, start );
		string s = line.substr ( start, end-start );
		start = end + 1;
		if ( end != string::npos || !s.empty () ) {
			if ( s [0] == '"' ) {
				if ( end == string::npos ) { // deals with comma at end of field in inverted commas at line end
					cols.push_back ( gen_strtrim2 ( s2 + "," ) );
				}
				else {
					s2 = s.substr ( 1 );
					append = true;
					if ( !s2.empty () && s2 [s2.length () - 1] == '"' ) {
						cols.push_back ( gen_strtrim2 ( s2.substr ( 0, s2.length () - 1 ) ) );
						append = false;
					}
				}
			}
			else if ( append ) {
				s2 += ',' + s;
				if ( s2 [s2.length () - 1] == '"' ) {
					cols.push_back ( gen_strtrim2 ( s2.substr ( 0, s2.length () - 1 ) ) );
					append = false;
				}
			}
			else
				cols.push_back ( gen_strtrim2 ( s ) );
		}
		if ( end == string::npos ) break;
	}
	return cols;
}
void MSViewerSearch::printBodyHTML ( ostream& os )
{
	if ( viewerParams.getSaveSettings () ) {
		printSaveSettingsHTML ( os );
	}
	else if ( viewerParams.getDeleteFlag () ) {
		printMSViewerRepositoryHMTL ( os );
	}
	else if ( viewerParams.getTabDelimitedText () ) {
		printMSViewerReportTabDelimited ( os, false, true );
	}
	else if ( viewerParams.getTabDelimitedTextWithURL () ) {
		printMSViewerReportTabDelimited ( os, true, true );
	}
	else if ( viewerParams.getFilteredViewerFiles () ) {
		printMSViewerReportTabDelimited ( os, false, false );
		printMSViewerFilteredPeakList ( os );
	}
	else if ( viewerParams.getViewerFiles () ) {
		printMSViewerReportViewerFiles ( os );
	}
	else {
		printMSViewerReportHMTL ( os );
	}
}
void MSViewerSearch::printSaveSettingsHTML ( ostream& os )
{
	string sKey = viewerParams.getKey ();
	os << "<p>" << endl;
	os << "The search key for the saved data set is <b>" << sKey << "</b><br />" << endl;
	os << "</p>" << endl;

	os << "<p>" << endl;
	os << "The data set can be accessed using the following URL:<br />" << endl;
	os << "</p>" << endl;
	ParameterList pList ( "mssearch", false, false, false, false, false );
	pList.addName ( "search_name", "msviewer" );
	pList.addName ( "report_title", "MS-Viewer" );
	pList.addName ( "search_key", sKey );
	string url = ParameterList::getServer () + pList.getURL ();
	os << "<p><a href=\"" << url << "\">" << url << "</a></p>" << endl;
	os << "<p>Click on the URL to display the report.</p>" << endl;
}
void MSViewerSearch::printMSViewerRepositoryHMTL ( ostream& os )
{
	string viewerDir = MSViewerParameters::getViewerRepositoryContainerPath ();
	StringVector sv = getRepositoryKeyList ( viewerDir, 10 );
	vector <ViewerDataSet> vvDataSet;
	for ( int i = 0 ; i < sv.size () ; i++ ) {
		string key = sv [i];
		string fname = viewerDir;
		fname += key [0];
		fname += SLASH;
		fname += key [1];
		fname += SLASH;
		fname += key;
		fname += SLASH;
		fname += "params.xml";
		ParameterList plist ( fname, false, false, false, false );
		vvDataSet.push_back ( ViewerDataSet ( key, plist.getStringValue ( "peak_list_filepath" ), plist.getStringValue ( "results_filepath" ), GenTime ( getXMLFileDate ( fname ) ) ) );
	}
	stable_sort ( vvDataSet.begin (), vvDataSet.end (), sortViewerDataSet () );
	startJavascript ( os );
	ViewerMSViewerLink vl;
	vl.printHTML ( os );
	endJavascript ( os );
	tableStart ( os, true );
		for ( int j = 0 ; j < vvDataSet.size () ; j++ ) {
			if	( j == 0 ) {
				tableRowStart ( os );
					tableHeader ( os, "Search Key" );
					tableHeader ( os, "Data Set" );
					tableHeader ( os, "Result File" );
					tableHeader ( os, "Creation Date" );
					tableHeader ( os, "Delete" );
				tableRowEnd ( os );
			}
			tableRowStart ( os );
				tableDataStart ( os );
					vl.write ( os, vvDataSet [j].searchKey );
				tableDataEnd ( os );
				tableCell ( os, vvDataSet [j].dataSet );
				tableCell ( os, vvDataSet [j].results );
				tableCell ( os, vvDataSet [j].date.getDateAndTime () );
				tableDataStart ( os );
					vl.write2 ( os, vvDataSet [j].searchKey );
				tableDataEnd ( os );
			tableRowEnd ( os );
		}
	tableEnd ( os );
}
void MSViewerSearch::printMSViewerReportHMTL ( ostream& os )
{
	int numRows = rows.size ();
	map <string, MSProductLink*> productLink;
	string linkSearchType = viewerParams.getLinkSearchType ();
	string bridgeFormula = viewerParams.getBridgeFormula ();
	if ( numRows ) {
		startJavascript ( os );
		for ( int b = 0 ; b < peakListCentroidFiles.size () ; b++ ) {
			string fname = peakListFPath + SLASH + peakListCentroidFiles [b];
			string fracName = peakListFractionNames [b];
			productLink [fracName] = new MSProductLink ( viewerParams.getInstrumentName (), fname, viewerParams.getParentTolerance (), viewerParams.getFragmentTolerance () );
			productLink [fracName]->printHTML ( os );
		}
		endJavascript ( os );
	}
	bool open = false;

	int numDataLines = numRows - extraTitleLines;
	int pageNumber = viewerParams.getPage ();
	string sRowsPerPage = viewerParams.getRowsPerPage ();
	int rowsPerPage;
	if ( sRowsPerPage == "All" ) {
		pageNumber = 1;
		const_cast <MSViewerParameters*> (&viewerParams)->setPage ( 1 );
		rowsPerPage = numDataLines;
	}
	else {
		rowsPerPage = atoi ( sRowsPerPage.c_str () );
		int numPages = numDataLines / rowsPerPage;
		if ( numDataLines % rowsPerPage ) numPages++;
		if ( pageNumber > numPages ) {
			const_cast <MSViewerParameters*> (&viewerParams)->setPage ( numPages );
			pageNumber = numPages;
		}
	}

	bool pFlag = viewerParams.getPListFlag ();
	MSViewerForm* msfv;
	ExpandableJavascriptBlock ejb ( "Display and Parameter Settings", open );
	ejb.printHeader ( os );
	const ParameterList* pList = viewerParams.getParameterList ();
	VectorConstParameterListPtr vcplp;
	vcplp.push_back ( pList );
	int nc = rows.empty () ? 0 : rows [0].size ();
	msfv = new MSViewerForm ( vcplp, pFlag, nc );
	bool saveOption = pFlag && !fractionName.empty () && !scanID.empty () && !peptide.empty () && !productLink.empty () && !charge.empty ();
	msfv->printHTMLForm ( os, saveOption, false );
	ejb.printFooter ( os );
	if ( numRows == 0 ) {
		os << "<p>" << endl;
		ErrorHandler::genError ()->message ( "No data to report.\n" );
		os << "</p>" << endl;
		return;
	}
	else {
		os << "<p>" << endl;
			os << "<b>Report contains " << numRows << " rows.</b>" << "<br />" << endl;
		os << "</p>" << endl;
	}

	if ( sRowsPerPage != "All" ) {
		string searchKey = viewerParams.getSearchKey ();
		if ( numDataLines > rowsPerPage && !searchKey.empty () ) {
			startJavascript ( os );
			ReportLinkProgram rlp ( "mssearch", "viewerLink2" );
			rlp.printHTML ( os, pList );
			endJavascript ( os );

			ReportLinks rLinks ( "", rowsPerPage, pageNumber, numDataLines );
			rLinks.printHTML ( os, rlp );
		}
	}
	int startRow = (pageNumber-1) * rowsPerPage;

	map <string, ViewerMSTagLink*> tagLink;
	if ( numRows && !fractionName.empty () && !scanID.empty () && !productLink.empty () && !charge.empty () ) {
		for ( int i = 0 ; i < peakListCentroidFiles.size () ; i++ ) {
			tagLink [peakListFractionNames [i]] = new ViewerMSTagLink ( peakListFPath + SLASH + peakListCentroidFiles [i], msfv );
			tagLink [peakListFractionNames [i]]->putDataFilename ( os );
		}
		for ( int j = extraTitleLines ; j < numRows ; j++ ) {
			int m = j - extraTitleLines;
			if ( m < startRow ) continue;
			if ( m >= startRow + rowsPerPage ) break;
			if ( !peptide [m].empty () ) {
				string tLinkName = "tag_link_" + gen_itoa ( m );
				if ( tagLink.find ( fractionName [m] ) != tagLink.end () ) {
					if ( spectrumIdentifier == SpecIDPPRT ) {
						tagLink [fractionName [m]]->putForm1 ( os, tLinkName, SpecID ( getSpecIDStringFromScanID ( scanID [m] ) ) );
					}
					else if ( spectrumIdentifier == SpecIDTitle )
						tagLink [fractionName [m]]->putForm2 ( os, tLinkName, scanID [m] );
					else if ( spectrumIdentifier == SpecIDSpecNo )
						tagLink [fractionName [m]]->putForm3 ( os, tLinkName, scanID [m] );
					else if ( spectrumIdentifier == SpecIDMZ )
						tagLink [fractionName [m]]->putForm4 ( os, tLinkName, scanID [m] );
					else if ( spectrumIdentifier == SpecIDScanNo )
						tagLink [fractionName [m]]->putForm5 ( os, tLinkName, scanID [m] );
					else if ( spectrumIdentifier == SpecIDScanNoMinusOne ) {
						int sc = atoi ( scanID [m].c_str () ) + 1;
						tagLink [fractionName [m]]->putForm5 ( os, tLinkName, gen_itoa ( sc ) );
					}
				}
			}
		}
	}
	if ( !titleLines.empty () ) os << "<p>" << endl;
	for ( int a = 0 ; a < titleLines.size () ; a++ ) {
		os << "<b>" << titleLines [a] << "</b>" << "<br />" << endl;
	}
	if ( !titleLines.empty () ) os << "</p>" << endl;

	BoolDeque removeColumn = viewerParams.getRemoveColumn ( numCols );
	tableStart ( os, true );
		tableRowStart ( os );
			for ( int i = 1 ; i <= numCols ; i++ ) {
				if ( !removeColumn [i-1] ) {
					if	( i != peptide2ColumnNumber ) {
						tableHeader ( os, gen_itoa ( i ) );
					}
				}
			}
		tableRowEnd ( os );
		if ( !headerLines.empty () ) {
			for ( int i = 0 ; i < headerLines.size () ; i++ ) {
				tableRowStart ( os );
					for ( int j = 0 ; j < headerLines [i].size () ; j++ ) {
						if ( !removeColumn [j] ) {
							if ( peptide2ColumnNumber ) {
								if ( j+1 == peptideColumnNumber ) tableHeader ( os, "Crosslinked Peptide" );
								else if ( j+1 != peptide2ColumnNumber ) tableHeader ( os, headerLines [i][j] );
							}
							else tableHeader ( os, headerLines [i][j] );
						}
					}
				tableRowEnd ( os );
			}
		}
		for ( int j = extraTitleLines ; j < numRows ; j++ ) {
			int m = j - extraTitleLines;
			if ( m < startRow ) continue;
			if ( m >= startRow + rowsPerPage ) break;
			tableRowStart ( os );
				for ( int k = 1 ; k <= numCols ; k++ ) {
					if ( !removeColumn [k-1] ) {
						if ( scanIDColumnNumber == k ) {
							if ( fractionName.empty () || scanID.empty () || peptide [m].empty () || productLink.empty () || charge.empty ()  )
								tableCell ( os, scanID [m], true );
							else {
								tableCellStart ( os, "", "", true );
									string s1 = scanID [m];
									string::size_type ss = s1.find ( "-" );
									if ( ss != string::npos ) {
										s1 = s1.substr ( 0, ss );
									}
									if ( tagLink.find ( fractionName [m] ) != tagLink.end () ) {
										ViewerMSTagLink::putFormLink ( os, "tag_link_" + gen_itoa ( m ), s1 );
									}
									else
										os << s1;
								tableCellEnd ( os );
							}
						}
						if	( peptideColumnNumber == k ) {
							if ( fractionName.empty () || scanID.empty () || peptide [m].empty () || productLink.empty () || charge.empty ()  )
								tableCell ( os, peptide [m], true );
							else {
								string mod;
								if ( allMods.empty () ) {
									string cMod = constantMods.empty () ? "" : constantMods [m];
									if ( cMod.empty () && !constMods.empty () ) cMod = getConstModsString ( peptide [m] );
									string vMod = variableMods.empty () ? "" : variableMods [m];
									mod = cMod.empty () ? vMod : ( vMod.empty () ? cMod : cMod + ';' + vMod );
								}
								else {
									mod = allMods [m];
								}
								tableCellStart ( os, "", "", true );
									if ( productLink.find ( fractionName [m] ) != productLink.end () ) {
										if ( spectrumIdentifier == SpecIDPPRT ) {
											SpecID specID ( getSpecIDStringFromScanID ( scanID [m] ) );
											if ( !peptide2.empty () ) {
												string lst = linkerColumnNumber ? linker [m] : linkSearchType;
												productLink [fractionName [m]]->write5 ( os, specID, peptide [m], peptide2 [m], charge [m], false, lst, bridgeFormula );
											}
											else
												productLink [fractionName [m]]->write5 ( os, specID, peptide [m], mod, charge [m], false );
										}
										else if ( spectrumIdentifier == SpecIDTitle )
											productLink [fractionName [m]]->write5 ( os, scanID [m], peptide [m], mod, charge [m], false );
										else if ( spectrumIdentifier == SpecIDSpecNo )
											productLink [fractionName [m]]->write6 ( os, scanID [m], peptide [m], mod, charge [m], false );
										else if ( spectrumIdentifier == SpecIDMZ )
											productLink [fractionName [m]]->write7 ( os, scanID [m], peptide [m], mod, charge [m], false );
										else if ( spectrumIdentifier == SpecIDScanNo )
											productLink [fractionName [m]]->write8 ( os, scanID [m], peptide [m], mod, charge [m], false );
										else if ( spectrumIdentifier == SpecIDScanNoMinusOne ) {
											int sc = atoi ( scanID [m].c_str () ) + 1;
											productLink [fractionName [m]]->write8 ( os, gen_itoa ( sc ), peptide [m], mod, charge [m], false );
										}
									}
									else
										os << peptide [m];
								tableCellEnd ( os );
							}
						}
						if	( peptideColumnNumber != k && peptide2ColumnNumber != k && scanIDColumnNumber != k )
							tableCell ( os, rows [j][k-1], true );
					}
				}
			tableRowEnd ( os );
		}
	tableEnd ( os );
	delete msfv;
}
void MSViewerSearch::printMSViewerReportTabDelimited ( ostream& os, bool url, bool printNumbers )
{
	PPTempFile pptf ( "", "" );
	string fullPath = pptf.getFullPath ();
	string filename = "report.txt";
	string actualPath = fullPath + SLASH + filename;
	string outputPath = pptf.getURL ();
	genCreateDirectory ( fullPath );
	GenOFStream ost ( actualPath, std::ios_base::out );

	if ( url ) {
		printMSViewerReportTabDelimitedTextWithURL ( ost, printNumbers );
	}
	else {
		printMSViewerReportTabDelimitedText ( ost, printNumbers );
	}

	ost.close ();
	os << "<a href=\"" << outputPath << "/" << filename << "\">";
	os << "Tab delimited results report. Click to view the report.";
	os << "</a><br />" << endl;
}
void MSViewerSearch::printMSViewerReportTabDelimitedText ( ostream& os, bool printNumbers )
{
	int numRows = rows.size ();

	if ( numRows == 0 ) return;

	int numDataLines = numRows - extraTitleLines;

	for ( int a = 0 ; a < titleLines.size () ; a++ ) {
		os << titleLines [a] << endl;
	}
	BoolDeque removeColumn = viewerParams.getRemoveColumn ( numCols );
	if ( printNumbers ) {
		delimitedRowStart ( os );
			for ( int i = 1 ; i <= numCols ; i++ ) {
				if ( !removeColumn [i-1] ) {
					if	( i != peptide2ColumnNumber ) delimitedHeader ( os, gen_itoa ( i ) );
				}
			}
		delimitedRowEnd ( os );
	}
	if ( !headerLines.empty () ) {
		for ( int i = 0 ; i < headerLines.size () ; i++ ) {
			delimitedRowStart ( os );
				for ( int j = 0 ; j < headerLines [i].size () ; j++ ) {
					if ( !removeColumn [j] ) {
						delimitedHeader ( os, headerLines [i][j] );
					}
				}
			delimitedRowEnd ( os );
		}
	}
	for ( int j = extraTitleLines ; j < numRows ; j++ ) {
		int m = j - extraTitleLines;
		if ( m >= numDataLines ) break;
		delimitedRowStart ( os );
			for ( int k = 1 ; k <= numCols ; k++ ) {
				if ( !removeColumn [k-1] ) {
					delimitedCell ( os, rows [j][k-1] );
				}
			}
		delimitedRowEnd ( os );
	}
}
void MSViewerSearch::printMSViewerReportTabDelimitedTextWithURL ( ostream& os, bool printNumbers )
{
	int numRows = rows.size ();

	if ( numRows == 0 ) return;

	map <string, MSProductLink*> productLink;
	for ( int b = 0 ; b < peakListCentroidFiles.size () ; b++ ) {
		string fname = peakListFPath + SLASH + peakListCentroidFiles [b];
		string fracName = peakListFractionNames [b];
		productLink [fracName] = new MSProductLink ( viewerParams.getInstrumentName (), fname, viewerParams.getParentTolerance (), viewerParams.getFragmentTolerance () );
	}

	int numDataLines = numRows - extraTitleLines;

	for ( int a = 0 ; a < titleLines.size () ; a++ ) {
		os << titleLines [a] << endl;
	}
	BoolDeque removeColumn = viewerParams.getRemoveColumn ( numCols );
	if ( printNumbers ) {
		delimitedRowStart ( os );
			for ( int i = 1 ; i <= numCols ; i++ ) {
				if ( !removeColumn [i-1] ) {
					if	( i != peptide2ColumnNumber ) delimitedHeader ( os, gen_itoa ( i ) );
				}
			}
			delimitedHeader ( os, gen_itoa ( numCols+1 ) );		// Extra column for URL
		delimitedRowEnd ( os );
	}
	if ( !headerLines.empty () ) {
		for ( int i = 0 ; i < headerLines.size () ; i++ ) {
			delimitedRowStart ( os );
				for ( int j = 0 ; j < headerLines [i].size () ; j++ ) {
					if ( !removeColumn [j] ) {
						int idx = j+1;
						delimitedHeader ( os, headerLines [i][j] );
						if ( i == 0 ) {
							if ( peptide2ColumnNumber ) {
								if ( idx == peptide2ColumnNumber )	delimitedHeader ( os, "URL" );
							}
							else {
								if ( idx == peptideColumnNumber )	delimitedHeader ( os, "URL" );
							}
						}
					}
				}
			delimitedRowEnd ( os );
		}
	}
	for ( int j = extraTitleLines ; j < numRows ; j++ ) {
		int m = j - extraTitleLines;
		if ( m >= numDataLines ) break;
		delimitedRowStart ( os );
			for ( int k = 1 ; k <= numCols ; k++ ) {
				if ( !removeColumn [k-1] ) {
					delimitedCell ( os, rows [j][k-1] );
					if ( peptide2ColumnNumber )	{
						if ( k == peptide2ColumnNumber )	urlDelimitedCell ( os, productLink, m );
					}
					else {
						if ( k == peptideColumnNumber )		urlDelimitedCell ( os, productLink, m );
					}
				}
			}
		delimitedRowEnd ( os );
	}
}
void MSViewerSearch::printMSViewerFilteredPeakList ( ostream& os )
{
	if ( spectrumIdentifier != SpecIDPPRT || extraTitleLines ) {
		ErrorHandler::genError ()->error ( "This option is not available for the current data set." );
	}
	int numRows = rows.size ();

	if ( numRows == 0 ) return;

	int numDataLines = numRows - extraTitleLines;

	IntVector sortLevel;
	if ( fractionColumnNumber ) sortLevel.push_back ( fractionColumnNumber );
	sortLevel.push_back ( scanIDColumnNumber );
	StringVector sortOrderType;
	if ( fractionColumnNumber ) sortOrderType.push_back ( "Alphabetic" );
	sortOrderType.push_back ( "Numeric" );
	StringVector sortOrderDirection;
	if ( fractionColumnNumber ) sortOrderDirection.push_back ( "Ascending" );
	sortOrderDirection.push_back ( "Ascending" );
	if ( !sortLevel.empty () && !rows.empty () ) {
		try {
			stable_sort ( rows.begin (), rows.end (), SortRowsByColumn ( numCols, sortLevel, sortOrderType, sortOrderDirection ) );
		}
		catch ( runtime_error e ) {
			ErrorHandler::genError ()->error ( e );
		}
	}

	PPTempFile pptf ( "", "" );
	string fullPath = pptf.getFullPath ();
	genCreateDirectory ( fullPath );

	FileList fileList ( peakListFPath, "", ".mgf", true );
	string fractionName;
	StringVector nameList = fileList.getNameList ();
	if ( !nameList.empty () ) {
		fractionName = nameList [0];		// Single fraction - fraction column may not be present
	}

	MSMSPeakListDataSetInfo* dsi = 0;
	string curFraction;
	string actualPath;
	StringVector actualPaths;
	for ( int j = extraTitleLines ; j < numRows ; j++ ) {
		int m = j - extraTitleLines;
		if ( m >= numDataLines ) break;
		if ( fractionColumnNumber ) {
			fractionName = rows [j][fractionColumnNumber-1];
		}
		string rt = rows [j][scanIDColumnNumber-1];
		if ( fractionName != curFraction ) {
			if ( !curFraction.empty () ) delete dsi;
			curFraction = fractionName;
			dsi = new MSMSPeakListDataSetInfo ( peakListFPath + "/" + fractionName + ".mgf" );
			actualPaths.push_back ( fullPath + SLASH + fractionName + ".mgf" );
		}
		GenOFStream ost ( actualPaths.back (), std::ios_base::out | std::ios_base::app );
		try {
			dsi->writePeakList ( ost, SpecID ( getSpecIDStringFromScanID ( rt ) ), "" );		// version assumed to not matter
		}
		catch ( runtime_error e ) {
			ErrorHandler::genError ()->error ( e );
		}
	}
	if ( dsi ) {
		delete dsi;
		string projectName = genFilenameFromPath ( peakListFPath );
		bool flag = gen7zaCreate ( fullPath + SLASH + projectName, fullPath + SLASH + "*.mgf", "zip" );
		os << "<br /><br />" << endl;
		if ( flag ) {
			printHTMLLink ( os, pptf.getURL (), projectName + ".zip", "Filtered Peak Lists" );
		}
		genUnlink ( actualPaths );
	}
}
void MSViewerSearch::urlDelimitedCell ( ostream& os, map <string, MSProductLink*>& productLink, int m )
{
	if ( fractionName.empty () || scanID.empty () || peptide [m].empty () || productLink.empty () || charge.empty ()  )
		delimitedEmptyCell ( os );
	else {
		string mod;
		if ( allMods.empty () ) {
			string cMod = constantMods.empty () ? "" : constantMods [m];
			if ( cMod.empty () && !constMods.empty () ) cMod = getConstModsString ( peptide [m] );
			string vMod = variableMods.empty () ? "" : variableMods [m];
			mod = cMod.empty () ? vMod : ( vMod.empty () ? cMod : cMod + ';' + vMod );
		}
		else {
			mod = allMods [m];
		}
		if ( productLink.find ( fractionName [m] ) != productLink.end () ) {
			ostringstream ost;
			if ( spectrumIdentifier == SpecIDPPRT ) {
				SpecID specID ( getSpecIDStringFromScanID ( scanID [m] ) );
				if ( !peptide2.empty () ) {
					string lst = linkerColumnNumber ? linker [m] : viewerParams.getLinkSearchType ();
					productLink [fractionName [m]]->write5 ( ost, specID, peptide [m], peptide2 [m], charge [m], true, lst, viewerParams.getBridgeFormula () );
				}
				else
					productLink [fractionName [m]]->write5 ( ost, specID, peptide [m], mod, charge [m], true );
			}
			else if ( spectrumIdentifier == SpecIDTitle )
				productLink [fractionName [m]]->write5 ( ost, scanID [m], peptide [m], mod, charge [m], true );
			else if ( spectrumIdentifier == SpecIDSpecNo )
				productLink [fractionName [m]]->write6 ( ost, scanID [m], peptide [m], mod, charge [m], true );
			else if ( spectrumIdentifier == SpecIDMZ )
				productLink [fractionName [m]]->write7 ( ost, scanID [m], peptide [m], mod, charge [m], true );
			else if ( spectrumIdentifier == SpecIDScanNo )
				productLink [fractionName [m]]->write8 ( ost, scanID [m], peptide [m], mod, charge [m], true );
			else if ( spectrumIdentifier == SpecIDScanNoMinusOne ) {
				int sc = atoi ( scanID [m].c_str () ) + 1;
				productLink [fractionName [m]]->write8 ( ost, gen_itoa ( sc ), peptide [m], mod, charge [m], true );
			}
			delimitedCell ( os, ost.str () );
		}
		else
			delimitedEmptyCell ( os );
	}
}
void MSViewerSearch::printMSViewerReportViewerFiles ( ostream& os )
{
	string sKey = viewerParams.getSearchKey ();
	string directoryPath = ParameterList::getVirtualDir ();
	directoryPath += "results/msviewer/";
	directoryPath += sKey [0];
	directoryPath += '/';
	directoryPath += sKey [1];
	directoryPath += '/';
	directoryPath += sKey;
	printHTMLLink ( os, directoryPath, genFilenameFromPath ( resultsFPath ), "Report" );
	os << "<br />" << endl;
	os << "<br />" << endl;
	directoryPath += '/';
	directoryPath += genFilenameFromPath ( peakListFPath );
	for ( int i = 0 ; i < peakListCentroidFiles.size () ; i++ ) {
		printHTMLLink ( os, directoryPath, peakListCentroidFiles [i], "Data Fraction" );
	}
}
string MSViewerSearch::getConstModsString ( const string& dbPeptide )
{
	string cModString;
	MapCharToStringConstIterator cur = constMods.find ( 'n' );
	if ( cur != constMods.end () ) {
		if ( !cModString.empty () ) cModString += ';';
		cModString += (*cur).second + "@N-term";
	}
	MapCharToStringConstIterator cur2 = constMods.find ( 'c' );
	if ( cur2 != constMods.end () ) {
		if ( !cModString.empty () ) cModString += ';';
		cModString += (*cur2).second + "@C-term";
	}
	for ( MapCharToStringConstIterator i = constMods.begin () ; i != constMods.end () ; i++ ) {
		char aa = (*i).first;
		string s = (*i).second;
		if ( aa != 'n' && aa != 'c' ) {
			int startInd = 0;
			for ( ; ; ) {
				int ind = dbPeptide.find ( aa, startInd );
				if ( ind == string::npos ) break;
				else {
					if ( !cModString.empty () ) cModString += ';';
					cModString += s + "@" + gen_itoa ( ind+1 );
					startInd = ind + 1;
				}
			}
		}
	}
	return cModString;
}
// Example tandem mod:  M [1651] 15.9949, T [1650] 203.079
string MSViewerSearch::convertTandemMod ( const string& tMod, int startAA )
{
	if ( genEmptyString ( tMod ) ) return "";
	string newMod;
	string::size_type start1 = 0;
	bool last = false;
	for ( ; ; ) {
		string curMod;
		string::size_type end1 = tMod.find ( ',', start1 );
		if ( end1 == string::npos ) {
			curMod = tMod.substr ( start1 );
			last = true;
		}
		else
			curMod = tMod.substr ( start1, end1-start1 );

		string::size_type start2 = curMod.find ( '[' ) + 1;
		string::size_type end2 = curMod.find ( ']' );
		int modRes = atoi ( curMod.substr ( start2, end2-start2 ).c_str () ) - startAA + 1;

		string::size_type start3 = end2 + 2;
		if ( !newMod.empty () ) newMod += ';';
		newMod += curMod.substr ( start3 ) + "@" + gen_itoa ( modRes );
		if ( last ) break;
		start1 = end1 + 1;
	}
	return newMod;
}
string MSViewerSearch::getSpecIDStringFromScanID ( const string& scanID )
{
	return scanID.find ( "-" ) == string::npos ? "1-" + scanID + "-1-1" : "1-" + scanID;
}
