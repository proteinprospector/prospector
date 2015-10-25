/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_viewer_form.cpp                                            *
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
#ifndef VIS_C
#include <stdexcept>
#endif
#include <lg_io.h>
#include <lu_getfil.h>
#include <lu_html.h>
#include <lu_table.h>
#include <lu_data_form.h>
#include <lu_viewer_form.h>

using std::string;
using std::ostream;
using std::endl;
using std::runtime_error;
using std::ostringstream;

//********************************************************************************
ViewerResultsFileConverter::ViewerResultsFileConverter ()
{
	init ();
}
ViewerResultsFileConverter& ViewerResultsFileConverter::instance ()
{
	static ViewerResultsFileConverter d;
	return d;
}
void ViewerResultsFileConverter::init ()
{
	GenCommentedIFStream ifs ( MsparamsDir::instance ().getParamPath ( "viewer_conv.txt" ) );
	string line;
	int phase = 1;
	string menuItem;
	ViewerFileConverter vfc;
	while ( ifs.getUncommentedLine ( line ) ) {
		if ( phase == 1 )		menuItem = line;
		else if ( phase == 2 )	vfc.script = line;
		else if ( phase == 3 )	vfc.numTitleLines = atoi ( line.c_str () );
		else if ( phase == 4 )	vfc.numHeaderLines = atoi ( line.c_str () );
		else if ( phase == 5 )	{
			if ( line == "TAB" )		vfc.delimiter = "\t";
			else if ( line == "CSV" )	vfc.delimiter = ",";
			else throw runtime_error ( "Illegal delimiter for entry " + menuItem + " in file viewer_conv.txt." );
		}
		else if ( phase == 6 )	vfc.spectrumIdentifier = line;
		else if ( phase == 7 )	vfc.fraction = line;
		else if ( phase == 8 )	vfc.scanID = line;
		else if ( phase == 9 )	vfc.peptide = line;
		else if ( phase == 10 )	vfc.charge = line;
		else if ( phase == 11 )	{
			vfc.modifications = line;
			converters [menuItem] = vfc;
		}
		phase = ( phase == 11 ) ? 1 : phase + 1;
	}
}
StringVector ViewerResultsFileConverter::getConverterList ()
{
	StringVector sv;
	sv.push_back ( "Protein Prospector Tab Delimited" );
	sv.push_back ( "Protein Prospector Crosslinked Peptides Tab Delimited" );
	sv.push_back ( "PRIDE XML" );
	sv.push_back ( "pepXML" );
	sv.push_back ( "BiblioSpec" );
	sv.push_back ( "Thermo MSF" );
	sv.push_back ( "NIST MSP" );
	sv.push_back ( "SpectraST sptxt" );
	for ( MapStringToViewerFileConverterConstIterator i = converters.begin () ; i != converters.end () ; i++ ) {
		sv.push_back ( (*i).first );
	}
	sv.push_back ( "Other" );
	return sv;
}
string ViewerResultsFileConverter::getScript ( const string& name ) const
{
	MapStringToViewerFileConverterConstIterator cur = converters.find ( name );
	if ( cur != converters.end () ) return (*cur).second.script;
	else {
		ErrorHandler::genError ()->message ( "Converter not found.\n" );
		return "";
	}
}
void ViewerResultsFileConverter::getScriptParameters ( const string& name, int& numTitleLines, int& numHeaderLines, string& delimiter ) const
{
	MapStringToViewerFileConverterConstIterator cur = converters.find ( name );
	if ( cur != converters.end () ) {
		numTitleLines	= (*cur).second.numTitleLines;
		numHeaderLines	= (*cur).second.numHeaderLines;
		delimiter		= (*cur).second.delimiter;
	}
	else ErrorHandler::genError ()->message ( "Converter not found.\n" );
}
void ViewerResultsFileConverter::getScriptParameters2 ( const string& name, string& spectrumIdentifier, string& fraction, string& scanID, string& peptide, string& charge, string& modifications ) const
{
	MapStringToViewerFileConverterConstIterator cur = converters.find ( name );
	if ( cur != converters.end () ) {
		spectrumIdentifier	= (*cur).second.spectrumIdentifier;
		fraction			= (*cur).second.fraction;
		scanID				= (*cur).second.scanID;
		peptide				= (*cur).second.peptide;
		charge				= (*cur).second.charge;
		modifications		= (*cur).second.modifications;
	}
	else ErrorHandler::genError ()->message ( "Converter not found.\n" );
}
StringVector ViewerResultsFileConverter::getNameList () const
{
	StringVector sv;
	for ( MapStringToViewerFileConverterConstIterator i = converters.begin () ; i != converters.end () ; i++ ) {
		sv.push_back ( (*i).first );
	}
	return sv;
}
//********************************************************************************
class FormItemSKey : public FormItemText {
public:
	FormItemSKey ( const string& searchKey = "" );
	static std::string getName () { return "search_key"; }
};
FormItemSKey::FormItemSKey ( const string& searchKey ) :
	FormItemText ( "Search Key", "", getName (), 12, 12, searchKey )
{
}
//********************************************************************************
class FormItemRowsPerPage : public FormItemSelect {
	static const char* options [];
public:
	FormItemRowsPerPage ();
	static string getName () { return "rows_per_page"; }
};
const char* FormItemRowsPerPage::options [] = {
	"20", "50", "100", "200", "500", "1000", "2000", "5000", "10000", "All", 0
};
FormItemRowsPerPage::FormItemRowsPerPage () :
	FormItemSelect ( "Rows Per Page", "", getName (), options, "20" )
{
}
//********************************************************************************
class FormItemPageNumber : public FormItemText {
public:
	FormItemPageNumber ( FormValidatingJavascript* fvj );
	static std::string getName () { return "page"; }
};
FormItemPageNumber::FormItemPageNumber ( FormValidatingJavascript* fvj ) :
	FormItemText ( "Page Number", "", getName (), 4, 4, "1", fvj->addPositiveNonZeroIntegerValidator ( getName (), "Page Number" ) )
{
}
//********************************************************************************
class FormItemNumTitleLines : public FormItemText {
public:
	FormItemNumTitleLines ( FormValidatingJavascript* fvj );
	static std::string getName () { return "num_title_lines"; }
};
FormItemNumTitleLines::FormItemNumTitleLines ( FormValidatingJavascript* fvj ) :
	FormItemText ( "Num Title Lines", "", getName (), 2, 2, "0", fvj->addPositiveIntegerValidator ( getName (), "Num Title Lines" ) )
{
}
//********************************************************************************
class FormItemNumHeaderLines : public FormItemText {
public:
	FormItemNumHeaderLines ( FormValidatingJavascript* fvj );
	static std::string getName () { return "num_header_lines"; }
};
FormItemNumHeaderLines::FormItemNumHeaderLines ( FormValidatingJavascript* fvj ) :
	FormItemText ( "Num Header Lines", "", getName (), 2, 2, "0", fvj->addPositiveIntegerValidator ( getName (), "Num Header Lines" ) )
{
}
//********************************************************************************
class FormItemResultsFileFormat : public FormItemSelect {
public:
	FormItemResultsFileFormat ();
	static string getName () { return "results_file_format"; }
};
FormItemResultsFileFormat::FormItemResultsFileFormat () :
	FormItemSelect ( "Results File Format", "", getName (), ViewerResultsFileConverter::instance ().getConverterList (), "Protein Prospector Tab Delimited", 1, "showResultsFileFormatItems( this.form." + getName () + " )" )
{
}
//********************************************************************************
class FormItemInstrumentFilter : public FormItemSelect {
	static const char* options [];
public:
	FormItemInstrumentFilter ();
	static string getName () { return "instrument_filter"; }
};
const char* FormItemInstrumentFilter::options [] = {
	"", "CID", "ETD", 0
};
FormItemInstrumentFilter::FormItemInstrumentFilter () :
	FormItemSelect ( "Instrument Filter", "", getName (), options, "" )
{
}
//********************************************************************************
class FormItemProbabilityLimit : public FormItemText {
public:
	FormItemProbabilityLimit ( FormValidatingJavascript* fvj ) :
		FormItemText ( "Probability Limit", "", getName (), 6, 6, "0.05", fvj->addPositiveFloatingPointValidator (getName (), "Probability Limit") ) {}
	static string getName () { return "probability_limit"; }
};
//********************************************************************************
class FormItemColumnSeparator : public FormItemSelect {
	static const char* options [];
public:
	FormItemColumnSeparator ();
	static string getName () { return "column_separator"; }
};
const char* FormItemColumnSeparator::options [] = {
	"Tab Delimited", "CSV", 0
};
FormItemColumnSeparator::FormItemColumnSeparator () :
	FormItemSelect ( "Column Separator", "", getName (), options, "Tab Delimited" )
{
}
//********************************************************************************
class FormItemSpectrumIdentifier : public FormItemSelect {
	static const char* options [];
public:
	FormItemSpectrumIdentifier ();
	static string getName () { return "spectrum_identifier"; }
};
const char* FormItemSpectrumIdentifier::options [] = {
	"Protein Prospector RT",
	"Scan Title (Mascot/X!Tandem)",
	"Spectrum Number",
	"m/z",
	"Scan Number",
	"Scan Number - 1", 0
};
FormItemSpectrumIdentifier::FormItemSpectrumIdentifier () :
	FormItemSelect ( "Spectrum Identifier", "", getName (), options, options [0] )
{
}
//********************************************************************************
class FormItemModifications : public FormItemSelect {
	static const char* options [];
public:
	FormItemModifications ();
	static string getName () { return "modifications"; }
};
const char* FormItemModifications::options [] = {
	"Variable Mods In Peptide",
	"All Mods In Peptide",
	"Variable Mods Column",
	"All Mods (1 Column)",
	"All Mods (2 Columns)", 0 
};
FormItemModifications::FormItemModifications () :
	FormItemSelect ( "Modification Reporting", "", getName (), options, options [2], 1, "showModificationItems( this.form." + getName () + " )" )
{
}
//********************************************************************************
class FormItemSortOrderType : public FormItemSelect {
	static const char* options [];
public:
	FormItemSortOrderType ( int i = 0 );
	static string getName ( int i )
	{
		if ( i == 0 )	return "sort_order_type";
		else			return "sort_order_type_" + gen_itoa ( i );
	}
};
const char* FormItemSortOrderType::options [] = {
	"Alphabetic", "Numeric", 0
};
FormItemSortOrderType::FormItemSortOrderType ( int i ) :
	FormItemSelect ( "", "", getName ( i ), options, "Alphabetic" )
{
}
//********************************************************************************
class FormItemSortOrderDirection : public FormItemSelect {
	static const char* options [];
public:
	FormItemSortOrderDirection ( int i = 0 );
	static string getName ( int i )
	{
		if ( i == 0 )	return "sort_order_direction";
		else			return "sort_order_direction_" + gen_itoa ( i );
	}
};
const char* FormItemSortOrderDirection::options [] = {
	"Ascending", "Descending", 0
};
FormItemSortOrderDirection::FormItemSortOrderDirection ( int i ) :
	FormItemSelect ( "", "", getName ( i ), options, "Ascending" )
{
}
//********************************************************************************

class FormItemColumnNumber : public FormItemSelect {
	StringVector getOptions ( int numCols );
public:
	FormItemColumnNumber ( const string& label, const string& field, int numCols = 100 );
	static string getName ( const string& field ) { return "column_num_" + field; }
};

FormItemColumnNumber::FormItemColumnNumber ( const string& label, const string& field, int numCols ) :
	FormItemSelect ( label, "", getName ( field ), getOptions (numCols), "Undefined" )
{
}
StringVector FormItemColumnNumber::getOptions ( int numCols )
{
	StringVector sv;
	sv.push_back ( "Undefined" );
	for ( int i = 1 ; i <= numCols ; i++ ) {
		sv.push_back ( gen_itoa ( i ) );
	}
	return sv;
}
//********************************************************************************

class FormItemRemoveColumn : public FormItemSelectMultiple {
	StringVector getOptions ( int numCols );
public:
	FormItemRemoveColumn ( int numCols );
	static std::string getName () { return "remove_column"; }
};
FormItemRemoveColumn::FormItemRemoveColumn ( int numCols ) :
	FormItemSelectMultiple ( "Remove<br />Column", "", getName (), getOptions (numCols), StringVector (), 8 )
{
}
StringVector FormItemRemoveColumn::getOptions ( int numCols )
{
	StringVector sv;
	for ( int i = 1 ; i <= numCols ; i++ ) {
		sv.push_back ( gen_itoa ( i ) );
	}
	return sv;
}
//********************************************************************************

class FormItemFilterType : public FormItemSelect {
	static const char* options [];
public:
	FormItemFilterType ( int i = 0 );
	static string getName ( int i )
	{
		if ( i == 0 )	return "filter_type";
		else			return "filter_type_" + gen_itoa ( i );
	}
};
const char* FormItemFilterType::options [] = {
	"Equals",
	"Not Equal To",
	"Greater Than Alphabetic",
	"Greater Than Numeric",
	"Less Than Alphabetic",
	"Less Than Numeric",
	"Contains",
	"Does Not Contain",
	"Prefix",
	"Suffix",
	0
};
FormItemFilterType::FormItemFilterType ( int i ) :
	FormItemSelect ( "", "", getName ( i ), options, "Equals" )
{
}
//********************************************************************************

class FormItemFilterValue : public FormItemTextArea {
public:
	FormItemFilterValue ( int i = 0 );
	static std::string getName ( int i )
	{
		if ( i == 0 )	return "filter_value";
		else			return "filter_value_" + gen_itoa ( i );
	}
};
FormItemFilterValue::FormItemFilterValue ( int i ) :
	FormItemTextArea ( "", "", getName ( i ), 2, 40, StringVector () )
{
}
//********************************************************************************

class FormItemUploadResults : public FormItemFile {
public:
	FormItemUploadResults () :
		FormItemFile ( "", "", getName (), 100, 256, "" ) {}
	static string getName () { return "upload_temp_results"; }	// temp in the name indicates that the upload is written to the temp directory
};
//********************************************************************************
class FormItemViewerOutputType : public FormItemSelect {
	static const char* options [];
public:
	FormItemViewerOutputType ();
	static string getName () { return "viewer_output_type"; }
};
FormItemViewerOutputType::FormItemViewerOutputType () :
	FormItemSelect ( "Output", "", getName (), options, "HTML" )
{
}
const char* FormItemViewerOutputType::options [] = {	"HTML",
														"Tab delimited text",
														"Tab delimited text with URL",
														"Filtered Viewer files",
														"Viewer files", 0 };

//********************************************************************************

const int MSViewerForm::maxSortLevels = 6;
const int MSViewerForm::maxFilterLevels = 2;
const int MSViewerForm::maxRemoveReplicateLevels = 5;

MSViewerForm::MSViewerForm ( const VectorConstParameterListPtr& params, bool createDataset, int numCols ) :
	createDataset ( createDataset ),
	ionTypeForm ( params ),
	searchForm ( params, "mstag", &fvj, 2 ),
	msmsToleranceForm ( &fvj, params, true, false, 2 ),
	variableModsForm ( params, "msms_", &fvj ),
	extraUserModsForm ( params, &fvj, false ),
	extraUserModsForm2 ( params, &fvj, true, 2 ),
	massModificationForm ( params, &fvj ),
	crosslinkingForm ( params, &fvj, true, 2 ),
	matrixModeForm ( params, "msms" ),
	filterPeaks ( true ),
	numCols ( numCols )
{
	create ( params );
}
void MSViewerForm::printJavascriptFunctions ( ostream& os )
{
	if ( createDataset ) {
		basicJavascriptFunctions ( os );
		startJavascript ( os );
		os << "function showModificationItems ( item ) {" << endl;
		os << "\t" << "var val = getSelectValue(item);" << endl;
		os << "\t" << "if ( val == 'Variable Mods In Peptide' ) {" << endl;
		os << "\t\t" << "showdiv ( 'div_cm' );" << endl;
		os << "\t\t" << "hidediv ( 'div_cm_col' );" << endl;
		os << "\t\t" << "hidediv ( 'div_vm_col' );" << endl;
		os << "\t\t" << "hidediv ( 'div_am_col' );" << endl;
		os << "\t" << "}" << endl;
		os << "\t" << "else if ( val == 'All Mods In Peptide' ) {" << endl;
		os << "\t\t" << "hidediv ( 'div_cm' );" << endl;
		os << "\t\t" << "hidediv ( 'div_cm_col' );" << endl;
		os << "\t\t" << "hidediv ( 'div_vm_col' );" << endl;
		os << "\t\t" << "hidediv ( 'div_am_col' );" << endl;
		os << "\t" << "}" << endl;
		os << "\t" << "else if ( val == 'Variable Mods Column' ) {" << endl;
		os << "\t\t" << "showdiv ( 'div_cm' );" << endl;
		os << "\t\t" << "hidediv ( 'div_cm_col' );" << endl;
		os << "\t\t" << "showdiv ( 'div_vm_col' );" << endl;
		os << "\t\t" << "hidediv ( 'div_am_col' );" << endl;
		os << "\t" << "}" << endl;
		os << "\t" << "else if ( val == 'All Mods (1 Column)' ) {" << endl;
		os << "\t\t" << "hidediv ( 'div_cm' );" << endl;
		os << "\t\t" << "hidediv ( 'div_cm_col' );" << endl;
		os << "\t\t" << "hidediv ( 'div_vm_col' );" << endl;
		os << "\t\t" << "showdiv ( 'div_am_col' );" << endl;
		os << "\t" << "}" << endl;
		os << "\t" << "else if ( val == 'All Mods (2 Columns)' ) {" << endl;
		os << "\t\t" << "hidediv ( 'div_cm' );" << endl;
		os << "\t\t" << "showdiv ( 'div_cm_col' );" << endl;
		os << "\t\t" << "showdiv ( 'div_vm_col' );" << endl;
		os << "\t\t" << "hidediv ( 'div_am_col' );" << endl;
		os << "\t" << "}" << endl;
		os << "}" << endl;
		os << "function showResultsFileFormatItems ( item ) {" << endl;
		os << "\t" << "var val = getSelectValue(item);" << endl;
		os << "\t" << "if ( val == 'Other' ) {" << endl;
		os << "\t\t" << "showdiv ( 'div_info_1' );" << endl;
		os << "\t\t" << "showdiv ( 'div_info_2' );" << endl;
		os << "\t\t" << "showdiv ( 'div_info_3' );" << endl;
		os << "\t\t" << "showdiv ( 'div_info_4' );" << endl;
		os << "\t\t" << "showdiv ( 'div_info_5' );" << endl;
		os << "\t\t" << "showdiv ( 'div_info_6' );" << endl;
		os << "\t\t" << "showdiv ( 'div_info_7' );" << endl;
		os << "\t\t" << "showdiv ( 'div_info_8' );" << endl;
		os << "\t\t" << "showdiv ( 'div_info_9' );" << endl;
		os << "\t\t" << "showdiv ( 'div_info_10' );" << endl;
		os << "\t\t" << "hidediv ( 'div_info_11' );" << endl;
		os << "\t" << "}" << endl;
		os << "\t" << "else if ( val == 'MaxQuant' ) {" << endl;
		os << "\t\t" << "hidediv ( 'div_info_1' );" << endl;
		os << "\t\t" << "hidediv ( 'div_info_2' );" << endl;
		os << "\t\t" << "hidediv ( 'div_info_3' );" << endl;
		os << "\t\t" << "hidediv ( 'div_info_4' );" << endl;
		os << "\t\t" << "hidediv ( 'div_info_5' );" << endl;
		os << "\t\t" << "hidediv ( 'div_info_6' );" << endl;
		os << "\t\t" << "hidediv ( 'div_info_7' );" << endl;
		os << "\t\t" << "hidediv ( 'div_info_8' );" << endl;
		os << "\t\t" << "hidediv ( 'div_info_9' );" << endl;
		os << "\t\t" << "hidediv ( 'div_info_10' );" << endl;
		os << "\t\t" << "showdiv ( 'div_info_11' );" << endl;
		os << "\t" << "}" << endl;
		os << "\t" << "else {" << endl;
		os << "\t\t" << "hidediv ( 'div_info_1' );" << endl;
		os << "\t\t" << "hidediv ( 'div_info_2' );" << endl;
		os << "\t\t" << "hidediv ( 'div_info_3' );" << endl;
		os << "\t\t" << "hidediv ( 'div_info_4' );" << endl;
		os << "\t\t" << "hidediv ( 'div_info_5' );" << endl;
		os << "\t\t" << "hidediv ( 'div_info_6' );" << endl;
		os << "\t\t" << "hidediv ( 'div_info_7' );" << endl;
		os << "\t\t" << "hidediv ( 'div_info_8' );" << endl;
		os << "\t\t" << "hidediv ( 'div_info_9' );" << endl;
		os << "\t\t" << "hidediv ( 'div_info_10' );" << endl;
		os << "\t\t" << "hidediv ( 'div_info_11' );" << endl;
		os << "\t" << "}" << endl;
		os << "}" << endl;
		os << "function showLinkSearchTypeItems ( item ) {" << endl;
		os << "\t" << "var val = getSelectValue(item);" << endl;
		os << "\t" << "if ( val == 'User Defined Link' ) {" << endl;
		os << "\t\t" << "showdiv ( 'div_ls' );" << endl;
		os << "\t" << "}" << endl;
		os << "\t" << "else {" << endl;
		os << "\t\t" << "hidediv ( 'div_ls' );" << endl;
		os << "\t" << "}" << endl;
		os << "}" << endl;
		endJavascript ( os );
	}
}
void MSViewerForm::createItems ()
{
	formItemMap [FormItemSearchName::getName ()]	= new FormItemSearchName ( "msviewer" );
	formItemMap [FormItemReportTitle::getName ()]	= new FormItemReportTitle ( "MS-Viewer" );
	formItemMap [FormItemVersion::getName ()]		= new FormItemVersion ();

	if ( createDataset ) {
		formItemMap [FormItemSaveParameters::getName ()]	= new FormItemSaveParameters ();

		formItemMap [FormItemUploadPeakList::getName ()]	= new FormItemUploadPeakList ();
		formItemMap [FormItemUploadResults::getName ()]		= new FormItemUploadResults ();

		formItemMap ["peak_list_filepath"]	= new FormItemText ( "", "", "peak_list_filepath", 20, 120, "" );
		formItemMap ["results_filepath"]	= new FormItemText ( "", "", "results_filepath", 20, 120, "" );

		formItemMap [FormItemResultsFileFormat::getName ()]			= new FormItemResultsFileFormat ();

		formItemMap [FormItemInstrumentFilter::getName ()]			= new FormItemInstrumentFilter ();
		formItemMap [FormItemProbabilityLimit::getName ()]			= new FormItemProbabilityLimit (&fvj);

		formItemMap [FormItemColumnSeparator::getName ()]			= new FormItemColumnSeparator ();
		formItemMap [FormItemNumTitleLines::getName ()]				= new FormItemNumTitleLines (&fvj);
		formItemMap [FormItemNumHeaderLines::getName ()]			= new FormItemNumHeaderLines (&fvj);

		formItemMap [FormItemColumnNumber::getName ("fraction")]	= new FormItemColumnNumber ( "Fraction Column (Required If Multiple Peaklists)", "fraction", numCols );
		formItemMap [FormItemSpectrumIdentifier::getName ()] = new FormItemSpectrumIdentifier ();
		formItemMap [FormItemColumnNumber::getName ("scan_id")]		= new FormItemColumnNumber ( "Column", "scan_id", numCols );

		formItemMap [FormItemColumnNumber::getName ("peptide")]		= new FormItemColumnNumber ( "Peptide Column", "peptide", numCols );
		formItemMap [FormItemColumnNumber::getName ("z")]			= new FormItemColumnNumber ( "Charge Column", "z", numCols );
		formItemMap [FormItemModifications::getName ()]				= new FormItemModifications ();
		formItemMap [FormItemColumnNumber::getName ("constant_mod")]= new FormItemColumnNumber ( "Constant Mods Column", "constant_mod", numCols );
		formItemMap [FormItemColumnNumber::getName ("variable_mod")]= new FormItemColumnNumber ( "Variable Mods Column", "variable_mod", numCols );
		formItemMap [FormItemColumnNumber::getName ("all_mod")]		= new FormItemColumnNumber ( "Mods Column", "all_mod", numCols );
		formItemMap [FormItemConstMod::getName ()]	= new FormItemConstMod ();
		formItemMap ["uploads_optional"]				= new FormItemCheckbox ( "", "", "uploads_optional", true );
	}
	formItemMap [FormItemParentMassConvert::getName ()]			= new FormItemParentMassConvert ( "MS" );
	formItemMap [FormItemFragmentMassesTolerance::getName ()]	= new FormItemFragmentMassesTolerance (&fvj);
	formItemMap [FormItemFragmentMassesToleranceUnits::getName ()]= new FormItemFragmentMassesToleranceUnits ();
	formItemMap [FormItemInstrumentName::getName ()]= new FormItemInstrumentName ();
	formItemMap [FormItemUseInstrumentIonTypes::getName ()]	= new FormItemUseInstrumentIonTypes ();
	formItemMap [FormItemRemoveColumn::getName ()]	= new FormItemRemoveColumn ( numCols );
	for ( int i = 1 ; i <= maxFilterLevels ; i++ ) {
		string num = gen_itoa (i);
		string label = "Filter " + num + " Column Number";
		string field = "filter_" + num;
		formItemMap [FormItemColumnNumber::getName (field)]	= new FormItemColumnNumber (label, field, numCols);
		formItemMap [FormItemFilterType::getName (i)]		= new FormItemFilterType (i);
		formItemMap [FormItemFilterValue::getName (i)]		= new FormItemFilterValue (i);
	}

	for ( int j = 1 ; j <= maxSortLevels ; j++ ) {
		string num = gen_itoa (j);
		string label = "Sort Level " + num + " Column Number";
		string field = "sort_level_" + num;
		formItemMap [FormItemColumnNumber::getName (field)]		= new FormItemColumnNumber (label, field, numCols);
		formItemMap [FormItemSortOrderType::getName (j)]		= new FormItemSortOrderType (j);
		formItemMap [FormItemSortOrderDirection::getName (j)]	= new FormItemSortOrderDirection (j);
	}
	for ( int k = 1 ; k <= maxRemoveReplicateLevels ; k++ ) {
		string num = gen_itoa (k);
		string label = ( k == 1 ) ? "Column Numbers" : "";
		string field = "replicate_test_" + num;
		formItemMap [FormItemColumnNumber::getName (field)]		= new FormItemColumnNumber (label, field, numCols);
	}
	formItemMap [FormItemSKey::getName ()]			= new FormItemSKey ();
	formItemMap [FormItemRowsPerPage::getName ()]	= new FormItemRowsPerPage ();
	formItemMap [FormItemPageNumber::getName ()]	= new FormItemPageNumber ( &fvj );
	formItemMap [FormItemMSMSPkFilter::getName ()]	= new FormItemMSMSPkFilter ( true );
	formItemMap [FormItemMaxMSMSPeaks::getName ()]	= new FormItemMaxMSMSPeaks ( &fvj, true );
	formItemMap [FormItemLinkSearchType::getName ()]= new FormItemLinkSearchType ( "No Link", "showLinkSearchTypeItems( this.form." + FormItemLinkSearchType::getName () + " )" );
	formItemMap [FormItemComposition::getName ("bridge")]	= new FormItemComposition ( "Bridge", "bridge", 1, "", &fvj );

	formItemMap [FormItemComment::getName()]						= new FormItemComment ();
	formItemMap [FormItemMaxReportedHits::getName("msms_")]			= new FormItemMaxReportedHits ( "msms_", "5",  &fvj, false  );
	formItemMap [FormItemMSMSPkFilter::getName (2)]					= new FormItemMSMSPkFilter ( false, 2 );
	formItemMap [FormItemMaxMSMSPeaks::getName (2)]					= new FormItemMaxMSMSPeaks ( &fvj, true, 2 );
	formItemMap [FormItemExpectationCalculationMethod::getName ()]	= new FormItemExpectationCalculationMethod ( "None" );
	formItemMap [FormItemMSMSMaxModifications::getName ()]			= new FormItemMSMSMaxModifications ( &fvj );
	formItemMap [FormItemMSMSMaxPeptidePermutations::getName ()]	= new FormItemMSMSMaxPeptidePermutations ( &fvj );
	formItemMap [FormItemUseInstrumentIonTypes::getName (2)]		= new FormItemUseInstrumentIonTypes ( true, 2 );

	formItemMap [FormItemInstrumentName::getName (2)] = new FormItemInstrumentName ( 2 );

	formItemMap [FormItemViewerOutputType::getName ()] = new FormItemViewerOutputType ();
}
void MSViewerForm::printHTML ( ostream& os )
{
	init_html ( os, "MS-Viewer", "viewerman.htm" );
	printHTMLForm ( os, false, true );
	printHTMLFooter ( os, "2011" );
}
void MSViewerForm::printHTMLForm ( ostream& os, bool saveParamsOption, bool initialForm )
{
	fvj.print ( os );
	if ( initialForm ) {
		os << "<p>" << endl;
			os << "<a href=\"http://prospector.ucsf.edu/prospector/html/misc/viewereg.htm\">";
			os << "Example MS-Viewer Datasets";
			os << "</a>";
			os << "<br />";
		os << "</p>" << endl;
		printHTMLFORMStart ( os, "post", "mssearch", true, true );
		tableStart ( os, true );
			tableRowStart ( os );
				tableHeaderStart ( os, "", "center", false );
					tableStart ( os );
						tableRowStart ( os );
							tableHeaderStart ( os, "", "center", true );
								printHTMLFORMSubmit ( os, "Get Existing Results" );
							tableHeaderEnd ( os );
						tableRowEnd ( os );
					tableEnd ( os );
				tableHeaderEnd ( os );
			tableRowEnd ( os );
			tableRowStart ( os );
				tableHeaderStart ( os, "", "center", false );
					tableStart ( os );
						tableRowStart ( os );
							tableHeaderStart ( os, "", "center", true );
								formItemMap [FormItemSKey::getName ()]->printHTML ( os );
							tableHeaderEnd ( os );
							tableHeaderStart ( os, "", "center", true );
								formItemMap [FormItemRowsPerPage::getName ()]->printHTML ( os );
							tableHeaderEnd ( os );
							tableHeaderStart ( os, "", "center", true );
								formItemMap [FormItemViewerOutputType::getName ()]->printHTML ( os );
							tableHeaderEnd ( os );
						tableRowEnd ( os );
					tableEnd ( os );
				tableHeaderEnd ( os );
			tableRowEnd ( os );
		tableEnd ( os );
		formItemMap [FormItemSearchName::getName ()]->printHTMLHidden ( os );
		formItemMap [FormItemReportTitle::getName ()]->printHTMLHidden ( os );
		printHTMLFORMEnd ( os );
		os << "<br />" << endl;
		os << "<br />" << endl;
	}
	printHTMLFORMStart ( os, "post", "mssearch", true, true );
	printJavascriptFunctions ( os );
	tableStart ( os, true );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true, 2 );
				tableStart ( os );
					tableRowStart ( os );
						tableHeaderStart ( os, "", "center", true );
							if ( initialForm )	printHTMLFORMSubmit ( os, "Upload New Results" );
							else				printHTMLFORMSubmit ( os, "Regenerate Report" );
						tableHeaderEnd ( os );
						if ( saveParamsOption ) {
							tableHeaderStart ( os, "", "center", true );
								formItemMap [FormItemSaveParameters::getName ()]->printHTML ( os );
							tableHeaderEnd ( os );
						}
					tableRowEnd ( os );
					if ( createDataset ) {
						if ( formItemMap ["peak_list_filepath"]->getValue ( 0 ).empty () ) {
							tableRowStart ( os );
								tableHeaderStart ( os, "", "left", true );
									os << "<div class=\"form_large_label\">Peak List File</div>" << endl;
									formItemMap [FormItemUploadPeakList::getName ()]->printHTML ( os );
								tableHeaderEnd ( os );
							tableRowEnd ( os );
						}
						if ( formItemMap ["results_filepath"]->getValue ( 0 ).empty () ) {
							tableRowStart ( os );
								tableHeaderStart ( os, "", "left", true );
									os << "<div class=\"form_large_label\">Results File</div>" << endl;
									formItemMap [FormItemUploadResults::getName ()]->printHTML ( os );
								tableHeaderEnd ( os );
							tableRowEnd ( os );
						}
					}
					tableRowStart ( os );
						tableHeaderStart ( os, "", "left", true );
							formItemMap [FormItemRowsPerPage::getName ()]->printHTML ( os );
							formItemMap [FormItemPageNumber::getName ()]->printHTML ( os );
							if ( !createDataset ) {
								formItemMap [FormItemViewerOutputType::getName ()]->printHTML ( os );
							}
						tableHeaderEnd ( os );
					tableRowEnd ( os );
					if ( createDataset ) {
						tableRowStart ( os );
							tableHeaderStart ( os, "", "left", true );
								formItemMap [FormItemResultsFileFormat::getName ()]->printHTML ( os );
							tableHeaderEnd ( os );
						tableRowEnd ( os );

						bool div_info_1, div_info_2, div_info_3, div_info_4, div_info_5, div_info_6, div_info_7, div_info_8, div_info_9, div_info_10, div_info_11;
						setResultsFileFormatVisualizationFlags ( formItemMap [FormItemResultsFileFormat::getName ()]->getValue ( 0 ), div_info_1, div_info_2, div_info_3, div_info_4, div_info_5, div_info_6, div_info_7, div_info_8, div_info_9, div_info_10, div_info_11 );
						tableRowStart ( os );
							tableHeaderStart ( os, "", "left", true );
								divStart ( os, "div_info_11", div_info_11 );
									formItemMap [FormItemInstrumentFilter::getName ()]->printHTML ( os );
									formItemMap [FormItemProbabilityLimit::getName ()]->printHTML ( os );
								divEnd ( os );
							tableHeaderEnd ( os );
						tableRowEnd ( os );
						tableRowStart ( os );
							tableHeaderStart ( os, "", "left", true );
								divStart ( os, "div_info_1", div_info_1 );
									formItemMap [FormItemColumnSeparator::getName ()]->printHTML ( os );
									formItemMap [FormItemNumTitleLines::getName ()]->printHTML ( os );
									formItemMap [FormItemNumHeaderLines::getName ()]->printHTML ( os );
								divEnd ( os );
							tableHeaderEnd ( os );
						tableRowEnd ( os );
						tableRowStart ( os );
							tableHeaderStart ( os, "", "left", true );
								divStart ( os, "div_info_2", div_info_2 );
									os << "<div class=\"form_large_label\">Spectrum Identification Columns</div>" << endl;
								divEnd ( os );
							tableHeaderEnd ( os );
						tableRowEnd ( os );
						tableRowStart ( os );
							tableHeaderStart ( os, "", "left", true );
								divStart ( os, "div_info_3", div_info_3 );
									formItemMap [FormItemColumnNumber::getName ("fraction")]->printHTML ( os );
								divEnd ( os );
							tableHeaderEnd ( os );
						tableRowEnd ( os );
						tableRowStart ( os );
							tableHeaderStart ( os, "", "left", true );
								divStart ( os, "div_info_4", div_info_4 );
									formItemMap [FormItemSpectrumIdentifier::getName ()]->printHTML ( os );
									formItemMap [FormItemColumnNumber::getName ("scan_id")]->printHTML ( os );
								divEnd ( os );
							tableHeaderEnd ( os );
						tableRowEnd ( os );
						tableRowStart ( os );
							tableHeaderStart ( os, "", "left", true );
								divStart ( os, "div_info_5", div_info_5 );
									os << "<div class=\"form_large_label\">Peptide Information Columns</div>" << endl;
								divEnd ( os );
							tableHeaderEnd ( os );
						tableRowEnd ( os );
						tableRowStart ( os );
							tableHeaderStart ( os, "", "left", true );
								divStart ( os, "div_info_6", div_info_6 );
									formItemMap [FormItemColumnNumber::getName ("peptide")]->printHTML ( os );
									formItemMap [FormItemColumnNumber::getName ("z")]->printHTML ( os );
								divEnd ( os );
							tableHeaderEnd ( os );
						tableRowEnd ( os );
						tableRowStart ( os );
							tableHeaderStart ( os, "", "left", true );
								divStart ( os, "div_info_7", div_info_7 );
									formItemMap [FormItemModifications::getName ()]->printHTML ( os );
								divEnd ( os );
							tableHeaderEnd ( os );
						tableRowEnd ( os );
						bool div_cm, div_cm_col, div_vm_col, div_am_col;
						setModVisualizationFlags ( formItemMap [FormItemModifications::getName ()]->getValue ( 0 ), div_cm, div_cm_col, div_vm_col, div_am_col );
						tableRowStart ( os );
							tableHeaderStart ( os, "", "left", true );
								divStart ( os, "div_info_8", div_info_8 );
								divStart ( os, "div_cm_col", div_cm_col );
									formItemMap [FormItemColumnNumber::getName ("constant_mod")]->printHTML ( os );
								divEnd ( os );
								divEnd ( os );
							tableHeaderEnd ( os );
						tableRowEnd ( os );
						tableRowStart ( os );
							tableHeaderStart ( os, "", "left", true );
								divStart ( os, "div_info_9", div_info_9 );
								divStart ( os, "div_vm_col", div_vm_col );
									formItemMap [FormItemColumnNumber::getName ("variable_mod")]->printHTML ( os );
								divEnd ( os );
								divEnd ( os );
							tableHeaderEnd ( os );
						tableRowEnd ( os );
						tableRowStart ( os );
							tableHeaderStart ( os, "", "left", true );
								divStart ( os, "div_info_10", div_info_10 );
								divStart ( os, "div_am_col", div_am_col );
									formItemMap [FormItemColumnNumber::getName ("all_mod")]->printHTML ( os );
								divEnd ( os );
								divEnd ( os );
							tableHeaderEnd ( os );
						tableRowEnd ( os );
						tableRowStart ( os );
							tableHeaderStart ( os, "", "left", true );
								divStart ( os, "div_cm", div_cm );
									formItemMap [FormItemConstMod::getName ()]->printHTML ( os );
								divEnd ( os );
							tableHeaderEnd ( os );
						tableRowEnd ( os );
					}
					if ( !initialForm ) {
						tableRowStart ( os );
							tableHeaderStart ( os, "", "left", true );
								formItemMap [FormItemRemoveColumn::getName ()]->printHTML ( os );
							tableHeaderEnd ( os );
						tableRowEnd ( os );
						tableRowStart ( os );
							tableHeaderStart ( os, "", "left", true );
								os << "<div class=\"form_large_label\">Report Filtering</div>" << endl;
							tableHeaderEnd ( os );
						tableRowEnd ( os );
						for ( int i = 1 ; i <= maxFilterLevels ; i++ ) {
							tableRowStart ( os );
								tableHeaderStart ( os, "", "left", true );
									string num = gen_itoa (i);
									string field = "filter_" + num;
									formItemMap [FormItemColumnNumber::getName (field)]->printHTML ( os );
									formItemMap [FormItemFilterType::getName (i)]->printHTML ( os );
									formItemMap [FormItemFilterValue::getName (i)]->printHTML ( os );
								tableHeaderEnd ( os );
							tableRowEnd ( os );
						}
						tableRowStart ( os );
							tableHeaderStart ( os, "", "left", true );
								os << "<div class=\"form_large_label\">Report Sorting</div>" << endl;
							tableHeaderEnd ( os );
						tableRowEnd ( os );
						for ( int j = 1 ; j <= maxSortLevels ; j++ ) {
							tableRowStart ( os );
								tableHeaderStart ( os, "", "left", true );
									formItemMap [FormItemColumnNumber::getName ("sort_level_" + gen_itoa (j))]->printHTML ( os );
									formItemMap [FormItemSortOrderType::getName (j)]->printHTML ( os );
									formItemMap [FormItemSortOrderDirection::getName (j)]->printHTML ( os );
								tableHeaderEnd ( os );
							tableRowEnd ( os );
						}
						tableRowStart ( os );
							tableHeaderStart ( os, "", "left", true );
								os << "<div class=\"form_large_label\">Remove Replicates</div>" << endl;
							tableHeaderEnd ( os );
						tableRowEnd ( os );
						tableRowStart ( os );
							tableHeaderStart ( os, "", "left", true );
							for ( int k = 1 ; k <= maxRemoveReplicateLevels ; k++ ) {
								formItemMap [FormItemColumnNumber::getName ("replicate_test_" + gen_itoa (k))]->printHTML ( os );
							}
							tableHeaderEnd ( os );
						tableRowEnd ( os );
					}
					tableRowStart ( os );
						tableHeaderStart ( os, "", "left", true );
							os << "<div class=\"form_large_label\">MS-Product Parameters</div>" << endl;
						tableHeaderEnd ( os );
					tableRowEnd ( os );
					tableRowStart ( os );
						tableHeaderStart ( os, "", "left", true );
							formItemMap [FormItemParentMassConvert::getName ()]->printHTML ( os );
							formItemMap [FormItemFragmentMassesTolerance::getName ()]->printHTML ( os );
							formItemMap [FormItemFragmentMassesToleranceUnits::getName ()]->printHTML ( os );
							formItemMap [FormItemInstrumentName::getName ()]->printHTML ( os );
						tableHeaderEnd ( os );
					tableRowEnd ( os );
					tableRowStart ( os );
						tableHeaderStart ( os, "", "left", true );
							formItemMap [FormItemMSMSPkFilter::getName ()]->printHTML ( os );
							formItemMap [FormItemMaxMSMSPeaks::getName ()]->printHTML ( os );
						tableHeaderEnd ( os );
					tableRowEnd ( os );
					bool div_ls;
					setLinkSearchVisualizationFlags ( formItemMap [FormItemLinkSearchType::getName ()]->getValue (0), div_ls );
					tableRowStart ( os );
						tableHeaderStart ( os, "", "left", true );
							tableStart ( os, false );
								tableRowStart ( os );
									tableHeaderStart ( os, "", "left", true );
										formItemMap [FormItemLinkSearchType::getName ()]->printHTML ( os );
									tableHeaderEnd ( os );
									tableHeaderStart ( os, "", "left", true );
										divStart ( os, "div_ls", div_ls );
											formItemMap [FormItemComposition::getName ("bridge")]->printHTML ( os );
										divEnd ( os );
									tableHeaderEnd ( os );
								tableRowEnd ( os );
							tableEnd ( os );
						tableHeaderEnd ( os );
					tableRowEnd ( os );
					tableRowStart ( os );
						tableHeaderStart ( os, "", "left", true );
							extraUserModsForm.printHTML ( os );
						tableHeaderEnd ( os );
					tableRowEnd ( os );
					tableRowStart ( os );
						tableHeaderStart ( os, "", "center", false );
							bool open = false;
							ExpandableJavascriptBlock ejb ( "Ion Types", open );
							ejb.printHeader ( os );
							formItemMap [FormItemUseInstrumentIonTypes::getName ()]->printHTML ( os );
							ionTypeForm.printHTML ( os );
							ejb.printFooter ( os );
						tableHeaderEnd ( os );
					tableRowEnd ( os );
					tableRowStart ( os );
						tableHeaderStart ( os, "", "left", true );
							os << "<div class=\"form_large_label\">MS-Tag Parameters</div>" << endl;
						tableHeaderEnd ( os );
					tableRowEnd ( os );
					tableRowStart ( os );
						tableHeaderStart ( os, "", "left", true );
							bool open2 = false;
							ExpandableJavascriptBlock ejb2 ( "", open2 );
							ejb2.printHeader ( os );
							tableStart ( os, true );
								searchForm.printHTML ( os );
								printHTMLFirstTable ( os );
								printHTMLSecondTable ( os );
							tableEnd ( os );
							ejb2.printFooter ( os );
						tableHeaderEnd ( os );
					tableRowEnd ( os );
				tableEnd ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
	tableEnd ( os );
	showHiddenItems ( os );
	printHTMLFORMEnd ( os );
}
void MSViewerForm::printHTMLFirstTable ( ostream& os )
{
	tableRowStart ( os );
		tableHeaderStart ( os, "", "left", true );
			formItemMap [FormItemComment::getName()]->printHTML ( os );
			os << "<br />" << endl;
			formItemMap [FormItemMaxReportedHits::getName("msms_")]->printHTML ( os );
			if ( filterPeaks ) {
				formItemMap [FormItemMSMSPkFilter::getName (2)]->printHTML ( os );
				formItemMap [FormItemMaxMSMSPeaks::getName (2)]->printHTML ( os );
			}
			os << "<br />" << endl;
			formItemMap [FormItemExpectationCalculationMethod::getName ()]->printHTML ( os );
		tableHeaderEnd ( os );
	tableRowEnd ( os );
}
void MSViewerForm::printHTMLSecondTable ( ostream& os )
{
	basicJavascriptFunctions ( os );
	massModCrosslinkingJavascriptFunctions ( os, 2 );
	tableRowStart ( os );
		tableHeaderStart ( os, "", "left", true );
			tableStart ( os );
				tableRowStart ( os );
					tableHeaderStart ( os, "", "left", true );
						variableModsForm.printHTML ( os );
					tableHeaderEnd ( os );
				tableRowEnd ( os );
				tableRowStart ( os );
					tableHeaderStart ( os, "", "left", true );
						extraUserModsForm2.printHTML ( os );
					tableHeaderEnd ( os );
				tableRowEnd ( os );
				tableRowStart ( os );
					tableHeaderStart ( os, "", "left", true );
						formItemMap [FormItemMSMSMaxModifications::getName ()]->printHTML ( os );
						formItemMap [FormItemMSMSMaxPeptidePermutations::getName ()]->printHTML ( os );
					tableHeaderEnd ( os );
				tableRowEnd ( os );
			tableEnd ( os );
			massModificationForm.printHTML ( os );
			os << "<br />" << endl;
			divStart ( os, "div_cl", massModificationForm.getDivCL () );
				crosslinkingForm.printHTML ( os );
				os << "<br />" << endl;
			divEnd ( os );
			matrixModeForm.printHTML ( os );
		tableHeaderEnd ( os );
		tableHeaderStart ( os, "", "left", true );
			msmsToleranceForm.printHTML ( os );
		tableHeaderEnd ( os );
	tableRowEnd ( os );
	tableRowStart ( os );
		tableHeaderStart ( os, "", "left", true );
			formItemMap [FormItemInstrumentName::getName (2)]->printHTML ( os );
		tableHeaderEnd ( os );
	tableRowEnd ( os );
}
void MSViewerForm::setModVisualizationFlags ( const string& val, bool& div_cm, bool& div_cm_col, bool& div_vm_col, bool& div_am_col ) const
{
	if ( val == "Variable Mods In Peptide" ) {
		div_cm = true;
		div_cm_col = false;
		div_vm_col = false;
		div_am_col = false;
	}
	else if ( val == "All Mods In Peptide" ) {
		div_cm = false;
		div_cm_col = false;
		div_vm_col = false;
		div_am_col = false;
	}
	else if ( val == "Variable Mods Column" ) {
		div_cm = true;
		div_cm_col = false;
		div_vm_col = true;
		div_am_col = false;
	}
	else if ( val == "All Mods (1 Column)" ) {
		div_cm = false;
		div_cm_col = false;
		div_vm_col = false;
		div_am_col = true;
	}
	else if ( val == "All Mods (2 Columns)" ) {
		div_cm = false;
		div_cm_col = true;
		div_vm_col = true;
		div_am_col = false;
	}
}
void MSViewerForm::setLinkSearchVisualizationFlags ( const string& val, bool& div_ls )
{
	if ( val == "User Defined Link" ) {
		div_ls = true;
	}
	else {
		div_ls = false;
	}
}
void MSViewerForm::setResultsFileFormatVisualizationFlags ( const string& val, bool& div_info_1, bool& div_info_2, bool& div_info_3, bool& div_info_4, bool& div_info_5, bool& div_info_6, bool& div_info_7, bool& div_info_8, bool& div_info_9, bool& div_info_10, bool& div_info_11 ) const
{
	if ( val == "Other" ) {
		div_info_1 = true;
		div_info_2 = true;
		div_info_3 = true;
		div_info_4 = true;
		div_info_5 = true;
		div_info_6 = true;
		div_info_7 = true;
		div_info_8 = true;
		div_info_9 = true;
		div_info_10 = true;
		div_info_11 = false;
	}
	else {
		div_info_1 = false;
		div_info_2 = false;
		div_info_3 = false;
		div_info_4 = false;
		div_info_5 = false;
		div_info_6 = false;
		div_info_7 = false;
		div_info_8 = false;
		div_info_9 = false;
		div_info_10 = false;
		if ( val == "MaxQuant" ) {
			div_info_11 = true;
		}
		else {
			div_info_11 = false;
		}
	}
}
void MSViewerForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];

		if ( createDataset ) {
			formItemMap ["peak_list_filepath"]->setValue ( p );
			formItemMap ["results_filepath"]->setValue ( p );

			formItemMap [FormItemResultsFileFormat::getName ()]->setValue ( p );
			formItemMap [FormItemInstrumentFilter::getName ()]->setValue ( p );
			formItemMap [FormItemProbabilityLimit::getName ()]->setValue ( p );
			formItemMap [FormItemColumnSeparator::getName ()]->setValue ( p );
			formItemMap [FormItemNumTitleLines::getName ()]->setValue ( p );
			formItemMap [FormItemNumHeaderLines::getName ()]->setValue ( p );

			formItemMap [FormItemColumnNumber::getName ("fraction")]->setValue ( p );
			formItemMap [FormItemSpectrumIdentifier::getName ()]->setValue ( p );
			formItemMap [FormItemColumnNumber::getName ("scan_id")]->setValue ( p );

			formItemMap [FormItemColumnNumber::getName ("peptide")]->setValue ( p );
			formItemMap [FormItemColumnNumber::getName ("z")]->setValue ( p );
			formItemMap [FormItemModifications::getName ()]->setValue ( p );
			formItemMap [FormItemColumnNumber::getName ("constant_mod")]->setValue ( p );
			formItemMap [FormItemColumnNumber::getName ("variable_mod")]->setValue ( p );
			formItemMap [FormItemColumnNumber::getName ("all_mod")]->setValue ( p );
			formItemMap [FormItemConstMod::getName ()]->setValue ( p );
		}
		else {
			formItemMap [FormItemSKey::getName ()]->setValue ( p );
		}
		formItemMap [FormItemRowsPerPage::getName ()]->setValue ( p );
		formItemMap [FormItemPageNumber::getName ()]->setValue ( p );

		formItemMap [FormItemParentMassConvert::getName ()]->setValue ( p );
		formItemMap [FormItemFragmentMassesTolerance::getName ()]->setValue ( p );
		formItemMap [FormItemFragmentMassesToleranceUnits::getName ()]->setValue ( p );
		formItemMap [FormItemInstrumentName::getName ()]->setValue ( p );
		formItemMap [FormItemUseInstrumentIonTypes::getName ()]->setValue ( p );

		for ( int i = 1 ; i <= maxFilterLevels ; i++ ) {
			string num = gen_itoa (i);
			string field = "filter_" + num;
			formItemMap [FormItemRemoveColumn::getName ()]->setValue ( p );
			formItemMap [FormItemColumnNumber::getName (field)]->setValue ( p );
			formItemMap [FormItemFilterType::getName (i)]->setValue ( p );
			formItemMap [FormItemFilterValue::getName (i)]->setValue ( p );
		}
		for ( int j = 1 ; j <= maxSortLevels ; j++ ) {
			string num = gen_itoa (j);
			string field = "sort_level_" + num;
			formItemMap [FormItemColumnNumber::getName (field)]->setValue ( p );
			formItemMap [FormItemSortOrderType::getName (j)]->setValue ( p );
			formItemMap [FormItemSortOrderDirection::getName (j)]->setValue ( p );
		}
		for ( int k = 1 ; k <= maxRemoveReplicateLevels ; k++ ) {
			string num = gen_itoa (k);
			string field = "replicate_test_" + num;
			formItemMap [FormItemColumnNumber::getName (field)]->setValue ( p );
		}
		formItemMap [FormItemMSMSPkFilter::getName ()]->setValue ( p );
		formItemMap [FormItemMaxMSMSPeaks::getName ()]->setValue ( p );
		formItemMap [FormItemLinkSearchType::getName ()]->setValue ( p );
		formItemMap [FormItemComposition::getName ("bridge")]->setValue ( p );

		formItemMap [FormItemComment::getName()]->setValue ( p );
		formItemMap [FormItemMaxReportedHits::getName("msms_")]->setValue ( p );
		formItemMap [FormItemMSMSPkFilter::getName (2)]->setValue ( p );
		formItemMap [FormItemMaxMSMSPeaks::getName (2)]->setValue ( p );
		formItemMap [FormItemExpectationCalculationMethod::getName ()]->setValue ( p );
		formItemMap [FormItemMSMSMaxModifications::getName ()]->setValue ( p );
		formItemMap [FormItemMSMSMaxPeptidePermutations::getName ()]->setValue ( p );
		formItemMap [FormItemUseInstrumentIonTypes::getName (2)]->setValue ( p );
		formItemMap [FormItemInstrumentName::getName (2)]->setValue ( p );
	}
}
void MSViewerForm::printTagParamsCGI ( ostream& os ) const
{
	ostringstream ostr;
	searchForm.printCGI ( ostr );
	msmsToleranceForm.printCGI ( ostr );
	variableModsForm.printCGI ( ostr );
	extraUserModsForm2.printCGI ( ostr );
	massModificationForm.printCGI ( ostr );
	crosslinkingForm.printCGI ( ostr );
	matrixModeForm.printCGI ( ostr );

	formItemMap.find ( FormItemComment::getName() )->second->printCGI ( ostr );
	formItemMap.find ( FormItemMaxReportedHits::getName("msms_") )->second->printCGI ( ostr );
	formItemMap.find ( FormItemMSMSPkFilter::getName (2) )->second->printCGI ( ostr );
	formItemMap.find ( FormItemMaxMSMSPeaks::getName (2) )->second->printCGI ( ostr );
	formItemMap.find ( FormItemExpectationCalculationMethod::getName () )->second->printCGI ( ostr );
	formItemMap.find ( FormItemMSMSMaxModifications::getName () )->second->printCGI ( ostr );
	formItemMap.find ( FormItemMSMSMaxPeptidePermutations::getName () )->second->printCGI ( ostr );
	formItemMap.find ( FormItemUseInstrumentIonTypes::getName (2) )->second->printCGI ( ostr );
	formItemMap.find ( FormItemInstrumentName::getName (2) )->second->printCGI ( ostr );
	string s = ostr.str ();
	s = genReplaceSubstrings ( s, "fragment_masses_tolerance2", "fragment_masses_tolerance" );
	s = genReplaceSubstrings ( s, "fragment_masses_tolerance_units2", "fragment_masses_tolerance_units" );
	s = genReplaceSubstrings ( s, "use_instrument_ion_types2", "use_instrument_ion_types" );
	s = genReplaceSubstrings ( s, "instrument_name2", "instrument_name" );
	s = genReplaceSubstrings ( s, "const_mod2", "const_mod" );
	s = genReplaceSubstrings ( s, "msms_pk_filter2", "msms_pk_filter" );
	s = genReplaceSubstrings ( s, "msms_max_peaks2", "msms_max_peaks" );
	s = genReplaceSubstrings ( s, "link_search_type2", "link_search_type" );
	s = genReplaceSubstrings ( s, "link_aa2", "link_aa" );
	s = genReplaceSubstrings ( s, "bridge_composition2", "bridge_composition" );
	s = genReplaceSubstrings ( s, "parent_mass_convert2", "parent_mass_convert" );
	extraUserModsForm2.replaceSubstrings ( s );
	os << s;
}
void MSViewerForm::printTagParamsHTMLJavascriptHidden ( ostream& os ) const
{
	ostringstream ostr;
	searchForm.printHTMLJavascriptHidden ( ostr );
	msmsToleranceForm.printHTMLJavascriptHidden ( ostr );
	variableModsForm.printHTMLJavascriptHidden ( ostr );
	extraUserModsForm2.printHTMLJavascriptHidden ( ostr );
	massModificationForm.printHTMLJavascriptHidden ( ostr );
	crosslinkingForm.printHTMLJavascriptHidden ( ostr );
	matrixModeForm.printHTMLJavascriptHidden ( ostr );

	formItemMap.find ( FormItemComment::getName() )->second->printHTMLJavascriptHidden ( ostr );
	formItemMap.find ( FormItemMaxReportedHits::getName("msms_") )->second->printHTMLJavascriptHidden ( ostr );
	formItemMap.find ( FormItemMSMSPkFilter::getName (2) )->second->printHTMLJavascriptHidden ( ostr );
	formItemMap.find ( FormItemMaxMSMSPeaks::getName (2) )->second->printHTMLJavascriptHidden ( ostr );
	formItemMap.find ( FormItemExpectationCalculationMethod::getName () )->second->printHTMLJavascriptHidden ( ostr );
	formItemMap.find ( FormItemMSMSMaxModifications::getName () )->second->printHTMLJavascriptHidden ( ostr );
	formItemMap.find ( FormItemMSMSMaxPeptidePermutations::getName () )->second->printHTMLJavascriptHidden ( ostr );
	formItemMap.find ( FormItemUseInstrumentIonTypes::getName (2) )->second->printHTMLJavascriptHidden ( ostr );
	formItemMap.find ( FormItemInstrumentName::getName (2) )->second->printHTMLJavascriptHidden ( ostr );
	string s = ostr.str ();
	s = genReplaceSubstrings ( s, "fragment_masses_tolerance2", "fragment_masses_tolerance" );
	s = genReplaceSubstrings ( s, "fragment_masses_tolerance_units2", "fragment_masses_tolerance_units" );
	s = genReplaceSubstrings ( s, "use_instrument_ion_types2", "use_instrument_ion_types" );
	s = genReplaceSubstrings ( s, "instrument_name2", "instrument_name" );
	s = genReplaceSubstrings ( s, "const_mod2", "const_mod" );
	s = genReplaceSubstrings ( s, "msms_pk_filter2", "msms_pk_filter" );
	s = genReplaceSubstrings ( s, "msms_max_peaks2", "msms_max_peaks" );
	s = genReplaceSubstrings ( s, "link_search_type2", "link_search_type" );
	s = genReplaceSubstrings ( s, "link_aa2", "link_aa" );
	s = genReplaceSubstrings ( s, "bridge_composition2", "bridge_composition" );
	s = genReplaceSubstrings ( s, "parent_mass_convert2", "parent_mass_convert" );
	extraUserModsForm2.replaceSubstrings ( s );
	os << s;
}
