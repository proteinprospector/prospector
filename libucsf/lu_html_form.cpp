/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_html_form.cpp                                              *
*                                                                             *
*  Created    : September 16th 2002                                           *
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
*  Copyright (2002-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lg_string.h>
#include <lg_time.h>
#include <lu_html_form.h>
#include <lu_prog.h>
#include <lu_table.h>
#include <lu_getfil.h>
#include <lu_pq_vector.h>
#include <lu_param_list.h>
using std::string;
using std::ostream;
using std::endl;

void printHTMLFORMJavascriptHidden2 ( ostream& os, const string& name, const string& value )
{
	os << "\tdocument.writeln ( \"";
	os << "<input";
	os << " ";
	os << "type=\\\"hidden\\\"";
	os << " ";
	os << "name=\\\"" << name << "\\\"";
	os << " ";
	os << "value=\\\"" << genEscapePath ( value ) << "\\\"";
	os << " />";
	os << "\" );" << endl;
}
void printHTMLFORMLabel ( ostream& os, const string& label, const string& manual )
{
	if ( !manual.empty () ) {
		os << "<a href=\"";
		os << "../html/instruct/";
		os << manual;
		os << "\">";
	}
	os << label;
	if ( !manual.empty () ) {
		os << "</a>";
	}
	os << endl;
}
void printHTMLFORMStart ( ostream& os, const string& method, const string& program, bool multipart, bool onSubmit, bool newWindow, const string& type )
{
	os << "<form";
	os << " ";
	os << "method=\"" << method << "\"";
	os << " ";
	os << "action=\"" << ProgramLink::getURLStart ( program, type ) << "\"";
	if ( multipart ) os << " enctype=\"multipart/form-data\"";
	if ( newWindow ) os << " target=\"_blank\"";
	if ( onSubmit ) {
		os << " onsubmit=\"return validateForm(this)\"";
	}
	os << ">" << endl;
}
void printHTMLFORMEnd ( ostream& os )
{
	os << "</form>" << endl << endl;
}
void printHTMLFORMSubmit ( ostream& os, const string& label, const string& color )
{
	os << "<table border=\"border\" cellspacing=\"1\" cellpadding=\"10\">" << endl;
	os << "\t" << "<tr valign=\"middle\">" << endl;
	os << "\t" << "<th bgcolor=\"#" << color << "\">" << endl;
	os << "\t\t" << "<input type=\"submit\" value=\"" << label << "\" />" << endl;
	os << "\t" << "</th>" << endl;
	os << "\t" << "</tr>" << endl;
	os << "</table>" << endl;
}
void printHTMLFORMCheckbox ( ostream& os, const string& label, const string& manual, const string& name, bool checked, const string& value, const string& onChangeFunction )
{
	if ( !label.empty () ) {
		printHTMLFORMLabel ( os, label, manual );
	}
	os << "<input";
	os << " ";
	os << "type=\"checkbox\"";
	os << " ";
	os << "name=\"" << name << "\"";
	os << " ";
	os << "value=\"" << value << "\"";
	if ( checked ) {
		os << " checked=\"checked\"";
	}
	if ( !onChangeFunction.empty () ) {
		os << " ";
		os << "onchange=\"" << onChangeFunction << "\"";
	}
	os << " />" << endl;
}
void printHTMLFORMText ( ostream& os, const string& label, const string& manual, const string& name, int size, int maxlength, const string& value, const string& onChangeFunction )
{
	if ( !label.empty () ) {
		printHTMLFORMLabel ( os, label, manual );
	}
	os << "<input";
	os << " ";
	os << "type=\"text\"";
	os << " ";
	os << "name=\"" << name << "\"";
	os << " ";
	os << "size=\"" << size << "\"";
	os << " ";
	os << "maxlength=\"" << maxlength << "\"";
	os << " ";
	os << "value=\"" << value << "\"";
	if ( !onChangeFunction.empty () ) {
		os << " ";
		os << "onchange=\"" << onChangeFunction << "\"";
	}
	os << " />" << endl;
}
void printHTMLFORMPassword ( ostream& os, const string& label, const string& manual, const string& name, int size, int maxlength, const string& value )
{
	if ( !label.empty () ) {
		printHTMLFORMLabel ( os, label, manual );
	}
	os << "<input";
	os << " ";
	os << "type=\"password\"";
	os << " ";
	os << "name=\"" << name << "\"";
	os << " ";
	os << "size=\"" << size << "\"";
	os << " ";
	os << "maxlength=\"" << maxlength << "\"";
	os << " ";
	os << "value=\"" << value << "\"";
	os << " />" << endl;
}
static void printHTMLFORMFile ( ostream& os, const string& label, const string& manual, const string& name, int size, int maxlength, const string& value )
{
	if ( !label.empty () ) {
		printHTMLFORMLabel ( os, label, manual );
	}
	os << "<input";
	os << " ";
	os << "type=\"file\"";
	os << " ";
	os << "name=\"" << name << "\"";
	os << " ";
	os << "size=\"" << size << "\"";
	os << " ";
	os << "maxlength=\"" << maxlength << "\"";
	os << " ";
	os << "value=\"" << value << "\"";
	os << " />" << endl;
}
void printHTMLFORMSelect ( ostream& os, const string& label, const string& manual, const string& name, const StringVector& value, const string& select, int size, const string& onChangeFunction, const StringVector& theClass )
{
	if ( !label.empty () ) {
		printHTMLFORMLabel ( os, label, manual );
	}
	os << "<select";
	if ( size != 1 ) {
		os << " size=\"" << size << "\"";
	}
	os << " name=\"" << name << "\"";
	if ( !onChangeFunction.empty () ) {
		os << " ";
		os << "onchange=\"" << onChangeFunction << "\"";
	}
	os << ">" << endl;

	bool selectFound = false;
	for ( StringVectorSizeType i = 0 ; i < value.size () ; i++ ) {
		os << "<option";
		if ( value [i] == select ) {
			os << " selected=\"selected\"";
			selectFound = true;
		}
		if ( !theClass.empty () && !theClass [i].empty () ) {
			os << " class=";
			os << "\"";
			os << theClass [i];
			os << "\"";
		}
		os << ">";
		os << value [i];
		os << "</option>";
		os << endl;
	}
	if ( !selectFound ) {		// This could be used to add an extra item to the menu
		os << "<option selected=\"selected\">" << select << "</option>" << endl;
	}
	os << "</select>" << endl;
}
void printHTMLFORMSelectOptgroup ( ostream& os, const string& label, const string& manual, const string& name, const StringVector& value, const StringVector& optgroup, const string& select, int size )
{
	if ( !label.empty () ) {
		printHTMLFORMLabel ( os, label, manual );
	}
	os << "<select";
	if ( size != 1 ) {
		os << " size=\"" << size << "\"";
	}
	os << " name=\"" << name << "\"";
	os << ">" << endl;

	string ogroup;
	bool first = true;
	for ( StringVectorSizeType i = 0 ; i < value.size () ; i++ ) {
		if ( !optgroup [i].empty () ) {
			if ( !first ) os << "</optgroup>" << endl;
			ogroup = optgroup [i]; 
			os << "<optgroup label=\"" << ogroup << "\">" << endl;
			first = false;
		}
		os << "<option";
		if ( value [i] == select ) {
			os << " selected=\"selected\"";
		}
		os << " value=\"" << ogroup + " " + value [i] << "\"";
		os << ">";
		os << value [i];
		os << "</option>";
		os << endl;
	}
	if ( !first ) os << "</optgroup>" << endl;
	os << "</select>" << endl;
}
void printHTMLFORMSelectMultiple ( ostream& os, const string& label, const string& manual, const string& name, const StringVector& value, const StringVector& select, int size, const string& onChangeFunction, const StringVector& theClass )
{
	tableStart ( os );
		tableRowStart ( os );
			if ( !label.empty () ) {
				tableHeaderStart ( os );
					printHTMLFORMLabel ( os, label, manual );
				tableHeaderEnd ( os );
			}
			tableHeaderStart ( os );
				os << "<select";
				os << " multiple=\"multiple\"";
				if ( size != 1 ) {
					os << " size=\"" << size << "\"";
				}
				os << " name=\"" << name << "\"";
				if ( !onChangeFunction.empty () ) {
					os << " ";
					os << "onchange=\"" << onChangeFunction << "\"";
				}
				os << ">" << endl;

				for ( StringVectorSizeType i = 0 ; i < value.size () ; i++ ) {
					os << "<option";
					for ( StringVectorSizeType j = 0 ; j < select.size () ; j++ ) {
						if ( select [j] == value [i] ) {
							os << " selected=\"selected\"";
							break;							// This is to prevent the case where the same value appears twice in the menu. This shouldn't happen 
						}									// but did once due to a bug.
					}
					if ( !theClass.empty () && !theClass [i].empty () ) {
						os << " class=";
						os << "\"";
						os << theClass [i];
						os << "\"";
					}
					os << ">";
					os << value [i];
					os << "</option>";
					os << endl;
				}
				os << "</select>" << endl;
			tableHeaderEnd ( os );
		tableRowEnd ( os );
	tableEnd ( os );
}
void printHTMLFORMTextArea ( ostream& os, const string& label, const string& manual, const string& name, int rows, int cols, const StringVector& value )
{
	tableStart ( os );
		tableRowStart ( os );
			if ( !label.empty () ) {
				tableHeaderStart ( os );
					printHTMLFORMLabel ( os, label, manual );
				tableHeaderEnd ( os );
			}
			tableHeaderStart ( os );
				os << "<textarea";
				os << " ";
				os << "name=\"" << name << "\"";
				os << " ";
				os << "rows=\"" << rows << "\"";
				os << " ";
				os << "cols=\"" << cols << "\"";
				os << ">";
				for ( StringVectorSizeType i = 0 ; i < value.size () ; i++ ) {
					os << value [i] << endl;
				}
				os << "</textarea>" << endl;
			tableHeaderEnd ( os );
		tableRowEnd ( os );
	tableEnd ( os );
}
void printHTMLFORMTextAreaSimple ( ostream& os, const string& label, const string& manual, const string& name, int rows, int cols, const StringVector& value )
{
	if ( !label.empty () ) {
		printHTMLFORMLabel ( os, label, manual );
		os << "<br />" << endl;
	}
	os << "<textarea";
	os << " ";
	os << "name=\"" << name << "\"";
	os << " ";
	os << "rows=\"" << rows << "\"";
	os << " ";
	os << "cols=\"" << cols << "\"";
	os << ">";
	for ( StringVectorSizeType i = 0 ; i < value.size () ; i++ ) {
		os << value [i] << endl;
	}
	os << "</textarea>" << endl;
}
void printHTMLFORMHiddenContainer ( ostream& os, const string& name, const StringVector& value )
{
	if ( value.size () ) {
		for ( StringVectorSizeType i = 0 ; i < value.size () ; i++ ) {
			os << "<input type=\"hidden\" name=\"" << name << "\" value=\"" << value [i] << "\" />" << endl;
		}
	}
}
void printHTMLFORMJavascriptHiddenContainer ( ostream& os, const string& name, const StringVector& value )
{
	if ( value.size () ) {
		os << "\tdocument.writeln ( \"";
		os << "<input";
		os << " ";
		os << "type=\\\"hidden\\\"";
		os << " ";
		os << "name=\\\"" << name << "\\\"";
		os << " ";
		os << "value=\\\"" << value [0] << "\" );" << endl;
		for ( StringVectorSizeType i = 1 ; i < value.size () ; i++ ) {
			os << "\tdocument.writeln ( \"" << value [i] << "\" );" << endl;
		}
		os << "\tdocument.writeln ( \"\\\" />\" );" << endl;
	}
}
void printHTMLFORMJavascriptHiddenContainer2 ( ostream& os, const string& name, const StringVector& value )
{
	for ( StringVectorSizeType i = 0 ; i < value.size () ; i++ ) {
		printHTMLFORMJavascriptHidden2 ( os, name, value [i] );
	}
}
void printHTMLFORMCGI ( ostream& os, const string& name, const string& value )
{
	os << name << "=" << escapeURL ( value ) << "&";
}
void printHTMLFooter ( ostream& os, const string& date )
{
	os << "</div>" << endl;		// centerbody
	printHTMLCopyright ( os, date );
	os << "</div>" << endl;		// content
	os << "</div>" << endl;		// prop
	os << "</body>" << endl;
	os << "</html>" << endl;
}
void printHTMLCopyright ( ostream& os, const string& startYear )
{
	static string copyrightHolders = "The Regents of the University of California";
	static string endYear = genCurrentYear ();
	string years = ( startYear == endYear ) ? endYear : startYear + "-" + endYear;
	os << "<div id=\"footer\">" << endl;
		os << "This program is the confidential and proprietary product of " + copyrightHolders + "." << endl;
		os << "Any unauthorized reproduction or transfer of this program is strictly prohibited." << endl;
		os << "<br />" << endl;
		os << "&#169; Copyright (" + years + ") " + copyrightHolders + "." << endl;
		os << "All rights reserved." << endl;
	os << "</div>" << endl;
}
void printHTMLFORMHidden ( ostream& os, const string& name, double value, int precision )
{
	os << "<input";
	os << " ";
	os << "type=\"hidden\"";
	os << " ";
	os << "name=\"" << name << "\"";
	os << " ";
	os << "value=\"";
	genPrint ( os, value, precision );
	os << "\"";
	os << " />" << endl;
}
void printHTMLFORMHiddenSigFig ( ostream& os, const string& name, double value, int precision )
{
	os << "<input";
	os << " ";
	os << "type=\"hidden\"";
	os << " ";
	os << "name=\"" << name << "\"";
	os << " ";
	os << "value=\"";
	genPrintSigFig ( os, value, precision );
	os << "\"";
	os << " />" << endl;
}
FormItem::FormItem ( const string& label, const string& manual, const string& name, const StringVector& value ) :
	label ( label ),
	manual ( manual ),
	_name ( name ),
	value ( value ),
	shown ( false )
{
}
FormItemCheckbox::FormItemCheckbox ( const string& label, const string& manual, const string& name, bool checked, const string& val, const string& changeFunction ) :
	FormItem ( label, manual, name, StringVector () ),
	checked ( checked ),
	changeFunction ( changeFunction )
{
	value.push_back ( val );
}
void FormItemCheckbox::printHTML ( ostream& os ) const
{
	shown = true;
	printHTMLFORMCheckbox ( os, label, manual, _name, checked, value [0], changeFunction );
}
void FormItemCheckbox::printHTMLHidden ( ostream& os ) const
{
	if ( checked ) printHTMLFORMHidden ( os, _name, value [0] );
}
void FormItemCheckbox::printCGI ( ostream& os ) const
{
	if ( checked ) printHTMLFORMCGI ( os, _name, value [0] );
}
void FormItemCheckbox::printHTMLJavascriptHidden ( ostream& os ) const
{
	if ( checked ) printHTMLFORMJavascriptHidden ( os, _name, value [0] );
}
void FormItemCheckbox::setValue ( const ParameterList* p, const string& n )
{
	checked = false;
	StringVector sv = p->getStringVectorValue ( ( n == "" ) ? _name : n );
	for ( StringVectorSizeType i = 0 ; i < sv.size () ; i++ ) {
		if ( sv [i] == value [0] ) {
			checked = true;
			break;
		}
	}
	//value [0] = p->getStringValue ( ( n == "" ) ? _name : n, value [0] );
}
void FormItemCheckbox::setValue ( bool c )
{
	checked = c;
}
FormItemText::FormItemText ( const string& label, const string& manual, const string& name, int size, int maxlength, const string& val, const string& validateFunction ) :
	FormItem ( label, manual, name, StringVector () ),
	size ( size ),
	maxlength ( maxlength ),
	validateFunction ( validateFunction )
{
	value.push_back ( val );
}
void FormItemText::printHTML ( ostream& os ) const
{
	shown = true;
	printHTMLFORMText ( os, label, manual, _name, size, maxlength, value [0], validateFunction );
}
void FormItemText::printHTMLHidden ( ostream& os ) const
{
	printHTMLFORMHidden ( os, _name, value [0] );
}
void FormItemText::printCGI ( ostream& os ) const
{
	printHTMLFORMCGI ( os, _name, value [0] );
}
void FormItemText::printHTMLJavascriptHidden ( ostream& os ) const
{
	printHTMLFORMJavascriptHidden2 ( os, _name, value [0] );
}
void FormItemText::setValue ( const ParameterList* p, const string& n )
{
	value [0] = p->getStringValue ( ( n == "" ) ? _name : n, value [0] );
}
FormItemPassword::FormItemPassword ( const string& label, const string& manual, const string& name, int size, int maxlength, const string& val ) :
	FormItem ( label, manual, name, StringVector () ),
	size ( size ),
	maxlength ( maxlength )
{
	value.push_back ( val );
}
void FormItemPassword::printHTML ( ostream& os ) const
{
	shown = true;
	printHTMLFORMPassword ( os, label, manual, _name, size, maxlength, value [0] );
}
void FormItemPassword::printHTMLHidden ( ostream& os ) const
{
	printHTMLFORMHidden ( os, _name, value [0] );
}
void FormItemPassword::printCGI ( ostream& os ) const
{
	printHTMLFORMCGI ( os, _name, value [0] );
}
void FormItemPassword::printHTMLJavascriptHidden ( ostream& os ) const
{
	printHTMLFORMJavascriptHidden2 ( os, _name, value [0] );
}
void FormItemPassword::setValue ( const ParameterList* p, const string& n )
{
	value [0] = p->getStringValue ( ( n == "" ) ? _name : n, value [0] );
}
FormItemFile::FormItemFile ( const string& label, const string& manual, const string& name, int size, int maxlength, const string& val ) :
	FormItem ( label, manual, name, StringVector () ),
	size ( size ),
	maxlength ( maxlength )
{
	value.push_back ( val );
}
void FormItemFile::printHTML ( ostream& os ) const
{
	shown = true;
	printHTMLFORMFile ( os, label, manual, _name, size, maxlength, value [0] );
}
void FormItemFile::printHTMLHidden ( ostream& os ) const
{
	printHTMLFORMHidden ( os, _name, value [0] );
}
void FormItemFile::printCGI ( ostream& os ) const
{
	printHTMLFORMCGI ( os, _name, value [0] );
}
void FormItemFile::printHTMLJavascriptHidden ( ostream& os ) const
{
	printHTMLFORMJavascriptHidden2 ( os, _name, value [0] );
}
void FormItemFile::setValue ( const ParameterList* p, const string& n )
{
	value [0] = p->getStringValue ( ( n == "" ) ? _name : n, value [0] );
}
FormItemMultiFile::FormItemMultiFile ( const string& name ) :
	FormItem ( "", "", name, StringVector () )
{
	value.push_back ( "" );
}
void FormItemMultiFile::printHTML ( ostream& os ) const
{
	shown = true;
	os << "<script src=\"" << HTMLDir::instance ().getVirtualHTMLDir () << "/js/multifile.js\"></script>" << endl;
	os << "<input type=\"file\" id=\"upload_file_element\" />" << endl;
	os << "<br />" << endl;
	os << "<div id=\"files_list\"></div>" << endl;
	os << "<script>" << endl;
	os << "var multi_selector = new MultiSelector( document.getElementById( 'files_list' ), '" + _name + "' );" << endl;
	os << "multi_selector.addElement( document.getElementById( 'upload_file_element' ) );" << endl;
	os << "</script>" << endl;
//	os << "<iframe";
//	os << " ";
//	os << "src=\"about:blank\"";
//	os << " ";
//	os << "name=\"upload_status\"";
//	os << " ";
//	os << "style=\"position:absolute;left:-9999px;\"";
//	os << "style=\"position:relative;left:0px;\"";
//	os << " ";
//	os << "width=\"250\"";
//	os << " ";
//	os << "height=\"100\"";
//	os << ">";
//	os << "</iframe>";
//	os << endl;
//display:none; to hide
}
void FormItemMultiFile::printHTMLHidden ( ostream& os ) const
{
	printHTMLFORMHidden ( os, _name, value [0] );
}
void FormItemMultiFile::printCGI ( ostream& os ) const
{
	printHTMLFORMCGI ( os, _name, value [0] );
}
void FormItemMultiFile::printHTMLJavascriptHidden ( ostream& os ) const
{
	printHTMLFORMJavascriptHidden2 ( os, _name, value [0] );
}
FormItemSelect::FormItemSelect ( const string& label, const string& manual, const string& name, const StringVector& value, const StringVector& optgroup, const string& select, int size, const string& changeFunction ) :
	FormItem ( label, manual, name, value ),
	optgroup ( optgroup ),
	select ( select ),
	size ( size ),
	changeFunction ( changeFunction )
{
}
FormItemSelect::FormItemSelect ( const string& label, const string& manual, const string& name, const StringVector& value, const string& select, int size, const string& changeFunction, const StringVector& theClass ) :
	FormItem ( label, manual, name, value ),
	select ( select ),
	size ( size ),
	changeFunction ( changeFunction ),
	theClass ( theClass )
{
}
FormItemSelect::FormItemSelect ( const string& label, const string& manual, const string& name, const char** options, const string& select, int size, const string& changeFunction ) :
	FormItem ( label, manual, name, getOptions ( options ) ),
	select ( select ),
	size ( size ),
	changeFunction ( changeFunction )
{
}
void FormItemSelect::printHTML ( ostream& os ) const
{
	shown = true;
	if ( optgroup.empty () )
		printHTMLFORMSelect ( os, label, manual, _name, value, select, size, changeFunction, theClass );
	else
		printHTMLFORMSelectOptgroup ( os, label, manual, _name, value, optgroup, select, size );
}
void FormItemSelect::printHTMLHidden ( ostream& os ) const
{
	printHTMLFORMHidden ( os, _name, select );
}
void FormItemSelect::printCGI ( ostream& os ) const
{
	printHTMLFORMCGI ( os, _name, select );
}
void FormItemSelect::printHTMLJavascriptHidden ( ostream& os ) const
{
	printHTMLFORMJavascriptHidden2 ( os, _name, select );
}
void FormItemSelect::setValue ( const ParameterList* p, const string& n )
{
	select = p->getStringValue ( ( n == "" ) ? _name : n, select );
}
void FormItemSelect::setValue ( const StringVector& sv )
{
	select = sv [0];
}
void FormItemSelect::setOptions ( const StringVector& sv )
{
	value = sv;
}
FormItemSelectMultiple::FormItemSelectMultiple ( const string& label, const string& manual, const string& name, const StringVector& value, const StringVector& select, int size, const string& changeFunction ) :
	FormItem ( label, manual, name, value ),
	select ( select ),
	size ( size ),
	changeFunction ( changeFunction )
{
}
FormItemSelectMultiple::FormItemSelectMultiple ( const string& label, const string& manual, const string& name, const StringVector& value, const StringVector& select, int size, const StringVector& theClass ) :
	FormItem ( label, manual, name, value ),
	select ( select ),
	size ( size ),
	theClass ( theClass )
{
}
FormItemSelectMultiple::FormItemSelectMultiple ( const string& label, const string& manual, const string& name, const char** options, const StringVector& select, int size ) :
	FormItem ( label, manual, name, getOptions ( options ) ),
	select ( select ),
	size ( size )
{
}
void FormItemSelectMultiple::printHTML ( ostream& os ) const
{
	shown = true;
	printHTMLFORMSelectMultiple ( os, label, manual, _name, value, select, size, changeFunction, theClass );
}
void FormItemSelectMultiple::printHTMLHidden ( ostream& os ) const
{
	printHTMLFORMHiddenContainer ( os, _name, select );
}
void FormItemSelectMultiple::printCGI ( ostream& os ) const
{
	if ( select.size () ) {
		for ( StringVectorSizeType i = 0 ; i < select.size () ; i++ ) {
			printHTMLFORMCGI ( os, _name, select [i] );
		}
	}
}
void FormItemSelectMultiple::printHTMLJavascriptHidden ( ostream& os ) const
{
	printHTMLFORMJavascriptHiddenContainer2 ( os, _name, select );
}
void FormItemSelectMultiple::setValue ( const ParameterList* p, const string& n )
{
	select = p->getUniqueStringVectorValue ( ( n == "" ) ? _name : n );
}
void FormItemSelectMultiple::setValue ( const StringVector& sv )
{
	select = sv;
}
void FormItemSelectMultiple::setOptions ( const ParameterList* p, const string& n )
{
	value = p->getUniqueStringVectorValue ( ( n == "" ) ? _name : n );
}
FormItemTextArea::FormItemTextArea ( const string& label, const string& manual, const string& name, int size, int maxlength, const StringVector& value ) :
	FormItem ( label, manual, name, value ),
	size ( size ),
	maxlength ( maxlength )
{
}
FormItemTextArea::FormItemTextArea ( const string& label, const string& manual, const string& name, int size, int maxlength, const char** v ) :
	FormItem ( label, manual, name, getOptions ( v ) ),
	size ( size ),
	maxlength ( maxlength )
{
}
void FormItemTextArea::printHTML ( ostream& os ) const
{
	shown = true;
	printHTMLFORMTextArea ( os, label, manual, _name, size, maxlength, value );
}
void FormItemTextArea::printHTMLHidden ( ostream& os ) const
{
	printHTMLFORMHiddenVector ( os, _name, value );
}
void FormItemTextArea::printCGI ( ostream& os ) const
{
	if ( value.size () ) {
		os << _name << "=";
		for ( StringVectorSizeType i = 0 ; i < value.size () ; i++ ) {
			os << value [i] << "%0D%0A";
		}
		os << "&";
	}
}
void FormItemTextArea::printHTMLJavascriptHidden ( ostream& os ) const
{
	printHTMLFORMJavascriptHiddenContainer ( os, _name, value );
}
void FormItemTextArea::setValue ( const ParameterList* p, const string& n )
{
	const char* v;
	if ( p->getValue ( ( n == "" ) ? _name : n, v ) ) {
		value.clear ();
		getPostQueryVector ( v, value, '\n' );
	}
}
FormItemTextAreaSimple::FormItemTextAreaSimple ( const string& label, const string& manual, const string& name, int size, int maxlength, const StringVector& value ) :
	FormItemTextArea ( label, manual, name, size, maxlength, value )
{
}
FormItemTextAreaSimple::FormItemTextAreaSimple ( const string& label, const string& manual, const string& name, int size, int maxlength, const char** v ) :
	FormItemTextArea ( label, manual, name, size, maxlength, v )
{
}
void FormItemTextAreaSimple::printHTML ( ostream& os ) const
{
	shown = true;
	printHTMLFORMTextAreaSimple ( os, label, manual, _name, size, maxlength, value );
}
