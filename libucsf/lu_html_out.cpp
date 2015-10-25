/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_html_out.cpp                                               *
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
*  Copyright (1996-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef VIS_C
#include <process.h>
#else
#include <sys/types.h>
#include <unistd.h>
#endif
#include <lg_new.h>
#include <lu_fas_ind.h>
#include <lu_html.h>
#include <lu_html_form.h>
#include <lu_srch_form.h>
#include <lu_table.h>
using std::string;
using std::ostream;
using std::endl;
using std::cout;
using std::showpos;
using std::noshowpos;

namespace {
void ucsfBanner ( ostream& os )
{
	if ( InfoParams::instance ().getBoolValue ( "ucsf_banner" ) ) {
		os << "<div id=\"headngnorm\">" << endl;
			os << "<img src=\"" << HTMLDir::instance ().getVirtualHTMLImagesDir () << "bannerblack.gif\" border=\"0\" usemap=\"#UCSF\" alt=\"\" />" << endl;
			os << "<map name=\"UCSF\" id=\"UCSF\">" << endl;
				os << "<area target=\"_blank\" shape=\"rect\" coords=\"89,6,291,23\" href=\"http://www.ucsf.edu\" alt=\"University of California, San Francisco\" />";
				os << "<area target=\"_blank\" shape=\"rect\" coords=\"303,6,371,23\" href=\"http://www.ucsf.edu/about_ucsf\" alt=\"About UCSF\" />";
				os << "<area target=\"_blank\" shape=\"rect\" coords=\"381,6,458,23\" href=\"http://www.google.com/u/ucsf\" alt=\"Search UCSF\" />";
				os << "<area target=\"_blank\" shape=\"rect\" coords=\"470,6,594,23\" href=\"http://www.ucsfhealth.org/\" alt=\"UCSF  medical center\" />";
			os << "</map>" << endl;
		os << "</div>" << endl;
	}
}
}

void init_html ( ostream& os, const string& title, const string& manual, bool reset )
{
	static bool called = false;
	if ( reset ) called = false;
	if ( !called ) {
		os << "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\"" << endl;
		os << "\"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">" << endl;
		os << "<html xmlns=\"http://www.w3.org/1999/xhtml\" xml:lang=\"en\" lang=\"en\">" << endl;
		os << "<head>" << endl;
		os << "<title>" << title << "</title>" << endl;
		os << "<link href=\"" << HTMLDir::instance ().getVirtualHTMLDir () << "/pform.css\" media=\"screen\" rel=\"stylesheet\" type=\"text/css\" />" << endl;
		os << "<link href=\"" << HTMLDir::instance ().getVirtualHTMLDir () << "/print.css\" media=\"print\" rel=\"stylesheet\" type=\"text/css\" />" << endl;
		os << "</head>" << endl;
		os << "<body>" << endl;
		os << "<div id=\"prop\">" << endl;
		os << "<div id=\"content\">" << endl;
		ucsfBanner ( os );
		OldResourceLinkTable ().printHTML ( os );
		os << "<div id=\"centerbody\">" << endl;
		if ( !title.empty () ) {
			if ( manual.empty () ) {
				os << "<div class=\"form_name\">" << title << "</div>" << endl;
			}
			else {
				os << "<div class=\"form_name\">" << endl;
				printHTMLFORMLabel ( os, title, manual );
				os << "</div>" << endl;
			}
		}
		called = true;
	}
	else {
		os << "<div id=\"error_header\">" << title << "</div>" << endl;
	}
}
static bool stopFlag = false;
void init_html_premature_stop ( const string& programName, bool stopButton, const StringVector& filesToDelete )
{
	string title = programName + " Search Results";
	init_html ( cout, title );

	if ( stopButton ) {
		if ( InfoParams::instance ().getBoolValue ( "kill_button", true ) ) {
			printHTMLFORMStart ( cout, "post", "kill" );

			cout << "<div id=\"abort\">" << endl;
			printHTMLFORMSubmit ( cout, "Abort Search", "FF0000" );
			cout << "<br />" << endl;
			cout << "</div>" << endl;

			printHTMLFORMHidden ( cout, "pid", getpid () );
			printHTMLFORMHiddenContainer ( cout, "delete", filesToDelete );

			cout << "</form>" << endl << endl;
			stopFlag = true;
		}
		else {
			cout << "<p>" << endl;
			cout << "Press stop on your browser if you wish to abort this ";
			cout << programName;
			cout << " search prematurely.";
			cout << "</p>" << endl;
		}
	}
}
void printAbortButton ( ostream& os, const string& searchJobKey )
{
	printHTMLFORMStart ( os, "post", "jobStatus" );

	os << "<div id=\"abort\">" << endl;
	printHTMLFORMSubmit ( os, "Abort Search", "FF0000" );
	os << "<br />" << endl;
	os << "</div>" << endl;

	printHTMLFORMHidden ( os, "search_key", searchJobKey );
	printHTMLFORMHidden ( os, "abort", 1 );

	os << "</form>" << endl << endl;
}
void printAbortFunctions ( ostream& os )
{
	if ( stopFlag ) {
		startJavascript ( os );
			os << "function getObject( id ) {" << endl;
			os << "\t" << "if ( document.all ) {" << endl;
			os << "\t\t" << "return document.all[id];" << endl;
			os << "\t" << "}" << endl;
			os << "\t" << "else {" << endl;
			os << "\t\t" << "return document.getElementById( id );" << endl;
			os << "\t" << "}" << endl;
			os << "}" << endl;

			os << "var section = getObject( 'abort' );" << endl;
			os << "section.style.display = 'none';" << endl;
		endJavascript ( os );
	}
}
void htmlPrintProteinSequence ( ostream& os, const string& protein, const IntVector& cleavageIndicies, const CharVector& coverageMap, bool compressedDisplay )
{
	static int aaPerRow = 80;
	static int aaPerRowPlusOne = aaPerRow+1;
	static int gapSpacing = 10;
	int i, j, k;
	bool showCleavages = !cleavageIndicies.empty ();
	int protLen = protein.length ();

	os << "<table cellspacing=\"3\">" << endl;
	for ( i = 0, k = 0 ; i < protLen ; i += aaPerRow ) {
		tableRowStart ( os );
		if ( compressedDisplay ) {
			tableDataStart (os );
				os << "<tt>" << i + 1 << "</tt>" << endl;
			tableDataEnd (os );
		}
		else {
			for ( j = 1 ; j < aaPerRowPlusOne ; j += gapSpacing ) {
				if ( i + j <= protLen ) {
					tableDataStart (os );
						os << "<tt>" << i + j << "</tt>" << endl;
					tableDataEnd (os );
				}
			}
			tableEmptyCell ( os );
			tableRowEnd ( os );
			tableRowStart ( os );
		}
		for ( j = 0 ; j < aaPerRow ; j++ ) {
			int posn = i + j;
			if ( posn % gapSpacing == 0 ) {
				tableDataStart (os );
				os << "<tt>";
			}
			if ( posn < protLen ) {
			if ( showCleavages && posn == cleavageIndicies [k] ) os << "<u>";
				if ( !coverageMap.empty () && coverageMap [posn] == 1 ) os << "<b><font color=\"#FF0000\">";
				if ( !coverageMap.empty () && coverageMap [posn] == 2 ) os << "<b><font color=\"#00AF00\">";
				os << protein [posn];
				if ( !coverageMap.empty () && coverageMap [posn] ) os << "</font></b>";
				if ( showCleavages && posn == cleavageIndicies [k] ) {
					os << "</u>";
					k++;
				}
			}
			if ( ( posn + 1 ) % gapSpacing == 0 ) {
				os << "</tt>" << endl;
				tableDataEnd (os );
			}
		}
		/*								This prints a number after the sequence
		if ( compressedDisplay ) {
			os << "</tt>" << endl;
			tableDataEnd (os );
			tableDataStart (os );
			os << "<tt>";
				os << ( ( i + j < protLen ) ? i + j : protLen );
			os << "</tt>" << endl;
			tableDataEnd (os );
			tableEmptyCell ( os );
			tableRowEnd ( os );
		}
		else
		*/
		tableEmptyCell ( os );
		tableRowEnd ( os );
	}
	os << "</table>" << endl;
	os << "<p />" << endl << endl;
}
void printHTMLLink ( ostream& os, const string& outputPath, const string& filename, const string& label )
{
	os << "<a href=\"";
	os << outputPath;
	os << "/";
	os << filename;
	os << "\">";
	os << label;
	os << " (" << filename << ")";
	os << "</a>";
	os << "<br />";
	os << endl;
}
void print_charge_superscript ( ostream& os, int charge )
{
	if ( charge > 1 || charge < 0 ) {
		os << "<sup>";
		os << showpos;
		os << charge;
		os << "</sup>";
	}
	else if ( charge == 1 ) return;
	else if ( charge == 0 ) os << "<sup>0</sup>";
	os << noshowpos;
}
void printNumHits ( ostream& os, const string& programName, int numHits, int maxReportedHits )
{
	os << programName << " search selects <b>" << numHits << "</b>";
	if ( numHits == 1 )	os << " entry";
	else				os << " entries";

	if ( numHits > maxReportedHits )
		os << " (results displayed for top <b>" << maxReportedHits << "</b> matches)";
	os << "." << endl;
	os << "<p />" << endl << endl;
}
void startJavascript ( ostream& os, const string& version )
{
	os << "<script type=\"text/javascript\" language=\"Javascript";
	if ( !version.empty () ) os << version;
	os << "\">" << endl;
	os << "<!--" << endl;
}
void endJavascript ( ostream& os )
{
	os << "//-->" << endl;
	os << "</script>" << endl << endl;
}
void htmlPreformattedTextFile ( ostream& os, const string& fileName )
{
	int numLines;
	char* info = getFileInfo ( fileName, '\n', 1, false, &numLines );

	os << "<pre>";
	for ( int i = 0 ; i < numLines ; i++ ) {
		char* line = ( i == 0 ) ? strtok ( info, "\n" ) : strtok ( NULL, "\n" );
		os << line << endl;
	}
	os << "</pre>";

	delete [] info;
}
HTMLErrorHandler::HTMLErrorHandler ()
{
	ErrorHandler::registration ( this );
}
void HTMLErrorHandler::errorDisplay ( const string& errorString, bool endProgram )
{
	init_html ( cout, "Error Message" );

	cout << "<font color=\"#FF0000\" size=\"+2\"><b>" << endl;
	messageDisplay ( errorString );
	cout << "</b></font>" << endl;

	if ( endProgram ) {
		cout << "</div>" << endl;		// centerbody
		cout << "<div id=\"footer\">" << endl;
		cout << "UCSF Mass Spectrometry Facility" << endl;
		cout << "</div>" << endl;
		cout << "</div>" << endl;		// content
		cout << "</div>" << endl;		// prop
		cout << "</body>" << endl;
		cout << "</html>" << endl;
	}
}
void HTMLErrorHandler::messageDisplay ( const string& messageString )
{
	char* tempMessage = gen_new_string ( messageString.c_str () );

	for ( int i = 0 ; ; i++ ) {
		char* messageLine = ( i == 0 ) ? strtok ( tempMessage, "\n" ) : strtok ( NULL, "\n" );
		if ( messageLine != NULL ) {
			cout << messageLine << "<br />" << endl;
		}
		else break;
	}
	delete [] tempMessage;
}
void basicJavascriptFunctions ( ostream& os )
{
	startJavascript ( os );

	os << "function showTableColumn(tableID, col, show) {" << endl;
	os << "\t" << "var stl;" << endl;
	os << "\t" << "if (show) stl = 'block'" << endl;
	os << "\t" << "else      stl = 'none';" << endl;
	os << "\t" << "var tbl  = document.getElementById(tableID);" << endl;
	os << "\t" << "var rows = tbl.getElementsByTagName('tr');" << endl;
	os << "\t" << "for (var row=0; row<rows.length;row++) {" << endl;
	os << "\t\t" << "var cells = rows[row].getElementsByTagName('td')" << endl;
	os << "\t\t" << "if ( cells.length ) cells[col].style.display=stl;" << endl;
	os << "\t\t" << "var cells2 = rows[row].getElementsByTagName('th')" << endl;
	os << "\t\t" << "if ( cells2.length ) cells2[col].style.display=stl;" << endl;
	os << "\t" << "}" << endl;
	os << "}" << endl;

	os << "function hidediv(id) {" << endl;
	os << "\t" << "if (document.getElementById) {" << endl;			// DOM3 = IE5, NS6
	os << "\t\t" << "document.getElementById(id).style.display = 'none';" << endl;
	os << "\t" << "}" << endl;
	os << "\t" << "else {" << endl;
	os << "\t\t" << "if (document.layers) {" << endl;				// Netscape 4
	os << "\t\t\t" << "document.id.display = 'none';" << endl;
	os << "\t\t" << "}" << endl;
	os << "\t\t" << "else {" << endl;								// IE 4
	os << "\t\t\t" << "document.all.id.style.display = 'none';" << endl;
	os << "\t\t" << "}" << endl;
	os << "\t" << "}" << endl;
	os << "}" << endl;

	os << "function showdiv(id) {" << endl;
	os << "\t" << "if (document.getElementById) {" << endl;			// DOM3 = IE5, NS6
	os << "\t\t" << "document.getElementById(id).style.display = 'block';" << endl;
	os << "\t" << "}" << endl;
	os << "\t" << "else {" << endl;
	os << "\t\t" << "if (document.layers) {" << endl;				// Netscape 4
	os << "\t\t\t" << "document.id.display = 'block';" << endl;
	os << "\t\t" << "}" << endl;
	os << "\t\t" << "else {" << endl;								// IE 4
	os << "\t\t\t" << "document.all.id.style.display = 'block';" << endl;
	os << "\t\t" << "}" << endl;
	os << "\t" << "}" << endl;
	os << "}" << endl;

	os << "function getSelectValue(name) {" << endl;
	os << "\t" << "var val;" << endl;
	os << "\t" << "if ( document.all ) {" << endl;
	os << "\t\t" << "val=name.options[name.selectedIndex].text;" << endl;
	os << "\t" << "}" << endl;
	os << "\t" << "else {" << endl;
	os << "\t\t" << "val=name.value;" << endl;
	os << "\t" << "}" << endl;
	os << "\t" << "return val;" << endl;
	os << "}" << endl;

	os << "function multiSelectSet ( field, value, flag ) {" << endl;
	os << "\t" << "for ( var i = 0; i < field.options.length ; i++ ) {" << endl;
	os << "\t\t" << "var opt = field.options [i];" << endl;
	os << "\t\t" << "if ( opt.text == value ) {" << endl;
	os << "\t\t\t" << "opt.selected = flag;" << endl;
	os << "\t\t" << "}" << endl;
	os << "\t" << "}" << endl;
	os << "}" << endl;

	os << "function addOptionToSelect( item, opt ) {" << endl;
	os << "\t" << "var op = document.createElement(\"option\");" << endl;
	os << "\t" << "op.text = opt" << endl;
	os << "\t" << "item.add(op);" << endl;
	os << "}" << endl;

	os << "function removeOptionsFromSelect( item ) {" << endl;
	os << "\t" << "var len1 = item.length;" << endl;
	os << "\t" << "for ( var i = len1 - 1 ; i >= 0 ; i-- ) {" << endl;
	os << "\t\t" << "item.remove(i);" << endl;
	os << "\t" << "}" << endl;
	os << "}" << endl;

	endJavascript ( os );
}
void updateJavascriptList ( const string& jsFilename, const StringVector& nameList, const string& selectedName )
{
	int numEntries;
	char* info = getFileInfo ( jsFilename, '\n', 1, false, &numEntries );

	GenOFStream ost ( jsFilename );

	bool print = true;
	for ( int i = 0 ; i < numEntries ; i++ ) {
		string line = ( i == 0 ) ? strtok ( info, "\n" ) : strtok ( NULL, "\n" );
		if ( print ) ost << line << '\n';
		if ( line [0] == '{' ) {
			print = false;
			bool selected = false;
			for ( StringVectorSizeType j = 0 ; j < nameList.size () ; j++ ) {
				if ( selected == false ) {
					if ( nameList [j] == selectedName ) {
						ost << "\tdocument.writeln ( \"<option selected=\\\"selected\\\">" << nameList [j] << "</option>\" );\n";
						selected = true;
						continue;
					}
				}
				ost << "\tdocument.writeln ( \"<option>" << nameList [j] << "</option>\" );\n";
			}
		}
		if ( line [0] == '}' ) {
			print = true;
			ost << line << '\n';
		}
	}
	delete [] info;
}

int ExpandableJavascriptBlock::blockNumber = 1;
string ExpandableJavascriptBlock::idNamePrefix = "ejb_";
bool ExpandableJavascriptBlock::functionsPrinted = false;

ExpandableJavascriptBlock::ExpandableJavascriptBlock ( const string& title, bool open ) :
	title ( title ),
	idName ( idNamePrefix + gen_itoa ( blockNumber++ ) ),
	linkName ( idName + "_link" ),
	defaultToggle ( open ? "+" : "-" )
{
}
void ExpandableJavascriptBlock::printHeader ( ostream& os )
{
	if ( !functionsPrinted ) {
		printFunctions ( os );
		functionsPrinted = true;
	}
	os << "[";
	os << "<a title=\"show/hide\"";
	os << " ";
	os << "id=\"" + linkName + "\"";
	os << " ";
	os << "href=\"javascript:void(0);\"";
	os << " ";
	os << "onclick=\"toggle(this, '" + idName + "');\"";
	os << ">";
	os << defaultToggle;
	os << "</a>";
	os << "]";
	os << " ";
	os << "<b>" << title << "</b>";
	os << endl;
	os << "<div id=\"" + idName + "\" style=\"display:none\">" << endl;
}
void ExpandableJavascriptBlock::printFooter ( ostream& os )
{
	os << "</div>" << endl;
	startJavascript ( os );
	os << "toggle(getObject('" + linkName + "'), '" + idName + "');" << endl;
	endJavascript ( os );
}
void ExpandableJavascriptBlock::printFunctions ( ostream& os )
{
	startJavascript ( os );
	os << "function onLoadFunction () {" << endl;
	os << "\t" << "return;" << endl;
	os << "}" << endl;
	os << "var ie4 = false;" << endl; 
	os << "if ( document.all ) {" << endl;
	os << "\t" << "ie4 = true;" << endl;
	os << "}" << endl;
	os << "function getObject ( id ){" << endl;
	os << "\t" << "if ( ie4 ){" << endl;
	os << "\t\t" << "return document.all[id];" << endl;
	os << "\t" << "}" << endl;
	os << "\t" << "else {" << endl;
	os << "\t\t" << "return document.getElementById ( id );" << endl;
	os << "\t" << "}" << endl;
	os << "}" << endl;
	os << "function toggle ( link, div_id ) {" << endl;
	os << "\t" << "var link_text = link.innerHTML;" << endl;
	os << "\t" << "var d = getObject ( div_id );" << endl;
	os << "\t" << "if ( link_text == '+' ) {" << endl;
	os << "\t\t" << "link.innerHTML = '&ndash;';" << endl;
	os << "\t\t" << "d.style.display = 'block';" << endl;
	os << "\t" << "}" << endl;
	os << "\t" << "else {" << endl;
	os << "\t\t" << "link.innerHTML = '+';" << endl;
	os << "\t\t" << "d.style.display = 'none';" << endl;
	os << "\t" << "}" << endl;
	os << "}" << endl;
	endJavascript ( os );
}
string UpdatingJavascriptMessage::outstandingSpanID;
UpdatingJavascriptMessage::UpdatingJavascriptMessage () :
	index ( 0 )
{
	static int instanceNumber = 0;
	instanceNumber++;
	id = "ujm_" + gen_itoa ( instanceNumber ) + "_";
}
void UpdatingJavascriptMessage::writeMessage ( ostream& os, int num )
{
	writeMessage ( os, gen_itoa ( num ) );
}
void UpdatingJavascriptMessage::startWriteMessage ( ostream& os )
{
	deletePreviousMessage ( os );
	index++;
	spanID = id + gen_itoa ( index );
	os << "<span id=\"" << spanID << "\">" << endl;
}
void UpdatingJavascriptMessage::addMessage ( ostream& os, const string& message ) const
{
	os << message << "<br />" << endl;
}
void UpdatingJavascriptMessage::endWriteMessage ( ostream& os )
{
	os << "</span>";
	os << endl;
	outstandingSpanID = spanID;
}
void UpdatingJavascriptMessage::writeMessage ( ostream& os, const string& message )
{
	deletePreviousMessage ( os );
	index++;
	spanID = id + gen_itoa ( index );
	os << "<span id=\"" << spanID << "\">";
	os << message;
	os << "</span>";
	os << endl;
	outstandingSpanID = spanID;
}
void UpdatingJavascriptMessage::deletePreviousMessage ( ostream& os ) const
{
	if ( index ) {
		os << "<script type=\"text/javascript\">";
		os << "document.getElementById(\"" << spanID << "\").style.display='none';";
		os << "</script>";
		outstandingSpanID = "";
	}
}
void UpdatingJavascriptMessage::deleteOutstandingSpanID ( ostream& os )
{
	if ( !outstandingSpanID.empty () ) {
		os << "<script type=\"text/javascript\">";
		os << "document.getElementById(\"" << outstandingSpanID << "\").style.display='none';";
		os << "</script>";
		outstandingSpanID = "";
	}
}
int CheckboxArraySettingJavascript::num = 0;
CheckboxArraySettingJavascript::CheckboxArraySettingJavascript ( ostream& os, const string& str, const string& buttonLabel, const BoolDeque& flagList )
{
	num++;
	printFunctions ( os, flagList );
	os << "<input type=\"button\" value=\"" << buttonLabel << " On\" onclick=\"" << "setCheckboxes" << num << "( " << str << ", true )\" />" << endl;
	os << "<input type=\"button\" value=\"" << buttonLabel << " Off\" onclick=\"" << "setCheckboxes" << num << "( " << str << ", false )\" />" << endl;
}
void CheckboxArraySettingJavascript::printFunctions ( ostream& os, const BoolDeque& flagList )
{
	startJavascript ( os );
	os << "function setCheckboxes" << num << " ( thisField, flag ) {" << endl;
	for ( int i = 0 ; i < flagList.size () ; i++ ) {
		os << "\t" << "thisField[" << i << "].checked = flag ? " << flagList [i] << " : " << !flagList [i] << ";" << endl;
	}
	os << "}" << endl;
	endJavascript ( os );
}
bool CheckboxSettingJavascript::functionsPrinted = false;
CheckboxSettingJavascript::CheckboxSettingJavascript ( const string& str )
{
	names.push_back ( str );
}
CheckboxSettingJavascript::CheckboxSettingJavascript ( const StringVector& names ) :
	names ( names )
{
}
void CheckboxSettingJavascript::print ( ostream& os ) const
{
	if ( !functionsPrinted ) {
		printFunctions ( os );
		functionsPrinted = true;
	}
	printButton ( os, "All On", true );
	printButton ( os, "All Off", false );
}
void CheckboxSettingJavascript::printButton ( ostream& os, const string& buttonLabel, bool flag ) const
{
	os << "<input";
	os << " ";
	os << "type=\"button\"";
	os << " ";
	os << "value=";
	os << "\"";
	os << buttonLabel;
	os << "\"";
	os << " ";
	os << "onclick=";
	os << "\"";
	for ( StringVectorSizeType i = 0 ; i < names.size () ; i++ ) {
		os << "setCheckboxes( ";
		os << names [i];
		os << ", ";
		if ( flag ) os << "true";
		else os << "false";
		os << " )";
		if ( i != names.size () - 1 ) os << "; ";
	}
	os << "\"";
	os << " />";
	os << endl;
}
void CheckboxSettingJavascript::printFunctions ( ostream& os )
{
	startJavascript ( os );
	os << "function setCheckboxes ( thisField, flag ) {" << endl;
	os << "\t" << "if ( thisField.length == undefined ) {" << endl;
	os << "\t\t" << "thisField.checked = flag;" << endl;
	os << "\t" << "}" << endl;
	os << "\t" << "else {" << endl;
	os << "\t\t" << "for ( i=0; i < thisField.length; i++ ) {" << endl;
	os << "\t\t\t" << "thisField[i].checked = flag;" << endl;
	os << "\t\t" << "}" << endl;
	os << "\t" << "}" << endl;
	os << "}" << endl;
	endJavascript ( os );
}
void refreshJavascript ( ostream& os, int timeout, const string& url, bool manualRedirectMessage, bool automaticUpdateMessage )
{
	if ( manualRedirectMessage ) {
		os << "<p>" << endl;
		os << "You should be redirected automatically. If not, click ";
		os << "<a href=\"";
		os << url;
		os << "\">";
		os << "here";
		os << "</a>.";
		os << endl;
		os << "</p>";
	}
	else if ( automaticUpdateMessage ) {
		os << "<p>" << endl;
		os << "If this page does not automatically update, click ";
		os << "<a href=\"";
		os << url;
		os << "\">";
		os << "here";
		os << "</a> to refresh.";
		os << endl;
		os << "</p>";
	}
	startJavascript ( os );

	os << "var sURL = ";
	if ( url.empty () ) os << "unescape(window.location.pathname);" << endl;
	else os << "\"" << url << "\";";
	os << endl;

	os << "function refresh()" << endl;
	os << "{" << endl;
    os << "\twindow.location.href = sURL;" << endl;
	os << "}" << endl;
	endJavascript ( os );

	startJavascript ( os, "1.1" );
	os << "function refresh ()" << endl;
	os << "{" << endl;
    os << "\twindow.location.replace( sURL );" << endl;
	os << "}" << endl;
	endJavascript ( os );

	startJavascript ( os, "1.2" );
	os << "function refresh ()" << endl;
	os << "{" << endl;
    if ( url.empty () ) os << "\twindow.location.reload ( true );" << endl;
	else os << "\twindow.location.replace( sURL );" << endl;
	os << "}" << endl;
	endJavascript ( os );

	startJavascript ( os );
    os << "setTimeout( \"refresh ()\", " << timeout << " );" << endl;
	endJavascript ( os );
}
void back2Javascript ( ostream& os )
{
	startJavascript ( os );
	os << "history.go(-2);" << endl;
	endJavascript ( os );
}
void deleteCookie ( ostream& os, const string& name )
{
	startJavascript ( os );

	os << "var d = new Date ();" << endl;
	os << "d.setTime ( d.getTime () - 1 );" << endl;
	os << "document.cookie = \'" << name << "=; expires = \' + d.toGMTString ();" << endl;

	endJavascript ( os );
}
bool setCookie ( ostream& os, const string& name, const string& value, bool persistant )
{
	if ( name.length () + value.length () > 4096 ) {	// Cookie too long
		return false;
	}
	else {
		startJavascript ( os );
		if ( persistant ) {
			os << "var nextYear = new Date ();" << endl;
			os << "nextYear.setFullYear ( nextYear.getFullYear () + 1 );" << endl;
			os << "document.cookie = \'" << name << "=" << value << "; expires = \' + nextYear.toGMTString ();" << endl;
		}
		else {
			os << "document.cookie = \'" << name << "=" << value << "\';" << endl;
		}
		endJavascript ( os );
		return true;
	}
}
void setValueFromCookie ( ostream& os, const string& name )
{
	startJavascript ( os );
		os << "var allcookies = document.cookie;" << endl;
		os << "var pos = allcookies.indexOf ( \"" << name << "=\" );" << endl;
		os << "if ( pos != -1 ) {" << endl;
		os << "\t" << "var start = pos + " << name.length () + 1 << ";" << endl;
		os << "\t" << "var end = allcookies.indexOf ( \";\" );" << endl;
		os << "\t" << "if (end == -1) end = allcookies.length;" << endl;
		os << "\t" << "document.forms[0]." << name << ".value = allcookies.substring(start,end);" << endl;
		os << "}" << endl;
	endJavascript ( os );
}
string getElementalFormulaHTML ( const string& formula )
{
	string s;
	bool sub = false;
	for ( StringSizeType i = 0 ; i < formula.length () ; i++ ) {
		char c = formula [i];
		if ( c != ' ' ) {
			if ( isdigit ( c ) ) {
				if ( !sub ) {
					s += "<sub>";
					sub = true;
				}
			}
			else {
				if ( sub ) {
					s += "</sub>";
					sub = false;
				}
			}
			s += c;
		}
	}
	if ( sub ) s += "</sub>";
	return s;
}
void compressedHTMLStyle ( ostream& os )
{
	static bool printed = false;
	if ( !printed ) {
		os << "<style>" << endl;
		os << ".compressed_style" << endl;
		os << "{" << endl;
		os << "background: #fff url(" << HTMLDir::instance ().getVirtualHTMLImagesDir () << "compressed.gif) no-repeat;" << endl;
		os << "padding-left: 18px;" << endl;
		os << "}" << endl;
		os << "</style>" << endl;
		printed = true;
	}
}
