/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_scsel_form.cpp                                             *
*                                                                             *
*  Created    : December 10th 2004                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2004-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef MYSQL_DATABASE
#include <algorithm>
#include <lu_scsel_form.h>
#include <lu_html.h>
#include <lu_cookie.h>
#include <lu_table.h>
#include <lu_html_form.h>
#include <lu_param_list.h>
#include <lu_srch_form.h>
#include <ld_init.h>
using std::sort;
using std::vector;
using std::string;
using std::ostream;
using std::copy;
using std::back_inserter;
using std::endl;

static string FORM	= "form";
static string USER	= "user";
static string PROJECT_NAME = "project_name";
static string RESULTS_FILE = "results_file";
static string COMPARE_FILES = "compare_files";
static string CALIBRATE		= "calibrate";
static string SELECT_RESULTS = "select_results";

SCSelProjectForm::SCSelProjectForm ( const ParameterList* paramList, bool setCookie ) :
	pdff ( new ProjectDateFilterForm ( paramList ) ),
	setCookie ( setCookie )
{
	user = paramList->getStringValue ( USER, user );
	VectorConstParameterListPtr pv;
	pv.push_back ( paramList );
	create ( pv );
}
void SCSelProjectForm::setOptions ()
{
	projects = MySQLPPSDDBase::instance ().getProjectList ( user, false, pdff->getProjectDateFilter (), true );
	sort ( projects.begin (), projects.end (), genStrcasecmpAscending () );
	for ( StringVectorSizeType i = 0 ; i < projects.size () ; i++ ) {
		bool comp = projects [i][projects [i].length ()-1] == '$';
		theClass.push_back ( comp ? "compressed_style" : "" );
		if ( comp ) projects [i] = projects [i].substr ( 0, projects [i].length () - 1 );
	}
}
void SCSelProjectForm::createItems ()
{
	formItemMap [USER] = new FormItemText ( "User", "", USER, 20, 60, user );
	formItemMap [FORM] = new FormItemText ( "", "", FORM, 20, 60, "search_compare" );
	formItemMap [SELECT_RESULTS] = new FormItemCheckbox ( "", "", SELECT_RESULTS, true );

	formItemMap [PROJECT_NAME] = new FormItemSelectMultiple ( "Projects", "", PROJECT_NAME, projects, StringVector (), 10, theClass );
}
void SCSelProjectForm::printHTML ( ostream& os )
{
	init_html ( os, "Search Compare Select Project(s)" );
	compressedHTMLStyle ( os );
	printJavascriptFunctions ( os );
	if ( setCookie ) setUserCookie ( os, user );
	printHTMLFORMStart ( os, "post", "msform" );
	tableStart ( os, true, "center" );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				if ( projects.empty () ) {
					printHTMLFORMLabel ( os, "No projects present", "" );
				}
				else {
					formItemMap [PROJECT_NAME]->printHTML ( os );
				}
			tableHeaderEnd ( os );
		tableRowEnd ( os );

		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				printHTMLFORMSubmit ( os, "Select" );
			tableHeaderEnd ( os );
		tableRowEnd ( os );

		pdff->printHTML ( os, true );
	tableEnd ( os );
	showHiddenItems ( os );
	printHTMLFORMEnd ( os );
	printHTMLFooter ( os, "2007" );
}
void SCSelProjectForm::printJavascriptFunctions ( ostream& os )
{
	basicJavascriptFunctions ( os );
	startJavascript ( os );
	pdff->printJavascriptFunctions ( os );
	endJavascript ( os );
}
SCSelResultsForm::SCSelResultsForm ( const ParameterList* paramList )
{
	user = paramList->getStringValue ( USER, user );
	UserInfo* userInfo = MySQLPPSDDBase::instance ().getUserInfo ( user );
	guest = userInfo->getIsGuest ();
	projects = paramList->getStringVectorValue ( PROJECT_NAME );
	form = paramList->getStringValue ( FORM, form );
	VectorConstParameterListPtr pv;
	pv.push_back ( paramList );
	create ( pv );
}
void SCSelResultsForm::setOptions ()
{
	for ( StringVectorSizeType i = 0 ; i < projects.size () ; i++ ) {
		string p = projects [i];
		string::size_type ind = p.find ( "/" );
		StringVector temp;
		if ( ind == string::npos ) {
			temp = MySQLPPSDDBase::instance ().getResultsList ( user, p );
		}
		else {
			temp = MySQLPPSDDBase::instance ().getResultsList ( p.substr ( 0, ind ), p.substr ( ind + 1 ) );
		}
		for ( StringVectorSizeType j = 0 ; j < temp.size () ; j++ ) {
			temp [j] = p + "/" + temp [j];
		}
		copy ( temp.begin (), temp.end (), back_inserter ( results ) );
	}
	sort ( results.begin (), results.end (), genStrcasecmpAscending () );
}
void SCSelResultsForm::createItems ()
{
	if ( !results.empty () ) {
		formItemMap [USER] = new FormItemText ( "User", "", USER, 20, 60, user );
		formItemMap [FORM] = new FormItemText ( "", "", FORM, 20, 60, form );
		formItemMap [RESULTS_FILE] = new FormItemSelect ( "Results File", "", RESULTS_FILE, results, results [0] );
		int resultsSize = genMin ( 20, (int)results.size () );
		formItemMap [COMPARE_FILES] = new FormItemSelectMultiple ( "", "", COMPARE_FILES, results, StringVector (), resultsSize );
		if ( !guest ) formItemMap [CALIBRATE] = new FormItemCheckbox ( "Calibrate", "", CALIBRATE, false );
	}
}
void SCSelResultsForm::printHTML ( ostream& os )
{
	if ( results.empty () ) {
		ErrorHandler::genError ()->error ( "No results selected.\n" );
	}
	else {
		init_html ( os, "Search Compare Select Results" );
		printHTMLFORMStart ( os, "post", "msform" );
		tableStart ( os, true, "center" );
			tableRowStart ( os );
				tableHeaderStart ( os, "", "left", true );
					formItemMap [RESULTS_FILE]->printHTML ( os );
					if ( !guest ) formItemMap [CALIBRATE]->printHTML ( os );
					if ( results.size () > 1 ) {
						os << "<br />" << endl;
						ExpandableJavascriptBlock ejb ( "Compare Files" );
						ejb.printHeader ( os );
						formItemMap [COMPARE_FILES]->printHTML ( os );
						ejb.printFooter ( os );
					}
				tableHeaderEnd ( os );
			tableRowEnd ( os );

			tableRowStart ( os );
				tableHeaderStart ( os, "", "center", true );
					printHTMLFORMSubmit ( os, "Run Search Compare" );
				tableHeaderEnd ( os );
			tableRowEnd ( os );
		tableEnd ( os);
		showHiddenItems ( os );
		printHTMLFORMEnd ( os );
		printHTMLFooter ( os, "2007" );
	}
}
#endif
