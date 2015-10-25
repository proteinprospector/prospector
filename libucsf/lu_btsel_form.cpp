/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_btsel_form.cpp                                             *
*                                                                             *
*  Created    : March 28th 2007                                               *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2007-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef MYSQL_DATABASE
#include <algorithm>
#include <lgen_file.h>
#include <lu_cookie.h>
#include <lu_getfil.h>
#include <lu_html.h>
#include <lu_btsel_form.h>
#include <lu_table.h>
#include <lu_param_list.h>
#include <lu_html_form.h>
#include <lu_srch_form.h>
#include <ld_init.h>

using std::sort;
using std::vector;
using std::string;
using std::ostream;
using std::copy;
using std::endl;
using std::back_inserter;

static string FORM = "form";
static string USER = "user";
static string PROJECT_NAME = "project_name";
static string RESULTS_FILE = "results_file";
static string SELECT_RESULTS = "select_results";
static string FILTER_PROJECT_LIST = "filter_project_list";

BTSelProjectForm::BTSelProjectForm ( const ParameterList* paramList, bool setCookie ) :
	pdff ( new ProjectDateFilterForm ( paramList ) ),
	setCookie ( setCookie )
{
	user = paramList->getStringValue ( USER, user );
	VectorConstParameterListPtr pv;
	pv.push_back ( paramList );
	create ( pv );
}
void BTSelProjectForm::setOptions ()
{
	bool instDirFlag = genFileExists ( MsparamsDir::instance ().getParamPath ( "inst_dir.txt" ) );
	bool repositoryFlag = genFileExists ( MsparamsDir::instance ().getParamPath ( "repository.xml" ) );
	if ( repositoryFlag || instDirFlag ) {
		projects.push_back ( "Make Project...." );
		projects.push_back ( "Import Project...." );
	}
	StringVector temp = MySQLPPSDDBase::instance ().getProjectList ( user, true, pdff->getProjectDateFilter (), true );
	if ( !temp.empty () ) {
		sort ( temp.begin (), temp.end (), genStrcasecmpAscending () );
		copy ( temp.begin (), temp.end (), back_inserter ( projects ) );
	}
	if ( projects.empty () ) {
		ErrorHandler::genError ()->error ( "No projects have been created yet. Please use Batch-Tag Web.\n" );
	}
	for ( StringVectorSizeType i = 0 ; i < projects.size () ; i++ ) {
		bool comp = projects [i][projects [i].length ()-1] == '$';
		theClass.push_back ( comp ? "compressed_style" : "" );
		if ( comp ) projects [i] = projects [i].substr ( 0, projects [i].length () - 1 );
	}
}
void BTSelProjectForm::createItems ()
{
	formItemMap [USER] = new FormItemText ( "User", "", USER, 20, 60, user );
	formItemMap [FORM] = new FormItemText ( "", "", FORM, 20, 60, "batchtag" );
	formItemMap [SELECT_RESULTS] = new FormItemCheckbox ( "", "", SELECT_RESULTS, true );
	formItemMap [FILTER_PROJECT_LIST] = new FormItemCheckbox ( "Filter Project List", "", FILTER_PROJECT_LIST, true, "1", "showFilterItems( this.form.filter_project_list )" );

	formItemMap [PROJECT_NAME] = new FormItemSelect ( "Project", "", PROJECT_NAME, projects, projects.empty () ? "" : projects [0], 1, "", theClass );
}
void BTSelProjectForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];
		formItemMap [FILTER_PROJECT_LIST]->setValue ( p );
	}
}
void BTSelProjectForm::printHTML ( ostream& os )
{
	init_html ( os, "Select Project" );
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
				printHTMLFORMSubmit ( os, "Select Project" );
			tableHeaderEnd ( os );
		tableRowEnd ( os );

		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				formItemMap [FILTER_PROJECT_LIST]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );

		pdff->printHTML ( os, formItemMap [FILTER_PROJECT_LIST]->getChecked () );
	tableEnd ( os);
	showHiddenItems ( os );
	printHTMLFORMEnd ( os );
	printHTMLFooter ( os, "2007" );
}
void BTSelProjectForm::printJavascriptFunctions ( ostream& os )
{
	basicJavascriptFunctions ( os );
	startJavascript ( os );
	os << "function showFilterItems ( item ) {" << endl;
	os << "\t" << "if ( item.checked ) {" << endl;
	os << "\t\t" << "showdiv ( 'show_projects' );" << endl;
	os << "\t" << "}" << endl;
	os << "\t" << "else {" << endl;
	os << "\t\t" << "hidediv ( 'show_projects' );" << endl;
	os << "\t" << "}" << endl;
	os << "}" << endl;
	pdff->printJavascriptFunctions ( os );
	endJavascript ( os );
}
BTSelResultsForm::BTSelResultsForm ( const ParameterList* paramList )
{
	user = paramList->getStringValue ( USER, user );
	project = paramList->getStringValue ( PROJECT_NAME, project );
	form = paramList->getStringValue ( FORM, form );
	VectorConstParameterListPtr pv;
	pv.push_back ( paramList );
	create ( pv );
}
void BTSelResultsForm::setOptions ()
{
	results.push_back ( "Default...." );
	StringVector temp = MySQLPPSDDBase::instance ().getResultsList ( user, project );
	sort ( temp.begin (), temp.end (), genStrcasecmpAscending () );
	copy ( temp.begin (), temp.end (), back_inserter ( results ) );
}
void BTSelResultsForm::createItems ()
{
	formItemMap [USER] = new FormItemText ( "User", "", USER, 20, 60, user );
	formItemMap [FORM] = new FormItemText ( "", "", FORM, 20, 60, form );
	formItemMap [PROJECT_NAME] = new FormItemText ( "", "", PROJECT_NAME, 20, 58, project );

	formItemMap [RESULTS_FILE] = new FormItemSelect ( "Results File", "", RESULTS_FILE, results, results [0] );
}
void BTSelResultsForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];
		form = p->getStringValue ( FORM, form );
		user = p->getStringValue ( USER, user );
	}
}
void BTSelResultsForm::printHTML ( ostream& os )
{
	init_html ( os, "Select Results" );
	printHTMLFORMStart ( os, "post", "msform" );
	tableStart ( os, true, "center" );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				formItemMap [RESULTS_FILE]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );

		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				printHTMLFORMSubmit ( os, "Select Results" );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
	tableEnd ( os);
	showHiddenItems ( os );
	printHTMLFORMEnd ( os );
	printHTMLFooter ( os, "2007" );
}
#endif
