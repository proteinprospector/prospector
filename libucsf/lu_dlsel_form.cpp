/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_dlsel_form.cpp                                             *
*                                                                             *
*  Created    : September 27th 2007                                           *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2007-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef MYSQL_DATABASE
#include <algorithm>
#include <lu_cookie.h>
#include <lu_html.h>
#include <lu_dlsel_form.h>
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
using std::back_inserter;
using std::endl;

static string FORM = "form";
static string USER = "user";
static string PROJECT_NAME = "project_name";
static string RESULTS_FILE = "results_file";
static string UPLOAD_TEMP = "upload_temp";
static string UPLOADS_OPTIONAL = "uploads_optional";

class FormItemAction : public FormItemSelect {
	static const char* options [];
public:
	FormItemAction ();
	static std::string getName () { return "action"; }
};
const char* FormItemAction::options [] = { "Delete", "Import", "Export", "Compress", "Uncompress", "Check", 0 };

FormItemAction::FormItemAction () :
	FormItemSelect ( "", "", getName (), options, "Delete", 1, "showResultsManagementItems( this.form.action )" )
{
}

class FormItemDataExport : public FormItemSelect {
	static const char* options [];
public:
	FormItemDataExport ();
	static std::string getName () { return "data_export"; }
};
const char* FormItemDataExport::options [] = { "All Data", "Peaklists Only", "No Data", 0 };

FormItemDataExport::FormItemDataExport () :
	FormItemSelect ( "Data To Export", "", getName (), options, options [0] )
{
}

class FormItemNewProjectName : public FormItemText {
public:
	FormItemNewProjectName ( FormValidatingJavascript& fvj );
	static std::string getName () { return "new_project_name"; }
};

FormItemNewProjectName::FormItemNewProjectName ( FormValidatingJavascript& fvj ) :
	FormItemText ( "New Project Name", "", getName (), 40, 58, "", fvj.addFilenameAllowBlankValidator ( getName (), "New Project Name" ) )
{
}

DeleteSelForm::DeleteSelForm ( const ParameterList* paramList, bool setCookie ) :
	pdff ( new ProjectDateFilterForm ( paramList ) ),
	fvj (),
	setCookie ( setCookie )
{
	user = paramList->getStringValue ( USER, user );
	VectorConstParameterListPtr pv;
	pv.push_back ( paramList );
	create ( pv );
}
void DeleteSelForm::setOptions ()
{
	projects = MySQLPPSDDBase::instance ().getProjectList ( user, true, pdff->getProjectDateFilter (), true );
	sort ( projects.begin (), projects.end (), genStrcasecmpAscending () );

	for ( StringVectorSizeType i = 0 ; i < projects.size () ; i++ ) {
		bool comp = projects [i][projects [i].length ()-1] == '$';
		theClass.push_back ( comp ? "compressed_style" : "" );
		if ( comp ) projects [i] = projects [i].substr ( 0, projects [i].length () - 1 );
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
void DeleteSelForm::createItems ()
{
	formItemMap [USER] = new FormItemText ( "User", "", USER, 20, 60, user );
	formItemMap [FORM] = new FormItemText ( "", "", FORM, 20, 60, "results_management" );

	int projectsSize = genMin ( 10, (int)projects.size () );
	formItemMap [PROJECT_NAME] = new FormItemSelectMultiple ( "Project", "", PROJECT_NAME, projects, StringVector (), projectsSize, theClass );
	
	int resultsSize = genMin ( 20, (int)results.size () );
	formItemMap [RESULTS_FILE] = new FormItemSelectMultiple ( "Results", "", RESULTS_FILE, results, StringVector (), resultsSize );

	formItemMap [FormItemAction::getName ()] = new FormItemAction ();

	formItemMap [UPLOAD_TEMP] = new FormItemFile ( "File", "", UPLOAD_TEMP, 100, 256, "" );
	formItemMap [UPLOADS_OPTIONAL] = new FormItemCheckbox ( "", "", UPLOADS_OPTIONAL, true );

	formItemMap [FormItemNewProjectName::getName ()] = new FormItemNewProjectName ( fvj );

	formItemMap [FormItemDataExport::getName ()] = new FormItemDataExport ();
}
void DeleteSelForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];

		formItemMap [FormItemAction::getName ()]->setValue ( p );
		formItemMap [FormItemNewProjectName::getName ()]->setValue ( p );
		formItemMap [FormItemDataExport::getName ()]->setValue ( p );
	}
}
void DeleteSelForm::printHTML ( ostream& os )
{
	init_html ( os, "Results Management" );
	compressedHTMLStyle ( os );
	if ( setCookie ) setUserCookie ( os, user );

	printHTMLFORMStart ( os, "post", "msform", true );
	printJavascriptFunctions ( os );
	tableStart ( os, true, "center" );
		bool projDiv, resDiv, showProjectsDiv, deDiv, upDiv, up2Div;
		setActionVisualizationFlags ( formItemMap [FormItemAction::getName ()]->getValue ( 0 ), projDiv, resDiv, showProjectsDiv, deDiv, upDiv, up2Div );
		tableRowStartDiv ( os, "proj", projDiv );
			tableHeaderStart ( os, "", "center", true );
				if ( projects.empty () ) {
					printHTMLFORMLabel ( os, "No projects present", "" );
				}
				else {
					formItemMap [PROJECT_NAME]->printHTML ( os );
				}
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStartDiv ( os, "res", resDiv );
			tableHeaderStart ( os, "", "center", true );
				if ( results.empty () ) {
					printHTMLFORMLabel ( os, "No results present", "" );
				}
				else {
					formItemMap [RESULTS_FILE]->printHTML ( os );
				}
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				tableStart ( os, false, "center" );
					tableRowStart ( os );
						tableHeaderStart ( os, "", "center", true );
							printHTMLFORMSubmit ( os, "Submit" );
						tableHeaderEnd ( os );
						tableHeaderStart ( os, "", "center", true );
							formItemMap [FormItemAction::getName ()]->printHTML ( os );
						tableHeaderEnd ( os );
					tableRowEnd ( os );
					tableRowStartDiv ( os, "de", deDiv );
						tableHeaderStart ( os, "", "center", true );
							formItemMap [FormItemDataExport::getName ()]->printHTML ( os );
						tableHeaderEnd ( os );
					tableRowEnd ( os );
				tableEnd ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStartDiv ( os, "up", upDiv );
			tableHeaderStart ( os, "", "center", true );
				formItemMap [UPLOAD_TEMP]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStartDiv ( os, "up2", up2Div );
			tableHeaderStart ( os, "", "center", true );
				formItemMap [FormItemNewProjectName::getName ()]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		pdff->printHTML ( os, showProjectsDiv );
	tableEnd ( os );

	formItemMap [USER]->printHTMLHidden ( os );
	formItemMap [FORM]->printHTMLHidden ( os );
	formItemMap [UPLOADS_OPTIONAL]->printHTMLHidden ( os );
	printHTMLFORMEnd ( os );

	printHTMLFooter ( os, "2007" );
}
void DeleteSelForm::setActionVisualizationFlags ( const string& val, bool& projDiv, bool& resDiv, bool& showProjectsDiv, bool& deDiv, bool& upDiv, bool& up2Div ) const
{
	if ( val == "Delete" ) {
		projDiv = true;
		resDiv = true;
		showProjectsDiv = true;
		deDiv = false;
		upDiv = false;
		up2Div = false;
	}
	else if ( val == "Import" ) {
		projDiv = false;
		resDiv = false;
		showProjectsDiv = false;
		deDiv = false;
		upDiv = true;
		up2Div = true;
	}
	else if ( val == "Compress" || val == "Uncompress" || val == "Check" ) {
		projDiv = true;
		resDiv = false;
		showProjectsDiv = true;
		deDiv = false;
		upDiv = false;
		up2Div = false;
	}
	else {								// Export
		projDiv = true;
		resDiv = false;
		showProjectsDiv = true;
		deDiv = true;
		upDiv = false;
		up2Div = false;
	}
}
void DeleteSelForm::printJavascriptFunctions ( ostream& os )
{
	basicJavascriptFunctions ( os );
	startJavascript ( os );
	os << "function showResultsManagementItems ( item ) {" << endl;
	os << "\t" << "var val = getSelectValue(item);" << endl;
	os << "\t" << "if ( val == 'Delete' ) {" << endl;
	os << "\t\t" << "showdiv ( 'proj' );" << endl;
	os << "\t\t" << "showdiv ( 'res' );" << endl;
	os << "\t\t" << "showdiv ( 'show_projects' );" << endl;
	os << "\t\t" << "hidediv ( 'up' );" << endl;
	os << "\t\t" << "hidediv ( 'up2' );" << endl;
	os << "\t\t" << "hidediv ( 'de' );" << endl;
	os << "\t" << "}" << endl;
	os << "\t" << "else if ( val == 'Import' ) {" << endl;
	os << "\t\t" << "hidediv ( 'proj' );" << endl;
	os << "\t\t" << "hidediv ( 'res' );" << endl;
	os << "\t\t" << "hidediv ( 'show_projects' );" << endl;
	os << "\t\t" << "hidediv ( 'de' );" << endl;
	os << "\t\t" << "showdiv ( 'up' );" << endl;
	os << "\t\t" << "showdiv ( 'up2' );" << endl;
	os << "\t" << "}" << endl;
	os << "\t" << "else if ( val == 'Compress' || val == 'Uncompress' || val == 'Check' ) {" << endl;
	os << "\t\t" << "showdiv ( 'proj' );" << endl;
	os << "\t\t" << "showdiv ( 'show_projects' );" << endl;
	os << "\t\t" << "hidediv ( 'de' );" << endl;
	os << "\t\t" << "hidediv ( 'res' );" << endl;
	os << "\t\t" << "hidediv ( 'up' );" << endl;
	os << "\t\t" << "hidediv ( 'up2' );" << endl;
	os << "\t" << "}" << endl;
	os << "\t" << "else {" << endl;							// Export
	os << "\t\t" << "showdiv ( 'proj' );" << endl;
	os << "\t\t" << "showdiv ( 'de' );" << endl;
	os << "\t\t" << "showdiv ( 'show_projects' );" << endl;
	os << "\t\t" << "hidediv ( 'res' );" << endl;
	os << "\t\t" << "hidediv ( 'up' );" << endl;
	os << "\t\t" << "hidediv ( 'up2' );" << endl;
	os << "\t" << "}" << endl;
	os << "}" << endl;
	pdff->printJavascriptFunctions ( os );
	endJavascript ( os );
}
#endif
