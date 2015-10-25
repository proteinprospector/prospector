/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_proj_form.cpp                                              *
*                                                                             *
*  Created    : February 11th 2005                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2005-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef MYSQL_DATABASE
#include <algorithm>
#include <lgen_file.h>
#include <lu_html.h>
#include <lu_table.h>
#include <lu_file_type.h>
#include <lu_proj_form.h>
#include <lu_getfil.h>
#include <lu_param_list.h>
using std::sort;
using std::string;
using std::vector;
using std::ostream;
using std::endl;
using std::copy;
using std::back_inserter;
using namespace FileTypes;

static string DATA_FILE_LIST = "data_file_list";
static string FORM = "form";
static string WRITE_FILE = "write_file";
static string PROJECT_NAME = "project_name";
static string USER = "user";

static void directoryStyle ( ostream& os )
{
	static bool printed = false;
	if ( !printed ) {
		os << "<style>" << endl;
		os << ".directory_style" << endl;
		os << "{" << endl;
		os << "background: #fff url(" << HTMLDir::instance ().getVirtualHTMLImagesDir () << "folder.gif) no-repeat;" << endl;
		os << "padding-left: 18px;" << endl;
		os << "}" << endl;
		os << "</style>" << endl;
		printed = true;
	}
}
MakeProjectForm::MakeProjectForm ( const ParameterList* paramList ) :
	finalForm ( false )
{
	user = paramList->getStringValue ( USER, user );
	dirList = paramList->getStringVectorValue ( "data_file_list" );
	if ( dirList.empty () ) dirList.push_back ( "" );
	VectorConstParameterListPtr pv;
	pv.push_back ( paramList );
	create ( pv );
}
void MakeProjectForm::createItems ()
{
	if ( user == "root" ) {
		ErrorHandler::genError ()->error ( "root user can't create projects.\n" );
	}
	string centroidDir = adjustPPOutputPath ( InfoParams::instance ().getCentroidDir () );
	StringVector fullFileList;
	StringVector selected;
	string dirClass = "directory_style";
	bool dirListEmpty = true;
	for ( StringVectorSizeType i = 0 ; i < dirList.size () ; i++ ) {
		if ( genIsDirectory ( centroidDir + dirList [i] ) ) {			// If the entry in the directory list is a directory
			StringVector fList = FileList ( centroidDir + dirList [i] ).getNameList ();		// Get all the directories in the directory
			for ( StringVectorSizeType j = 0 ; j < fList.size () ; j++ ) {
				string s = dirList [i];
				if ( s.length () ) s += "/";
				s += fList [j];
				fList [j] = s;				// Append the directory to the path
			}
			if ( !fList.empty () ) dirListEmpty = false;	// There are more sub directories
			copy ( fList.begin (), fList.end (), back_inserter ( fullFileList ) );
		}
	}
	sort ( fullFileList.begin (), fullFileList.end (), genStrcasecmpAscending () );		// Sort the directories into order
	StringVector theClass ( fullFileList.size (), dirClass );
	for ( StringVectorSizeType i2 = 0 ; i2 < dirList.size () ; i2++ ) {		// Select the files that have already been selected
		if ( !genIsDirectory ( centroidDir + dirList [i2] ) ) {
			fullFileList.push_back ( dirList [i2] );
			selected.push_back ( dirList [i2] );
		}
	}
	StringVector fileOnlyList;
	for ( StringVectorSizeType ii = 0 ; ii < dirList.size () ; ii++ ) {
		if ( genIsDirectory ( centroidDir + dirList [ii] ) ) {			// If the entry in the directory list is a directory
			string d = centroidDir + dirList [ii];
			FileList f1 ( d, "", "", false );
			StringVector fList;
			for ( int x = 0 ; x < f1.size () ; x++ ) {
				string f = f1 [x];
				string fp = d + SLASH + f;
				if ( isFileType ( f, MGF ) )		fList.push_back ( f );
				else if ( isFileType ( f, MZXML ) )	fList.push_back ( f );
				else if ( isFileType ( f, TXT ) )	fList.push_back ( f );
				else if ( isFileType ( f, MZDATA ) )fList.push_back ( f );
				else if ( isFileType ( f, XML ) )	fList.push_back ( f );
			}
			for ( StringVectorSizeType j = 0 ; j < fList.size () ; j++ ) {
				string s = dirList [ii];
				if ( s.length () ) s += "/";
				s += fList [j];
				fList [j] = s;
			}
			copy ( fList.begin (), fList.end (), back_inserter ( fileOnlyList ) );
		}
	}
	sort ( fileOnlyList.begin (), fileOnlyList.end (), genStrcasecmpAscending () );
	copy ( fileOnlyList.begin (), fileOnlyList.end (), back_inserter ( fullFileList ) );
	theClass.resize ( fullFileList.size () );
	if ( dirListEmpty ) {		// No more subdirectories
		formItemMap [WRITE_FILE] = new FormItemCheckbox ( "", "", WRITE_FILE, true );
		formItemMap [PROJECT_NAME] = new FormItemText ( "Project Name", "", PROJECT_NAME, 40, 58, "", fvj.addFilenameValidator ( PROJECT_NAME, "Project Name" ) );
		finalForm = true;
	}
	int size = genMin ( (int)fullFileList.size (), 20 );
	if ( fullFileList.size () == 1 && selected.empty () ) {
		selected.push_back ( fullFileList [0] ); // If there is only one item select it.
	}
	formItemMap [DATA_FILE_LIST] = new FormItemSelectMultiple ( "", "", DATA_FILE_LIST, fullFileList, selected, size, theClass );
	formItemMap [FORM] = new FormItemText ( "FORM", "", FORM, 12, 20, "makeproject" );
	formItemMap [USER] = new FormItemText ( "User", "", USER, 20, 60, user );
}
void MakeProjectForm::printHTML ( ostream& os )
{
	init_html ( os, "Make Project" );
	directoryStyle ( os );
	if ( finalForm ) {
		fvj.print ( os );
		printHTMLFORMStart ( os, "post", "msform", true, true );
	}
	else printHTMLFORMStart ( os, "post", "msform" );
	tableStart ( os, true, "center" );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				if ( finalForm )
					printHTMLFORMSubmit ( os, "Make Project" );
				else
					printHTMLFORMSubmit ( os, "Select" );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				formItemMap [DATA_FILE_LIST]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		if ( finalForm ) {
			tableRowStart ( os );
				tableHeaderStart ( os, "", "left", true );
					formItemMap [PROJECT_NAME]->printHTML ( os );
				tableHeaderEnd ( os );
			tableRowEnd ( os );
		}
	tableEnd ( os );
	showHiddenItems ( os );
	printHTMLFORMEnd ( os );
	printHTMLFooter ( os, "2005" );
}

class FormItemUploadProject : public FormItemFile {
public:
	FormItemUploadProject () :
		FormItemFile ( "", "", getName (), 100, 256, "" ) {}
	static string getName () { return "upload_project"; }
};

ImportProjectForm::ImportProjectForm ( const ParameterList* paramList )
{
	user = paramList->getStringValue ( USER, user );
	VectorConstParameterListPtr pv;
	pv.push_back ( paramList );
	create ( pv );
}
void ImportProjectForm::createItems ()
{
	if ( user == "root" ) {
		ErrorHandler::genError ()->error ( "root user can't create projects.\n" );
	}
	formItemMap [FormItemUploadProject::getName ()]	= new FormItemUploadProject ();
	formItemMap [PROJECT_NAME] = new FormItemText ( "Project Name", "", PROJECT_NAME, 40, 58, "", fvj.addFilenameAllowBlankValidator ( PROJECT_NAME, "Project Name" ) );
	formItemMap [FORM] = new FormItemText ( "FORM", "", FORM, 12, 20, "importproject" );
	formItemMap [USER] = new FormItemText ( "User", "", USER, 20, 60, user );
}
void ImportProjectForm::printHTML ( ostream& os )
{
	init_html ( os, "Import Project" );
	fvj.print ( os );
	printHTMLFORMStart ( os, "post", "msform", true, true );
	tableStart ( os, true, "center" );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				printHTMLFORMSubmit ( os, "Import Project" );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				os << "<div class=\"form_large_label\">Upload Project From File</div>" << endl;
				formItemMap [FormItemUploadProject::getName ()]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true );
				formItemMap [PROJECT_NAME]->printHTML ( os );
				os << "(Leave Blank to Keep the Same Name)" << endl;
			tableHeaderEnd ( os );
		tableRowEnd ( os );
	tableEnd ( os );
	showHiddenItems ( os );
	printHTMLFORMEnd ( os );
	printHTMLFooter ( os, "2008" );
}
#endif
