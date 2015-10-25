/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_add_user.cpp                                               *
*                                                                             *
*  Created    : June 28th 2007                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2007-2011) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef MYSQL_DATABASE
#include <lu_html.h>
#include <lu_add_user.h>
#include <lu_table.h>
#include <lu_html_form.h>

using std::string;
using std::vector;
using std::ostream;
using std::endl;

static string USER	= "user";
static string EMAIL	= "email";
static string PASSWORD = "password";
static string RETYPE_PASSWORD = "retype_password";
static string FIRST_NAME = "first_name";
static string LAST_NAME = "last_name";

static string FORM = "form";
static string SUBMIT = "submit";

AddUserForm::AddUserForm ()
{
	create ( VectorConstParameterListPtr () );
}
void AddUserForm::createItems ()
{
	formItemMap [USER] = new FormItemText ( "User Name", "", USER, 10, 10, "" );
	formItemMap [EMAIL] = new FormItemText ( "Email", "", EMAIL, 50, 60, "" );
	formItemMap [PASSWORD] = new FormItemPassword ( "Password", "", PASSWORD, 16, 60, "" );
	formItemMap [RETYPE_PASSWORD] = new FormItemPassword ( "Retype Password", "", RETYPE_PASSWORD, 16, 60, "" );
	formItemMap [FIRST_NAME] = new FormItemText ( "First Name (optional)", "", FIRST_NAME, 24, 60, "" );
	formItemMap [LAST_NAME] = new FormItemText ( "Surname (optional)", "", LAST_NAME, 24, 60, "" );

	formItemMap [FORM] = new FormItemText ( "", "", FORM, 20, 60, "add_user" );
	formItemMap [SUBMIT] = new FormItemCheckbox ( "", "", SUBMIT, true );
}
void AddUserForm::printHTML ( ostream& os )
{
	init_html ( os, "Add User" );
	printHTMLFORMStart ( os, "post", "msform" );
	tableStart ( os, true, "center" );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true );
				formItemMap [USER]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [EMAIL]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [PASSWORD]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [RETYPE_PASSWORD]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [FIRST_NAME]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [LAST_NAME]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );

		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				printHTMLFORMSubmit ( os, "Add User" );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
	tableEnd ( os);
	showHiddenItems ( os );
	printHTMLFORMEnd ( os );
	printHTMLFooter ( os, "2007" );
}
#endif
