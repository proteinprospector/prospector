/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_login_form.cpp                                             *
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
*  Copyright (2007-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef MYSQL_DATABASE
#include <lu_html.h>
#include <lu_login_form.h>
#include <lu_table.h>
#include <lu_param_list.h>
#include <lu_html_form.h>
using std::string;
using std::vector;
using std::ostream;
using std::endl;

namespace {
	string FORM	= "form";
	string USER	= "user";
	string PASSWORD	= "password";
	string NEW_PASSWORD	= "new_password";
	string RETYPE_NEW_PASSWORD	= "retype_new_password";
	string SELECT_PROJECT = "select_project";
}

LoginForm::LoginForm ( const ParameterList* paramList )
{
	formValue = paramList->getStringValue ( FORM, formValue );
	VectorConstParameterListPtr pv;
	pv.push_back ( paramList );
	create ( pv );
}
void LoginForm::createItems ()
{
	formItemMap [USER] = new FormItemText ( "User", "", USER, 20, 60, ParameterList::getCookieValue ( USER ) );
	formItemMap [PASSWORD] = new FormItemPassword ( "Password", "", PASSWORD, 16, 60, "" );
	formItemMap [NEW_PASSWORD] = new FormItemPassword ( "New Password", "", NEW_PASSWORD, 16, 60, "" );
	formItemMap [RETYPE_NEW_PASSWORD] = new FormItemPassword ( "Retype New Password", "", RETYPE_NEW_PASSWORD, 16, 60, "" );
	formItemMap [FORM] = new FormItemText ( "", "", FORM, 20, 60, formValue );
	if ( formValue != "search_table" && formValue != "batchtagweb" ) {
		formItemMap [SELECT_PROJECT] = new FormItemCheckbox ( "", "", SELECT_PROJECT, true );
	}
}
void LoginForm::printHTML ( ostream& os )
{
	string title;
	init_html ( os, "Login" );
	printHTMLFORMStart ( os, "post", "msform" );
	tableStart ( os, true, "center" );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true );
				formItemMap [USER]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );

		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true );
				formItemMap [PASSWORD]->printHTML ( os );
				os << "<br />" << endl;
				ExpandableJavascriptBlock ejb ( "Change Password" );
				ejb.printHeader ( os );
				formItemMap [NEW_PASSWORD]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [RETYPE_NEW_PASSWORD]->printHTML ( os );
				ejb.printFooter ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );

		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				ParameterList pList ( "msform", false, false, false, false, false );
				pList.addName ( "form", "add_user" );
				os << "<a href=\"" << pList.getURL () << "\">Add User</a>" << endl;
			tableHeaderEnd ( os );
		tableRowEnd ( os );

		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				printHTMLFORMSubmit ( os, "Login" );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
	tableEnd ( os);
	showHiddenItems ( os );
	printHTMLFORMEnd ( os );
	printHTMLFooter ( os, "2007" );
}
#endif
