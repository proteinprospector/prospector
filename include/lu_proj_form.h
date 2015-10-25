/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_proj_form.h                                                *
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
*  Copyright (2005-2008) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_proj_form_h
#define __lu_proj_form_h

#include <ostream>
#include <lu_form_valid.h>
#include <lu_html_form.h>
#include <lu_pros_form.h>

class MakeProjectForm : public ProspectorForm {
	FormValidatingJavascript fvj;
	StringVector dirList;
	bool finalForm;
	std::string user;
protected:
	void createItems ();
public:
	MakeProjectForm ( const ParameterList* paramList );
	virtual void printHTML ( std::ostream& os );
};

class ImportProjectForm : public ProspectorForm {
	FormValidatingJavascript fvj;
	std::string user;
protected:
	void createItems ();
public:
	ImportProjectForm ( const ParameterList* paramList );
	virtual void printHTML ( std::ostream& os );
};

#endif /* ! __lu_proj_form_h */
