/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_scsel_form.h                                               *
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

#ifndef __lu_scsel_form_h
#define __lu_scsel_form_h

#include <lgen_define.h>
#include <lu_pros_form.h>

class ProjectDateFilterForm;

class SCSelProjectForm : public ProspectorForm {
	ProjectDateFilterForm* pdff;
protected:
	bool setCookie;
	std::string user;
	StringVector projects;
	StringVector theClass;

	void setOptions ();
	void createItems ();
	void printJavascriptFunctions ( std::ostream& os );
public:
	SCSelProjectForm ( const ParameterList* paramList, bool setCookie );
	virtual void printHTML ( std::ostream& os );
};

class SCSelResultsForm : public ProspectorForm {
protected:
	std::string form;
	bool guest;
	std::string user;
	StringVector projects;
	StringVector results;

	void setOptions ();
	void createItems ();
public:
	SCSelResultsForm ( const ParameterList* paramList );
	virtual void printHTML ( std::ostream& os );
};

#endif /* ! __lu_scsel_form_h */
