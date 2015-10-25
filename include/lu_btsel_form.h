/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_btsel_form.h                                               *
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

#ifndef __lu_btsel_form_h
#define __lu_btsel_form_h

#include <lu_pros_form.h>

class ProjectDateFilterForm;

class BTSelProjectForm : public ProspectorForm {
	ProjectDateFilterForm* pdff;
protected:
	bool setCookie;
	std::string user;
	StringVector projects;
	StringVector theClass;

	void setOptions ();
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
	void printJavascriptFunctions ( std::ostream& os );
public:
	BTSelProjectForm ( const ParameterList* paramList, bool setCookie );
	virtual void printHTML ( std::ostream& os );
};

class BTSelResultsForm : public ProspectorForm {
protected:
	std::string form;
	std::string user;
	std::string project;
	StringVector results;

	void setOptions ();
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	BTSelResultsForm ( const ParameterList* paramList );
	virtual void printHTML ( std::ostream& os );
};

#endif /* ! __lu_btsel_form_h */
