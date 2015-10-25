/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_dlsel_form.h                                               *
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
*  Copyright (2007-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_dlsel_form_h
#define __lu_dlsel_form_h

#include <lu_form_valid.h>
#include <lu_pros_form.h>

class ProjectDateFilterForm;

class DeleteSelForm : public ProspectorForm {
	ProjectDateFilterForm* pdff;
	FormValidatingJavascript fvj;
protected:
	bool setCookie;
	std::string user;
	StringVector projects;
	StringVector results;
	StringVector theClass;

	void setOptions ();
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
	void printJavascriptFunctions ( std::ostream& os );
	void setActionVisualizationFlags ( const std::string& val, bool& projDiv, bool& resDiv, bool& showProjectsDiv, bool& deDiv, bool& upDiv, bool& up2Div ) const;
public:
	DeleteSelForm ( const ParameterList* paramList, bool setCookie );
	virtual void printHTML ( std::ostream& os );
};

#endif /* ! __lu_dlsel_form_h */
