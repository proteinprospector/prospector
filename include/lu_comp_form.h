/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_comp_form.h                                                *
*                                                                             *
*  Created    : January 11th 2005                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2005-2009) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_comp_form_h
#define __lu_comp_form_h

#include <ostream>
#include <lu_srch_form.h>
#include <lu_form_valid.h>

class MSCompForm : public ProspectorForm {
	FormValidatingJavascript fvj;
	UserAAForm userAAForm;
	ImmoniumIonForm immoniumIonForm;
	MSCompExtraIonForm extraIonForm;
	MSCompAbsentIonForm absentIonForm;
	MSCompModifiedIonForm modifiedIonForm;
	SaveResultsForm saveResultsForm;
protected:
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	MSCompForm ( const VectorConstParameterListPtr& params );
	virtual void printHTML ( std::ostream& os );
};

#endif /* ! __lu_comp_form_h */
