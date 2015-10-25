/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_faind_form.h                                               *
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

#ifndef __lu_faind_form_h
#define __lu_faind_form_h

#include <ostream>
#include <lu_form_valid.h>
#include <lu_srch_form.h>

class FAIndexForm : public ProspectorForm {
	FormValidatingJavascript fvj;
	ProteinMWForm proteinMWForm;
	ProteinPIForm proteinPIForm;
	SearchResultsForm searchResultsForm;
protected:
	void setValues ( const VectorConstParameterListPtr& params ) {}
	void createItems ();
public:
	FAIndexForm ( const VectorConstParameterListPtr& params );
	virtual void printHTML ( std::ostream& os );
};

#endif /* ! __lu_faind_form_h */
