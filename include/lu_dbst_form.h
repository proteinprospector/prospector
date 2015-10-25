/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_dbst_form.h                                                *
*                                                                             *
*  Created    : January 3rd 2005                                              *
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

#ifndef __lu_dbst_form_h
#define __lu_dbst_form_h

#include <ostream>
#include <lu_srch_form.h>
#include <lu_data_form.h>
#include <lu_form_valid.h>

class DBStatForm : public ProspectorForm {
protected:
	FormValidatingJavascript fvj;
	HeaderForm headerForm;
	SearchForm searchForm;

	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	DBStatForm ( const VectorConstParameterListPtr& params );
	virtual void printHTML ( std::ostream& os );
};

#endif /* ! __lu_dbst_form_h */
