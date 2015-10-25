/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_patt_form.h                                                *
*                                                                             *
*  Created    : January 4th 2005                                              *
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

#ifndef __lu_patt_form_h
#define __lu_patt_form_h

#include <ostream>
#include <lu_form_valid.h>
#include <lu_srch_form.h>

class MSPatternForm : public ProspectorForm {
protected:
	FormValidatingJavascript fvj;
	HeaderForm headerForm;
	SearchForm searchForm;

	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	MSPatternForm ( const VectorConstParameterListPtr& params );
	virtual void printHTML ( std::ostream& os );
};

#endif /* ! __lu_patt_form_h */
