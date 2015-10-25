/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_filter_form.h                                              *
*                                                                             *
*  Created    : September 3rd 2012                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2012-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_filter_form_h
#define __lu_filter_form_h

#include <lu_form_valid.h>
#include <lu_srch_form.h>

class MSFilterForm : public ProspectorForm {
	FormValidatingJavascript fvj;
	MPlusHForm mPlusHForm;
protected:
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	MSFilterForm ( const VectorConstParameterListPtr& params );
	virtual void printHTML ( std::ostream& os );
	void printHTMLForm ( std::ostream& os );
};

#endif /* ! __lu_filter_form_h */
