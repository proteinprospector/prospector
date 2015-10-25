/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_login_form.h                                               *
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
*  Copyright (2007-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_login_form_h
#define __lu_login_form_h

#include <lu_pros_form.h>

class LoginForm : public ProspectorForm {
protected:
	std::string formValue;

	void createItems ();
public:
	LoginForm ( const ParameterList* paramList );
	virtual void printHTML ( std::ostream& os );
};

#endif /* ! __lu_login_form_h */
