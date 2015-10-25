/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_add_user.h                                                 *
*                                                                             *
*  Created    : June 28th 2007                                                *
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

#ifndef __lu_add_user_h
#define __lu_add_user_h

#include <lu_pros_form.h>

class AddUserForm : public ProspectorForm {
protected:
	void createItems ();
public:
	AddUserForm ();
	virtual void printHTML ( std::ostream& os );
};

#endif /* ! __lu_add_user_h */
