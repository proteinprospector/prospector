/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lg_email.h                                                    *
*                                                                             *
*  Created    : June 29th 2007                                                *
*                                                                             *
*  Purpose    : Checks email format.                   .                      *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2007-2009) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lg_email_h
#define __lg_email_h

#include <string>

bool isValidEmailAddress ( const char* address );
bool isValidName ( const std::string& name );
bool isValidUser ( const std::string& user, std::string::size_type min, std::string::size_type max );

#endif /* ! __lg_email_h */
