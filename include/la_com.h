/******************************************************************************
*                                                                             *
*  Library    : libanalyst                                                    *
*                                                                             *
*  Filename   : la_com.h                                                      *
*                                                                             *
*  Created    : June 3rd 2005                                                 *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2005-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __la_com_h
#define __la_com_h

#include <string>
#include <comdef.h>

class ComSingleton {
	ComSingleton ();
public:
	~ComSingleton ();
	static void instance ();
};

class _UBSTR {
	BSTR m_bstr;
public:
	_UBSTR ( const char* psz );
	_UBSTR ( const std::string& psz );
	_UBSTR ( const wchar_t* pwsz );
	~_UBSTR ();
	BSTR getBSTR () const;
};
std::string toString ( BSTR& bstr );

#endif /* ! __la_com_h */
