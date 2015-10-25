/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lgen_net.cpp                                                  *
*                                                                             *
*  Created    : February 26th 2008                                            *
*                                                                             *
*  Purpose    :                                        .                      *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2008-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifndef VIS_C
#include <unistd.h>
#endif
#ifdef VIS_C
#include <windows.h>
#endif
#include <lgen_net.h>

using std::string;

Hostname::Hostname ()
{
#ifdef VIS_C
	WSADATA WSAData;
	WSAStartup ( MAKEWORD ( 1, 0 ), &WSAData );	// Initialize winsock dll
#endif
	char h [128] = "";
	gethostname ( h, sizeof(h) - 1 );
#ifdef VIS_C
	WSACleanup ();
#endif
	hostname = string ( h );
}
Hostname& Hostname::instance ()
{
	static Hostname h;
	return h;
}
