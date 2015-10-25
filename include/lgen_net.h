/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lgen_net.h                                                    *
*                                                                             *
*  Created    : February 26th 2008                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2008-2008) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lgen_net_h
#define __lgen_net_h

#include <string>

class Hostname {
	std::string hostname;
	Hostname ();
public:
	static Hostname& instance ();
	std::string getHostname () const { return hostname; }
};

#endif /* ! __lgen_net_h */
