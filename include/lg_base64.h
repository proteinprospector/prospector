/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lg_base64.h                                                   *
*                                                                             *
*  Created    : April 14th 2005                                               *
*                                                                             *
*  Purpose    : Base 64 encoding and decoding          .                      *
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

#ifndef __lg_base64_h
#define __lg_base64_h

#include <string>

std::string base64_encode ( unsigned char const* bytes_to_encode, unsigned int in_len );
void base64Decode ( const std::string& inStr, std::string& outStr );

#endif /* ! __lg_base64_h */
