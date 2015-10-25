/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lgen_uncompress.h                                             *
*                                                                             *
*  Created    : March 2nd 2007                                                *
*                                                                             *
*  Purpose    : Uncompression functions.                                      *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2007-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lgen_uncompress_h
#define __lgen_uncompress_h

#include <string>

std::string genUnzip ( const std::string& filename, bool deleteFile, bool createSubdirectory );
void gen7zaUncompress ( const std::string& filename );
bool gen7zaCreate ( const std::string& filename );
bool gen7zaCreate ( const std::string& archiveName, const std::string& filename, const std::string& type );
std::string genPreprocessFile ( const std::string& filename );

#endif /* ! __lgen_uncompress_h */
