/******************************************************************************
*                                                                             *
*  Library    : libanalyst                                                    *
*                                                                             *
*  Filename   : la_wiff_cent.h                                                *
*                                                                             *
*  Created    : June 2nd 2005                                                 *
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

#ifndef __la_wiff_cent_h
#define __la_wiff_cent_h

#include <string>

class WiffCentroidData {
public:
	WiffCentroidData ( const std::string& rawPath, const std::string& centroidPath );
};

#endif /* ! __la_wiff_cent_h */
