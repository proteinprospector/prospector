/******************************************************************************
*                                                                             *
*  Library    : libanalyst                                                    *
*                                                                             *
*  Filename   : la_wiff_xic.h                                                 *
*                                                                             *
*  Created    : March 31st 2009                                               *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2009-2009) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __la_wiff_xic_h
#define __la_wiff_xic_h

#include <string>
#include <lgen_define.h>

class WiffExploreInstance2;

class WiffFile2 {
protected:
	WiffExploreInstance2* inst;
	DoubleVector times;
	DoubleVector intensities;
	void initCache (  const std::string& file, const std::string& file2 );
public:
	WiffFile2 ( const std::string& file, const std::string& file2 );
	~WiffFile2 ();
};

#endif /* ! __la_wiff_xic_h */
