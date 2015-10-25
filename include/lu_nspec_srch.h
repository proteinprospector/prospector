/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_nspec_srch.h                                               *
*                                                                             *
*  Created    : February 12th 2002                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2002-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_nspec_srch_h
#define __lu_nspec_srch_h

#include <lu_spep_srch.h>
#include <lu_dig_srch.h>

class MSNonSpecificParameters;

class MSNonSpecificSearch : public MSSingleSearch {
	MSNonSpecificParameters& nonSpecificParams;
	PeakContainer parentPeaks;
	std::vector <NonSpecificSearch*> nonSpecificSearch;
public:
	MSNonSpecificSearch ( MSNonSpecificParameters& nonSpecificParams );
	~MSNonSpecificSearch ();
	void printBodyHTML ( std::ostream& os );
	void printBodyXML ( std::ostream& os );
	void printParamsBodyHTML ( std::ostream& os ) const;
	void printProteinHTML ( std::ostream& os, int searchIndex, int proteinIndex ) const;
	void printResultsHTML ( std::ostream& os, int i );
	void printResultsXML ( std::ostream& os, int i );
};

#endif /* ! __lu_nspec_srch_h */
