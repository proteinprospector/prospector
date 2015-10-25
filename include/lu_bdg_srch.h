/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_bdg_srch.h                                                 *
*                                                                             *
*  Created    : September 27th 2001                                           *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2008) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_bdg_srch_h
#define __lu_bdg_srch_h

#include <lu_links.h>
#include <lu_get_link.h>
#include <lu_dig_srch.h>

class MSBridgeParameters;

class MSBridgeSearch : public MSSingleSearch {
	const MSBridgeParameters& bridgeParams;
	PeakContainer parentPeaks;
	LinkInfo* linkInfo;
	std::vector <LinksSearch*> linksSearch;

	PeakContainer getParentPeaks ( MSBridgeParameters& bridgeParams );
public:
	MSBridgeSearch ( MSBridgeParameters& bridgeParams );
	void printParamsBodyHTML ( std::ostream& os ) const;
	void printProteinHTML ( std::ostream& os, int searchIndex, int proteinIndex ) const;
	void printResultsHTML ( std::ostream& os, int i );
	void printResultsXML ( std::ostream& os, int i );
};

MSSingleSearch* getMSBridgeSearch ( MSBridgeParameters& bridgeParams );

#endif /* ! __lu_bdg_srch_h */
