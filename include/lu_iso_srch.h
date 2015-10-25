/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_iso_srch.h                                                 *
*                                                                             *
*  Created    : October 22nd 2001                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2010) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_iso_srch_h
#define __lu_iso_srch_h

#include <lu_program.h>

class MSIsotopeParameters;
class IsotopePeakStats;
class GraphData;

class MSIsotopeSearch : public MSProgram {
	const MSIsotopeParameters& isoParams;
	GraphData* graphData;
	void printSummaryReportHTML ( std::ostream& os, const IsotopePeakStats* ips );
	void printSummaryReportXML ( std::ostream& os, const IsotopePeakStats* ips );
	void printSummaryReportTabDelimitedText ( std::ostream& os, const IsotopePeakStats* ips );
	void printDetailedReportHTML ( std::ostream& os, const IsotopePeakStats* ips );
	void printDetailedReportXML ( std::ostream& os, const IsotopePeakStats* ips );
	void printDetailedReportTabDelimitedText ( std::ostream& os, const IsotopePeakStats* ips );
	void setGraphData ();
public:
	MSIsotopeSearch ( const MSIsotopeParameters& params );
	~MSIsotopeSearch ();
	void printBodyHTML ( std::ostream& os );
	void printBodyXML ( std::ostream& os );
	void printBodyTabDelimitedText ( std::ostream& os );
	void printParamsBodyHTML ( std::ostream& os ) const;
};

#endif /* ! __lu_iso_srch_h */
