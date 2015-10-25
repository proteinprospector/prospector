/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_sing_fit.h                                                 *
*                                                                             *
*  Created    : October 29th 2001                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_sing_fit_h
#define __lu_sing_fit_h

#include <lu_charge.h>
#include <lu_prod_par.h>
#include <lu_iso_par.h>
#include <lu_coverage.h>

class SingleFitSearch {
protected:
	const PeakContainer& peaks;
	MSProductLink productLink;
	MSIsotopeLink isotopeLink;
	PeakMatchContext peakMatchContext;
	BoolDeque peakUsed;
	std::vector <CoverageMap> coverageMap;
	DoubleVector errors;
	IntVector coverCount;
	double aveError;
	double sdevError;
	double percentTic;
protected:
	virtual bool areHits () const = 0;
	int numUnique;
public:
	SingleFitSearch ( const PeakContainer& peaks );
	virtual ~SingleFitSearch ();
	BoolDeque getPeakUsed () const { return peakUsed; }
	int getNumUnique () const { return numUnique; }
	double getPercentCoverage ( int i ) const { return coverageMap [i].getPercentCoverage (); }
	double getMeanError () const { return aveError; }
	double getTolerance () const { return 2 * sdevError; }
	CoverageMap getCoverageMap ( int i ) const { return coverageMap [i]; }
	CharVector getAACovered ( int i ) const { return coverageMap [i].getAACovered (); }
	virtual void printHTML ( std::ostream& os, bool hideLinks );
	virtual void printHTMLBody ( std::ostream& os, bool hideLinks ) = 0;
	virtual void printXML ( std::ostream& os );
	virtual void printXMLBody ( std::ostream& os ) = 0;
	void calculateTIC ( const PeakContainer& peaks );
	void calculateStats ();
	virtual void calculateNumUnique () = 0;
	void printCoverageHTML ( std::ostream& os ) const;
	void printCoverCountHTML ( std::ostream& os, int i ) const;
	void printStatsHeaderHTML ( std::ostream& os ) const;
	void printStatsHTML ( std::ostream& os, int i ) const;
	void printStatsXML ( std::ostream& os ) const;
};

#endif /* ! __lu_sing_fit_h */
