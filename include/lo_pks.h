/******************************************************************************
*                                                                             *
*  Library    : liboracle                                                     *
*                                                                             *
*  Filename   : lo_pks.h                                                      *
*                                                                             *
*  Created    :                                                               *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2002-2009) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifndef __lo_pks_h
#define __lo_pks_h

#include <vector>
#include <oci.h>
#include <lo_init.h>
class JobRunItems;

class Peak {
	double mass;
	double intensity;
	double area;
	bool monoisotopic;
	friend class Spectrum;
public:
	Peak ( double mass, double intensity, double area, bool monoisotopic ) :
		mass ( mass ), intensity ( intensity ),
		area ( area ), monoisotopic ( monoisotopic ) {}
	double getMass () const { return mass; }
};

class Spectrum {
	std::vector <Peak> peaks;
	unsigned int runNumber;
	unsigned int jobRunItem;
	double precursorMass;
	bool msmsFlag;
	bool msFlag;
	std::string spotLabel;

	static void writeCommentLine ( std::ostream& os, const std::string& spotLabel, bool msFlag, int currentRunNumber, int jobRunItem );
	static bool retainIsotopes;
	static double minimumArea;
	friend class SpotSetPeaks;
public:
	Spectrum ( unsigned int runNumber, unsigned int jobRunItem, bool msFlag, double precursorMass, const std::string& spotLabel );
	void write ( std::ostream& os ) const;
	unsigned int getJobRunItem () const { return jobRunItem; }
	std::string getSpotLabel () const { return spotLabel; }
	bool getMsmsFlag () const { return msmsFlag; }
	bool getMsFlag () const { return msFlag; }
	double getPrecursorMass () const { return precursorMass; }
	static void setRetainIsotopes ( bool r ) { retainIsotopes = r; }
	static void setMinimumArea ( double m ) { minimumArea = m; }
};

class SpotSetPeaks : public OracleStatement {
	std::vector <Spectrum*> spectra;
	int currentSpectrum;
	unsigned int currentRunItemID;
	IntVector runNumbers;
	bool saveSpectrum;
	void addPeaks ( const JobRunItems* jobRunItems, const double* mass, const double* intensity, const double* area, const double* isotopeFlag, const double* runItemID, const unsigned int nPeaks );
	static char intensityType;
	static const unsigned int blockSize;
public:
	SpotSetPeaks ( unsigned int nSpotSetID, OracleConnection& oc, const JobRunItems* jobRunItems, const IntVector& runNumber );
	virtual ~SpotSetPeaks ();
	static void setIntensityType ( char i ) { intensityType = i; }
	void writeSpectraSingleFile ( std::ostream& os ) const;
};

#endif /* ! __lo_pks_h */
