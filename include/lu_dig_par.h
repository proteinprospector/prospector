/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_dig_par.h                                                  *
*                                                                             *
*  Created    : June 13th 2001                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_dig_par_h
#define __lu_dig_par_h

#include <lg_io.h>
#include <lu_coverage.h>
#include <lu_prog_par.h>
#include <lu_fas_enz.h>
#include <lu_sim_ent.h>
#include <lu_mod_frag.h>
#include <lu_usermod.h>
#ifdef CHEM_SCORE
#include <lu_chem_sc.h>
#endif

class MSDigestParameters : public MSProgramParameters {
protected:
	AAInitInfo aaInitInfo;
private:
	SingleEntryParameters singleEntryParameters;
	EnzymeParameters enzymeParameters;
	bool hideProteinSequence;
	bool reportMultCharge;
	bool hideHtmlLinks;
	DigestFragmentParameters digestFragmentParameters;
	std::string instrumentName;
	MultipleModification2 multipleModification;
#ifdef CHEM_SCORE
	ChemScoreParameters chemScoreParams;
#endif
	bool separateProteinsFlag;
	bool bullBreeseFlag;
	bool HPLCFlag;
	CharVector coverageMap;
	int maxHits;
public:
	MSDigestParameters ( const ParameterList* params );

	AAInitInfo getAAInitInfo () const { return aaInitInfo; }
	SingleEntryParameters getSingleEntryParameters () const { return singleEntryParameters; }
	EnzymeParameters getEnzymeParameters () const { return enzymeParameters; }
	bool getHideProteinSequence () const { return hideProteinSequence; }
	bool getReportMultCharge () const { return reportMultCharge; }
	bool getHideHTMLLinks () const { return hideHtmlLinks; }
	DigestFragmentParameters getDigestFragmentParameters () const { return digestFragmentParameters; }
	std::string getInstrumentName () const { return instrumentName; }
	bool getSeparateProteinsFlag () const { return separateProteinsFlag; }
	CharVector getCoverageMap () const { return coverageMap; }
	int getMaxHits () const { return maxHits; }

	virtual void printHTML ( std::ostream& os ) const;
};

class MSDigestLink {
	static bool init;
	static std::string urlParams;
	static std::string linkName;
	static StringVector database;
	static int num;
	int ind;
	void initLinks ( const std::string& programName );
public:
	MSDigestLink ( const std::string& programName, const std::string& db = "" );
	void printHTML ( std::ostream& os ) const;
	static void write ( std::ostream& os, int indexNumber, int dnaReadingFrame, int openReadingFrame, int num );
	static void write ( std::ostream& os, int indexNumber, int dnaReadingFrame, int openReadingFrame, const CoverageMap& coverageMap, int num );
};

class MSDigestLinkNameValueStream : public GenCommentedIFStream {
	MapStringToStringVector params;
public:
	MSDigestLinkNameValueStream ( bool process );
	std::string getStringValue ( const std::string& name, const std::string& defaultVal = "" ) const;
	StringVector getNameList () const;
};
#endif /* ! __lu_dig_par_h */
