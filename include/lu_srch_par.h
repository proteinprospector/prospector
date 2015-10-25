/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_srch_par.h                                                 *
*                                                                             *
*  Created    : June 16th 2001                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_srch_par_h
#define __lu_srch_par_h

#include <lu_fas_enz.h>
#include <lu_pre_srch.h>
#include <lu_fasta.h>
#include <lu_mass.h>
#include <lu_prog_par.h>
#include <lu_sim_ent.h>

class MSSearchParameters : public MSProgramParameters {

	std::string prefix;
protected:
	AAInitInfo aaInitInfo;
private:
	StringVector database;
protected:
	EnzymeParameters enzymeParameters;
private:
	bool tempOverride;
	int maxReportedHits;
	int maxHits;
	std::string comment;
	PreSearchInfo preSearchInfo;
	bool preSearchOnly;
	bool detailedReport;
	int dnaFrameTranslation;
	int maxNTermAA;
	StringVector dbPath;

	SingleEntryParameters singleEntryParameters;
	mutable bool preSearchDone;
public:
	MSSearchParameters ( const ParameterList* params, const std::string& prefix = "", const StringVector& db = StringVector () );
	virtual ~MSSearchParameters ();

	AAInitInfo getAAInitInfo () const { return aaInitInfo; }
	StringVector getDatabase () const { return database; }
	EnzymeParameters getEnzymeParameters () const { return enzymeParameters; }
	std::string getEnzyme () const { return enzymeParameters.getEnzyme (); }
	bool getNoEnzyme () const { return enzymeParameters.getNoEnzyme (); }
	std::string getAllowNonSpecific () const { return enzymeParameters.getAllowNonSpecific (); }
	int getMissedCleavages () const { return enzymeParameters.getMissedCleavages (); }
	bool getTempOverride () const { return tempOverride; }
	int getMaxReportedHits () const { return maxReportedHits; }
	int getMaxHits () const { return maxHits; }
	std::string getComment () const { return comment; }
	bool getPreSearchOnly () const { return preSearchOnly; }
	bool getDetailedReport () const { return detailedReport; }
	int getDNAFrameTranslation () const { return dnaFrameTranslation; }
	int getMaxNTermAA () const { return maxNTermAA; }

	const IntVector& getIndicies ( int num ) const { return preSearchInfo.getIndicies ( num ); }
	int getNumIndicies ( int num ) const { return preSearchInfo.getNumIndicies ( num ); }
	void doSearch ( const FastaServerPtrVector fs ) const;
	void printPreSearchHTML ( std::ostream& os ) const { preSearchInfo.printHTML ( os ); }
	void printPreSearchBodyXML ( std::ostream& os ) const { preSearchInfo.printBodyXML ( os ); }
	bool getPreSearchInfoFromFile () const { return preSearchInfo.getFromFile (); }

	void printHTML ( std::ostream& os ) const;
	void putCGI ( std::ostream& os, bool noEnzyme = false ) const;
	void putHiddenFormJavascriptEntry ( std::ostream& os, bool noEnzyme = false ) const;

	PairStringBool createTemporaryDatabase ( bool randomFlag, bool reverseFlag ) const { return singleEntryParameters.createTemporaryDatabase ( randomFlag, reverseFlag ); }
};

#endif /* ! __lu_srch_par_h */
