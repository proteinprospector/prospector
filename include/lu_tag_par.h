/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_tag_par.h                                                  *
*                                                                             *
*  Created    : July 3rd 2001                                                 *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_tag_par_h
#define __lu_tag_par_h

#include <lu_composit.h>
#include <lu_srch_par.h>
#include <lu_mass_frag.h>
#include <lu_comb_perm.h>
#include <lu_sim_ent.h>
#include <lu_spec_id.h>
#include <lu_mut_mtrx.h>
#include <lu_df_info.h>
#include <lu_pk_filter.h>
#include <lu_tol.h>

class LinkInfo;

class MSTagParameters : public MSSearchParameters {

	MSMSPeakFilterOptions* msmsPeakFilterOptions;

	double parentMass;
	int parentCharge;
	SpecID specID;
	MSMSDataSetInfo* dataSetInfo;

	ToleranceInfo parentMassTolerance;
	ToleranceInfo productMassTolerance;
	double systematicError;
	MassInfo massInfo;

	ConsideredAA consideredAA; 
	CompositionSearchParameters compSearchParams;
	std::string regularExpression;

	ModificationParameters modificationParameters;
	LinkInfo* linkInfo;
	int numSavedCrosslinkHits;

	BiemannParameters biemannParams;

	std::string expectationMethod;

	bool scoreHistogramOnly;
public:
	MSTagParameters ( const ParameterList* params, bool noFilesOK = false );
	~MSTagParameters ();

	MSMSPeakFilterOptions* getMSMSPeakFilterOptions () const { return msmsPeakFilterOptions; }

	MSMSDataSetInfo* getDataSetInfo () { return dataSetInfo; }

	Tolerance* getParentMassTolerance () const { return parentMassTolerance.getTolerance (); }
	Tolerance* getProductMassTolerance () const { return productMassTolerance.getTolerance (); }
	double getParentMassSystematicError () const { return systematicError; }
	bool getMonoisotopicFlag () const { return massInfo.getMonoisotopicFlag (); }
	bool getMonoParentAverageFragments () const { return massInfo.getMonoParentAverageFragments (); }
	bool getAverageParentMonoFragments () const { return massInfo.getAverageParentMonoFragments (); }

	std::string getAAList () const { return consideredAA.getAAList (); }
	CompositionSearchParameters getCompSearchParams () const { return compSearchParams; }
	std::string getRegularExpression () const { return regularExpression; }
	bool getRegExpSearch () const { return regularExpression != "."; }

	ModificationParameters getModificationParameters () const { return modificationParameters; }
	double getMassModTolerance () const { return ( 1.0 + modificationParameters.getDefect () ) / 2.0; }
	bool getAllowErrors () const { return modificationParameters.getAllowErrors (); }
	const BiemannParameters& getBiemannParameters () const { return biemannParams; }
	std::string getIonSeries () const { return biemannParams.getIonSeries (); }

	unsigned int getLinkAminoAcid ( int num ) const;
	std::string getBridgeFormula () const;
	LinkInfo* getLinkInfo () const { return linkInfo; }
	bool isCrosslinking () const;

	std::string getExpectationMethod () const { return expectationMethod; }

	void printHTML ( std::ostream& os ) const;
	char getAllowNonSpecificType () const;

	bool getScoreHistogramOnly () const { return scoreHistogramOnly; }

	int getNumSavedHits () const;

	int getNumSavedCrosslinkHits () const { return numSavedCrosslinkHits; }

	void setDataSetInfo ( const ParameterList* params, bool noFilesOK = false );
};

#endif /* ! __lu_tag_par_h */
