/******************************************************************************
*                                                                             *
*  Program    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_prod_par.h                                                 *
*                                                                             *
*  Created    : June 18th 2001                                                *
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

#ifndef __lu_prod_par_h
#define __lu_prod_par_h

#include <lg_string.h>
#include <lu_mass.h>
#include <lu_prog.h>
#include <lu_prog_par.h>
#include <lu_spec_id.h>
#include <lu_pk_filter.h>
#include <lu_tol.h>
#include <lu_df_info.h>

class BiemannParameters;
class LinkInfo;

class MSProductParameters : public MSProgramParameters  {

	MSMSPeakFilterOptions* msmsPeakFilterOptions;
	SpecID specID;
	MSMSDataSetInfo dataSetInfo;
	ToleranceInfo parentMassTolerance;
	ToleranceInfo productMassTolerance;
	MassInfo* massInfo;

	std::vector <ElementalFormula> userAAElemForms;

	int maxCharge;
	bool calibrate;
	double calTolerance;
	std::string dataPlotted;
	std::string instrumentName;

	BiemannParameters* biemannParams;
	StringVectorVector sequence;
	StringVector nTermName;
	StringVector cTermName;
	StringVector nLossName;
	BoolDeque sFlag;
	std::vector <MapStringConstModPtr> constMods;
	bool alternative;
	bool discriminating;
	LinkInfo* linkInfo;

	static std::string MAX_CHARGE;
	static std::string MAX_CHARGE_DEFAULT;
	static void printPSILink ( std::ostream& os, const std::string& s, bool term, const std::string& mainURL, const std::string& startRange, const std::string& endRange );
	void printSequence ( std::ostream& os, const StringVector& seq, const std::string& nName, const std::string& cName, const std::string& lName, int charge, int colour ) const;
	static bool getSequenceElements ( ParameterList* params, const std::string& num, std::string& sequence, std::string& nTermName, std::string& cTermName, std::string& neutralLossName );
public:
	MSProductParameters ( ParameterList* params );
	~MSProductParameters ();
	int initMaxCharge ( const std::string& maxChargeString );

	StringVector getSequence ( StringVectorVectorSizeType idx ) const
	{
		return idx < sequence.size () ? sequence [idx] : StringVector ();
	}

	MSMSPeakFilterOptions* getMSMSPeakFilterOptions () const { return msmsPeakFilterOptions; }
	MSMSDataSetInfo* getDataSetInfo () { return &dataSetInfo; }
	Tolerance* getParentMassTolerance () const { return parentMassTolerance.getTolerance (); }
	Tolerance* getProductMassTolerance () const { return productMassTolerance.getTolerance (); }
	bool getMonoisotopicFlag () const { return massInfo->getMonoisotopicFlag (); }
	bool getMonoParentAverageFragments () const { return massInfo->getMonoParentAverageFragments (); }
	bool getAverageParentMonoFragments () const { return massInfo->getAverageParentMonoFragments (); }

	const BiemannParameters* getBiemannParameters () const { return biemannParams; }
	int getMaxCharge () const { return maxCharge; }
	bool getCalibrate () const { return calibrate; }
	double getCalTolerance () const { return calTolerance; }
	bool getRawData () const { return dataPlotted == "Raw" || dataPlotted == "Raw and Centroid"; }
	bool getPlotCentroids () const { return dataPlotted == "Raw and Centroid"; }
	bool getAlternative () const { return alternative; }
	bool getDiscriminating () const { return discriminating; }
	LinkInfo* getLinkInfo () const { return linkInfo; }
	bool getShow ( int i ) const { return sFlag.empty () || i >= sFlag.size () ? true : sFlag [i]; }

	MapStringConstModPtr getConstMods ( int idx ) const
	{
		return constMods [idx];
	}

	SpecID getSpecID () { return specID; }
	int getNumSequences () const { return sequence.size (); }
	void printHTML ( std::ostream& os );
};

class PeptideSequence;

class MSProductLink {
	static int num;
	int index;
	std::string getLinkName () const
	{
		if ( num == 0 ) return "msprodLink";
		else			return "msprodLink" + gen_itoa ( index );
	}

	const MSMSDataSetInfo* dataSetInfo;
	int dataSet;

	std::string dataFilename;

	std::string instrumentName;
	std::string searchKey;
	SpecID specID;
	const Tolerance* pTol;
	const Tolerance* tol;

	void start ( std::ostream& os, bool url ) const;
	void end ( std::ostream& os, const std::string& sequence, bool url ) const;
	void end ( std::ostream& os, bool url, const std::string& sequence, const std::string& sequence2 = "" ) const;
	void putCGI ( std::ostream& os ) const;
	void writeSequenceParams ( std::ostream& os, const std::string& peptide ) const;
	void writeLink ( std::ostream& os, const std::string& sequence, const std::string& mods, int maxCharge, bool url ) const;
public:
	MSProductLink ( const std::string& programName, const MSMSDataSetInfo* dataSetInfo = 0, int dataSet = -1 );
	MSProductLink ( const std::string& instrumentName, const std::string& projectName, const SpecID& specID, const Tolerance* pTol, const Tolerance* tol );
	MSProductLink ( const std::string& instrumentName, const std::string& dataFilename, const Tolerance* pTol, const Tolerance* tol );
	void write1 ( std::ostream& os, const PeptideSequence& ps, bool hideLinks, const int maxCharge, char aa ) const;
	void write2 ( std::ostream& os, const std::string& sequence, const PeptideSequence& ps, bool hideLinks, const int maxCharge = 0, double err = 0.0 ) const;
	void write2 ( std::ostream& os, const CharVector& previousAA, const StringVector& sequence, const CharVector& nextAA, const std::vector <PeptideSequence>& ps, bool hideLinks, const int maxCharge, const LinkInfo* linkInfo ) const;
	void write3 ( std::ostream& os, const std::string& label, int maxCharge ) const;
	void write4 ( std::ostream& os, const std::string& sequence, const std::string& mods, int maxCharge ) const;
	void write4 ( std::ostream& os, const std::string& dbPeptide, const std::string& sequence, const std::string& nTermName, const std::string& cTermName, const std::string& neutralLossName, const std::string& mods, int maxCharge ) const;
	void write4 ( std::ostream& os, const StringVector& previousAA, const StringVector& sequence, const StringVector& nTermName, const StringVector& cTermName, const StringVector& neutralLossName, const StringVector& nextAA, bool hideLinks, const int maxCharge, const LinkInfo* linkInfo ) const;
	void write5 ( std::ostream& os, const SpecID& spID, const std::string& sequence, const std::string& mods, int maxCharge, bool url ) const;				// Viewer
	void write5 ( std::ostream& os, const SpecID& spID, const std::string& sequence1, const std::string& sequence2, int maxCharge, bool url, const std::string& linkSearchType, const std::string& bridgeComposition ) const;		// Viewer
	void write5 ( std::ostream& os, const std::string& title, const std::string& sequence, const std::string& mods, int maxCharge, bool url ) const;		// Viewer
	void write6 ( std::ostream& os, const std::string& index, const std::string& sequence, const std::string& mods, int maxCharge, bool url ) const;		// Viewer
	void write7 ( std::ostream& os, const std::string& mOverZ, const std::string& sequence, const std::string& mods, int maxCharge, bool url ) const;		// Viewer
	void write8 ( std::ostream& os, const std::string& scanNumber, const std::string& sequence, const std::string& mods, int maxCharge, bool url ) const;	// Viewer
	void printHTML ( std::ostream& os ) const;
};

#endif /* ! __lu_prod_par_h */
