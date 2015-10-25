/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_tag_frag.h                                                 *
*                                                                             *
*  Created    : July 21st 1996                                                *
*                                                                             *
*  Purpose    : Contains most of the search algorithms for MS-Tag.            *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1996-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_tag_frag_h
#define __lu_tag_frag_h

#include <lgen_define.h>
#include <lu_mut_mtrx.h>
#include <lu_charge.h>
#include <lu_histogram.h>

class BiemannParameters;
class UnmatchedCompositionSearch;
class MSTagParameters;
class CompositionSearch;

class TagMatch {
	Peak* parentPeak;
	int unmatched;
	ScoreType score;
	PeptideSequence mutatedSequence;
public:
	TagMatch ( Peak* parentPeak, int unmatched, ScoreType score, const PeptideSequence& mutatedSequence ) :
		parentPeak ( parentPeak ),
		unmatched ( unmatched ),
		score ( score ),
		mutatedSequence ( mutatedSequence ) {}
	TagMatch ( Peak* parentPeak, int unmatched, ScoreType score, const std::string& peptide ) :
		parentPeak ( parentPeak ),
		unmatched ( unmatched ),
		score ( score ),
		mutatedSequence ( peptide ) {}

	Peak* getParentPeak () const { return parentPeak; }
	int getUnmatched () const { return unmatched; }
	ScoreType getScore () const { return score; }
	PeptideSequence getMutatedSequence () const { return mutatedSequence; }
	SetPairIntString getModIndicies () const { return mutatedSequence.getModIndicies (); }
};
typedef std::vector <TagMatch> TagMatchVector;
typedef TagMatchVector::size_type TagMatchVectorSizeType;

class MSMSSearch {
protected:
	double parentTolerance;
	Peak parentPeak;
	PeakContainer peaks;
	double parentMass;
	double parentMassPlusNegTolerance;
	double parentMassMinusPosTolerance;
	CompositionSearch* compositionSearch;
	UnmatchedCompositionSearch* unmatchedCompositionSearch;
	double getParentTolerance ( MSMSDataPoint* file );
public:
	MSMSSearch ( MSMSDataPoint* file, MSTagParameters& params );
	virtual ~MSMSSearch ();
	const Peak* getParentPeak () const { return &parentPeak; }
	const PeakContainer* getPeaks () const { return ( &peaks );}
	double getParentMass () const { return ( parentMass );}
	virtual bool getSpectrumRetained () const { return peaks.getSpectrumRetained (); }
	virtual bool checkOffset ( double moleWt ) const { return true;}
	virtual bool checkMatch ( double molWt ) const { return true;}
	virtual double getParentTolerance () const { return ( parentTolerance );}
	double getParentMassPlusNegTolerance () const { return ( parentMassPlusNegTolerance );}
	double getParentMassMinusPosTolerance () const { return ( parentMassMinusPosTolerance );}
	virtual CompositionSearch* getCompositionSearch () const { return ( compositionSearch );}
	virtual bool doMatch ( const std::string& peptide, bool nTermPeptide, bool cTermPeptide, double mol_wt, TagMatchVector& tagMatch, const ScoreType& minScore );
	virtual void printExpectationHTML ( std::ostream& os, double score, int numSavedSpectra ) {}
	virtual void printExpectationXML ( std::ostream& os, int numSavedSpectra ) {}
};
typedef std::vector<MSMSSearch*>::iterator MSMSSearchIterator;
typedef std::vector<MSMSSearch*>::const_iterator MSMSSearchConstIterator;

class MSSeqSearch : public MSMSSearch {
protected:
	double tagMass1;
	double tagMass1Tolerance;
	DoubleVector seqTag;
	DoubleVector seqTagTolerance;
	int numMasses;
	int numPeaks;
	double tagMass3;
	double tagMass3Tolerance;
public:
	MSSeqSearch ( MSMSDataPoint* file, MSTagParameters& params );
	virtual ~MSSeqSearch ();
	bool doMatch ( const std::string& peptide, bool nTermPeptide, bool cTermPeptide, double mol_wt, TagMatchVector& tagMatch, const ScoreType& minScore );
	virtual bool doSeqMatch ( const std::string& peptide, double mol_wt, TagMatchVector& tagMatch );
	double getParentMassPlusNegTolerance () const { return ( parentMass + parentTolerance + 135.0 ); }
	double getParentMassMinusPosTolerance () const { return ( parentMass - ( parentTolerance + 135.0 ) ); }
};

class MSSeqMass1Mass3Search : public MSSeqSearch {
public:
	MSSeqMass1Mass3Search ( MSMSDataPoint* file, MSTagParameters& params ) :
		MSSeqSearch ( file, params ) {}
	virtual ~MSSeqMass1Mass3Search ();
	virtual bool doSeqMatch ( const std::string& peptide, double mol_wt, TagMatchVector& tagMatch );
	double getParentTolerance () const { return ( parentTolerance + 135.0 );}
};

class MSSeqMass1SeqSearch : public MSSeqSearch {
public:
	MSSeqMass1SeqSearch ( MSMSDataPoint* file, MSTagParameters& params ) :
		MSSeqSearch ( file, params ) {}
	virtual ~MSSeqMass1SeqSearch ();
	virtual bool doSeqMatch ( const std::string& peptide, double mol_wt, TagMatchVector& tagMatch );
	double getParentTolerance () const { return ( parentTolerance + 135.0 );}
};

class MSSeqMass3SeqSearch : public MSSeqSearch {
public:
	MSSeqMass3SeqSearch ( MSMSDataPoint* file, MSTagParameters& params ) :
		MSSeqSearch ( file, params ) {}
	virtual ~MSSeqMass3SeqSearch ();
	virtual bool doSeqMatch ( const std::string& peptide, double mol_wt, TagMatchVector& tagMatch );
	double getParentTolerance () const { return ( parentTolerance + 135.0 );}
};

class LinkInfo;

class MSTagSearch : public MSMSSearch {
	friend class UnmatchedCompositionSearch;
protected:
	MassTypeVectorVector nIonFragTag;
	MassTypeVectorVector cIonFragTag;

	MassTypeVector maxNFragTag;
	MassTypeVector maxCFragTag;
	MassType maxNFragTagInData;
	MassType maxCFragTagInData;
	MassTypeVector minNFragTag;
	MassTypeVector minCFragTag;
	MassType minNFragTagInData;
	MassType minCFragTagInData;
	int numPeaks;
	MassTypeVector fragTolerance;

	ScoreTypeVector ionMatched;

	ScoreType previousScore;
	int previousNumUnmatchedIons;
	MassType previousNTerminusWt;
	MassType previousCTerminusWt;
	std::string previousSequence;

	bool doublyChargedIons;
	bool triplyChargedIons;
	int precursorCharge;
	bool spectrumRetained;
	SurvivalHistogram survHist;

	static void initXLVariables ( const LinkInfo* linkInfo );
	static void initFragTagFlags ( const BiemannParameters& bp );
	void initFragTags ();
	void xLinkMatchCID ( const PeptideSequence& ps, double molWt, double diff );
	void xLinkMatchCID2 ( const PeptideSequence& ps );
	void xLinkMatchETD ( const PeptideSequence& ps, double molWt, double diff );
	void xLinkMatchETD2 ( const PeptideSequence& ps );
	void fragmentMatch ( const std::string& sequence, int& numUnmatchedIons, ScoreType& score, double molWt = 0.0, double diff = 0.0 );
	void fragmentMatch ( const PeptideSequence& ps, int& numUnmatchedIons, ScoreType& score, double molWt = 0.0, double diff = 0.0 );
	void doNIons ( const std::string& sequence );
	void doNIonsESI_TRAP_CID_low_res ( const std::string& sequence );
	void doNIonsESI_Q_CID ( const std::string& sequence );
	void doNIonsESI_ETD_low_res ( const std::string& sequence );
	void doNIonsESI_ETD_high_res ( const std::string& sequence );
	void doCIons ( const std::string& sequence );
	void doCIonsESI_TRAP_CID_low_res ( const std::string& sequence );
	void doCIonsESI_Q_CID ( const std::string& sequence );
	void doCIonsESI_ETD_low_res ( const std::string& sequence );
	void doCIonsESI_ETD_high_res ( const std::string& sequence );
	void doInternalIons ( const std::string& sequence );

	static bool ESI_TRAP_CID_low_res;
	static bool ESI_Q_CID;
	static bool ESI_ETD_low_res;
	static bool ESI_ETD_high_res;
	static bool etd;

	static MassType internalNTerminusWt;
	static MassType nTerminusWt;
	static MassType cTerminusWt;
	static bool flagsSet;
	static bool a_nh3_flag;
	static bool a_h2o_flag;
	static bool a_flag;
	static bool a_h3po4_flag;
	static bool b_h2o_flag;
	static bool b_plus_h2o_flag;
	static bool c_ladder_flag;
	static bool b_nh3_flag;
	static bool b_soch4_flag;
	static bool b_h3po4_flag;
	static bool b_flag;
	static bool cPlus2DaFlag;
	static bool cPlus1DaFlag;
	static bool c_flag;
	static bool cp2_flag;
	static bool cMinus1DaFlag;
	static bool d_flag;

	static bool v_flag;
	static bool w_flag;
	static bool x_flag;
	static bool y_h2o_flag;
	static bool y_nh3_flag;
	static bool y_soch4_flag;
	static bool y_h3po4_flag;
	static bool y_flag;
	static bool Y_flag;
	static bool z_flag;
	static bool zp2_flag;
	static bool zPlus1Dap2_flag;
	static bool zPlus1DaFlag;
	static bool zPlus2DaFlag;
	static bool zPlus3DaFlag;
	static bool n_ladder_flag;

	static bool bp2_flag;
	static bool bp3_flag;
	static bool bp2_nh3_flag;
	static bool bp2_h2o_flag;
	static bool bp2_soch4_flag;
	static bool bp2_h3po4_flag;
	static bool yp2_flag;
	static bool yp3_flag;
	static bool yp2_nh3_flag;
	static bool yp2_h2o_flag;
	static bool yp2_soch4_flag;
	static bool yp2_h3po4_flag;

	static ScoreType a_nh3_score;
	static ScoreType a_h2o_score;
	static ScoreType a_score;
	static ScoreType a_h3po4_score;
	static ScoreType b_h2o_score;
	static ScoreType b_plus_h2o_score;
	static ScoreType c_ladder_score;
	static ScoreType b_nh3_score;
	static ScoreType b_soch4_score;
	static ScoreType b_h3po4_score;
	static ScoreType b_score;
	static ScoreType cPlus2DaScore;
	static ScoreType cPlus1DaScore;
	static ScoreType c_score;
	static ScoreType cMinus1DaScore;
	static ScoreType d_score;

	static ScoreType v_score;
	static ScoreType w_score;
	static ScoreType x_score;
	static ScoreType y_h2o_score;
	static ScoreType y_nh3_score;
	static ScoreType y_soch4_score;
	static ScoreType y_h3po4_score;
	static ScoreType y_score;
	static ScoreType Y_score;
	static ScoreType z_score;
	static ScoreType zPlus1DaScore;
	static ScoreType zPlus2DaScore;
	static ScoreType zPlus3DaScore;
	static ScoreType n_ladder_score;

	static ScoreType bp2_score;
	static ScoreType bp3_score;
	static ScoreType bp2_nh3_score;
	static ScoreType bp2_h2o_score;
	static ScoreType bp2_soch4_score;
	static ScoreType bp2_h3po4_score;
	static ScoreType cp2_score;

	static ScoreType yp2_score;
	static ScoreType yp3_score;
	static ScoreType yp2_nh3_score;
	static ScoreType yp2_h2o_score;
	static ScoreType yp2_soch4_score;
	static ScoreType yp2_h3po4_score;
	static ScoreType zp2_score;
	static ScoreType zPlus1Dap2Score;

	static ScoreType internal_a_score;
	static ScoreType internal_b_score;
	static ScoreType internal_loss_score;

	static ScoreType unmatched_score;

	static ScoreType ESI_ETD_low_res_b_scoreZ2NC;
	static ScoreType ESI_ETD_low_res_c_scoreZ2NC;
	static ScoreType ESI_ETD_low_res_cMinus1DaZ2NC;
	static ScoreType ESI_ETD_low_res_y_scoreZ2NC;
	static ScoreType ESI_ETD_low_res_z_scoreZ2NC;
	static ScoreType ESI_ETD_low_res_zPlus1DaZ2NC;

	static ScoreType ESI_ETD_low_res_b_scoreZ2N;
	static ScoreType ESI_ETD_low_res_c_scoreZ2N;
	static ScoreType ESI_ETD_low_res_cMinus1DaZ2N;
	static ScoreType ESI_ETD_low_res_y_scoreZ2N;
	static ScoreType ESI_ETD_low_res_z_scoreZ2N;
	static ScoreType ESI_ETD_low_res_zPlus1DaZ2N;

	static ScoreType ESI_ETD_low_res_b_scoreZ2C;
	static ScoreType ESI_ETD_low_res_c_scoreZ2C;
	static ScoreType ESI_ETD_low_res_cMinus1DaZ2C;
	static ScoreType ESI_ETD_low_res_y_scoreZ2C;
	static ScoreType ESI_ETD_low_res_z_scoreZ2C;
	static ScoreType ESI_ETD_low_res_zPlus1DaZ2C;

	static ScoreType ESI_ETD_low_res_b_scoreZ2;
	static ScoreType ESI_ETD_low_res_c_scoreZ2;
	static ScoreType ESI_ETD_low_res_cMinus1DaZ2;
	static ScoreType ESI_ETD_low_res_y_scoreZ2;
	static ScoreType ESI_ETD_low_res_z_scoreZ2;
	static ScoreType ESI_ETD_low_res_zPlus1DaZ2;

	static ScoreType ESI_ETD_low_res_c_scoreZ3NC;
	static ScoreType ESI_ETD_low_res_cp2_scoreZ3NC;
	static ScoreType ESI_ETD_low_res_cMinus1DaZ3NC;
	static ScoreType ESI_ETD_low_res_y_scoreZ3NC;
	static ScoreType ESI_ETD_low_res_z_scoreZ3NC;
	static ScoreType ESI_ETD_low_res_zPlus1DaZ3NC;
	static ScoreType ESI_ETD_low_res_zp2_scoreZ3NC;
	static ScoreType ESI_ETD_low_res_zPlus1Dap2Z3NC;

	static ScoreType ESI_ETD_low_res_c_scoreZ3N;
	static ScoreType ESI_ETD_low_res_cp2_scoreZ3N;
	static ScoreType ESI_ETD_low_res_cMinus1DaZ3N;
	static ScoreType ESI_ETD_low_res_y_scoreZ3N;
	static ScoreType ESI_ETD_low_res_z_scoreZ3N;
	static ScoreType ESI_ETD_low_res_zPlus1DaZ3N;
	static ScoreType ESI_ETD_low_res_zp2_scoreZ3N;
	static ScoreType ESI_ETD_low_res_zPlus1Dap2Z3N;

	static ScoreType ESI_ETD_low_res_c_scoreZ3C;
	static ScoreType ESI_ETD_low_res_cp2_scoreZ3C;
	static ScoreType ESI_ETD_low_res_cMinus1DaZ3C;
	static ScoreType ESI_ETD_low_res_y_scoreZ3C;
	static ScoreType ESI_ETD_low_res_z_scoreZ3C;
	static ScoreType ESI_ETD_low_res_zPlus1DaZ3C;
	static ScoreType ESI_ETD_low_res_zp2_scoreZ3C;
	static ScoreType ESI_ETD_low_res_zPlus1Dap2Z3C;

	static ScoreType ESI_ETD_low_res_c_scoreZ3;
	static ScoreType ESI_ETD_low_res_cp2_scoreZ3;
	static ScoreType ESI_ETD_low_res_cMinus1DaZ3;
	static ScoreType ESI_ETD_low_res_y_scoreZ3;
	static ScoreType ESI_ETD_low_res_z_scoreZ3;
	static ScoreType ESI_ETD_low_res_zPlus1DaZ3;
	static ScoreType ESI_ETD_low_res_zp2_scoreZ3;
	static ScoreType ESI_ETD_low_res_zPlus1Dap2Z3;

	static ScoreType ESI_ETD_low_res_c_scoreZ4NC;
	static ScoreType ESI_ETD_low_res_cp2_scoreZ4NC;
	static ScoreType ESI_ETD_low_res_y_scoreZ4NC;
	static ScoreType ESI_ETD_low_res_z_scoreZ4NC;
	static ScoreType ESI_ETD_low_res_zPlus1DaZ4NC;
	static ScoreType ESI_ETD_low_res_zp2_scoreZ4NC;
	static ScoreType ESI_ETD_low_res_zPlus1Dap2Z4NC;

	static ScoreType ESI_ETD_low_res_c_scoreZ4N;
	static ScoreType ESI_ETD_low_res_cp2_scoreZ4N;
	static ScoreType ESI_ETD_low_res_y_scoreZ4N;
	static ScoreType ESI_ETD_low_res_z_scoreZ4N;
	static ScoreType ESI_ETD_low_res_zPlus1DaZ4N;
	static ScoreType ESI_ETD_low_res_zp2_scoreZ4N;
	static ScoreType ESI_ETD_low_res_zPlus1Dap2Z4N;

	static ScoreType ESI_ETD_low_res_c_scoreZ4C;
	static ScoreType ESI_ETD_low_res_cp2_scoreZ4C;
	static ScoreType ESI_ETD_low_res_y_scoreZ4C;
	static ScoreType ESI_ETD_low_res_z_scoreZ4C;
	static ScoreType ESI_ETD_low_res_zPlus1DaZ4C;
	static ScoreType ESI_ETD_low_res_zp2_scoreZ4C;
	static ScoreType ESI_ETD_low_res_zPlus1Dap2Z4C;

	static ScoreType ESI_ETD_low_res_c_scoreZ4;
	static ScoreType ESI_ETD_low_res_cp2_scoreZ4;
	static ScoreType ESI_ETD_low_res_y_scoreZ4;
	static ScoreType ESI_ETD_low_res_z_scoreZ4;
	static ScoreType ESI_ETD_low_res_zPlus1DaZ4;
	static ScoreType ESI_ETD_low_res_zp2_scoreZ4;
	static ScoreType ESI_ETD_low_res_zPlus1Dap2Z4;

	static ScoreType ESI_ETD_low_res_c_scoreZ5NC;
	static ScoreType ESI_ETD_low_res_cp2_scoreZ5NC;
	static ScoreType ESI_ETD_low_res_y_scoreZ5NC;
	static ScoreType ESI_ETD_low_res_z_scoreZ5NC;
	static ScoreType ESI_ETD_low_res_zPlus1DaZ5NC;
	static ScoreType ESI_ETD_low_res_zp2_scoreZ5NC;
	static ScoreType ESI_ETD_low_res_zPlus1Dap2Z5NC;

	static ScoreType ESI_ETD_low_res_c_scoreZ5N;
	static ScoreType ESI_ETD_low_res_cp2_scoreZ5N;
	static ScoreType ESI_ETD_low_res_y_scoreZ5N;
	static ScoreType ESI_ETD_low_res_z_scoreZ5N;
	static ScoreType ESI_ETD_low_res_zPlus1DaZ5N;
	static ScoreType ESI_ETD_low_res_zp2_scoreZ5N;
	static ScoreType ESI_ETD_low_res_zPlus1Dap2Z5N;

	static ScoreType ESI_ETD_low_res_c_scoreZ5C;
	static ScoreType ESI_ETD_low_res_cp2_scoreZ5C;
	static ScoreType ESI_ETD_low_res_y_scoreZ5C;
	static ScoreType ESI_ETD_low_res_z_scoreZ5C;
	static ScoreType ESI_ETD_low_res_zPlus1DaZ5C;
	static ScoreType ESI_ETD_low_res_zp2_scoreZ5C;
	static ScoreType ESI_ETD_low_res_zPlus1Dap2Z5C;

	static ScoreType ESI_ETD_low_res_c_scoreZ5;
	static ScoreType ESI_ETD_low_res_cp2_scoreZ5;
	static ScoreType ESI_ETD_low_res_y_scoreZ5;
	static ScoreType ESI_ETD_low_res_z_scoreZ5;
	static ScoreType ESI_ETD_low_res_zPlus1DaZ5;
	static ScoreType ESI_ETD_low_res_zp2_scoreZ5;
	static ScoreType ESI_ETD_low_res_zPlus1Dap2Z5;

	static ScoreType ESI_ETD_hi_res_b_scoreZ2NC;
	static ScoreType ESI_ETD_hi_res_c_scoreZ2NC;
	static ScoreType ESI_ETD_hi_res_cMinus1DaZ2NC;
	static ScoreType ESI_ETD_hi_res_y_scoreZ2NC;
	static ScoreType ESI_ETD_hi_res_z_scoreZ2NC;
	static ScoreType ESI_ETD_hi_res_zPlus1DaZ2NC;

	static ScoreType ESI_ETD_hi_res_b_scoreZ2N;
	static ScoreType ESI_ETD_hi_res_c_scoreZ2N;
	static ScoreType ESI_ETD_hi_res_cMinus1DaZ2N;
	static ScoreType ESI_ETD_hi_res_y_scoreZ2N;
	static ScoreType ESI_ETD_hi_res_z_scoreZ2N;
	static ScoreType ESI_ETD_hi_res_zPlus1DaZ2N;

	static ScoreType ESI_ETD_hi_res_b_scoreZ2C;
	static ScoreType ESI_ETD_hi_res_c_scoreZ2C;
	static ScoreType ESI_ETD_hi_res_cMinus1DaZ2C;
	static ScoreType ESI_ETD_hi_res_y_scoreZ2C;
	static ScoreType ESI_ETD_hi_res_z_scoreZ2C;
	static ScoreType ESI_ETD_hi_res_zPlus1DaZ2C;

	static ScoreType ESI_ETD_hi_res_b_scoreZ2;
	static ScoreType ESI_ETD_hi_res_c_scoreZ2;
	static ScoreType ESI_ETD_hi_res_cMinus1DaZ2;
	static ScoreType ESI_ETD_hi_res_y_scoreZ2;
	static ScoreType ESI_ETD_hi_res_z_scoreZ2;
	static ScoreType ESI_ETD_hi_res_zPlus1DaZ2;

	static ScoreType ESI_ETD_hi_res_c_scoreZ3NC;
	static ScoreType ESI_ETD_hi_res_cMinus1DaZ3NC;
	static ScoreType ESI_ETD_hi_res_y_scoreZ3NC;
	static ScoreType ESI_ETD_hi_res_z_scoreZ3NC;
	static ScoreType ESI_ETD_hi_res_zPlus1DaZ3NC;

	static ScoreType ESI_ETD_hi_res_c_scoreZ3N;
	static ScoreType ESI_ETD_hi_res_cMinus1DaZ3N;
	static ScoreType ESI_ETD_hi_res_y_scoreZ3N;
	static ScoreType ESI_ETD_hi_res_z_scoreZ3N;
	static ScoreType ESI_ETD_hi_res_zPlus1DaZ3N;

	static ScoreType ESI_ETD_hi_res_c_scoreZ3C;
	static ScoreType ESI_ETD_hi_res_cMinus1DaZ3C;
	static ScoreType ESI_ETD_hi_res_y_scoreZ3C;
	static ScoreType ESI_ETD_hi_res_z_scoreZ3C;
	static ScoreType ESI_ETD_hi_res_zPlus1DaZ3C;

	static ScoreType ESI_ETD_hi_res_c_scoreZ3;
	static ScoreType ESI_ETD_hi_res_cMinus1DaZ3;
	static ScoreType ESI_ETD_hi_res_y_scoreZ3;
	static ScoreType ESI_ETD_hi_res_z_scoreZ3;
	static ScoreType ESI_ETD_hi_res_zPlus1DaZ3;

	static ScoreType ESI_ETD_hi_res_c_scoreZ4NC;
	static ScoreType ESI_ETD_hi_res_y_scoreZ4NC;
	static ScoreType ESI_ETD_hi_res_z_scoreZ4NC;
	static ScoreType ESI_ETD_hi_res_zPlus1DaZ4NC;

	static ScoreType ESI_ETD_hi_res_c_scoreZ4N;
	static ScoreType ESI_ETD_hi_res_y_scoreZ4N;
	static ScoreType ESI_ETD_hi_res_z_scoreZ4N;
	static ScoreType ESI_ETD_hi_res_zPlus1DaZ4N;

	static ScoreType ESI_ETD_hi_res_c_scoreZ4C;
	static ScoreType ESI_ETD_hi_res_y_scoreZ4C;
	static ScoreType ESI_ETD_hi_res_z_scoreZ4C;
	static ScoreType ESI_ETD_hi_res_zPlus1DaZ4C;

	static ScoreType ESI_ETD_hi_res_c_scoreZ4;
	static ScoreType ESI_ETD_hi_res_y_scoreZ4;
	static ScoreType ESI_ETD_hi_res_z_scoreZ4;
	static ScoreType ESI_ETD_hi_res_zPlus1DaZ4;

	static ScoreType ESI_ETD_hi_res_c_scoreZ5NC;
	static ScoreType ESI_ETD_hi_res_y_scoreZ5NC;
	static ScoreType ESI_ETD_hi_res_z_scoreZ5NC;
	static ScoreType ESI_ETD_hi_res_zPlus1DaZ5NC;

	static ScoreType ESI_ETD_hi_res_c_scoreZ5N;
	static ScoreType ESI_ETD_hi_res_y_scoreZ5N;
	static ScoreType ESI_ETD_hi_res_z_scoreZ5N;
	static ScoreType ESI_ETD_hi_res_zPlus1DaZ5N;

	static ScoreType ESI_ETD_hi_res_c_scoreZ5C;
	static ScoreType ESI_ETD_hi_res_y_scoreZ5C;
	static ScoreType ESI_ETD_hi_res_z_scoreZ5C;
	static ScoreType ESI_ETD_hi_res_zPlus1DaZ5C;

	static ScoreType ESI_ETD_hi_res_c_scoreZ5;
	static ScoreType ESI_ETD_hi_res_y_scoreZ5;
	static ScoreType ESI_ETD_hi_res_z_scoreZ5;
	static ScoreType ESI_ETD_hi_res_zPlus1DaZ5;

	static ScoreType ESI_Q_CID_a_scoreZ1NC;
	static ScoreType ESI_Q_CID_a_nh3_scoreZ1NC;
	static ScoreType ESI_Q_CID_a_h2o_scoreZ1NC;
	static ScoreType ESI_Q_CID_a_h3po4_scoreZ1NC;
	static ScoreType ESI_Q_CID_b_scoreZ1NC;
	static ScoreType ESI_Q_CID_b_nh3_scoreZ1NC;
	static ScoreType ESI_Q_CID_b_h2o_scoreZ1NC;
	static ScoreType ESI_Q_CID_b_soch4_scoreZ1NC;
	static ScoreType ESI_Q_CID_b_h3po4_scoreZ1NC;

	static ScoreType ESI_Q_CID_a_scoreZ1N;
	static ScoreType ESI_Q_CID_a_nh3_scoreZ1N;
	static ScoreType ESI_Q_CID_a_h2o_scoreZ1N;
	static ScoreType ESI_Q_CID_a_h3po4_scoreZ1N;
	static ScoreType ESI_Q_CID_b_scoreZ1N;
	static ScoreType ESI_Q_CID_b_nh3_scoreZ1N;
	static ScoreType ESI_Q_CID_b_h2o_scoreZ1N;
	static ScoreType ESI_Q_CID_b_soch4_scoreZ1N;
	static ScoreType ESI_Q_CID_b_h3po4_scoreZ1N;

	static ScoreType ESI_Q_CID_a_scoreZ1C;
	static ScoreType ESI_Q_CID_a_nh3_scoreZ1C;
	static ScoreType ESI_Q_CID_a_h2o_scoreZ1C;
	static ScoreType ESI_Q_CID_a_h3po4_scoreZ1C;
	static ScoreType ESI_Q_CID_b_scoreZ1C;
	static ScoreType ESI_Q_CID_b_nh3_scoreZ1C;
	static ScoreType ESI_Q_CID_b_h2o_scoreZ1C;
	static ScoreType ESI_Q_CID_b_soch4_scoreZ1C;
	static ScoreType ESI_Q_CID_b_h3po4_scoreZ1C;

	static ScoreType ESI_Q_CID_a_scoreZ1;
	static ScoreType ESI_Q_CID_a_nh3_scoreZ1;
	static ScoreType ESI_Q_CID_a_h2o_scoreZ1;
	static ScoreType ESI_Q_CID_a_h3po4_scoreZ1;
	static ScoreType ESI_Q_CID_b_scoreZ1;
	static ScoreType ESI_Q_CID_b_nh3_scoreZ1;
	static ScoreType ESI_Q_CID_b_h2o_scoreZ1;
	static ScoreType ESI_Q_CID_b_soch4_scoreZ1;
	static ScoreType ESI_Q_CID_b_h3po4_scoreZ1;

	static ScoreType ESI_Q_CID_a_scoreZ2NC;
	static ScoreType ESI_Q_CID_a_nh3_scoreZ2NC;
	static ScoreType ESI_Q_CID_a_h2o_scoreZ2NC;
	static ScoreType ESI_Q_CID_a_h3po4_scoreZ2NC;
	static ScoreType ESI_Q_CID_b_scoreZ2NC;
	static ScoreType ESI_Q_CID_b_nh3_scoreZ2NC;
	static ScoreType ESI_Q_CID_b_h2o_scoreZ2NC;
	static ScoreType ESI_Q_CID_b_soch4_scoreZ2NC;
	static ScoreType ESI_Q_CID_b_h3po4_scoreZ2NC;

	static ScoreType ESI_Q_CID_a_scoreZ2N;
	static ScoreType ESI_Q_CID_a_nh3_scoreZ2N;
	static ScoreType ESI_Q_CID_a_h2o_scoreZ2N;
	static ScoreType ESI_Q_CID_a_h3po4_scoreZ2N;
	static ScoreType ESI_Q_CID_b_scoreZ2N;
	static ScoreType ESI_Q_CID_b_nh3_scoreZ2N;
	static ScoreType ESI_Q_CID_b_h2o_scoreZ2N;
	static ScoreType ESI_Q_CID_b_soch4_scoreZ2N;
	static ScoreType ESI_Q_CID_b_h3po4_scoreZ2N;

	static ScoreType ESI_Q_CID_a_scoreZ2C;
	static ScoreType ESI_Q_CID_a_nh3_scoreZ2C;
	static ScoreType ESI_Q_CID_a_h2o_scoreZ2C;
	static ScoreType ESI_Q_CID_a_h3po4_scoreZ2C;
	static ScoreType ESI_Q_CID_b_scoreZ2C;
	static ScoreType ESI_Q_CID_b_nh3_scoreZ2C;
	static ScoreType ESI_Q_CID_b_h2o_scoreZ2C;
	static ScoreType ESI_Q_CID_b_soch4_scoreZ2C;
	static ScoreType ESI_Q_CID_b_h3po4_scoreZ2C;

	static ScoreType ESI_Q_CID_a_scoreZ2;
	static ScoreType ESI_Q_CID_a_nh3_scoreZ2;
	static ScoreType ESI_Q_CID_a_h2o_scoreZ2;
	static ScoreType ESI_Q_CID_a_h3po4_scoreZ2;
	static ScoreType ESI_Q_CID_b_scoreZ2;
	static ScoreType ESI_Q_CID_b_nh3_scoreZ2;
	static ScoreType ESI_Q_CID_b_h2o_scoreZ2;
	static ScoreType ESI_Q_CID_b_soch4_scoreZ2;
	static ScoreType ESI_Q_CID_b_h3po4_scoreZ2;

	static ScoreType ESI_Q_CID_a_scoreZ3NC;
	static ScoreType ESI_Q_CID_a_nh3_scoreZ3NC;
	static ScoreType ESI_Q_CID_a_h2o_scoreZ3NC;
	static ScoreType ESI_Q_CID_a_h3po4_scoreZ3NC;
	static ScoreType ESI_Q_CID_b_scoreZ3NC;
	static ScoreType ESI_Q_CID_b_nh3_scoreZ3NC;
	static ScoreType ESI_Q_CID_b_h2o_scoreZ3NC;
	static ScoreType ESI_Q_CID_b_soch4_scoreZ3NC;
	static ScoreType ESI_Q_CID_b_h3po4_scoreZ3NC;

	static ScoreType ESI_Q_CID_a_scoreZ3N;
	static ScoreType ESI_Q_CID_a_nh3_scoreZ3N;
	static ScoreType ESI_Q_CID_a_h2o_scoreZ3N;
	static ScoreType ESI_Q_CID_a_h3po4_scoreZ3N;
	static ScoreType ESI_Q_CID_b_scoreZ3N;
	static ScoreType ESI_Q_CID_b_nh3_scoreZ3N;
	static ScoreType ESI_Q_CID_b_h2o_scoreZ3N;
	static ScoreType ESI_Q_CID_b_soch4_scoreZ3N;
	static ScoreType ESI_Q_CID_b_h3po4_scoreZ3N;

	static ScoreType ESI_Q_CID_a_scoreZ3C;
	static ScoreType ESI_Q_CID_a_nh3_scoreZ3C;
	static ScoreType ESI_Q_CID_a_h2o_scoreZ3C;
	static ScoreType ESI_Q_CID_a_h3po4_scoreZ3C;
	static ScoreType ESI_Q_CID_b_scoreZ3C;
	static ScoreType ESI_Q_CID_b_nh3_scoreZ3C;
	static ScoreType ESI_Q_CID_b_h2o_scoreZ3C;
	static ScoreType ESI_Q_CID_b_soch4_scoreZ3C;
	static ScoreType ESI_Q_CID_b_h3po4_scoreZ3C;

	static ScoreType ESI_Q_CID_a_scoreZ3;
	static ScoreType ESI_Q_CID_a_nh3_scoreZ3;
	static ScoreType ESI_Q_CID_a_h2o_scoreZ3;
	static ScoreType ESI_Q_CID_a_h3po4_scoreZ3;
	static ScoreType ESI_Q_CID_b_scoreZ3;
	static ScoreType ESI_Q_CID_b_nh3_scoreZ3;
	static ScoreType ESI_Q_CID_b_h2o_scoreZ3;
	static ScoreType ESI_Q_CID_b_soch4_scoreZ3;
	static ScoreType ESI_Q_CID_b_h3po4_scoreZ3;

	static ScoreType ESI_Q_CID_a_scoreZ4NC;
	static ScoreType ESI_Q_CID_a_nh3_scoreZ4NC;
	static ScoreType ESI_Q_CID_a_h2o_scoreZ4NC;
	static ScoreType ESI_Q_CID_a_h3po4_scoreZ4NC;
	static ScoreType ESI_Q_CID_b_scoreZ4NC;
	static ScoreType ESI_Q_CID_b_nh3_scoreZ4NC;
	static ScoreType ESI_Q_CID_b_h2o_scoreZ4NC;
	static ScoreType ESI_Q_CID_b_soch4_scoreZ4NC;
	static ScoreType ESI_Q_CID_b_h3po4_scoreZ4NC;

	static ScoreType ESI_Q_CID_a_scoreZ4N;
	static ScoreType ESI_Q_CID_a_nh3_scoreZ4N;
	static ScoreType ESI_Q_CID_a_h2o_scoreZ4N;
	static ScoreType ESI_Q_CID_a_h3po4_scoreZ4N;
	static ScoreType ESI_Q_CID_b_scoreZ4N;
	static ScoreType ESI_Q_CID_b_nh3_scoreZ4N;
	static ScoreType ESI_Q_CID_b_h2o_scoreZ4N;
	static ScoreType ESI_Q_CID_b_soch4_scoreZ4N;
	static ScoreType ESI_Q_CID_b_h3po4_scoreZ4N;

	static ScoreType ESI_Q_CID_a_scoreZ4C;
	static ScoreType ESI_Q_CID_a_nh3_scoreZ4C;
	static ScoreType ESI_Q_CID_a_h2o_scoreZ4C;
	static ScoreType ESI_Q_CID_a_h3po4_scoreZ4C;
	static ScoreType ESI_Q_CID_b_scoreZ4C;
	static ScoreType ESI_Q_CID_b_nh3_scoreZ4C;
	static ScoreType ESI_Q_CID_b_h2o_scoreZ4C;
	static ScoreType ESI_Q_CID_b_soch4_scoreZ4C;
	static ScoreType ESI_Q_CID_b_h3po4_scoreZ4C;

	static ScoreType ESI_Q_CID_a_scoreZ4;
	static ScoreType ESI_Q_CID_a_nh3_scoreZ4;
	static ScoreType ESI_Q_CID_a_h2o_scoreZ4;
	static ScoreType ESI_Q_CID_a_h3po4_scoreZ4;
	static ScoreType ESI_Q_CID_b_scoreZ4;
	static ScoreType ESI_Q_CID_b_nh3_scoreZ4;
	static ScoreType ESI_Q_CID_b_h2o_scoreZ4;
	static ScoreType ESI_Q_CID_b_soch4_scoreZ4;
	static ScoreType ESI_Q_CID_b_h3po4_scoreZ4;

	static ScoreType ESI_Q_CID_a_scoreZ5NC;
	static ScoreType ESI_Q_CID_a_nh3_scoreZ5NC;
	static ScoreType ESI_Q_CID_a_h2o_scoreZ5NC;
	static ScoreType ESI_Q_CID_a_h3po4_scoreZ5NC;
	static ScoreType ESI_Q_CID_b_scoreZ5NC;
	static ScoreType ESI_Q_CID_b_nh3_scoreZ5NC;
	static ScoreType ESI_Q_CID_b_h2o_scoreZ5NC;
	static ScoreType ESI_Q_CID_b_soch4_scoreZ5NC;
	static ScoreType ESI_Q_CID_b_h3po4_scoreZ5NC;

	static ScoreType ESI_Q_CID_a_scoreZ5N;
	static ScoreType ESI_Q_CID_a_nh3_scoreZ5N;
	static ScoreType ESI_Q_CID_a_h2o_scoreZ5N;
	static ScoreType ESI_Q_CID_a_h3po4_scoreZ5N;
	static ScoreType ESI_Q_CID_b_scoreZ5N;
	static ScoreType ESI_Q_CID_b_nh3_scoreZ5N;
	static ScoreType ESI_Q_CID_b_h2o_scoreZ5N;
	static ScoreType ESI_Q_CID_b_soch4_scoreZ5N;
	static ScoreType ESI_Q_CID_b_h3po4_scoreZ5N;

	static ScoreType ESI_Q_CID_a_scoreZ5C;
	static ScoreType ESI_Q_CID_a_nh3_scoreZ5C;
	static ScoreType ESI_Q_CID_a_h2o_scoreZ5C;
	static ScoreType ESI_Q_CID_a_h3po4_scoreZ5C;
	static ScoreType ESI_Q_CID_b_scoreZ5C;
	static ScoreType ESI_Q_CID_b_nh3_scoreZ5C;
	static ScoreType ESI_Q_CID_b_h2o_scoreZ5C;
	static ScoreType ESI_Q_CID_b_soch4_scoreZ5C;
	static ScoreType ESI_Q_CID_b_h3po4_scoreZ5C;

	static ScoreType ESI_Q_CID_a_scoreZ5;
	static ScoreType ESI_Q_CID_a_nh3_scoreZ5;
	static ScoreType ESI_Q_CID_a_h2o_scoreZ5;
	static ScoreType ESI_Q_CID_a_h3po4_scoreZ5;
	static ScoreType ESI_Q_CID_b_scoreZ5;
	static ScoreType ESI_Q_CID_b_nh3_scoreZ5;
	static ScoreType ESI_Q_CID_b_h2o_scoreZ5;
	static ScoreType ESI_Q_CID_b_soch4_scoreZ5;
	static ScoreType ESI_Q_CID_b_h3po4_scoreZ5;

	static ScoreType ESI_Q_CID_y_scoreZ1NC;
	static ScoreType ESI_Q_CID_y_nh3_scoreZ1NC;
	static ScoreType ESI_Q_CID_y_h2o_scoreZ1NC;
	static ScoreType ESI_Q_CID_y_soch4_scoreZ1NC;
	static ScoreType ESI_Q_CID_y_h3po4_scoreZ1NC;

	static ScoreType ESI_Q_CID_y_scoreZ1N;
	static ScoreType ESI_Q_CID_y_nh3_scoreZ1N;
	static ScoreType ESI_Q_CID_y_h2o_scoreZ1N;
	static ScoreType ESI_Q_CID_y_soch4_scoreZ1N;
	static ScoreType ESI_Q_CID_y_h3po4_scoreZ1N;

	static ScoreType ESI_Q_CID_y_scoreZ1C;
	static ScoreType ESI_Q_CID_y_nh3_scoreZ1C;
	static ScoreType ESI_Q_CID_y_h2o_scoreZ1C;
	static ScoreType ESI_Q_CID_y_soch4_scoreZ1C;
	static ScoreType ESI_Q_CID_y_h3po4_scoreZ1C;

	static ScoreType ESI_Q_CID_y_scoreZ1;
	static ScoreType ESI_Q_CID_y_nh3_scoreZ1;
	static ScoreType ESI_Q_CID_y_h2o_scoreZ1;
	static ScoreType ESI_Q_CID_y_soch4_scoreZ1;
	static ScoreType ESI_Q_CID_y_h3po4_scoreZ1;

	static ScoreType ESI_Q_CID_y_scoreZ2NC;
	static ScoreType ESI_Q_CID_y_nh3_scoreZ2NC;
	static ScoreType ESI_Q_CID_y_h2o_scoreZ2NC;
	static ScoreType ESI_Q_CID_y_soch4_scoreZ2NC;
	static ScoreType ESI_Q_CID_y_h3po4_scoreZ2NC;

	static ScoreType ESI_Q_CID_y_scoreZ2N;
	static ScoreType ESI_Q_CID_y_nh3_scoreZ2N;
	static ScoreType ESI_Q_CID_y_h2o_scoreZ2N;
	static ScoreType ESI_Q_CID_y_soch4_scoreZ2N;
	static ScoreType ESI_Q_CID_y_h3po4_scoreZ2N;

	static ScoreType ESI_Q_CID_y_scoreZ2C;
	static ScoreType ESI_Q_CID_y_nh3_scoreZ2C;
	static ScoreType ESI_Q_CID_y_h2o_scoreZ2C;
	static ScoreType ESI_Q_CID_y_soch4_scoreZ2C;
	static ScoreType ESI_Q_CID_y_h3po4_scoreZ2C;

	static ScoreType ESI_Q_CID_y_scoreZ2;
	static ScoreType ESI_Q_CID_y_nh3_scoreZ2;
	static ScoreType ESI_Q_CID_y_h2o_scoreZ2;
	static ScoreType ESI_Q_CID_y_soch4_scoreZ2;
	static ScoreType ESI_Q_CID_y_h3po4_scoreZ2;

	static ScoreType ESI_Q_CID_y_scoreZ3NC;
	static ScoreType ESI_Q_CID_y_nh3_scoreZ3NC;
	static ScoreType ESI_Q_CID_y_h2o_scoreZ3NC;
	static ScoreType ESI_Q_CID_y_soch4_scoreZ3NC;
	static ScoreType ESI_Q_CID_y_h3po4_scoreZ3NC;

	static ScoreType ESI_Q_CID_y_scoreZ3N;
	static ScoreType ESI_Q_CID_y_nh3_scoreZ3N;
	static ScoreType ESI_Q_CID_y_h2o_scoreZ3N;
	static ScoreType ESI_Q_CID_y_soch4_scoreZ3N;
	static ScoreType ESI_Q_CID_y_h3po4_scoreZ3N;

	static ScoreType ESI_Q_CID_y_scoreZ3C;
	static ScoreType ESI_Q_CID_y_nh3_scoreZ3C;
	static ScoreType ESI_Q_CID_y_h2o_scoreZ3C;
	static ScoreType ESI_Q_CID_y_soch4_scoreZ3C;
	static ScoreType ESI_Q_CID_y_h3po4_scoreZ3C;

	static ScoreType ESI_Q_CID_y_scoreZ3;
	static ScoreType ESI_Q_CID_y_nh3_scoreZ3;
	static ScoreType ESI_Q_CID_y_h2o_scoreZ3;
	static ScoreType ESI_Q_CID_y_soch4_scoreZ3;
	static ScoreType ESI_Q_CID_y_h3po4_scoreZ3;

	static ScoreType ESI_Q_CID_y_scoreZ4NC;
	static ScoreType ESI_Q_CID_y_nh3_scoreZ4NC;
	static ScoreType ESI_Q_CID_y_h2o_scoreZ4NC;
	static ScoreType ESI_Q_CID_y_soch4_scoreZ4NC;
	static ScoreType ESI_Q_CID_y_h3po4_scoreZ4NC;

	static ScoreType ESI_Q_CID_y_scoreZ4N;
	static ScoreType ESI_Q_CID_y_nh3_scoreZ4N;
	static ScoreType ESI_Q_CID_y_h2o_scoreZ4N;
	static ScoreType ESI_Q_CID_y_soch4_scoreZ4N;
	static ScoreType ESI_Q_CID_y_h3po4_scoreZ4N;

	static ScoreType ESI_Q_CID_y_scoreZ4C;
	static ScoreType ESI_Q_CID_y_nh3_scoreZ4C;
	static ScoreType ESI_Q_CID_y_h2o_scoreZ4C;
	static ScoreType ESI_Q_CID_y_soch4_scoreZ4C;
	static ScoreType ESI_Q_CID_y_h3po4_scoreZ4C;

	static ScoreType ESI_Q_CID_y_scoreZ4;
	static ScoreType ESI_Q_CID_y_nh3_scoreZ4;
	static ScoreType ESI_Q_CID_y_h2o_scoreZ4;
	static ScoreType ESI_Q_CID_y_soch4_scoreZ4;
	static ScoreType ESI_Q_CID_y_h3po4_scoreZ4;

	static ScoreType ESI_Q_CID_y_scoreZ5NC;
	static ScoreType ESI_Q_CID_y_nh3_scoreZ5NC;
	static ScoreType ESI_Q_CID_y_h2o_scoreZ5NC;
	static ScoreType ESI_Q_CID_y_soch4_scoreZ5NC;
	static ScoreType ESI_Q_CID_y_h3po4_scoreZ5NC;

	static ScoreType ESI_Q_CID_y_scoreZ5N;
	static ScoreType ESI_Q_CID_y_nh3_scoreZ5N;
	static ScoreType ESI_Q_CID_y_h2o_scoreZ5N;
	static ScoreType ESI_Q_CID_y_soch4_scoreZ5N;
	static ScoreType ESI_Q_CID_y_h3po4_scoreZ5N;

	static ScoreType ESI_Q_CID_y_scoreZ5C;
	static ScoreType ESI_Q_CID_y_nh3_scoreZ5C;
	static ScoreType ESI_Q_CID_y_h2o_scoreZ5C;
	static ScoreType ESI_Q_CID_y_soch4_scoreZ5C;
	static ScoreType ESI_Q_CID_y_h3po4_scoreZ5C;

	static ScoreType ESI_Q_CID_y_scoreZ5;
	static ScoreType ESI_Q_CID_y_nh3_scoreZ5;
	static ScoreType ESI_Q_CID_y_h2o_scoreZ5;
	static ScoreType ESI_Q_CID_y_soch4_scoreZ5;
	static ScoreType ESI_Q_CID_y_h3po4_scoreZ5;

	static ScoreType ESI_TRAP_CID_low_res_a_scoreZ1NC;
	static ScoreType ESI_TRAP_CID_low_res_b_scoreZ1NC;
	static ScoreType ESI_TRAP_CID_low_res_b_nh3_scoreZ1NC;
	static ScoreType ESI_TRAP_CID_low_res_b_h2o_scoreZ1NC;
	static ScoreType ESI_TRAP_CID_low_res_b_soch4_scoreZ1NC;
	static ScoreType ESI_TRAP_CID_low_res_b_h3po4_scoreZ1NC;

	static ScoreType ESI_TRAP_CID_low_res_a_scoreZ1N;
	static ScoreType ESI_TRAP_CID_low_res_b_scoreZ1N;
	static ScoreType ESI_TRAP_CID_low_res_b_nh3_scoreZ1N;
	static ScoreType ESI_TRAP_CID_low_res_b_h2o_scoreZ1N;
	static ScoreType ESI_TRAP_CID_low_res_b_soch4_scoreZ1N;
	static ScoreType ESI_TRAP_CID_low_res_b_h3po4_scoreZ1N;

	static ScoreType ESI_TRAP_CID_low_res_a_scoreZ1C;
	static ScoreType ESI_TRAP_CID_low_res_b_scoreZ1C;
	static ScoreType ESI_TRAP_CID_low_res_b_nh3_scoreZ1C;
	static ScoreType ESI_TRAP_CID_low_res_b_h2o_scoreZ1C;
	static ScoreType ESI_TRAP_CID_low_res_b_soch4_scoreZ1C;
	static ScoreType ESI_TRAP_CID_low_res_b_h3po4_scoreZ1C;

	static ScoreType ESI_TRAP_CID_low_res_a_scoreZ1;
	static ScoreType ESI_TRAP_CID_low_res_b_scoreZ1;
	static ScoreType ESI_TRAP_CID_low_res_b_nh3_scoreZ1;
	static ScoreType ESI_TRAP_CID_low_res_b_h2o_scoreZ1;
	static ScoreType ESI_TRAP_CID_low_res_b_soch4_scoreZ1;
	static ScoreType ESI_TRAP_CID_low_res_b_h3po4_scoreZ1;

	static ScoreType ESI_TRAP_CID_low_res_a_scoreZ2NC;
	static ScoreType ESI_TRAP_CID_low_res_b_scoreZ2NC;
	static ScoreType ESI_TRAP_CID_low_res_bp2_scoreZ2NC;
	static ScoreType ESI_TRAP_CID_low_res_b_nh3_scoreZ2NC;
	static ScoreType ESI_TRAP_CID_low_res_bp2_nh3_scoreZ2NC;
	static ScoreType ESI_TRAP_CID_low_res_b_h2o_scoreZ2NC;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h2o_scoreZ2NC;
	static ScoreType ESI_TRAP_CID_low_res_b_soch4_scoreZ2NC;
	static ScoreType ESI_TRAP_CID_low_res_bp2_soch4_scoreZ2NC;
	static ScoreType ESI_TRAP_CID_low_res_b_h3po4_scoreZ2NC;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ2NC;

	static ScoreType ESI_TRAP_CID_low_res_a_scoreZ2N;
	static ScoreType ESI_TRAP_CID_low_res_b_scoreZ2N;
	static ScoreType ESI_TRAP_CID_low_res_bp2_scoreZ2N;
	static ScoreType ESI_TRAP_CID_low_res_b_nh3_scoreZ2N;
	static ScoreType ESI_TRAP_CID_low_res_bp2_nh3_scoreZ2N;
	static ScoreType ESI_TRAP_CID_low_res_b_h2o_scoreZ2N;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h2o_scoreZ2N;
	static ScoreType ESI_TRAP_CID_low_res_b_soch4_scoreZ2N;
	static ScoreType ESI_TRAP_CID_low_res_bp2_soch4_scoreZ2N;
	static ScoreType ESI_TRAP_CID_low_res_b_h3po4_scoreZ2N;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ2N;

	static ScoreType ESI_TRAP_CID_low_res_a_scoreZ2C;
	static ScoreType ESI_TRAP_CID_low_res_b_scoreZ2C;
	static ScoreType ESI_TRAP_CID_low_res_bp2_scoreZ2C;
	static ScoreType ESI_TRAP_CID_low_res_b_nh3_scoreZ2C;
	static ScoreType ESI_TRAP_CID_low_res_bp2_nh3_scoreZ2C;
	static ScoreType ESI_TRAP_CID_low_res_b_h2o_scoreZ2C;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h2o_scoreZ2C;
	static ScoreType ESI_TRAP_CID_low_res_b_soch4_scoreZ2C;
	static ScoreType ESI_TRAP_CID_low_res_bp2_soch4_scoreZ2C;
	static ScoreType ESI_TRAP_CID_low_res_b_h3po4_scoreZ2C;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ2C;

	static ScoreType ESI_TRAP_CID_low_res_a_scoreZ2;
	static ScoreType ESI_TRAP_CID_low_res_b_scoreZ2;
	static ScoreType ESI_TRAP_CID_low_res_bp2_scoreZ2;
	static ScoreType ESI_TRAP_CID_low_res_b_nh3_scoreZ2;
	static ScoreType ESI_TRAP_CID_low_res_bp2_nh3_scoreZ2;
	static ScoreType ESI_TRAP_CID_low_res_b_h2o_scoreZ2;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h2o_scoreZ2;
	static ScoreType ESI_TRAP_CID_low_res_b_soch4_scoreZ2;
	static ScoreType ESI_TRAP_CID_low_res_bp2_soch4_scoreZ2;
	static ScoreType ESI_TRAP_CID_low_res_b_h3po4_scoreZ2;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ2;

	static ScoreType ESI_TRAP_CID_low_res_a_scoreZ3NC;
	static ScoreType ESI_TRAP_CID_low_res_b_scoreZ3NC;
	static ScoreType ESI_TRAP_CID_low_res_bp2_scoreZ3NC;
	static ScoreType ESI_TRAP_CID_low_res_bp3_scoreZ3NC;
	static ScoreType ESI_TRAP_CID_low_res_b_nh3_scoreZ3NC;
	static ScoreType ESI_TRAP_CID_low_res_bp2_nh3_scoreZ3NC;
	static ScoreType ESI_TRAP_CID_low_res_b_h2o_scoreZ3NC;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h2o_scoreZ3NC;
	static ScoreType ESI_TRAP_CID_low_res_b_soch4_scoreZ3NC;
	static ScoreType ESI_TRAP_CID_low_res_bp2_soch4_scoreZ3NC;
	static ScoreType ESI_TRAP_CID_low_res_b_h3po4_scoreZ3NC;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ3NC;

	static ScoreType ESI_TRAP_CID_low_res_a_scoreZ3N;
	static ScoreType ESI_TRAP_CID_low_res_b_scoreZ3N;
	static ScoreType ESI_TRAP_CID_low_res_bp2_scoreZ3N;
	static ScoreType ESI_TRAP_CID_low_res_bp3_scoreZ3N;
	static ScoreType ESI_TRAP_CID_low_res_b_nh3_scoreZ3N;
	static ScoreType ESI_TRAP_CID_low_res_bp2_nh3_scoreZ3N;
	static ScoreType ESI_TRAP_CID_low_res_b_h2o_scoreZ3N;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h2o_scoreZ3N;
	static ScoreType ESI_TRAP_CID_low_res_b_soch4_scoreZ3N;
	static ScoreType ESI_TRAP_CID_low_res_bp2_soch4_scoreZ3N;
	static ScoreType ESI_TRAP_CID_low_res_b_h3po4_scoreZ3N;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ3N;

	static ScoreType ESI_TRAP_CID_low_res_a_scoreZ3C;
	static ScoreType ESI_TRAP_CID_low_res_b_scoreZ3C;
	static ScoreType ESI_TRAP_CID_low_res_bp2_scoreZ3C;
	static ScoreType ESI_TRAP_CID_low_res_bp3_scoreZ3C;
	static ScoreType ESI_TRAP_CID_low_res_b_nh3_scoreZ3C;
	static ScoreType ESI_TRAP_CID_low_res_bp2_nh3_scoreZ3C;
	static ScoreType ESI_TRAP_CID_low_res_b_h2o_scoreZ3C;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h2o_scoreZ3C;
	static ScoreType ESI_TRAP_CID_low_res_b_soch4_scoreZ3C;
	static ScoreType ESI_TRAP_CID_low_res_bp2_soch4_scoreZ3C;
	static ScoreType ESI_TRAP_CID_low_res_b_h3po4_scoreZ3C;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ3C;

	static ScoreType ESI_TRAP_CID_low_res_a_scoreZ3;
	static ScoreType ESI_TRAP_CID_low_res_b_scoreZ3;
	static ScoreType ESI_TRAP_CID_low_res_bp2_scoreZ3;
	static ScoreType ESI_TRAP_CID_low_res_bp3_scoreZ3;
	static ScoreType ESI_TRAP_CID_low_res_b_nh3_scoreZ3;
	static ScoreType ESI_TRAP_CID_low_res_bp2_nh3_scoreZ3;
	static ScoreType ESI_TRAP_CID_low_res_b_h2o_scoreZ3;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h2o_scoreZ3;
	static ScoreType ESI_TRAP_CID_low_res_b_soch4_scoreZ3;
	static ScoreType ESI_TRAP_CID_low_res_bp2_soch4_scoreZ3;
	static ScoreType ESI_TRAP_CID_low_res_b_h3po4_scoreZ3;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ3;

	static ScoreType ESI_TRAP_CID_low_res_a_scoreZ4NC;
	static ScoreType ESI_TRAP_CID_low_res_b_scoreZ4NC;
	static ScoreType ESI_TRAP_CID_low_res_bp2_scoreZ4NC;
	static ScoreType ESI_TRAP_CID_low_res_bp3_scoreZ4NC;
	static ScoreType ESI_TRAP_CID_low_res_b_nh3_scoreZ4NC;
	static ScoreType ESI_TRAP_CID_low_res_bp2_nh3_scoreZ4NC;
	static ScoreType ESI_TRAP_CID_low_res_b_h2o_scoreZ4NC;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h2o_scoreZ4NC;
	static ScoreType ESI_TRAP_CID_low_res_b_soch4_scoreZ4NC;
	static ScoreType ESI_TRAP_CID_low_res_bp2_soch4_scoreZ4NC;
	static ScoreType ESI_TRAP_CID_low_res_b_h3po4_scoreZ4NC;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ4NC;

	static ScoreType ESI_TRAP_CID_low_res_a_scoreZ4N;
	static ScoreType ESI_TRAP_CID_low_res_b_scoreZ4N;
	static ScoreType ESI_TRAP_CID_low_res_bp2_scoreZ4N;
	static ScoreType ESI_TRAP_CID_low_res_bp3_scoreZ4N;
	static ScoreType ESI_TRAP_CID_low_res_b_nh3_scoreZ4N;
	static ScoreType ESI_TRAP_CID_low_res_bp2_nh3_scoreZ4N;
	static ScoreType ESI_TRAP_CID_low_res_b_h2o_scoreZ4N;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h2o_scoreZ4N;
	static ScoreType ESI_TRAP_CID_low_res_b_soch4_scoreZ4N;
	static ScoreType ESI_TRAP_CID_low_res_bp2_soch4_scoreZ4N;
	static ScoreType ESI_TRAP_CID_low_res_b_h3po4_scoreZ4N;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ4N;

	static ScoreType ESI_TRAP_CID_low_res_a_scoreZ4C;
	static ScoreType ESI_TRAP_CID_low_res_b_scoreZ4C;
	static ScoreType ESI_TRAP_CID_low_res_bp2_scoreZ4C;
	static ScoreType ESI_TRAP_CID_low_res_bp3_scoreZ4C;
	static ScoreType ESI_TRAP_CID_low_res_b_nh3_scoreZ4C;
	static ScoreType ESI_TRAP_CID_low_res_bp2_nh3_scoreZ4C;
	static ScoreType ESI_TRAP_CID_low_res_b_h2o_scoreZ4C;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h2o_scoreZ4C;
	static ScoreType ESI_TRAP_CID_low_res_b_soch4_scoreZ4C;
	static ScoreType ESI_TRAP_CID_low_res_bp2_soch4_scoreZ4C;
	static ScoreType ESI_TRAP_CID_low_res_b_h3po4_scoreZ4C;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ4C;

	static ScoreType ESI_TRAP_CID_low_res_a_scoreZ4;
	static ScoreType ESI_TRAP_CID_low_res_b_scoreZ4;
	static ScoreType ESI_TRAP_CID_low_res_bp2_scoreZ4;
	static ScoreType ESI_TRAP_CID_low_res_bp3_scoreZ4;
	static ScoreType ESI_TRAP_CID_low_res_b_nh3_scoreZ4;
	static ScoreType ESI_TRAP_CID_low_res_bp2_nh3_scoreZ4;
	static ScoreType ESI_TRAP_CID_low_res_b_h2o_scoreZ4;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h2o_scoreZ4;
	static ScoreType ESI_TRAP_CID_low_res_b_soch4_scoreZ4;
	static ScoreType ESI_TRAP_CID_low_res_bp2_soch4_scoreZ4;
	static ScoreType ESI_TRAP_CID_low_res_b_h3po4_scoreZ4;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ4;

	static ScoreType ESI_TRAP_CID_low_res_a_scoreZ5NC;
	static ScoreType ESI_TRAP_CID_low_res_b_scoreZ5NC;
	static ScoreType ESI_TRAP_CID_low_res_bp2_scoreZ5NC;
	static ScoreType ESI_TRAP_CID_low_res_bp3_scoreZ5NC;
	static ScoreType ESI_TRAP_CID_low_res_b_nh3_scoreZ5NC;
	static ScoreType ESI_TRAP_CID_low_res_bp2_nh3_scoreZ5NC;
	static ScoreType ESI_TRAP_CID_low_res_b_h2o_scoreZ5NC;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h2o_scoreZ5NC;
	static ScoreType ESI_TRAP_CID_low_res_b_soch4_scoreZ5NC;
	static ScoreType ESI_TRAP_CID_low_res_bp2_soch4_scoreZ5NC;
	static ScoreType ESI_TRAP_CID_low_res_b_h3po4_scoreZ5NC;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ5NC;

	static ScoreType ESI_TRAP_CID_low_res_a_scoreZ5N;
	static ScoreType ESI_TRAP_CID_low_res_b_scoreZ5N;
	static ScoreType ESI_TRAP_CID_low_res_bp2_scoreZ5N;
	static ScoreType ESI_TRAP_CID_low_res_bp3_scoreZ5N;
	static ScoreType ESI_TRAP_CID_low_res_b_nh3_scoreZ5N;
	static ScoreType ESI_TRAP_CID_low_res_bp2_nh3_scoreZ5N;
	static ScoreType ESI_TRAP_CID_low_res_b_h2o_scoreZ5N;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h2o_scoreZ5N;
	static ScoreType ESI_TRAP_CID_low_res_b_soch4_scoreZ5N;
	static ScoreType ESI_TRAP_CID_low_res_bp2_soch4_scoreZ5N;
	static ScoreType ESI_TRAP_CID_low_res_b_h3po4_scoreZ5N;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ5N;

	static ScoreType ESI_TRAP_CID_low_res_a_scoreZ5C;
	static ScoreType ESI_TRAP_CID_low_res_b_scoreZ5C;
	static ScoreType ESI_TRAP_CID_low_res_bp2_scoreZ5C;
	static ScoreType ESI_TRAP_CID_low_res_bp3_scoreZ5C;
	static ScoreType ESI_TRAP_CID_low_res_b_nh3_scoreZ5C;
	static ScoreType ESI_TRAP_CID_low_res_bp2_nh3_scoreZ5C;
	static ScoreType ESI_TRAP_CID_low_res_b_h2o_scoreZ5C;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h2o_scoreZ5C;
	static ScoreType ESI_TRAP_CID_low_res_b_soch4_scoreZ5C;
	static ScoreType ESI_TRAP_CID_low_res_bp2_soch4_scoreZ5C;
	static ScoreType ESI_TRAP_CID_low_res_b_h3po4_scoreZ5C;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ5C;

	static ScoreType ESI_TRAP_CID_low_res_a_scoreZ5;
	static ScoreType ESI_TRAP_CID_low_res_b_scoreZ5;
	static ScoreType ESI_TRAP_CID_low_res_bp2_scoreZ5;
	static ScoreType ESI_TRAP_CID_low_res_bp3_scoreZ5;
	static ScoreType ESI_TRAP_CID_low_res_b_nh3_scoreZ5;
	static ScoreType ESI_TRAP_CID_low_res_bp2_nh3_scoreZ5;
	static ScoreType ESI_TRAP_CID_low_res_b_h2o_scoreZ5;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h2o_scoreZ5;
	static ScoreType ESI_TRAP_CID_low_res_b_soch4_scoreZ5;
	static ScoreType ESI_TRAP_CID_low_res_bp2_soch4_scoreZ5;
	static ScoreType ESI_TRAP_CID_low_res_b_h3po4_scoreZ5;
	static ScoreType ESI_TRAP_CID_low_res_bp2_h3po4_scoreZ5;

	static ScoreType ESI_TRAP_CID_low_res_y_scoreZ1NC;
	static ScoreType ESI_TRAP_CID_low_res_y_nh3_scoreZ1NC;
	static ScoreType ESI_TRAP_CID_low_res_y_h2o_scoreZ1NC;
	static ScoreType ESI_TRAP_CID_low_res_y_soch4_scoreZ1NC;
	static ScoreType ESI_TRAP_CID_low_res_y_h3po4_scoreZ1NC;

	static ScoreType ESI_TRAP_CID_low_res_y_scoreZ1N;
	static ScoreType ESI_TRAP_CID_low_res_y_nh3_scoreZ1N;
	static ScoreType ESI_TRAP_CID_low_res_y_h2o_scoreZ1N;
	static ScoreType ESI_TRAP_CID_low_res_y_soch4_scoreZ1N;
	static ScoreType ESI_TRAP_CID_low_res_y_h3po4_scoreZ1N;

	static ScoreType ESI_TRAP_CID_low_res_y_scoreZ1C;
	static ScoreType ESI_TRAP_CID_low_res_y_nh3_scoreZ1C;
	static ScoreType ESI_TRAP_CID_low_res_y_h2o_scoreZ1C;
	static ScoreType ESI_TRAP_CID_low_res_y_soch4_scoreZ1C;
	static ScoreType ESI_TRAP_CID_low_res_y_h3po4_scoreZ1C;

	static ScoreType ESI_TRAP_CID_low_res_y_scoreZ1;
	static ScoreType ESI_TRAP_CID_low_res_y_nh3_scoreZ1;
	static ScoreType ESI_TRAP_CID_low_res_y_h2o_scoreZ1;
	static ScoreType ESI_TRAP_CID_low_res_y_soch4_scoreZ1;
	static ScoreType ESI_TRAP_CID_low_res_y_h3po4_scoreZ1;

	static ScoreType ESI_TRAP_CID_low_res_y_scoreZ2NC;
	static ScoreType ESI_TRAP_CID_low_res_yp2_scoreZ2NC;
	static ScoreType ESI_TRAP_CID_low_res_y_nh3_scoreZ2NC;
	static ScoreType ESI_TRAP_CID_low_res_yp2_nh3_scoreZ2NC;
	static ScoreType ESI_TRAP_CID_low_res_y_h2o_scoreZ2NC;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h2o_scoreZ2NC;
	static ScoreType ESI_TRAP_CID_low_res_y_soch4_scoreZ2NC;
	static ScoreType ESI_TRAP_CID_low_res_yp2_soch4_scoreZ2NC;
	static ScoreType ESI_TRAP_CID_low_res_y_h3po4_scoreZ2NC;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ2NC;

	static ScoreType ESI_TRAP_CID_low_res_y_scoreZ2N;
	static ScoreType ESI_TRAP_CID_low_res_yp2_scoreZ2N;
	static ScoreType ESI_TRAP_CID_low_res_y_nh3_scoreZ2N;
	static ScoreType ESI_TRAP_CID_low_res_yp2_nh3_scoreZ2N;
	static ScoreType ESI_TRAP_CID_low_res_y_h2o_scoreZ2N;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h2o_scoreZ2N;
	static ScoreType ESI_TRAP_CID_low_res_y_soch4_scoreZ2N;
	static ScoreType ESI_TRAP_CID_low_res_yp2_soch4_scoreZ2N;
	static ScoreType ESI_TRAP_CID_low_res_y_h3po4_scoreZ2N;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ2N;

	static ScoreType ESI_TRAP_CID_low_res_y_scoreZ2C;
	static ScoreType ESI_TRAP_CID_low_res_yp2_scoreZ2C;
	static ScoreType ESI_TRAP_CID_low_res_y_nh3_scoreZ2C;
	static ScoreType ESI_TRAP_CID_low_res_yp2_nh3_scoreZ2C;
	static ScoreType ESI_TRAP_CID_low_res_y_h2o_scoreZ2C;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h2o_scoreZ2C;
	static ScoreType ESI_TRAP_CID_low_res_y_soch4_scoreZ2C;
	static ScoreType ESI_TRAP_CID_low_res_yp2_soch4_scoreZ2C;
	static ScoreType ESI_TRAP_CID_low_res_y_h3po4_scoreZ2C;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ2C;

	static ScoreType ESI_TRAP_CID_low_res_y_scoreZ2;
	static ScoreType ESI_TRAP_CID_low_res_yp2_scoreZ2;
	static ScoreType ESI_TRAP_CID_low_res_y_nh3_scoreZ2;
	static ScoreType ESI_TRAP_CID_low_res_yp2_nh3_scoreZ2;
	static ScoreType ESI_TRAP_CID_low_res_y_h2o_scoreZ2;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h2o_scoreZ2;
	static ScoreType ESI_TRAP_CID_low_res_y_soch4_scoreZ2;
	static ScoreType ESI_TRAP_CID_low_res_yp2_soch4_scoreZ2;
	static ScoreType ESI_TRAP_CID_low_res_y_h3po4_scoreZ2;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ2;

	static ScoreType ESI_TRAP_CID_low_res_y_scoreZ3NC;
	static ScoreType ESI_TRAP_CID_low_res_yp2_scoreZ3NC;
	static ScoreType ESI_TRAP_CID_low_res_yp3_scoreZ3NC;
	static ScoreType ESI_TRAP_CID_low_res_y_nh3_scoreZ3NC;
	static ScoreType ESI_TRAP_CID_low_res_yp2_nh3_scoreZ3NC;
	static ScoreType ESI_TRAP_CID_low_res_y_h2o_scoreZ3NC;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h2o_scoreZ3NC;
	static ScoreType ESI_TRAP_CID_low_res_y_soch4_scoreZ3NC;
	static ScoreType ESI_TRAP_CID_low_res_yp2_soch4_scoreZ3NC;
	static ScoreType ESI_TRAP_CID_low_res_y_h3po4_scoreZ3NC;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ3NC;

	static ScoreType ESI_TRAP_CID_low_res_y_scoreZ3N;
	static ScoreType ESI_TRAP_CID_low_res_yp2_scoreZ3N;
	static ScoreType ESI_TRAP_CID_low_res_yp3_scoreZ3N;
	static ScoreType ESI_TRAP_CID_low_res_y_nh3_scoreZ3N;
	static ScoreType ESI_TRAP_CID_low_res_yp2_nh3_scoreZ3N;
	static ScoreType ESI_TRAP_CID_low_res_y_h2o_scoreZ3N;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h2o_scoreZ3N;
	static ScoreType ESI_TRAP_CID_low_res_y_soch4_scoreZ3N;
	static ScoreType ESI_TRAP_CID_low_res_yp2_soch4_scoreZ3N;
	static ScoreType ESI_TRAP_CID_low_res_y_h3po4_scoreZ3N;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ3N;

	static ScoreType ESI_TRAP_CID_low_res_y_scoreZ3C;
	static ScoreType ESI_TRAP_CID_low_res_yp2_scoreZ3C;
	static ScoreType ESI_TRAP_CID_low_res_yp3_scoreZ3C;
	static ScoreType ESI_TRAP_CID_low_res_y_nh3_scoreZ3C;
	static ScoreType ESI_TRAP_CID_low_res_yp2_nh3_scoreZ3C;
	static ScoreType ESI_TRAP_CID_low_res_y_h2o_scoreZ3C;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h2o_scoreZ3C;
	static ScoreType ESI_TRAP_CID_low_res_y_soch4_scoreZ3C;
	static ScoreType ESI_TRAP_CID_low_res_yp2_soch4_scoreZ3C;
	static ScoreType ESI_TRAP_CID_low_res_y_h3po4_scoreZ3C;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ3C;

	static ScoreType ESI_TRAP_CID_low_res_y_scoreZ3;
	static ScoreType ESI_TRAP_CID_low_res_yp2_scoreZ3;
	static ScoreType ESI_TRAP_CID_low_res_yp3_scoreZ3;
	static ScoreType ESI_TRAP_CID_low_res_y_nh3_scoreZ3;
	static ScoreType ESI_TRAP_CID_low_res_yp2_nh3_scoreZ3;
	static ScoreType ESI_TRAP_CID_low_res_y_h2o_scoreZ3;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h2o_scoreZ3;
	static ScoreType ESI_TRAP_CID_low_res_y_soch4_scoreZ3;
	static ScoreType ESI_TRAP_CID_low_res_yp2_soch4_scoreZ3;
	static ScoreType ESI_TRAP_CID_low_res_y_h3po4_scoreZ3;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ3;

	static ScoreType ESI_TRAP_CID_low_res_y_scoreZ4NC;
	static ScoreType ESI_TRAP_CID_low_res_yp2_scoreZ4NC;
	static ScoreType ESI_TRAP_CID_low_res_yp3_scoreZ4NC;
	static ScoreType ESI_TRAP_CID_low_res_y_nh3_scoreZ4NC;
	static ScoreType ESI_TRAP_CID_low_res_yp2_nh3_scoreZ4NC;
	static ScoreType ESI_TRAP_CID_low_res_y_h2o_scoreZ4NC;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h2o_scoreZ4NC;
	static ScoreType ESI_TRAP_CID_low_res_y_soch4_scoreZ4NC;
	static ScoreType ESI_TRAP_CID_low_res_yp2_soch4_scoreZ4NC;
	static ScoreType ESI_TRAP_CID_low_res_y_h3po4_scoreZ4NC;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ4NC;

	static ScoreType ESI_TRAP_CID_low_res_y_scoreZ4N;
	static ScoreType ESI_TRAP_CID_low_res_yp2_scoreZ4N;
	static ScoreType ESI_TRAP_CID_low_res_yp3_scoreZ4N;
	static ScoreType ESI_TRAP_CID_low_res_y_nh3_scoreZ4N;
	static ScoreType ESI_TRAP_CID_low_res_yp2_nh3_scoreZ4N;
	static ScoreType ESI_TRAP_CID_low_res_y_h2o_scoreZ4N;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h2o_scoreZ4N;
	static ScoreType ESI_TRAP_CID_low_res_y_soch4_scoreZ4N;
	static ScoreType ESI_TRAP_CID_low_res_yp2_soch4_scoreZ4N;
	static ScoreType ESI_TRAP_CID_low_res_y_h3po4_scoreZ4N;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ4N;

	static ScoreType ESI_TRAP_CID_low_res_y_scoreZ4C;
	static ScoreType ESI_TRAP_CID_low_res_yp2_scoreZ4C;
	static ScoreType ESI_TRAP_CID_low_res_yp3_scoreZ4C;
	static ScoreType ESI_TRAP_CID_low_res_y_nh3_scoreZ4C;
	static ScoreType ESI_TRAP_CID_low_res_yp2_nh3_scoreZ4C;
	static ScoreType ESI_TRAP_CID_low_res_y_h2o_scoreZ4C;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h2o_scoreZ4C;
	static ScoreType ESI_TRAP_CID_low_res_y_soch4_scoreZ4C;
	static ScoreType ESI_TRAP_CID_low_res_yp2_soch4_scoreZ4C;
	static ScoreType ESI_TRAP_CID_low_res_y_h3po4_scoreZ4C;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ4C;

	static ScoreType ESI_TRAP_CID_low_res_y_scoreZ4;
	static ScoreType ESI_TRAP_CID_low_res_yp2_scoreZ4;
	static ScoreType ESI_TRAP_CID_low_res_yp3_scoreZ4;
	static ScoreType ESI_TRAP_CID_low_res_y_nh3_scoreZ4;
	static ScoreType ESI_TRAP_CID_low_res_yp2_nh3_scoreZ4;
	static ScoreType ESI_TRAP_CID_low_res_y_h2o_scoreZ4;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h2o_scoreZ4;
	static ScoreType ESI_TRAP_CID_low_res_y_soch4_scoreZ4;
	static ScoreType ESI_TRAP_CID_low_res_yp2_soch4_scoreZ4;
	static ScoreType ESI_TRAP_CID_low_res_y_h3po4_scoreZ4;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ4;

	static ScoreType ESI_TRAP_CID_low_res_y_scoreZ5NC;
	static ScoreType ESI_TRAP_CID_low_res_yp2_scoreZ5NC;
	static ScoreType ESI_TRAP_CID_low_res_yp3_scoreZ5NC;
	static ScoreType ESI_TRAP_CID_low_res_y_nh3_scoreZ5NC;
	static ScoreType ESI_TRAP_CID_low_res_yp2_nh3_scoreZ5NC;
	static ScoreType ESI_TRAP_CID_low_res_y_h2o_scoreZ5NC;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h2o_scoreZ5NC;
	static ScoreType ESI_TRAP_CID_low_res_y_soch4_scoreZ5NC;
	static ScoreType ESI_TRAP_CID_low_res_yp2_soch4_scoreZ5NC;
	static ScoreType ESI_TRAP_CID_low_res_y_h3po4_scoreZ5NC;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ5NC;

	static ScoreType ESI_TRAP_CID_low_res_y_scoreZ5N;
	static ScoreType ESI_TRAP_CID_low_res_yp2_scoreZ5N;
	static ScoreType ESI_TRAP_CID_low_res_yp3_scoreZ5N;
	static ScoreType ESI_TRAP_CID_low_res_y_nh3_scoreZ5N;
	static ScoreType ESI_TRAP_CID_low_res_yp2_nh3_scoreZ5N;
	static ScoreType ESI_TRAP_CID_low_res_y_h2o_scoreZ5N;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h2o_scoreZ5N;
	static ScoreType ESI_TRAP_CID_low_res_y_soch4_scoreZ5N;
	static ScoreType ESI_TRAP_CID_low_res_yp2_soch4_scoreZ5N;
	static ScoreType ESI_TRAP_CID_low_res_y_h3po4_scoreZ5N;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ5N;

	static ScoreType ESI_TRAP_CID_low_res_y_scoreZ5C;
	static ScoreType ESI_TRAP_CID_low_res_yp2_scoreZ5C;
	static ScoreType ESI_TRAP_CID_low_res_yp3_scoreZ5C;
	static ScoreType ESI_TRAP_CID_low_res_y_nh3_scoreZ5C;
	static ScoreType ESI_TRAP_CID_low_res_yp2_nh3_scoreZ5C;
	static ScoreType ESI_TRAP_CID_low_res_y_h2o_scoreZ5C;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h2o_scoreZ5C;
	static ScoreType ESI_TRAP_CID_low_res_y_soch4_scoreZ5C;
	static ScoreType ESI_TRAP_CID_low_res_yp2_soch4_scoreZ5C;
	static ScoreType ESI_TRAP_CID_low_res_y_h3po4_scoreZ5C;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ5C;

	static ScoreType ESI_TRAP_CID_low_res_y_scoreZ5;
	static ScoreType ESI_TRAP_CID_low_res_yp2_scoreZ5;
	static ScoreType ESI_TRAP_CID_low_res_yp3_scoreZ5;
	static ScoreType ESI_TRAP_CID_low_res_y_nh3_scoreZ5;
	static ScoreType ESI_TRAP_CID_low_res_yp2_nh3_scoreZ5;
	static ScoreType ESI_TRAP_CID_low_res_y_h2o_scoreZ5;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h2o_scoreZ5;
	static ScoreType ESI_TRAP_CID_low_res_y_soch4_scoreZ5;
	static ScoreType ESI_TRAP_CID_low_res_yp2_soch4_scoreZ5;
	static ScoreType ESI_TRAP_CID_low_res_y_h3po4_scoreZ5;
	static ScoreType ESI_TRAP_CID_low_res_yp2_h3po4_scoreZ5;

	static int a_nh3_index;
	static int a_h2o_index;
	static int a_index;
	static int a_h3po4_index;
	static int b_h2o_index;
	static int b_plus_h2o_index;
	static int c_ladder_index;
	static int b_nh3_index;
	static int b_soch4_index;
	static int b_h3po4_index;
	static int b_index;
	static int cPlus2DaIndex;
	static int cPlus1DaIndex;
	static int c_index;
	static int cp2_index;
	static int cMinus1DaIndex;

	static int v_index;
	static int x_index;
	static int Y_index;
	static int y_h2o_index;
	static int y_nh3_index;
	static int y_soch4_index;
	static int y_h3po4_index;
	static int y_index;
	static int z_index;
	static int zp2_index;
	static int zPlus1Dap2_index;
	static int zPlus1DaIndex;
	static int zPlus2DaIndex;
	static int zPlus3DaIndex;

	static int bp2_index;
	static int bp3_index;
	static int bp2_nh3_index;
	static int bp2_h2o_index;
	static int bp2_soch4_index;
	static int bp2_h3po4_index;
	static int yp2_index;
	static int yp3_index;
	static int yp2_nh3_index;
	static int yp2_h2o_index;
	static int yp2_soch4_index;
	static int yp2_h3po4_index;

	static int numNIonTypes;
	static int numCIonTypes;
	static int numInternalIonTypes;

	static bool checkInternalIons;
	static bool checkNIons;
	static bool checkCIons;
	static bool checkSatelliteIons;
	static bool checkNIonAmmoniaLoss;
	static bool checkNIonPosChargeBearing;
	static bool checkNIonPhosphorylation;
	static bool checkNIonM;
	static bool checkNIonWaterLoss;
	static bool checkCIonPosChargeBearing;
	static bool checkCIonAmmoniaLoss;
	static bool checkCIonPhosphorylation;
	static bool checkCIonM;
	static bool checkCIonWaterLoss;

	static unsigned int dIonExcludeMask;
	static unsigned int vIonExcludeMask;
	static unsigned int wIonExcludeMask;
	static unsigned int ammoniaLossMask;
	static unsigned int posChargeBearingMask;
	static unsigned int waterLossMask;

	static MassType maxInternalIonMass;
	static size_t getHistogramLimit ( const MSTagParameters& params );
	void fMatch ( const std::string& sequence, int& numUnmatchedIons, ScoreType& score );
	static double xLinkMass;
	static double xLinkImmoniumMass;
	static ScoreType maxPScore;
	static ScoreType pCIDScore;
	static ScoreType pCIDXLScore;
	static ScoreType pCIDImmScore;
	static ScoreType pCIDH2OScore;
	static ScoreType pCIDNH3Score;
	static ScoreType pCID2H2OScore;
	static ScoreType pCIDXLH2OScore;
	static ScoreType pCIDXLNH3Score;
	static ScoreType pCIDImmH2OScore;
	static ScoreType pCIDImmNH3Score;
	static ScoreType pETDScore;
	static ScoreType pETDXLScore;
	static bool crosslinking;
	static int initialPosCharge;
public:
	MSTagSearch ( MSMSDataPoint* file, MSTagParameters& params );
	virtual ~MSTagSearch ();
	bool getSpectrumRetained () const { return spectrumRetained && peaks.getSpectrumRetained (); }
	void setSpectrumRetained ( bool flag ) { spectrumRetained = flag; }
	int getHistogramSize () const { return survHist.getSize (); }
	bool doMatch ( const std::string& peptide, bool nTermPeptide, bool cTermPeptide, double molWt, TagMatchVector& tagMatch, const ScoreType& minScore );
	void printExpectationHTML ( std::ostream& os, double score, int numSavedSpectra ) { survHist.printHTML ( os, score, numSavedSpectra ); }
	double getEvalue ( double score ) const { return survHist.getEValue ( score ); }
	bool getEValueFlag () const { return survHist.getEValueFlag (); }
	void printExpectationXML ( std::ostream& os, int numSavedSpectra ) { survHist.printXML ( os, numSavedSpectra ); }
	void fragmentMatch2 ( const PeptideSequence& ps, int& numUnmatchedIons, ScoreType& score, bool reset );
	static void setNTerminusWt ( const MassType nt ) { nTerminusWt = nt; }
	static void setCTerminusWt ( const MassType ct ) { cTerminusWt = ct; }
};

class MSTagSearchAllowErrors : public MSTagSearch {
	double precursorTolerance;
	static ModificationTable* modificationTable;
	static double maxNegParentError;
	static double maxPosParentError;
	static double maxParentError;
	static MassType defaultNTerminusWT;
	static MassType defaultCTerminusWT;
	static bool massMods;
	static void resetMSTagSearchAllowErrors ( const ModificationParameters& modificationParameters );
public:
	MSTagSearchAllowErrors ( MSMSDataPoint* file, MSTagParameters& params );
	virtual ~MSTagSearchAllowErrors ();
	bool doMatch ( const std::string& peptide, bool nTermPeptide, bool cTermPeptide, double mol_wt, TagMatchVector& tagMatch, const ScoreType& minScore );
	bool checkOffset ( double moleWt ) const
	{
		return modificationTable->checkOffset ( parentMass - moleWt - parentTolerance, parentMass - moleWt + parentTolerance );
	}
	bool checkMatch ( double molWt ) const
	{
		if ( genAbsDiff ( molWt, parentMass ) < parentTolerance ) {
			return true;
		}
		return modificationTable->checkMatch ( parentMass - molWt + parentTolerance );
	}
	double getParentTolerance () const { return ( parentTolerance + maxParentError ); }
	static void resetNextMods ();
	static ModificationTable* getModificationTable () { return modificationTable; }
};

MSMSSearch* getMSTagSearch ( MSMSDataPoint* file, MSTagParameters& params );

typedef std::vector <MSTagSearch> MSTagSearchContainer;
typedef MSTagSearchContainer::iterator MSTagSearchIterator;

class UnmatchedCompositionSearch {
	UIntVector ionMask;
	int numIons;
	static ScoreType immonium_score;
public:
	UnmatchedCompositionSearch ( const PeakContainer& peaks );
	int doCompositionSearch ( const std::string& peptide ) const;
	void doCompositionSearch ( const std::string& peptide, ScoreTypeVector& ionMatched ) const;
	static void setImmoniumScore ( ScoreType score ) { immonium_score = score; }
};

#endif /* ! __lu_tag_frag_h */
