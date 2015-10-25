/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_frag_mtch.cpp                                              *
*                                                                             *
*  Created    : July 21st 1996                                                *
*                                                                             *
*  Purpose    : msfit style search algorithms.                                *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1996-2009) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lgen_error.h>
#include <lu_frag_mtch.h>
#include <lu_getfil.h>
#include <lu_fit_srch.h>
#include <lu_mass_elem.h>
using std::fill;
using std::vector;

MSFitSearch* getMSFitSearch ( MSDataPoint* dataSet, const MSFitParameters& params )
{
	if ( params.getModificationParameters ().getAllowErrors () )
		return new MSFitAllowErrorsSearch ( dataSet, params );
	else {
		if ( params.getMultipleModification ().getModificationsAllowed () )
			return new MSFitModifiedSearch ( dataSet, params );
		else
			return new MSFitSearch ( dataSet, params );
	}
}

const char MowseScore::SCORE_MATCH = 1;

const double MowseScore::MAX_MOWSE_PROTEIN_MASS = 100000.0;
const double MowseScore::MOWSE_PROTEIN_AMU_BIN = 10000.0;
const int MowseScore::MAX_MOWSE_PROTEIN_BIN_INDEX = 10;		// MAX_MOWSE_PROTEIN_MASS / MOWSE_PROTEIN_AMU_BIN
const int MowseScore::NUM_MOWSE_PROTEIN_BINS = 11;			// MAX_MOWSE_PROTEIN_BIN_INDEX + 1

const double MowseScore::MAX_MOWSE_PEPTIDE_MASS = 3000.0;
const double MowseScore::MOWSE_PEPTIDE_AMU_BIN = 100.0;
const int MowseScore::NUM_MOWSE_PEPTIDE_BINS = 30;			// MAX_MOWSE_PEPTIDE_MASS / MOWSE_PEPTIDE_AMU_BIN

MowseScore::MowseScore ( double mowsePFactor ) :
	mowsePFactor ( mowsePFactor )
{
	mowseScore.resize ( NUM_MOWSE_PROTEIN_BINS );
	for ( int i = 0 ; i < NUM_MOWSE_PROTEIN_BINS ; i++ ) {
		mowseScore [i].resize ( NUM_MOWSE_PEPTIDE_BINS );
		fill ( mowseScore [i].begin (), mowseScore [i].end (), 0.0 );
	}
	statsNeedUpdating = true;
}
void MowseScore::setMowseArray ( double proteinMW )
{
	mowse = &(mowseScore [genMin ((int)(proteinMW/MOWSE_PROTEIN_AMU_BIN), MAX_MOWSE_PROTEIN_BIN_INDEX)]);
}
void MowseScore::accumulateMowseScore ( double fragmentMass, bool singleCleavage )
{
	if ( fragmentMass < MAX_MOWSE_PEPTIDE_MASS ) {
		((*mowse) [(int)(fragmentMass/MOWSE_PEPTIDE_AMU_BIN)]) += ( singleCleavage ) ? 1.0 : mowsePFactor;
		statsNeedUpdating = true;
	}
}
void MowseScore::assembleMowseStats ()
{
	double maxMowseScore = 0.0;
	int i, j;
	double total;

	for ( i = 0 ; i < NUM_MOWSE_PROTEIN_BINS ; i++ ) {
		for ( j = 0, total = 0.0 ; j < NUM_MOWSE_PEPTIDE_BINS ; j++ ) {
			total += mowseScore [i][j];
		}
		for ( j = 0 ; j < NUM_MOWSE_PEPTIDE_BINS ; j++ ) {
			if ( total != 0.0 ) mowseScore [i][j] /= total;
			maxMowseScore = genMax( maxMowseScore, mowseScore [i][j] );
		}
		for ( j = 0 ; j < NUM_MOWSE_PEPTIDE_BINS ; j++ ) {
			if ( maxMowseScore != 0.0 ) mowseScore [i][j] /= maxMowseScore;
		}
	}
	statsNeedUpdating = false;
/* Potentially useful diagnostics. These are the mowse normalised cell frequencies. */
/*
	cout << "<pre>" << endl;
	for ( i = 0 ; i < NUM_MOWSE_PROTEIN_BINS ; i++ ) {
		for ( j = 0 ; j < NUM_MOWSE_PEPTIDE_BINS ; j++ ) {
			if ( j != 0 ) cout << "\t";
			cout << setprecision ( 4 ) << mowseScore [i][j];
		}
		cout << endl;
	}
	cout << "</pre>" << endl;
*/
}
void MowseScore::calculateMowseScores ( HitStats* hs, double proteinMW, const DoubleVector& peakMass )
{
	if ( hs->mowseScore == 0.0 ) {
		if ( statsNeedUpdating ) {
			assembleMowseStats ();
		}
		double pn = 1.0;

		int numPeaks = peakMass.size ();
		for ( int i = 0 ; i < numPeaks ; i++ ) {
			if ( hs->massMatched [i] == SCORE_MATCH && peakMass [i] < MAX_MOWSE_PEPTIDE_MASS ) {
				pn *= mowseScore [genMin ((int)(proteinMW/MOWSE_PROTEIN_AMU_BIN), MAX_MOWSE_PROTEIN_BIN_INDEX)][(int)(peakMass[i]/MOWSE_PEPTIDE_AMU_BIN)];
				if ( hs->mowseMissedCleavages [i] != 0 ) pn /= mowsePFactor;
			}
		}
		if ( pn )											// Prevent division by zero
			hs->mowseScore = 50000.0 / ( pn * proteinMW );
		else
			hs->mowseScore = 0.0;
	}
}
const char MSFitSearch::SCORE_MATCH = 1;
const char MSFitSearch::NO_SCORE_MATCH = 2;

int MSFitSearch::maxFitPeaks = InfoParams::instance ().getIntValue ( "max_msfit_peaks", 1000 );
MSFitSearch::MSFitSearch ( MSDataPoint* dataSet, const MSFitParameters& params ) :
	peaks ( dataSet, params.getMSPeakFilterOptions (), params.getPeakContainerInfo () )
{
	if ( peaks.size () > maxFitPeaks ) {
		std::string err;
		err += "For MS-Fit the maximum number of peaks in a single spectrum after filtering is ";
		err += gen_itoa ( maxFitPeaks );
		err += ". \n";
		err += "Spectrum ";
		err +=  dataSet->getSpecID ();
		err +=  " has ";
		err += gen_itoa ( peaks.size () );
		err += " peaks.\n";
		err += "Make sure that the data has the isotope peaks removed.\n";
		err += "The charge needs to be specified if it is not 1.\n";
		err += "Singly charged peaks less than around 800 Da are not generally very useful for this type of search.\n";
		ErrorHandler::genError ()->error ( err );
	}
	if ( peaks.size () < params.getMinMatches () ) {
		ErrorHandler::genError ()->message ( "Spectrum " + dataSet->getSpecID () + " has less peaks than the minimum number required to match.\n" );
	}
	init ( params );
}
MSFitSearch::~MSFitSearch ()
{
	delete mowseScore;
}
void MSFitSearch::init ( const MSFitParameters& params )
{
	numPeaks = peaks.size ();
	peakMassLowerBound.resize ( numPeaks );
	tolerance.resize ( numPeaks );
	peakMass.resize ( numPeaks );
	massMatched.resize ( numPeaks );
	for ( int i = 0 ; i < numPeaks ; i++ ) {
		peakMass [i] = peaks [i]->getMass ();
		tolerance [i] = peaks [i]->getTolerance ();
		peakMassLowerBound [i] = peakMass [i] - tolerance [i];
	}
	const MowseInfo& mowseInfo = params.getMowseInfo ();
	if ( numPeaks > 0 ) {
		lowMass = peakMass [0] - tolerance [0];
		highMass = peakMass [numPeaks-1] + tolerance [numPeaks-1];
	}
	if ( mowseInfo.getMowseOn () ) {
		mowseMissedCleavages.resize ( numPeaks );
		mowseScore = new MowseScore ( mowseInfo.getMowsePFactor () );
	}
	else mowseScore = 0;
	const EnzymeParameters& enzymeParameters = params.getEnzymeParameters ();
	compMask = enzymeParameters.getCompMask ();
	compMaskTypeAnd = enzymeParameters.getCompMaskType () == "AND";
	compMaskTypeOr = enzymeParameters.getCompMaskType () == "OR";
	missedCleavages = enzymeParameters.getMissedCleavages ();
}
MSFitModifiedSearch::MSFitModifiedSearch ( MSDataPoint* dataSet, const MSFitParameters& params ) :
	MSFitSearch ( dataSet, params )
{
	init ( params );
}
MSFitModifiedSearch::~MSFitModifiedSearch ()
{
}
void MSFitModifiedSearch::init ( const MSFitParameters& params )
{
	const MultipleModification& multiMod = params.getMultipleModification ();
	bool monoisotopicFlag = params.getMonoisotopicFlag ();

	pyroglutamicAcidFlag	= multiMod.getPyroglutamicAcidFlag ();
	user1Flag				= multiMod.getUser1Flag ();
	user2Flag				= multiMod.getUser2Flag ();
	user3Flag				= multiMod.getUser3Flag ();
	user4Flag				= multiMod.getUser4Flag ();
	oxidationFlag			= multiMod.getOxidationFlag ();
	acetylationFlag			= multiMod.getAcetylationFlag ();
	incompleteCysFlag		= multiMod.getIncompleteCysFlag ();

	PossFragMods::initialiseMasks ( multiMod );
	acetylationMask		= PossFragMods::getAcetylationMask ();
	pyroglutamicAcidMask= PossFragMods::getPyroglutamicAcidMask ();
	oxidationMask		= PossFragMods::getOxidationMask ();
	cysMask				= PossFragMods::getCysMask ();
	user1Mask			= PossFragMods::getUser1Mask ();
	user2Mask			= PossFragMods::getUser2Mask ();
	user3Mask			= PossFragMods::getUser3Mask ();
	user4Mask			= PossFragMods::getUser4Mask ();

	ElementalFormula pyroglutamicAcidModElemForm ( "N-1 H-3" );
	ElementalFormula acetylationModElemForm ( "C-3 H-7 N-1 S-1" );
	ElementalFormula oxidationModElemForm ( "O" );
	ElementalFormula incompleteCysModElemForm ( "C6 H10 N2 O2 S1" );

	const vector <Usermod*>& userMod = multiMod.getUserMods ();
	if ( monoisotopicFlag ) {
		pyroglutamicAcidMod	= formula_to_monoisotopic_mass ( pyroglutamicAcidModElemForm );
		acetylationMod		= formula_to_monoisotopic_mass ( acetylationModElemForm );
		oxidationMod		= formula_to_monoisotopic_mass ( oxidationModElemForm );
		incompleteCysMod	= formula_to_monoisotopic_mass ( incompleteCysModElemForm ) - amino_acid_wt ['C'];	// The formula is for acrylamide modified cys
		user1Mod			= userMod.empty () ? 0.0 : formula_to_monoisotopic_mass ( userMod [0]->getElementalFormulaString ().c_str () );
		user2Mod			= ( userMod.size () < 2 ) ? 0.0 : formula_to_monoisotopic_mass ( userMod [1]->getElementalFormulaString ().c_str () );
		user3Mod			= ( userMod.size () < 3 ) ? 0.0 : formula_to_monoisotopic_mass ( userMod [2]->getElementalFormulaString ().c_str () );
		user4Mod			= ( userMod.size () < 4 ) ? 0.0 : formula_to_monoisotopic_mass ( userMod [3]->getElementalFormulaString ().c_str () );
	}
	else {
		pyroglutamicAcidMod	= formula_to_average_mass ( pyroglutamicAcidModElemForm );
		acetylationMod		= formula_to_average_mass ( acetylationModElemForm );
		oxidationMod		= formula_to_average_mass ( oxidationModElemForm );
		incompleteCysMod	= formula_to_average_mass ( incompleteCysModElemForm ) - amino_acid_wt ['C'];	// The formula is for acrylamide modified cys
		user1Mod			= userMod.empty () ? 0.0 : formula_to_average_mass ( userMod [0]->getElementalFormulaString ().c_str () );
		user2Mod			= ( userMod.size () < 2 ) ? 0.0 : formula_to_average_mass ( userMod [1]->getElementalFormulaString ().c_str () );
		user3Mod			= ( userMod.size () < 3 ) ? 0.0 : formula_to_average_mass ( userMod [2]->getElementalFormulaString ().c_str () );
		user4Mod			= ( userMod.size () < 4 ) ? 0.0 : formula_to_average_mass ( userMod [3]->getElementalFormulaString ().c_str () );
	}
	highModifiedMass = highMass - acetylationMod;	/* This is the highest possible mass with the current
													modifications as acetylation is the highest negative
													mass shift and with the other negative mass shift modification
													(pyro glu) only one modification is allowed per peptide */
}
ModificationTable* MSFitAllowErrorsSearch::modificationTable = 0;
double MSFitAllowErrorsSearch::maxParentError = 0.0;
MSFitAllowErrorsSearch::MSFitAllowErrorsSearch ( MSDataPoint* dataSet, const MSFitParameters& params ) :
	MSFitSearch ( dataSet, params )
{
	init ( params );
}
MSFitAllowErrorsSearch::~MSFitAllowErrorsSearch ()
{
}
void MSFitAllowErrorsSearch::init ( const MSFitParameters& params )
{
	minParentIonMatches = params.getMinParentIonMatches ();
	enzymeTerminalSpecificity = get_enzyme_terminal_specificity ();
	if ( modificationTable == 0 ) {
		modificationTable = new ModificationTable ( params.getModificationParameters () );
		maxParentError = modificationTable->getMaxParentError ();
	}
	lowMass += modificationTable->getMostNegMassShift ();
	highMass += modificationTable->getMostPosMassShift ();
	double max_tolerance = 0.0;
	for ( int i = 0 ; i < numPeaks ; i++ ) {
		peakMassLowerBound [i] -= modificationTable->getMostNegMassShift ();
		max_tolerance = genMax( max_tolerance, tolerance [i] );
	}
}
int MSFitSearch::matchFragments ( char* protein, const IntVector& cleavageIndex )
{
	int numFragments = cleavageIndex.size ();
	int j, k;

	if ( mowseScore ) {
		ProteinMW pmw ( protein );
		mowseScore->setMowseArray ( pmw.getMass () );
		fill ( mowseMissedCleavages.begin (), mowseMissedCleavages.end (), 0 );
	}
	numScoreMatches = 0;
	fill ( massMatched.begin (), massMatched.end (), 0 );
	DoubleVector& enzymeFragmentMassArray = get_cleaved_masses ( protein, cleavageIndex );
	int missedCleavageLimit = missedCleavages;
	for ( int i = 0 ; i < numFragments ; i++ ) {
		char* startPep = ( i == 0 ) ? protein : protein + cleavageIndex[i-1] + 1;
		double fragmentMass = terminal_wt;
		for ( k = 0, j = i ; j < numFragments ; j++ ) {
			if ( j > missedCleavageLimit ) break;
			char* endPep = protein + cleavageIndex[j] + 1;
			fragmentMass += enzymeFragmentMassArray [j];
			if ( cnbr_digest && j == missedCleavageLimit && protein [cleavageIndex [j]] == 'M' ) {
				fragmentMass += cnbr_homoserine_lactone_mod;
			}
			if ( mowseScore ) mowseScore->accumulateMowseScore ( fragmentMass, j == i );
			if ( fragmentMass > highMass ) break;
			if ( fragmentMass > lowMass ) {
				for ( ; k < numPeaks && fragmentMass > peakMassLowerBound [k] ; k++ ) {
					if ( !massMatched [k] ) {
						if ( genAbsDiff ( fragmentMass, peakMass [k] ) < tolerance [k] ) {
							if ( !compMask || checkComposition ( startPep, endPep - startPep ) ) {
								numScoreMatches++;
								massMatched [k] = SCORE_MATCH;
								if ( mowseScore ) mowseMissedCleavages [k] = j - i;
							}
						}
					}
				}
			}
		}
		missedCleavageLimit++;
	}
	return numScoreMatches;
}
int MSFitModifiedSearch::matchFragments ( char* protein, const IntVector& cleavageIndex )
{
	int numFragments = cleavageIndex.size ();
	int j, k;
	int numNoScoreMatches = 0;

	if ( mowseScore ) {
		ProteinMW pmw ( protein );
		mowseScore->setMowseArray ( pmw.getMass () );
		fill ( mowseMissedCleavages.begin (), mowseMissedCleavages.end (), 0 );
	}
	numScoreMatches = 0;
	fill ( massMatched.begin (), massMatched.end (), 0 );
	DoubleVector& enzymeFragmentMassArray = get_cleaved_masses ( protein, cleavageIndex );
	int missedCleavageLimit = missedCleavages;
	for ( int i = 0 ; i < numFragments ; i++ ) {
		char* startPep = ( i == 0 ) ? protein : protein + cleavageIndex[i-1] + 1;
		char* pep = startPep;
		unsigned int mask = aa_composition_mask [*pep];
		double fragmentMass = terminal_wt;
		int user1 = 0;
		int user2 = 0;
		int user3 = 0;
		int user4 = 0;
		int oxidation = 0;
		int pyroglutamic_acid = pyroglutamicAcidFlag && ( mask & pyroglutamicAcidMask );
		int acetylation = ( i == 0 ) ? acetylationFlag && ( mask & acetylationMask ) : 0;
		int incompleteCys = 0;
		for ( k = 0, j = i ; j < numFragments ; j++ ) {
			bool unmodified_match = false;
			if ( j > missedCleavageLimit ) break;
			fragmentMass += enzymeFragmentMassArray [j];
			if ( cnbr_digest && j == missedCleavageLimit && protein [cleavageIndex [j]] == 'M' ) {
				fragmentMass += cnbr_homoserine_lactone_mod;
				if ( oxidationFlag ) oxidation--;
			}
			if ( mowseScore ) mowseScore->accumulateMowseScore ( fragmentMass, j == i );
			if ( fragmentMass > highModifiedMass ) break;		/* Don't bother totting up modification sites if mass too high */
			char* endPep = protein + cleavageIndex[j] + 1;
			for ( ; pep < endPep ; pep++ ) {
				mask = aa_composition_mask [*pep];
				if ( user1Flag )		if ( mask & user1Mask )		user1++;
				if ( user2Flag )		if ( mask & user2Mask )		user2++;
				if ( user3Flag )		if ( mask & user3Mask )		user3++;
				if ( user4Flag )		if ( mask & user4Mask )		user4++;
				if ( oxidationFlag )	if ( mask & oxidationMask )	oxidation++;
				if ( incompleteCysFlag )if ( mask & cysMask )		incompleteCys++;
			}
			if ( oxidation || pyroglutamic_acid || acetylation || user1 || user2 || user3 || user4 || incompleteCys ) {
				for ( int a = 0 ; a <= pyroglutamic_acid ; a++ ) {
					for ( int b = 0 ; b <= user1 ; b++ ) {
						for ( int c = 0 ; c <= user2 ; c++ ) {
							for ( int d = 0 ; d <= user3 ; d++ ) {
								for ( int e = 0 ; e <= user4 ; e++ ) {
									for ( int f = 0 ; f <= acetylation ; f++ ) {
										int oxidationSites = ( f == 1 ) ? oxidation - 1 : oxidation;
										for ( int g = 0 ; g <= oxidationSites ; g++ ) {
											for ( int h = 0 ; h <= incompleteCys ; h++ ) {
												double fragmentMass2 = fragmentMass + b * user1Mod + c * user2Mod + d * user3Mod + e * user4Mod + g * oxidationMod + h * incompleteCysMod;
												if ( a ) fragmentMass2 += pyroglutamicAcidMod;
												if ( f ) fragmentMass2 += acetylationMod;
												if ( fragmentMass2 > lowMass && fragmentMass2 < highMass ) {
													for ( k = 0 ; k < numPeaks && fragmentMass2 > peakMassLowerBound [k] ; k++ ) {
														if ( !massMatched [k] ) {
															if ( genAbsDiff ( fragmentMass2, peakMass [k] ) < tolerance [k] ) {
																if ( !compMask || checkComposition ( startPep, endPep - startPep ) ) {
																	if ( b || f ) {
																		massMatched [k] = NO_SCORE_MATCH;
																		numNoScoreMatches++;
																	}
																	else {
																		if ( !a && !g && !h ) unmodified_match = true;
																		if ( ( a || g ) && unmodified_match == false ) {
																			massMatched [k] = NO_SCORE_MATCH;
																			numNoScoreMatches++;
																		}
																		else {
																			massMatched [k] = SCORE_MATCH;
																			numScoreMatches++;
																		}
																	}
																	if ( mowseScore ) mowseMissedCleavages [k] = j - i;
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
			else {
				if ( fragmentMass > highMass ) break;
				if ( fragmentMass > lowMass ) {
					for ( ; k < numPeaks && fragmentMass > peakMassLowerBound [k] ; k++ ) {
						if ( !massMatched [k] ) {
							if ( genAbsDiff ( fragmentMass, peakMass [k] ) < tolerance [k] ) {
								if ( !compMask || checkComposition ( startPep, endPep - startPep ) ) {
									massMatched [k] = SCORE_MATCH;
									numScoreMatches++;
									if ( mowseScore ) mowseMissedCleavages [k] = j - i;
								}
							}
						}
					}
				}
			}
		}
		missedCleavageLimit++;
	}
	return numScoreMatches + numNoScoreMatches;
}
int MSFitAllowErrorsSearch::matchFragments ( char* protein, const IntVector& cleavageIndex )
{
	int numFragments = cleavageIndex.size ();
	int j, k, m;
	int numNoScoreMatches = 0;

	if ( mowseScore ) {
		ProteinMW pmw ( protein );
		mowseScore->setMowseArray ( pmw.getMass () );
		fill ( mowseMissedCleavages.begin (), mowseMissedCleavages.end (), 0 );
	}
	fill ( massMatched.begin (), massMatched.end (), 0 );
	numScoreMatches = 0;
	DoubleVector& enzymeFragmentMassArray = get_cleaved_masses ( protein, cleavageIndex );
	int missedCleavageLimit = missedCleavages;
	for ( int i = 0 ; i < numFragments ; i++ ) {
		char* startPep = ( i == 0 ) ? protein : protein + cleavageIndex[i-1] + 1;
		double fragmentMass = terminal_wt;
		for ( j = i ; j < numFragments ; j++ ) {
			if ( j > missedCleavageLimit ) break;
			char* endPep = protein + cleavageIndex[j] + 1;
			fragmentMass += enzymeFragmentMassArray [j];
			if ( cnbr_digest && j == missedCleavageLimit && protein [cleavageIndex [j]] == 'M' ) {
				fragmentMass += cnbr_homoserine_lactone_mod;
			}
			if ( mowseScore ) mowseScore->accumulateMowseScore ( fragmentMass, j == i );
			if ( fragmentMass > highMass ) break;
			if ( fragmentMass > lowMass ) {
				for ( k = 0 ; k < numPeaks ; k++ ) {
					if ( massMatched [k] == 0 || numScoreMatches < minParentIonMatches ) {
						if ( genAbsDiff ( fragmentMass, peakMass [k] ) < tolerance [k] ) {
							if ( !compMask || checkComposition ( startPep, endPep - startPep ) ) {
								if ( massMatched [k] == NO_SCORE_MATCH ) numNoScoreMatches--;
								massMatched [k] = SCORE_MATCH;
								numScoreMatches++;
								if ( mowseScore ) mowseMissedCleavages [k] = j - i;
							}
						}
						else if ( genAbsDiff ( fragmentMass, peakMass [k] ) < tolerance [k] + maxParentError ) {
							if ( !compMask || checkComposition ( startPep, endPep - startPep ) ) {
								double diff = peakMass [k] - fragmentMass;

								ModificationPair modPair = modificationTable->getPossibleModifications ( diff - tolerance [k], diff + tolerance [k] );
								if ( modPair.second ) {
									for ( ModificationVectorConstIterator l = modPair.first ; l < modPair.first + modPair.second ; l++ ) {
										int startM = ( i == 0 ) ? 0 : cleavageIndex [i-1] + 1;
										int endM = cleavageIndex [j];
										char terminalSpec = (*l)->getTerminalSpecificity ();
										char expectedResidue = (*l)->getExpectedResidue ();
										if ( terminalSpec == 'C' ) startM = endM;
										if ( terminalSpec ) endM = startM;
										if ( enzymeTerminalSpecificity == 'N' && i != 0 && terminalSpec != 'C' ) startM++;
										if ( enzymeTerminalSpecificity == 'C' && j != numFragments - 1 && terminalSpec != 'N' ) endM--;
										for ( m = startM ; m <= endM ; m++ ) {
											if ( protein [m] == expectedResidue ) {
												if ( massMatched [k] == 0 ) numNoScoreMatches++;
												massMatched [k] = NO_SCORE_MATCH;
												if ( mowseScore ) mowseMissedCleavages [k] = j - i;
												goto label;
											}
										}
									}
									label:;
								}
							}
						}
					}
				}
			}
		}
		missedCleavageLimit++;
	}
	if ( numScoreMatches < minParentIonMatches ) return ( 0 );
	else return numNoScoreMatches + numScoreMatches;
}
HitStats MSFitSearch::getProteinFragmentStats ()
{
	HitStats hs;
	hs.massMatched = massMatched;
	if ( mowseScore ) {
		hs.mowseMissedCleavages = mowseMissedCleavages;
	}
	return hs;
}
void MSFitSearch::calculateMowseScores ( HitStats* hs, double proteinMW )
{
	mowseScore->calculateMowseScores ( hs, proteinMW, peakMass );
}
bool MSFitSearch::checkComposition ( const char* fragment, int len )
{
	if ( compMaskTypeAnd ) {
		if ( ( compMask & string_to_mask ( fragment, len ) ) == compMask ) {
			return true;
		}
		return false;
	}
	if ( compMaskTypeOr ) {
		if ( ( compMask & string_to_mask ( fragment, len ) ) != 0 ) {
			return true;
		}
		return false;
	}
	return true;
}
