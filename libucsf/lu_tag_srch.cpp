/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_tag_srch.cpp                                               *
*                                                                             *
*  Created    : July 21st 1996                                                *
*                                                                             *
*  Purpose    : Interface code to the MS-Seq and MS-Tag search algorithms.    *
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
#include <queue>
#include <lg_io.h>
#include <lgen_reg_exp.h>
#include <lu_aa_calc.h>
#include <lu_mat_score.h>
#include <lu_html.h>
#include <lp_frame.h>
#include <lu_tag_srch.h>
#include <lu_delim.h>
#include <lu_tag_par.h>
#include <lu_table.h>
using std::vector;
using std::ostringstream;
using std::ostream;
using std::string;
using std::endl;
using std::lower_bound;
using std::random_shuffle;
using std::stable_sort;
using std::set_intersection;
using std::inserter;

static Histogram hist;

TagSearch* getTagSearch ( MSTagParameters& tagParams )
{
	if ( tagParams.getDataSetInfo ()->getNoDataFlag () ) {
		ErrorHandler::genError ()->error ( "No spectra in data file.\n" );
		return 0;
	}
	if ( tagParams.getNoEnzyme () ) {
		if ( tagParams.getAllowErrors () )
			return new NoEnzymeSearch ( tagParams );
		else
			return new NoEnzymeIntSearch ( tagParams );
	}
	else {
		if ( tagParams.getAllowNonSpecificType () == '0' )
			return new YesEnzymeSearch ( tagParams );
		else {
			if ( tagParams.getAllowErrors () )
				return new NoEnzymeSearch ( tagParams );
			else
				return new NoEnzymeIntSearch ( tagParams );
		}
	}
}
TagHits::TagHits ( const vector <MSMSSearch*>& msMSSearch, MSTagParameters& params ) :
	DatabaseHits (),
	numSearches ( msMSSearch.size () )
{
	for ( int i = 0 ; i < numSearches ; i++ )
		hits.push_back ( new TagHitsContainer ( msMSSearch [i], params ) );
}
TagHits::~TagHits ()
{
	for ( int i = 0 ; i < numSearches ; i++ )
		delete hits [i];
}
void TagHits::checkRegularExpression ( RegularExpression& regExp )
{
	for ( int i = 0 ; i < numSearches ; i++ ) {
		hits [i]->checkRegularExpression ( regExp );
	}
}
void TagHits::sortAndRank ()
{
	for ( int i = 0 ; i < numSearches ; i++ ) {
		hits [i]->rankHits ();
	}
}
void TagHits::prune ()
{
	for ( int i = 0 ; i < numSearches ; i++ ) {
		hits [i]->pruneHits ();
	}
}
void TagHits::addHit ( const TagHit& hit, int dataSet )
{
	hits [dataSet]->push_back ( hit );
}
class sortMSMSSearchByParentMass {
public:
	int operator () ( const MSMSSearch* a, const double b ) const
	{
		return ( a->getParentMassPlusNegTolerance () < b );
	}
};
char TagSearch::allowNonSpecificType;
int TagSearch::missedCleavages;
TagSearch::TagSearch ( MSTagParameters& params ) :
	DatabaseSearch ( params ),
	tagParams ( params ),
	maxTagMatches ( params.getMaxHits () ),
	randomSearch ( params.isRandomSearch () )
{
	SurvivalHistogram::setExpectationMethod ( params.getExpectationMethod () );
	MSMSDataSetInfo* dsi = params.getDataSetInfo ();
	for ( int i = 0 ; i < dsi->getNumDataSets () ; i++ ) {
		msMSSearch.push_back ( getMSTagSearch ( dsi->getDataSet ( i ), params ) );
	}
	tagHits = new TagHits ( msMSSearch, params );
	numSearches = msMSSearch.size ();
	EnzymeParameters enzymeParameters = params.getEnzymeParameters ();
	compMask = enzymeParameters.getCompMask ();
	compMaskTypeAnd = enzymeParameters.getCompMaskType () == "AND";
	compMaskTypeOr = enzymeParameters.getCompMaskType () == "OR";
	initNonSpecific ( params );
}
TagSearch::~TagSearch ()
{
	delete tagHits;
	for ( int i = 0 ; i < msMSSearch.size () ; i++ ) {
		delete msMSSearch [i];
	}
}
void TagSearch::initNonSpecific ( const MSTagParameters& params )
{
	allowNonSpecificType = params.getAllowNonSpecificType ();
	if ( allowNonSpecificType != '0' && params.getEnzyme ().find ( "CNBr" ) != string::npos ) {
		ErrorHandler::genError ()->error ( "Allow non specific cleavages option not supported for CNBr digest type.\n" );
	}
	missedCleavages = params.getMissedCleavages ();
	if ( allowNonSpecificType != 'E' ) {
		init_fasta_enzyme_function ( params.getEnzyme () );
	}
	checkD = false;
	switch ( allowNonSpecificType ) {
		case 'C':
			nonSpecificTermini = &nonSpecificCTermini;
			break;
		case 'N':
			nonSpecificTermini = &nonSpecificNTermini;
			break;
		case 'D':
			checkD = true;
			nonSpecificTermini = &nonSpecificNTermini;
			break;
		case '1':
			nonSpecificTermini = &nonSpecific1Termini;
			break;
		case '2':
			nonSpecificTermini = &nonSpecific2Termini;
			break;
		case 'E':
			nonSpecificTermini = &noEnzyme;
			break;
	}
}
bool TagSearch::nonSpecificCTermini ( const IntVector& cleavageIndex, int start, int end )
{
	bool flag = false;
	IntVectorConstIterator ii;
	if ( start == 0 ) {
		ii = cleavageIndex.begin ();
		flag = true;
	}
	else {
		int s = start-1;
		ii = lower_bound ( cleavageIndex.begin (), cleavageIndex.end (), s );
		flag = ( *ii == s );
		ii++;
	}
	if ( flag ) {	// Check number of missed cleavages
		int e = end-1;
		IntVectorConstIterator jj = lower_bound ( cleavageIndex.begin (), cleavageIndex.end (), e );
		flag = ( jj - ii <= missedCleavages );
	}
	return flag;
}
bool TagSearch::nonSpecificNTermini ( const IntVector& cleavageIndex, int start, int end )
{
	int e = end-1;
	IntVectorConstIterator jj = lower_bound ( cleavageIndex.begin (), cleavageIndex.end (), e );
	bool flag = ( *jj == e );
	if ( flag ) {	// Check number of missed cleavages
		IntVectorConstIterator ii;
		if ( start < cleavageIndex [0] ) {
			ii = cleavageIndex.begin ();
		}
		else {
			int s = start-1;
			ii = lower_bound ( cleavageIndex.begin (), cleavageIndex.end (), s );
			if ( *ii == s ) ii++;
		}
		flag = ( jj - ii <= missedCleavages );
	}
	return flag;
}
bool TagSearch::nonSpecific1Termini ( const IntVector& cleavageIndex, int start, int end )
{
	return nonSpecificCTermini ( cleavageIndex, start, end ) || nonSpecificNTermini ( cleavageIndex, start, end );
}
bool TagSearch::nonSpecific2Termini ( const IntVector& cleavageIndex, int start, int end )
{	// Just check missed cleavages
	IntVectorConstIterator ii;
	if ( start < cleavageIndex [0] ) {
		ii = cleavageIndex.begin ();
	}
	else {
		int s = start-1;
		ii = lower_bound ( cleavageIndex.begin (), cleavageIndex.end (), s );
		if ( *ii == s ) ii++;
	}
	int e = end-1;
	IntVectorConstIterator jj = lower_bound ( cleavageIndex.begin (), cleavageIndex.end (), e );
	return ( jj - ii <= missedCleavages );
}
void TagSearch::doSearch ()
{
	int pruneInterval = getPruneInterval ( numSearches );
	if ( randomSearch ) {
		char* readingFrame;
		ProteinHit::addFS ( fs [0], 0 );
		srand ( static_cast<unsigned int> ( time ( 0 ) ) );
		int aver = 0;
		do {
			FrameIterator* fi = new FrameIterator ( fs [0], params.getIndicies ( 0 ), dnaFrameTranslationPairVector [0], params.getTempOverride () );
			for ( int i = 1 ; ( readingFrame = fi->getNextFrame () ) != NULL ; i++ ) {
				random_shuffle ( readingFrame, readingFrame + strlen ( readingFrame ) );
				tag_search ( *fi, readingFrame );
				if ( i % ( getActualPruneInterval ( i, pruneInterval ) ) == 0 ) tagHits->prune ();
			}
			int minProcessed;
			if ( aver == 0 ) {		// First time through
				for ( int j = 0 ; j < numSearches ; j++ ) {
					aver += static_cast <MSTagSearch*> (msMSSearch [j])->getHistogramSize ();
				}
				aver /= numSearches;
				minProcessed = aver / 5;
				for ( int k = 0 ; k < numSearches ; k++ ) {
					if ( static_cast <MSTagSearch*> (msMSSearch [k])->getHistogramSize () < minProcessed ) {
						static_cast <MSTagSearch*> (msMSSearch [k])->setSpectrumRetained ( false );
					}
				}
			}
			delete fi;
			FrameIterator::decrementNum ();
		} while ( continueSearch () );
		FrameIterator::incrementNum ();
	}
	else {
		for ( int i = 0 ; i < fs.size () ; i++ ) {
			doNormalSearch ( fs [i], i );
			FrameIterator::resetElapsedTime ( 1 );
		}
	}
	FrameIterator::searchFinished ();

	if ( tagParams.getRegExpSearch () ) {
		RegularExpression regExp ( tagParams.getRegularExpression () );
		tagHits->checkRegularExpression ( regExp );
	}
	tagHits->sortAndRank ();
	numHits = tagHits->size ();
	databaseHits = tagHits;
}
void TagSearch::doNormalSearch ( FastaServer* fsPtr, int num )
{
	char* readingFrame;
	ProteinHit::addFS ( fsPtr, num );
	int pruneInterval = getPruneInterval ( numSearches );
	FrameIterator fi ( fsPtr, params.getIndicies ( num ), dnaFrameTranslationPairVector [num], params.getTempOverride () );
	for ( int i = 1 ; ( readingFrame = fi.getNextFrame () ) != NULL ; i++ ) {
		tag_search ( fi, readingFrame );
		if ( i % ( getActualPruneInterval ( i, pruneInterval ) ) == 0 ) tagHits->prune ();
	}
}
bool TagSearch::continueSearch () const
{
	int n = 0;
	bool flag = false;
	for ( int i = 0 ; i < numSearches ; i++ ) {
		if ( msMSSearch [i]->getSpectrumRetained () ) {
			n++;
			flag = true;
		}
	}
	if ( ((double)n / numSearches) < 0.1 ) return false;
	return flag;
}
void TagSearch::addTagHit ( int searchNumber, const TagMatchVector& tagMatch, const string& sequence, const FrameIterator& fi, int previousAA, int nextAA, int startAA )
{
	for ( TagMatchVectorSizeType i = 0 ; i < tagMatch.size () ; i++ ) {
		if ( tagHits->size ( searchNumber ) > maxTagMatches ) {
			ostringstream err;
			err << params.getReportTitle ();
			err << " has matched the data to too many entries (>" << maxTagMatches << ").\n";
			err << "Please increase the specificity of the search parameters.\n";
			ErrorHandler::genError ()->error ( err.str () );
		}
		tagHits->addHit ( TagHit ( fi, sequence, previousAA, nextAA, startAA, tagMatch [i] ), searchNumber );
	}
}
void TagSearch::printParamsBodyHTML ( ostream& os ) const
{
	tagParams.printHTML ( os );
}
void TagSearch::printXMLHits ( ostream& os ) const
{
	for ( int i = 0 ; i < numSearches ; i++ ) {
		os << "<d>\n";
		delimitedRowStart ( os );
			tagParams.getDataSetInfo ()->printSpecInfoDelimited ( os, i );
			delimitedCell ( os, tagHits->getNumPeaks ( i ) );
			tagHits->getParentPeak (i)->printDelimited ( os, instInf->getParentPeakPrecision () );
		delimitedRowEnd ( os );
		msMSSearch [i]->printExpectationXML ( os, tagHits->getNumSavedHits ( i ) );
		if ( !randomSearch ) tagHits->printDelimited ( os, i );
		os << "</d>\n";
	}
}
void TagSearch::printBodyXML ( ostream& os, bool showPreSearch )
{
	if ( showPreSearch ) tagParams.printPreSearchBodyXML ( os );
	printXMLHits ( os );
}
void TagSearch::printHTMLHitsReport ( ostream& os )
{
	if ( tagParams.getScoreHistogramOnly () ) hist.drawGraph ( os );
	else {
		if ( numSearches == 1 ) {
			os << "<div class=\"results_header\">Results</div>" << endl;
			if ( tagHits->getNumSavedHits ( 0 ) ) msMSSearch [0]->printExpectationHTML ( os, tagHits->getScore ( 0, 0 ), tagHits->getNumSavedHits ( 0 ) );
			tagHits->printHTMLReport ( os, 0 );
		}
		else {
			for ( int i = 0 ; i < numSearches ; i++ ) {
				os << "<div class=\"results_header\">Results for Spot ID " << tagParams.getDataSetInfo ()->getSpot ( i ) << "</div>" << endl;
				if ( tagHits->getNumSavedHits ( i ) ) msMSSearch [i]->printExpectationHTML ( os, tagHits->getScore ( i, 0 ), tagHits->getNumSavedHits ( i ) );
				tagHits->printHTMLReport ( os, i );
				os << "<br /><br />" << endl;
			}
		}
	}
}
bool TagSearch::checkComposition ( const string& fragment )
{
	if ( compMaskTypeAnd ) {
		if ( ( compMask & string_to_mask ( fragment ) ) == compMask ) {
			return true;
		}
		return false;
	}
	if ( compMaskTypeOr ) {
		if ( ( compMask & string_to_mask ( fragment ) ) != 0 ) {
			return true;
		}
		return false;
	}
	return true;
}
NoEnzymeSearch::NoEnzymeSearch ( MSTagParameters& params ) :
	TagSearch ( params )
{
	for ( int i = 0 ; i < numSearches ; i++ ) {
		startMasses.push_back ( msMSSearch[i]->getParentMassMinusPosTolerance () );
		endMasses.push_back ( msMSSearch[i]->getParentMassPlusNegTolerance () );
	}
	doSearch ();
}
NoEnzymeSearch::~NoEnzymeSearch () {}
void NoEnzymeSearch::tag_search ( const FrameIterator& fi, const char* frame )
{
	const IntVector& cleavageIndex = ( allowNonSpecificType == 'E' ) ? IntVector () : enzyme_fragmenter ( frame );
	int num_aa = strlen ( frame );
	double* startProtein = get_protein_double_mass_array ( frame );
	double* endProtein = startProtein + num_aa;
	for ( int i = 0 ; i < numSearches ; i++ ) {
		if ( msMSSearch[i]->getSpectrumRetained () ) {
			double startMass = startMasses [i];
			double endMass = endMasses [i];
			double* start = startProtein;
			double* end = startProtein;
			double sum = terminal_wt;
			for ( ; end <= endProtein || sum >= startMass ; ) {
				if ( sum < startMass ) {
					sum += *end++;
				}
				else if ( sum > endMass ) {
					sum -= *start++;
					while ( sum > startMass && end - start > 1 ) {
						sum -= *--end;
					}
				}
				else {
					if ( end - start > 0 ) {
						int s = start-startProtein;
						int e = end-startProtein;
						if ( !checkD || ( s != 0 && frame [s-1] == 'D' ) ) {
							if ( nonSpecificTermini ( cleavageIndex, s, e ) ) {
								string possMatch ( &frame [s], end - start );
								if ( msMSSearch [i]->doMatch ( possMatch, ( start == startProtein ), ( end == endProtein ), sum, tagMatch, tagHits->getMinScore ( i ) ) ) {
									if ( !compMask || checkComposition ( possMatch ) ) {
										addTagHit ( i, tagMatch, possMatch, fi, ( start == startProtein ) ? '-' : frame [s-1], ( end == endProtein ) ? '-' : frame [e], s+1 );
									}
								}
							}
						}
						if ( end == endProtein ) {
							sum -= *start++;
							while ( sum > startMass && end - start > 1 ) {
								sum -= *--end;
							}
						}
						else {
							sum += *end++;
						}
					}
					else {
						if ( end != endProtein ) sum += *end++;
						else break;
					}
				}
			}
		}
	}
}
NoEnzymeIntSearch::NoEnzymeIntSearch ( MSTagParameters& params ) :
	TagSearch ( params )
{
	for ( int i = 0 ; i < numSearches ; i++ ) {
		double nom_parent_mass = msMSSearch[i]->getParentMass () / nominal_mass_offset;
		int integer_tolerance = (int)(msMSSearch[i]->getParentTolerance ()) + 1;
		startMasses.push_back ( (int)nom_parent_mass - integer_tolerance );
		endMasses.push_back ( (int)nom_parent_mass + integer_tolerance );
	}
	doSearch ();
}
NoEnzymeIntSearch::~NoEnzymeIntSearch () {}
/* Bug fixed version - needs checking
void NoEnzymeIntSearch::tag_search ( const FrameIterator& fi, const char* frame )
{
	const IntVector& cleavageIndex = ( allowNonSpecificType == 'E' ) ? IntVector () : enzyme_fragmenter ( frame );
	int num_aa = strlen ( frame );
	int* startProtein = get_protein_int_mass_array ( frame );
	int* endProtein = startProtein + num_aa;
	for ( int i = 0 ; i < numSearches ; i++ ) {
		if ( msMSSearch[i]->getSpectrumRetained () ) {
			int startMass	= startMasses [i];
			int endMass	= endMasses [i];
			double parentMass = msMSSearch[i]->getParentMass ();
			double parentTolerance = msMSSearch[i]->getParentTolerance ();
			int* start = startProtein;
			int* end = startProtein;
			int sum = (int)terminal_wt;
			for ( ; end <= endProtein || sum >= startMass ; ) {
				if ( sum < startMass ) {
					sum += *end++;
				}
				else if ( sum > endMass ) {
					sum -= *start++;

					while ( sum > startMass && end - start > 1 ) {
						sum -= *--end;
					}
				}
				else {
					if ( end - start > 0 ) {
						int s = start-startProtein;
						int e = end-startProtein;
						if ( !checkD || ( s != 0 && frame [s-1] == 'D' ) ) {
							if ( nonSpecificTermini ( cleavageIndex, s, e ) ) {
								string possMatch ( &frame [s], end - start );
								double mol_wt = peptide_formula_to_molecular_weight ( possMatch.c_str () );
								if ( genAbsDiff ( mol_wt, parentMass ) < parentTolerance ) {
									if ( msMSSearch [i]->doMatch ( possMatch, ( start == startProtein ), ( end == endProtein ), mol_wt, tagMatch, tagHits->getMinScore ( i ) ) ) {
										if ( !compMask || checkComposition ( possMatch ) ) {
											addTagHit ( i, tagMatch, possMatch, fi, ( start == startProtein ) ? '-' : frame [s-1], ( end == endProtein ) ? '-' : frame [e], s+1 );
										}
									}
								}
							}
						}
						if ( end == endProtein ) {
							sum -= *start++;
							while ( sum > startMass && end - start > 1 ) {
								sum -= *--end;
							}
						}
						else {
							sum += *end++;
						}
					}
					else {
						if ( end != endProtein ) sum += *end++;
						else break;
					}
				}
			}
		}
	}
}
*/
void NoEnzymeIntSearch::tag_search ( const FrameIterator& fi, const char* frame )
{
	const IntVector& cleavageIndex = ( allowNonSpecificType == 'E' ) ? IntVector () : enzyme_fragmenter ( frame );
	int num_aa = strlen ( frame );
	int* startProtein = get_protein_int_mass_array ( frame );
	int* endProtein = startProtein + num_aa;
	for ( int i = 0 ; i < numSearches ; i++ ) {
		if ( msMSSearch[i]->getSpectrumRetained () ) {
			int startMass	= startMasses [i];
			int endMass	= endMasses [i];
			double parentMass = msMSSearch[i]->getParentMass ();
			double parentTolerance = msMSSearch[i]->getParentTolerance ();
			int* start = startProtein;
			int* end = startProtein;
			int sum = (int)terminal_wt;
			for ( ; end <= endProtein ; ) {
				if ( sum < startMass ) {
					sum += *end++;
				}
				else if ( sum > endMass ) {
					sum -= *start++;
				}
				else {
					int s = start-startProtein;
					int e = end-startProtein;
					if ( !checkD || ( s != 0 && frame [s-1] == 'D' ) ) {
						if ( nonSpecificTermini ( cleavageIndex, s, e ) ) {
							string possMatch ( &frame [s], end - start );
							double mol_wt = peptide_formula_to_molecular_weight ( possMatch.c_str () );
							if ( genAbsDiff ( mol_wt, parentMass ) < parentTolerance ) {
								if ( msMSSearch [i]->doMatch ( possMatch, ( start == startProtein ), ( end == endProtein ), mol_wt, tagMatch, tagHits->getMinScore ( i ) ) ) {
									if ( !compMask || checkComposition ( possMatch ) ) {
										addTagHit ( i, tagMatch, possMatch, fi, ( start == startProtein ) ? '-' : frame [s-1], ( end == endProtein ) ? '-' : frame [e], s+1 );
									}
								}
							}
						}
					}
					sum -= *start++;
				}
			}
		}
	}
}
YesEnzymeSearch::YesEnzymeSearch ( MSTagParameters& params ) :
	TagSearch ( params )
{
	startLimit = msMSSearch.front ()->getParentMassMinusPosTolerance ();
	endLimit = msMSSearch.back ()->getParentMassPlusNegTolerance ();
	cleavedLimit = ( cnbr_digest ) ? endLimit - cnbr_homoserine_lactone_mod : endLimit;
	doSearch ();
}
YesEnzymeSearch::~YesEnzymeSearch () {}
void YesEnzymeSearch::tag_search ( const FrameIterator& fi, const char* frame )
{
	IntVector& cleavage_index = enzyme_fragmenter ( frame );
	int numEnzymeFragments = cleavage_index.size ();
	DoubleVector& enzyme_fragment_mass_array = get_cleaved_masses_to_limit ( frame, cleavage_index, cleavedLimit );
	int missedCleavageLimit = missedCleavages;
	bool first;
	int hitLength;
	string possMatch;
	possMatch.reserve ( 80 );

	if ( cnbr_digest ) {
		for ( int i = 0 ; i < numEnzymeFragments ; i++ ) {
			bool nTermProt = ( i == 0 );
			double mol_wt = terminal_wt;
			int previousCleavageIndex = ( i == 0 ) ? -1 : cleavage_index [i-1];
			for ( int j = i ; j <= missedCleavageLimit ; j++ ) {
				if ( j >= numEnzymeFragments ) break;
				mol_wt += enzyme_fragment_mass_array [j];
				if ( frame [cleavage_index [j]] == 'M' ) mol_wt += cnbr_homoserine_lactone_mod;	//CNBr only							//CNBr only
				if ( mol_wt > endLimit ) break;
				if ( mol_wt > startLimit ) {
					MSTagSearchAllowErrors::resetNextMods ();
					int ind = lower_bound ( msMSSearch.begin (), msMSSearch.end (), mol_wt, sortMSMSSearchByParentMass () ) - msMSSearch.begin ();
					first = true;
					for ( ; ind < msMSSearch.size () ; ind++ ) {
						MSMSSearch* msms = msMSSearch [ind];
						if ( msms->getSpectrumRetained () ) {
							if ( msms->getParentMassMinusPosTolerance () > mol_wt ) break;
							if ( first ) {
								hitLength = cleavage_index [j] - previousCleavageIndex;
								possMatch = string ( frame + previousCleavageIndex + 1, hitLength );
								first = false;
							}
							if ( !compMask || checkComposition ( possMatch ) ) {
								bool cTermProt = ( j == numEnzymeFragments-1 );
								if ( possMatch [hitLength - 1] == 'M' ) possMatch[hitLength-1] = 'h';	//CNBr only
								if (  msms->doMatch ( possMatch, nTermProt, cTermProt, mol_wt, tagMatch, tagHits->getMinScore ( ind ) ) ) {
									addTagHit ( ind, tagMatch, possMatch, fi, nTermProt ? '-' : frame[previousCleavageIndex], cTermProt ? '-' : frame [cleavage_index[j] + 1], previousCleavageIndex + 2 );
								}
							}
						}
					}
				}
				if ( frame [cleavage_index [j]] == 'M' ) mol_wt -= cnbr_homoserine_lactone_mod;	//CNBr only							//CNBr only
			}
			missedCleavageLimit++;
		}
	}
	else {
		for ( int i = 0 ; i < numEnzymeFragments ; i++ ) {
			bool nTermProt = ( i == 0 );
			double mol_wt = terminal_wt;
			int previousCleavageIndex = ( i == 0 ) ? -1 : cleavage_index [i-1];
			for ( int j = i ; j <= missedCleavageLimit ; j++ ) {
				if ( j >= numEnzymeFragments ) break;
				mol_wt += enzyme_fragment_mass_array [j];
				if ( mol_wt > endLimit ) break;
				if ( mol_wt > startLimit ) {
					MSTagSearchAllowErrors::resetNextMods ();
					int ind = lower_bound ( msMSSearch.begin (), msMSSearch.end (), mol_wt, sortMSMSSearchByParentMass () ) - msMSSearch.begin ();
					first = true;
					for ( ; ind < msMSSearch.size () ; ind++ ) {
						MSMSSearch* msms = msMSSearch [ind];
						if ( msms->getSpectrumRetained () ) {
							if ( msms->getParentMassMinusPosTolerance () > mol_wt ) break;
							if ( msms->checkMatch ( mol_wt ) ) {
								if ( first ) {
									hitLength = cleavage_index [j] - previousCleavageIndex;
									possMatch = string ( frame + previousCleavageIndex + 1, hitLength );
									first = false;
								}
								if ( !compMask || checkComposition ( possMatch ) ) {
									bool cTermProt = ( j == numEnzymeFragments-1 );
									if ( msms->doMatch ( possMatch, nTermProt, cTermProt, mol_wt, tagMatch, tagHits->getMinScore ( ind ) ) ) {
										addTagHit ( ind, tagMatch, possMatch, fi, nTermProt ? '-' : frame[previousCleavageIndex], cTermProt ? '-' : frame [cleavage_index[j] + 1], previousCleavageIndex + 2 );
									}
								}
							}
						}
					}
				}
			}
			missedCleavageLimit++;
		}
	}
}
const ScoreType TagHit::SCORE_DIFF_NOT_CALCULATED = -999999;
bool TagHit::eValFlag = false;
TagHit::TagHit ( const string& sequence, char previousAA, char nextAA, int startAA, const TagMatch& tagMatch ) :
	sequence ( sequence ),
	previousAA ( previousAA ),
	nextAA ( nextAA ),
	startAA ( startAA ),
	unmatched ( tagMatch.getUnmatched () ),
	hitSequence ( tagMatch.getMutatedSequence () ),
	peakMatch ( getPeakMatch ( tagMatch.getParentPeak () ) ),
	score ( tagMatch.getScore () ),
	scoreDiff ( SCORE_DIFF_NOT_CALCULATED ),
	rank ( -1 )
{
}
TagHit::TagHit ( const FrameIterator& fi, const string& sequence, char previousAA, char nextAA, int startAA, const TagMatch& tagMatch ) :
	ProteinHit ( fi.getFS (), fi.getEntry (), fi.getFrameTranslation (), fi.getFrame () ),
	sequence ( sequence ),
	previousAA ( previousAA ),
	nextAA ( nextAA ),
	startAA ( startAA ),
	unmatched ( tagMatch.getUnmatched () ),
	hitSequence ( tagMatch.getMutatedSequence () ),
	peakMatch ( getPeakMatch ( tagMatch.getParentPeak () ) ),
	score ( tagMatch.getScore () ),
	scoreDiff ( SCORE_DIFF_NOT_CALCULATED )
{
}
PeakMatch TagHit::getPeakMatch ( Peak* parentPeak )
{
	double peptideMW = hitSequence.getMW ();
	if ( parentPeak->getAverageFlag () ) {
		peptideMW *= MONOISOTOPIC_TO_AVERAGE_FACTOR;
		parentPeak->setMassAverage ();
	}
	return PeakMatch ( parentPeak, peptideMW );
}
void TagHit::printHTMLProteinHeader ( ostream& os ) const
{
	ProteinHit::printHTMLHeader2 ( os );
}
void TagHit::printHTMLHeader ( ostream& os, const string& parentToleranceUnits ) const
{
	tableRowStart ( os );
		tableHeader ( os, "Rank" );
		tableHeader ( os, "#<br />Unmatched<br />Ions" );
		tableHeader ( os, "Sequence" );
		tableHeader ( os, "Score" );
		if ( eValFlag ) tableHeader ( os, "Expect" );
		tableHeaderStart ( os );
			os << mh_plus_html << "<br />Calculated<br />(Da)" << endl;
		tableHeaderEnd ( os );
		tableHeaderStart ( os );
			os << "Error<br />(" << parentToleranceUnits << ")" << endl;
		tableHeaderEnd ( os );
		ProteinHit::printHTMLHeader2 ( os );
	tableRowEnd ( os );
}
void TagHit::printDelimited ( ostream& os, const PeakMatchContext& peakMatchContext ) const
{
	delimitedRowStart ( os );
		if ( rank == -1 ) {
			delimitedCell ( os, "M" );
			delimitedCell ( os, unmatched );
			hitSequence.printDelimited ( os, sequence, peakMatch.getMassDiff () );
			delimitedCell ( os, getScoreOutput (), 1 );
		}
		else {
			delimitedCell ( os, rank );
			delimitedCell ( os, unmatched );
			if ( hitSequence.getMmod () == 1 )
				delimitedCell ( os, 0.0, 1 );
			else
				delimitedCellSigFig ( os, peakMatch.getError ( peakMatchContext ), 5 );
			hitSequence.printDelimited ( os, sequence, peakMatch.getMassDiff () );
			delimitedCell ( os, startAA );
			delimitedCell ( os, getScoreOutput (), 1 );
			if ( scoreDiff != SCORE_DIFF_NOT_CALCULATED )
				delimitedCell ( os, getScoreDiffOutput (), 1 );
			else
				delimitedCell ( os, "---" );
			printDelimitedAccNum ( os );
		}
	delimitedRowEnd ( os );
}
void TagHit::printDelimitedContinue ( ostream& os ) const
{
	if ( rank == -1 ) return;
	delimitedRowStart ( os );
		delimitedCell ( os, string ( "+" ) + gen_itoa ( startAA ) );
		printDelimitedAccNum ( os );
	delimitedRowEnd ( os );
}
void TagHit::printHTML ( ostream& os, const PeakMatchContext& peakMatchContext, const MSProductLink& productLink, double eValue ) const
{
	if ( rank == -1 ) return;
	tableRowStart ( os );
		tableCell ( os, rank );
		tableHeaderStart ( os, "", "center" );
			os << unmatched << endl;
		tableHeaderEnd ( os );

		printHTMLSequence ( os, productLink );

		tableHeaderStart ( os, "", "center" );
			genPrint ( os, getScoreOutput (), 1 );
			os << endl;
		tableHeaderEnd ( os );
		if ( eValFlag ) {
			tableHeaderStart ( os, "", "center" );
				genPrintSigFig ( os, eValue, 2 );
				os << endl;
			tableHeaderEnd ( os );
		}
		peakMatch.printTagHitHTML ( os, peakMatchContext, hitSequence.getMmod () == 1 );
		printHTMLHit ( os );
	tableRowEnd ( os );
}
void TagHit::printHTMLSequence ( ostream& os, const MSProductLink& productLink ) const
{
	tableHeaderStart ( os, "", "", true );
		os << "(" << previousAA << ")";
		productLink.write2 ( os, sequence, hitSequence, false, peakMatch.getCharge (), peakMatch.getMassDiff () );
		os << "(" << nextAA << ")";
		os << endl;
	tableHeaderEnd ( os );
}
TagHitsContainer::TagHitsContainer ( const MSMSSearch* msMSSearch, MSTagParameters& params ) :
	msMSSearch ( msMSSearch ),
	parentPeak ( msMSSearch->getParentPeak () ),
	peaks ( msMSSearch->getPeaks () ),
	params ( params ),
	minScore ( -999999 )
{
}
void TagHitsContainer::checkRegularExpression ( RegularExpression& regExp )
{
	int count = 0;
	for ( TagHitVectorSizeType i = 0 ; i < tHits.size () ; i++ ) {
		if ( regExp.isPresent ( tHits [i].getHitSequence ().getSequence ().c_str () ) ) {
			tHits [count++] = tHits [i];
		}
	}
	tHits.erase ( tHits.begin () + count, tHits.end () );
}
void TagHitsContainer::pruneHits ()
{
	sortHits ();
	numHits = tHits.size ();
	numSavedHits = getNumHitsToSave ( genMin( numHits, params.getNumSavedHits () ) );
	numRetainedHits = genMin( numHits, numSavedHits + 1 );	// Retain one extra hit for score difference calculation
	tHits.erase ( tHits.begin () + numRetainedHits, tHits.end () );
	if ( numRetainedHits == numSavedHits + 1 ) minScore = tHits.back ().getScore ();
}
int TagHitsContainer::getNumHitsToSave ( int numToSave )
{
	int a;
	int maxNumSavedHits = numToSave * 10;		// Never save more than 10 times the number of saved hits
	int numUnique = 0;
	for ( a = 0 ; a < numHits ; a++ ) {
		if ( a > maxNumSavedHits ) break;
		if ( a == 0 || !matrixMatch ( tHits [a].getSequence (), tHits [a-1].getSequence () ) ) {
			if ( numUnique == numToSave ) break;
			numUnique++;
		}
	}
	return a;
}
void TagHitsContainer::rankHits ()
{
	pruneHits ();
	bool scoreHistogramOnly = params.getScoreHistogramOnly ();
	bool crosslinking = params.isCrosslinking ();
	int maxReportedHits = params.getMaxReportedHits ();
	if ( scoreHistogramOnly ) {
		for ( int x = 0 ; x < numRetainedHits-1 ; x++ ) {
			if ( x == 0 || !matrixMatch ( tHits [x].getSequence (), tHits [x-1].getSequence () ) ) {
				double s = tHits [x].getScore () / SCORE_TYPE_MULTIPLIER;
				hist.add ( s );
			}
		}
	}
	else {
		if ( crosslinking ) calculateCrossLinks ();
		numSavedHits = getNumHitsToSave ( genMin( numHits, maxReportedHits ) );
		numRetainedHits = genMin( numHits, numSavedHits + 1 );	// Retain one extra hit for score difference calculation
		tHits.erase ( tHits.begin () + numRetainedHits, tHits.end () );
		for ( int x = 0 ; x < numRetainedHits-1 ; x++ ) {
			ScoreType scoreDiff = tHits [x].getScore () - tHits [numRetainedHits-1].getScore ();
			tHits [x].setScoreDiff ( scoreDiff );
		}
	}
	for ( int i = 0, rank = 1 ; i < numSavedHits ; i++ ) {
		if ( i != 0 ) {
			if ( tHits [i].getScore () < tHits [i-1].getScore () ) rank++;
		}
		tHits [i].rank = rank;
	}
	tHits.erase ( tHits.begin () + numSavedHits, tHits.end () );	// Get rid of extra retained hit
	if ( !crosslinking && !scoreHistogramOnly && maxReportedHits <= 5 ) calculateModScores ();
}
void TagHitsContainer::push_back ( const T& hit )
{
	static int pruneLimit = params.isCrosslinking () ? params.getNumSavedCrosslinkHits () * 10 : 1000;
	static int pruneLevel = genMax( pruneLimit, params.getNumSavedHits ()+1 );
	tHits.push_back ( hit );
	if ( tHits.size () > pruneLevel ) pruneHits ();
}
void TagHitsContainer::printHTMLReport ( ostream& os, int dataSet )
{
	printNumHits ( os, "MS-Tag", numHits, params.getMaxReportedHits () );

	parentPeak->printHTML ( os, "Parent mass", instInf->getParentPeakPrecision () );

	if ( params.getMonoParentAverageFragments () ) {
		modified_mass_convert ( false );
		peaks->setMassAverage ();
	}
	peaks->printPeaksHTML ( os, instInf->getFragmentPeakPrecision () );

	if ( msMSSearch->getCompositionSearch () ) msMSSearch->getCompositionSearch ()->printHTML ( os );

	MSMSDataSetInfo* dataSetInfo = params.getDataSetInfo ();
	startJavascript ( os );
	MSProductLink productLink ( "mstag", dataSetInfo, dataSet );
	productLink.printHTML ( os );
	endJavascript ( os );

	if ( numHits < 1 )
		os << "\n <p /><b>The data could not be matched to any peptide sequences. You might try decreasing the specificity of the search parameters.</b><p />";
	else
		printHTML ( os, productLink );
}
void TagHitsContainer::printDelimited ( ostream& os )
{
	PeakMatchContext pmc ( params.getParentMassTolerance (), instInf->getParentPeakPrecision (), parentPeak->getNonUnitChargeData () );
	for ( int i = 0 ; i < tHits.size () ; i++ ) {
		if ( i != 0 && tHits [i].getHitSequence () == tHits [i-1].getHitSequence () )
			tHits [i].printDelimitedContinue ( os );
		else
			tHits [i].printDelimited ( os, pmc );
	}
	if ( !combinedTagHits.empty () ) {
		int num = 1;
		for ( CombinedTagHitSetConstIterator i = combinedTagHits.begin () ; i != combinedTagHits.end () ; i++ ) {
			if ( num > params.getMaxReportedHits () ) break;
			(*i)->printDelimited ( os, num++ );
		}
	}
}
void TagHitsContainer::printHTML ( ostream& os, const MSProductLink& productLink ) const
{
	TagHit::setEValFlag ( static_cast <const MSTagSearch*> (msMSSearch)->getEValueFlag () );
	printCrossLinksHTML ( os, productLink );
	string parentToleranceUnits = params.getParentMassTolerance ()->getUnitsString ();
	os << "<table cellspacing=\"3\">" << endl;
	PeakMatchContext peakMatchContext ( params.getParentMassTolerance (), instInf->getParentPeakPrecision (), parentPeak->getNonUnitChargeData () );
	for ( int i = 0 ; i < tHits.size () ; i++ ) {
		const TagHit& hit = tHits [i];
		if ( i == 0 ) hit.printHTMLHeader ( os, parentToleranceUnits );
		hit.printHTML ( os, peakMatchContext, productLink, static_cast <const MSTagSearch*> (msMSSearch)->getEvalue ( hit.getScore () ) );
	}
	os << "</table>" << endl;
}
class sortTagMatchByScore {
public:
	int operator () ( const TagMatch& a, const TagMatch& b ) const
	{
		return ( a.getScore () > b.getScore () );
	}
};
void TagHitsContainer::calculateModScores ()
{
	if ( numSavedHits == 0 ) return;
	SetString ss;
	ScoreType minScore = tHits [numSavedHits-1].getScore ();
	for ( int i = 0 ; i < numSavedHits ; i++ ) {
		const TagHit& hit = tHits [i];
		const string& seq = hit.getSequence ();
		std::pair <SetStringIterator, bool> flag = ss.insert ( seq );
		if ( flag.second ) {							// If this is a new sequence
			SetPairIntString si = hit.getModIndicies ();// Set of residue/mod
			if ( !si.empty () ) {
				TagMatchVector tagMatch;
				double mol_wt = peptide_formula_to_molecular_weight ( seq.c_str () );
				char previousAA = hit.getPreviousAA ();
				char nextAA = hit.getNextAA ();
				int startAA = hit.getStartAA ();
				bool nTermFlag = previousAA == '-';
				bool cTermFlag = nextAA == '-';
				const_cast <MSTagSearchAllowErrors*> (static_cast <const MSTagSearchAllowErrors*> (msMSSearch))->doMatch ( seq, nTermFlag, cTermFlag, mol_wt, tagMatch, -9999999 );	// Re-search this peptide
				stable_sort ( tagMatch.begin (), tagMatch.end (), sortTagMatchByScore () );
				for ( int j = 0 ; j < tagMatch.size () ; j++ ) {
					SetPairIntString si2 = tagMatch [j].getModIndicies ();
					SetPairIntString si3;
					set_intersection ( si.begin (), si.end (), si2.begin (), si2.end (), inserter ( si3, si3.begin () ) );
					if ( si3.size () < si.size () ) {
						si = si3;
						if ( tagMatch [j].getScore () <= minScore ) {	// Add any new hits - these should represent other mod positions
							tHits.push_back ( TagHit ( seq, previousAA, nextAA, startAA, tagMatch [j] ) );
						}
					}
				}
			}
		}
	}
}
typedef std::map <string, CombinedTagHit*> CombinedTagHitMap;
typedef CombinedTagHitMap::const_iterator CombinedTagHitMapConstIterator;
void TagHitsContainer::calculateCrossLinks ()
{
	Tolerance* cTol = params.getParentMassTolerance ();
	PeakMatchContext pmc ( cTol, instInf->getParentPeakPrecision (), parentPeak->getNonUnitChargeData () );
	CombinedTagHitMap repeats;
	int maxRepHits = params.getMaxReportedHits ();
	std::priority_queue<ScoreType, vector <ScoreType>, std::greater <ScoreType> > bestScores;
	for ( int a = 0 ; a < 1 ; a++ ) {
		double linkMass = massConvert ( params.getBridgeFormula ().c_str () );
		unsigned int linkAA1 = params.getLinkAminoAcid ( 1 );
		unsigned int linkAA2 = params.getLinkAminoAcid ( 2 );
		for ( int i = 0 ; i < numSavedHits ; i++ ) {
			const TagHit& hit1 = tHits [i];
			double accurateMass = linkMass;
			double nLossMass1 = hit1.getNLossMass ();
			accurateMass += hit1.getMassModRemainder ();
			accurateMass -= cation_wt;
			accurateMass -= nLossMass1;
			for ( int j = i ; j < numSavedHits ; j++ ) {
				const TagHit& hit2 = tHits [j];
				double nLossMass2 = hit2.getNLossMass ();
				if ( nLossMass1 != nLossMass2 ) continue;
				double am = accurateMass + hit2.getMassModRemainder ();
				if ( genAbsDiff ( am, parentPeak->getMass () ) < cTol->getTolerance ( parentPeak->getMOverZ (), parentPeak->getCharge () ) ) {
					unsigned int mask1 = hit1.getMassModAAMask ();
					unsigned int mask2 = hit2.getMassModAAMask ();
					if ( (( mask1 & linkAA1 ) && ( mask2 & linkAA2 )) || (( mask2 & linkAA1 ) && ( mask1 & linkAA2 )) ) {
						double error = pmc.getError ( parentPeak->getMass (), am, parentPeak->getCharge () );
						const string& a1 = hit1.getAccessionNumber ();
						int startAA1 = hit1.getStartAA ();
						const string& a2 = hit2.getAccessionNumber ();
						int startAA2 = hit2.getStartAA ();
						string totalSeq = hit1.getHitSequenceString () + " " + hit2.getHitSequenceString ();
						CombinedTagHitMapConstIterator cur = repeats.find ( totalSeq );

						if ( cur == repeats.end () ) {
							CombinedTagHit* cth = new CombinedTagHit ( hit1, i+1, hit2, j+1, linkMass, error, msMSSearch );
							ScoreType sc = cth->getScore ();
							if ( bestScores.size () < maxRepHits ) {
								bestScores.push ( sc );
								repeats [totalSeq] = cth;
								combinedTagHits.insert ( cth );
							}
							else if ( sc > bestScores.top () ) {
								bestScores.push ( sc );
								bestScores.pop ();
								repeats [totalSeq] = cth;
								combinedTagHits.insert ( cth );
							}
							else {
								repeats [totalSeq] = 0;
								delete cth;
							}
						}
						else {
							CombinedTagHit* cth = (*cur).second;
							if ( cth ) cth->add ( hit1, hit2 );
						}
					}
				}
			}
		}
	}
}
void TagHitsContainer::printCrossLinksHTML ( ostream& os, const MSProductLink& productLink ) const
{
	if ( !combinedTagHits.empty () ) {
		os << "<div class=\"results_section_header\">Crosslinked Peptide Hits</div>" << endl;
		os << "<table cellspacing=\"3\">" << endl;
		Tolerance* cTol = params.getParentMassTolerance ();
		PeakMatchContext pmc ( cTol, instInf->getParentPeakPrecision (), parentPeak->getNonUnitChargeData () );
		int num = 1;
		for ( CombinedTagHitSetConstIterator i = combinedTagHits.begin () ; i != combinedTagHits.end () ; i++ ) {
			if ( num > params.getMaxReportedHits () ) break;
			if ( num == 1 ) (*i)->printHTMLHeader ( os, params.getParentMassTolerance ()->getUnitsString () );
			(*i)->printHTML ( os, num++, pmc, productLink, params.getLinkInfo () );
		}
		os << "</table>" << endl;
	}
}
CombinedTagHit::CombinedTagHit ( const TagHit& h1, int rank1, const TagHit& h2, int rank2, double linkMass, double error, const MSMSSearch* msMSSearch ) :
	hit1 ( h1 ),
	hit2 ( h2 ),
	error ( error )
{
	hit1.rank = rank1;
	hit2.rank = rank2;
	double mm1 = linkMass + hit2.getMassModRemainder () - hit1.getNLossMass ();	// Must be calculated before being set
	double mm2 = linkMass + hit1.getMassModRemainder () - hit2.getNLossMass ();
	mm1 -= cation_wt;
	mm2 -= cation_wt;
	hit1.setMassMod ( mm1 );
	hit2.setMassMod ( mm2 );

	PeptideSequence ps1 = ( rank1 < rank2 ) ? hit1.getHitSequence () : hit2.getHitSequence ();
	MassType nTerminusWt = static_cast <MassType> (n_terminus_wt * MASS_TYPE_MULTIPLIER);
	MassType cTerminusWt = static_cast <MassType> (c_terminus_wt * MASS_TYPE_MULTIPLIER);
	ps1.applyModifications ( nTerminusWt, cTerminusWt );
	MSTagSearch::setNTerminusWt ( nTerminusWt );
	MSTagSearch::setCTerminusWt ( cTerminusWt );
	numUnmatchedIons = 0;
	score = 0;
	const_cast <MSTagSearchAllowErrors*> (static_cast <const MSTagSearchAllowErrors*> (msMSSearch))->fragmentMatch2 ( ps1, numUnmatchedIons, score, true );

	firstPeptideScore = score;
	PeptideSequence ps2 = ( rank1 < rank2 ) ? hit2.getHitSequence () : hit1.getHitSequence ();
	nTerminusWt = static_cast <MassType> (n_terminus_wt * MASS_TYPE_MULTIPLIER);
	cTerminusWt = static_cast <MassType> (c_terminus_wt * MASS_TYPE_MULTIPLIER);
	ps2.applyModifications ( nTerminusWt, cTerminusWt );
	MSTagSearch::setNTerminusWt ( nTerminusWt );
	MSTagSearch::setCTerminusWt ( cTerminusWt );
	numUnmatchedIons = 0;
	score = 0;
	const_cast <MSTagSearchAllowErrors*> (static_cast <const MSTagSearchAllowErrors*> (msMSSearch))->fragmentMatch2 ( ps2, numUnmatchedIons, score, false );
	schani.insert ( CombinedHitAccNumInfo ( h1, h2 ) );
}
void CombinedTagHit::printHTMLHeader ( ostream& os, const string& parentToleranceUnits ) const
{
	tableRowStart ( os );
		tableHeader ( os, "Rank" );
		tableHeader ( os, "#<br />Unmatched<br />Ions" );
		tableHeader ( os, "Sequence" );
		tableHeader ( os, "Score" );
		tableHeaderStart ( os );
			os << "Error<br />(" << parentToleranceUnits << ")" << endl;
		tableHeaderEnd ( os );
		tableHeader ( os, "MS-Tag<br />Score" );
		tableHeader ( os, "Rank" );
		tableHeader ( os, "Low<br />Score" );
		tableHeader ( os, "XLink<br />AA" );
		hit1.printHTMLProteinHeader ( os );
	tableRowEnd ( os );
}
void CombinedTagHit::printHTML ( ostream& os, int rank, const PeakMatchContext& pmc, const MSProductLink& productLink, const LinkInfo* linkInfo ) const
{
	tableEmptyRow ( os );
	for ( std::set <CombinedHitAccNumInfo, SortCombinedTagAscending2>::const_iterator i = schani.begin () ; i != schani.end () ; i++ ) {
		tableRowStart ( os );
			if ( i == schani.begin () ) {
				tableCell ( os, rank );
				tableHeaderStart ( os, "", "center" );
					os << numUnmatchedIons << endl;
				tableHeaderEnd ( os );
				os << "<th nowrap=\"nowrap\" rowspan=\"2\">" << endl;
					printHTMLSequence ( os, productLink, linkInfo );
				os << "</th>" << endl;
				tableCell ( os, (double) score / (double) SCORE_TYPE_MULTIPLIER, 1, true );
				tableCellSigFig ( os, error, pmc.getErrorSigFig (), true );
				tableCell ( os, (double) hit1.getScore () / (double) SCORE_TYPE_MULTIPLIER, 1, true );
				tableCell ( os, hit1.rank, true );
				tableCell ( os, (double) (score - getFirstPeptideScore ()) / (double) SCORE_TYPE_MULTIPLIER, 1, true );
			}
			else {
				tableEmptyNCells ( os, 8 );
			}
			(*i).printHTMLHit1 ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableEmptyNCells ( os, 4 );
			if ( i == schani.begin () ) {
				tableCell ( os, (double) hit2.getScore () / (double) SCORE_TYPE_MULTIPLIER, 1, true );
				tableCell ( os, hit2.rank, true );
				tableEmptyNCells ( os, 1 );
			}
			else {
				tableEmptyNCells ( os, 4 );
			}
			(*i).printHTMLHit2 ( os );
		tableEmptyRow ( os );
		tableRowEnd ( os );
	}
}
void CombinedTagHit::printHTMLSequence ( ostream& os, const MSProductLink& productLink, const LinkInfo* linkInfo ) const
{
	CharVector previousAA;
	CharVector nextAA;
	vector <PeptideSequence> hitSequence;
	StringVector sequence;
	previousAA.push_back ( hit1.getPreviousAA () );
	previousAA.push_back ( hit2.getPreviousAA () );
	nextAA.push_back ( hit1.getNextAA () );
	nextAA.push_back ( hit2.getNextAA () );
	hitSequence.push_back ( hit1.getHitSequence () );
	hitSequence.push_back ( hit2.getHitSequence () );
	sequence.push_back ( hit1.getSequence () );
	sequence.push_back ( hit2.getSequence () );
	productLink.write2 ( os, previousAA, sequence, nextAA, hitSequence, false, hit1.getCharge (), linkInfo );
	os << endl;
}
void CombinedTagHit::printDelimited ( ostream& os, int rank ) const
{
	for ( std::set <CombinedHitAccNumInfo, SortCombinedTagAscending2>::const_iterator i = schani.begin () ; i != schani.end () ; i++ ) {
		delimitedRowStart ( os );
		if ( i == schani.begin () ) {
				delimitedCell ( os, "X" + gen_itoa ( rank ) );
				delimitedCell ( os, numUnmatchedIons );
				delimitedCellSigFig ( os, error, 5 );

				hit1.getHitSequence ().printDelimited ( os, hit1.getSequence (), 0.0 );
				(*i).printDelimitedHit1 ( os );

				hit2.getHitSequence ().printDelimited ( os, hit2.getSequence (), 0.0 );
				(*i).printDelimitedHit2 ( os );

				delimitedCell ( os, (double) score / (double) SCORE_TYPE_MULTIPLIER, 1 );
		// Extra info
				delimitedCell ( os, (double) hit1.getScore () / (double) SCORE_TYPE_MULTIPLIER, 1 );
				delimitedCell ( os, hit1.rank );
				delimitedCell ( os, (double) hit2.getScore () / (double) SCORE_TYPE_MULTIPLIER, 1 );
				delimitedCell ( os, hit2.rank );
				delimitedCell ( os, (double) getFirstPeptideScore () / (double) SCORE_TYPE_MULTIPLIER, 1 );
				delimitedCell ( os, 0 );		// This field contained the number of crosslink ions. After v5.14.1 is just set to zero.
		//
		}
		else {
			(*i).printDelimitedHitContinue ( os );
		}
		delimitedRowEnd ( os );
	}
}
void CombinedTagHit::add ( const TagHit& h1, const TagHit& h2 )
{
	schani.insert ( CombinedHitAccNumInfo ( h1, h2 ) );
}
CombinedHitAccNumInfo::CombinedHitAccNumInfo ( const TagHit& hit1, const TagHit& hit2 ) :
	hit1 ( hit1 ),
	hit2 ( hit2 )
{
}
void CombinedHitAccNumInfo::printHTMLHit1 ( ostream& os ) const
{
	tableCell ( os, hit1.getMmodProteinIdx (), true );
	hit1.printHTMLHit ( os );
}
void CombinedHitAccNumInfo::printHTMLHit2 ( ostream& os ) const
{
	tableCell ( os, hit2.getMmodProteinIdx (), true );
	hit2.printHTMLHit ( os );
}
void CombinedHitAccNumInfo::printDelimitedHit1 ( ostream& os ) const
{
	delimitedCell ( os, hit1.getStartAA () );
	hit1.printDelimitedAccNum ( os );
}
void CombinedHitAccNumInfo::printDelimitedHit2 ( ostream& os ) const
{
	delimitedCell ( os, hit2.getStartAA () );
	hit2.printDelimitedAccNum ( os );
}
void CombinedHitAccNumInfo::printDelimitedHitContinue ( ostream& os ) const
{
	delimitedCell ( os, string ( "+" ) + gen_itoa ( hit1.getStartAA () ) );
	hit1.printDelimitedAccNum ( os );
	delimitedCell ( os, hit2.getStartAA () );
	hit2.printDelimitedAccNum ( os );
}
