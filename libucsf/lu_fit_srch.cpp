/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_fit_srch.cpp                                               *
*                                                                             *
*  Created    : July 13th 2001                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lu_acc_link.h>
#include <lp_frame.h>
#include <lu_html.h>
#include <lu_fit_srch.h>
#include <lu_mat_score.h>
#include <lu_spep_srch.h>
#include <lu_aa_calc.h>
#include <lu_frag_mtch.h>
#include <lu_param_list.h>
#include <lu_table.h>
using std::ostream;
using std::string;
using std::vector;
using std::endl;
using std::sort;
using std::count;
using std::cout;
using std::make_pair;

class sort_hits_by_pi {
	const vector <FitHit>& hits;
public:
	sort_hits_by_pi ( const vector <FitHit>& hits ) :
		hits ( hits )
	{
	}
	int operator () ( const PairIntIntVector& a, const PairIntIntVector& b ) const
	{
		return ( hits [a.first].getProteinPI () < hits [b.first].getProteinPI () );
	}
};

class sort_hits_by_protein_mw {
	const vector <FitHit>& hits;
public:
	sort_hits_by_protein_mw ( const vector <FitHit>& hits ) :
		hits ( hits )
	{
	}
	int operator () ( const PairIntIntVector& a, const PairIntIntVector& b ) const
	{
		return ( hits [a.first].getProteinMW () < hits [b.first].getProteinMW () );
	}
};

class sort_fit_hits {
public:
	int operator () ( const FitHit& a, const FitHit& b ) const
	{
		if ( a.numMatches == b.numMatches )
			return a.getIndex () < b.getIndex ();
		else
			return b.numMatches < a.numMatches;
	}
};

class sort_mowse_hits {
public:
	int operator () ( const FitHit& a, const FitHit& b ) const
	{
		if ( a.hitStats.mowseScore == b.hitStats.mowseScore ) {
			if ( a.numMatches == b.numMatches )
				return a.getIndex () < b.getIndex ();
			else
				return b.numMatches < a.numMatches;
		}
		else
			return b.hitStats.mowseScore < a.hitStats.mowseScore;
	}
};

FitHits::FitHits ( const vector <MSFitSearch*>& msFitSearch, MSFitParameters& params ) :
	DatabaseHits (),
	numSearches ( msFitSearch.size () )
{
	for ( int i = 0 ; i < numSearches ; i++ ) {
		hits.push_back ( new FitHitsContainer ( msFitSearch [i], params, params.getDataSetInfo ()->getSpotID ( i ), i ) );
	}
}
void FitHits::sortAndRank ()
{
	for ( int i = 0 ; i < numSearches ; i++ ) {
		hits [i]->sortHits ();
		hits [i]->eraseHits ();
	}
}
FitSearch::FitSearch ( MSFitParameters& params ) :
	DatabaseSearch ( params ),
	fitParams ( params ),
	maxHits ( params.getMaxHits () )
{
	MSDataSetInfo* dsi = params.getDataSetInfo ();
	for ( int i = 0 ; i < dsi->getNumDataSets () ; i++ ) {
		fitSearch.push_back ( getMSFitSearch ( dsi->getDataSet ( i ), params ) );
	}
	fitHits = new FitHits ( fitSearch, params );
	numSearches = fitSearch.size ();
	doSearch ();
	UpdatingJavascriptMessage::deleteOutstandingSpanID ( cout );
	UpdatingJavascriptMessage ujm;
	ujm.writeMessage ( cout, "Calculating detailed results please wait. This can take some time." );
	for ( int m = 0 ; m < numSearches ; m++ ) {
		FitHitsContainer* fhc = fitHits->getHit ( m );
		(*fhc).setSingleFitSearch ();
		(*fhc).calculateMainAndSupplementaryHits ();
		(*fhc).sortHitsByPIorMW ();
	}
	numHits = fitHits->size ();
	databaseHits = fitHits;
	ujm.deletePreviousMessage ( cout );
}
FitSearch::~FitSearch ()
{
	delete fitHits;
	for ( int i = 0 ; i < fitSearch.size () ; i++ ) {
		delete fitSearch [i];
	}
}
void FitSearch::doSearch ()
{
	init_fasta_enzyme_function ( params.getEnzyme () );
	for ( int i = 0 ; i < fs.size () ; i++ ) {
		doSearch ( fs [i], i );
		FrameIterator::resetElapsedTime ( 1 );
	}
	if ( fitParams.getMowseFlag () ) calculateMowseScores ();
	fitHits->sortAndRank ();
}
void FitSearch::doSearch ( FastaServer* fsPtr, int num )
{
	int minMatches = fitParams.getMinMatches ();
	ProteinHit::addFS ( fsPtr, num );
	FrameIterator fi ( fsPtr, params.getIndicies ( num ), dnaFrameTranslationPairVector [num], params.getTempOverride () );
	char* frame;
	while ( ( frame = fi.getNextFrame () ) != NULL ) {
		const IntVector& cleavageIndex = enzyme_fragmenter ( frame );
		for ( MSFitSearchPtrVectorSizeType i = 0 ; i < fitSearch.size () ; i++ ) {
			int numMatches = fitSearch [i]->matchFragments ( frame, cleavageIndex );

			if ( numMatches >= minMatches ) {
				fitHits->addHit ( FitHit ( fsPtr, fi.getEntry (), fi.getFrameTranslation (), fi.getFrame (), numMatches, fitSearch [i]->getProteinFragmentStats () ), i );
			}
		}
	}
}
void FitSearch::calculateMowseScores ()
{
	for ( MSFitSearchPtrVectorSizeType i = 0 ; i < fitSearch.size () ; i++ ) {
		FitHitsContainer* fhc = fitHits->getHit ( i );
		for ( int j = 0 ; j < fhc->size () ; j++ ) {
			FitHit& fh = (*fhc)[j];
			fitSearch [i]->calculateMowseScores ( &(fh.hitStats), fh.getProteinMW () );
		}
	}
}
void FitSearch::printHTMLHitsJavascript ( ostream& os ) const
{
	startJavascript ( os );
	StringVector dbPaths;
	for ( int i = 0 ; i < fs.size () ; i++ ) {
		string fPath = fs [i]->getFilePath ();
		if ( fPath.empty () ) fPath = fs [i]->getFileName ();
		dbPaths.push_back ( fPath );
		AccessionNumberLinkInfo anli;
		anli.printHTML ( os, fs [i]->getFileName () );
		MSDigestLink digestLink ( params.getSearchName (), dbPaths.back () );
		digestLink.printHTML ( os );
		MSNonSpecificLink nSpecLink ( fPath );
		nSpecLink.putHidden ( os );
		MSBridgeLink bLink ( fPath );
		bLink.putHidden ( os );
	}
	MSFitLink fLink;
	fLink.putHidden ( os );
	endJavascript ( os );
}
void FitSearch::printHTMLHitsReport ( ostream& os )
{
	for ( int i = 0 ; i < numSearches ; i++ ) {
//		os << "<div class=\"results_header\">Data Set " << i+1 << " Results</div>" << endl;
		fitHits->printHTMLReport ( os, i );
	}
}
void FitSearch::printXMLHits ( ostream& os ) const
{
	os << "<" << params.getSearchName () << "_hits>" << endl;
	for ( int i = 0 ; i < numSearches ; i++ ) {
		os << "<data_set>" << endl;
		fitHits->printXMLReport ( os, i );
		os << "</data_set>" << endl;
	}
	os << "</" << params.getSearchName () << "_hits>" << endl;
}
FitHit::FitHit ( FastaServer* fs, int ind, int drf, int orf, int numMatches, const HitStats& hitStats ) :
	ProteinHit ( fs, ind, drf, orf ),
	numMatches ( numMatches ),
	hitStats ( hitStats ),
	numHomology ( -1 )
{
}
void FitHit::setSingleFitSearch ( const MSFitParameters& params, const PeakContainer& peaks )
{
	char* frame = fs->getProtein ( getDatabaseEntry () );
	EnzymeFragmentContainer enzFrags ( frame, params.getEnzymeParameters () );
	vector <EnzymeFragmentContainer> vEnzFrags;
	vEnzFrags.push_back ( enzFrags );
	static AACalculator aaCalc ( params.getMonoisotopicFlag (), params.getAAInitInfo ().getConstMods () );
	if ( params.getModificationParameters ().getAllowErrors () ) {
		sfs = new ModificationSearch ( enzFrags, ModificationTable ( params.getModificationParameters () ), peaks );
	}
	else {
		LinkInfo* li = new LinkInfo;
		sfs = new LinksSearch ( vEnzFrags, peaks, li, aaCalc );
		delete li;
	}
	numUnique = sfs->getNumUnique ();
}
void FitHit::printHTMLHeader ( ostream& os, int numPeaks, bool scoreFlag, bool homFlag ) const
{
	tableRowStart ( os );
		if ( numPeaks > 0 ) {
			tableHeader ( os, "Protein<br />Hit<br />Number" );
			if ( scoreFlag ) tableHeader ( os, "MOWSE<br />Score" );
			tableHeaderStart ( os );
				os << "#&nbsp;pep";
				os << "<br />";
				os << "#&nbsp;mat";
				os << "<br />";
				os << "%&nbsp;mat";
				os << "<br />";
				os << numPeaks;
				os << "&nbsp;pks";
				os << endl;
			tableHeaderEnd ( os );
			sfs->printStatsHeaderHTML ( os );
			if ( numHomology != -1 ) {
				if ( homFlag ) tableHeader ( os, "# Hom<br />Pep" );
				else tableHeader ( os, "# Hom<br />Prot" );
			}
		}
		ProteinHit::printHTMLHeader ( os );
	tableRowEnd ( os );
}
void FitHit::printHTML ( ostream& os, int numPeaks, bool scoreFlag, bool homFlag ) const
{
	tableRowStart ( os );
	if ( numPeaks > 0 ) {
		tableCell ( os, hitNumber, true );
		if ( scoreFlag ) {
			tableHeaderStart ( os );
				if ( hitStats.mowseScore == 0.0 )
					os << "---";
				else
					genPrintSigFig ( os, hitStats.mowseScore , 3 );
			os << endl;
			tableHeaderEnd ( os );
		}
		tableHeaderStart ( os, "", "", true );
			os << numUnique;
			os << '/';
			os << numMatches;
			os << '/';
			genPrint ( os, getPercentMatches ( numPeaks ), 0 );
			os << endl;
		tableHeaderEnd ( os );
		sfs->printStatsHTML ( os, 0 );
		if ( numHomology != -1 ) {
			tableDataStart ( os, "", "center" );
				if ( homFlag ) os << numHomology << endl;
				else {
					if ( numHomology ) os << "<a href=\"#HM" << hitNumber << "\">" << numHomology << "</a>";
					else os << "No" << endl;
				}
			tableDataEnd ( os );
		}
	}
	ProteinHit::printHTMLHit2 ( os );
	tableRowEnd ( os );
}
void FitHit::printCoverageHTML ( ostream& os ) const
{
	sfs->printCoverageHTML ( os );
}
void FitHit::printXML ( ostream& os, int numPeaks, bool scoreFlag ) const
{
	if ( numPeaks > 0 ) {
		if ( scoreFlag ) {
			ParameterList::printDoubleXMLSigFig ( os, "mowse_score", hitStats.mowseScore, 3 );
		}
		ParameterList::printXML ( os, "num_peptides", numUnique );
		ParameterList::printXML ( os, "num_matches", numMatches );
		ParameterList::printDoubleXMLFixed ( os, "percent_matches", getPercentMatches ( numPeaks ), 0 );
	}
	ProteinHit::printXMLHit ( os );
}
void FitHit::printHTMLDetailSummary ( ostream& os, int i, int numPeaks ) const
{
	os << "<font size=\"+2\">";
	os << "<a name=\"" << i << "\"></a>";
	os << "</font>";
	os << " \n";
	os << "<b>" << i+1 << ".</b>";
	os << " \n";
	os << numMatches << "/" << numPeaks << " matches";
	os << " \n";
	os << "(" << 100 * numMatches / numPeaks << "%).";
	os << " \n";

	os << "<br />" << endl;
	ProteinHit::printHit ( os );
}
void FitHit::printHTMLDetail ( ostream& os, const MSFitParameters& params, const PeakContainer& peaks ) const
{
	std::map <const FastaServer*, int>::const_iterator cur = idxMap.find ( fs );

	sfs->printHTML ( os, false );
	const BoolDeque& peakUsed = sfs->getPeakUsed ();

	int numUnmatchedMasses = count ( peakUsed.begin (), peakUsed.end (), false );
	ParameterList::printHTML ( os, "Num Unmatched Masses", numUnmatchedMasses );

	static int idx = 1;
	string idxStr = gen_itoa ( idx );
	peaks.putHiddenFormJavascriptEntry ( os, "ms_parent_mass", peakUsed, idxStr );	// Write out the data

	MSBridgeLink::write ( os, getIndex (), getDNAReadingFrame (), getOpenReadingFrame (), params.getMultipleModification2 (), idxStr, (*cur).second );

	MSNonSpecificLink::write ( os, getIndex (), getDNAReadingFrame (), getOpenReadingFrame (), idxStr, (*cur).second );

	MSFitLink::write ( os, params.getMinMatches (), peaks, peakUsed, idxStr );
	idx++;

	os << "<br />" << endl;
	os << "<br />" << endl;

	sfs->printCoverCountHTML ( os, 0 );

	os << "Coverage Map for This Hit (MS-Digest index #): " << endl;
	MSDigestLink::write ( os, getIndex (), getDNAReadingFrame (), getOpenReadingFrame (), sfs->getCoverageMap ( 0 ), (*cur).second );
	os << "<p />" << endl;
}
void FitHit::printXMLDetail ( ostream& os, const MSFitParameters& params, const PeakContainer& peaks ) const
{
	os << "<protein_match>" << endl;
	printXML ( os, peaks.size (), params.getMowseFlag () );
	sfs->printXML ( os );
	os << "</protein_match>" << endl;

	//const BoolDeque& peakUsed = linksSearch.getPeakUsed ();

	//linksSearch.printCoverCountHTML ( os );

	//const BoolDeque& aaCovered = linksSearch.getAACovered ();

	//linksSearch.printStatsHTML ( os );
}
FitHitsContainer::FitHitsContainer ( const MSFitSearch* fitSearch, MSFitParameters& params, const string& spot, int searchNumber ) :
	fitSearch ( fitSearch ),
	peaks ( fitSearch->getPeaks () ),
	params ( params ),
	spot ( spot ),
	searchNumber ( searchNumber )
{
}
void FitHitsContainer::sortHits ()
{
	if ( params.getMowseFlag () ) sort ( hits.begin (), hits.end (), sort_mowse_hits () );
	else sort ( hits.begin (), hits.end (), sort_fit_hits () );
}
void FitHitsContainer::eraseHits ()
{
	numHits = hits.size ();
	hits.erase ( hits.begin () + genMin( numHits, params.getMaxReportedHits () ), hits.end () );
	for ( int i = 0 ; i < numHits ; i++ ) {
		hits [i].setHitNumber ( i+1 );
	}
}
void FitHitsContainer::sortHitsByPIorMW ()
{
	if ( params.getSortType () == "pI Sort" ) sort ( mainHits.begin (), mainHits.end (), sort_hits_by_pi ( hits ) );
	if ( params.getSortType () == "MW Sort" ) sort ( mainHits.begin (), mainHits.end (), sort_hits_by_protein_mw ( hits ) );
}
void FitHitsContainer::setSingleFitSearch ()
{
	for ( FitHitVectorSizeType i = 0 ; i < hits.size () ; i++ ) {
		hits [i].setSingleFitSearch ( params, peaks );
	}
}
void FitHitsContainer::calculateMainAndSupplementaryHits ()
{
	if ( params.getModificationParameters ().getAllowErrors () )
		calculateMainAndSupplementaryHits ( calculateNumHomologyMatchesMS );
	else
		calculateMainAndSupplementaryHits ( calculateNumHomologyMatchesLS );
}
void FitHitsContainer::calculateMainAndSupplementaryHits ( int (*numHomologyMassesCalculator) ( const FitHit&, const FitHit& ) )
{
	string reportHomologousProteins = params.getReportHomologousProteins ();
	set_score_matrix ( "BLOSUM62" );
	IntVector currentHomologyList;
	for ( FitHitVectorSizeType i = 0 ; i < hits.size () ; i++ ) {
		bool newHomologyMatch = true;
		if ( i != 0 ) {
			for ( IntVectorSizeType j = 0 ; j < currentHomologyList.size () ; j++ ) {
				int numMatches = numHomologyMassesCalculator ( hits [i], hits [currentHomologyList [j]] );
				if ( numMatches ) {
					newHomologyMatch = false;
					if ( reportHomologousProteins == "All" || ( reportHomologousProteins == "Interesting" && hits [i].getNumUnique () != numMatches ) ) {
						hits [i].setNumHomology ( numMatches );
						mainHits [j].second.push_back ( i );
					}
				}
			}
		}
		if ( newHomologyMatch ) {
			currentHomologyList.push_back ( i );
			mainHits.push_back ( make_pair ( i, IntVector () ) );
		}
	}
	for ( VectorPairIntIntVectorSizeType k = 0 ; k < mainHits.size () ; k++ ) {
		if ( reportHomologousProteins != "None" ) {
			hits [mainHits[k].first].setNumHomology ( mainHits[k].second.size () );
		}
	}
}
int FitHitsContainer::calculateNumHomologyMatchesMS ( const FitHit& hit, const FitHit& homHit )
{
	const ModificationSearch* ls = static_cast <ModificationSearch*> (hit.getSfs ());
	const ModificationSearch* ls2 = static_cast <ModificationSearch*> (homHit.getSfs ());
	return ls->calculateNumHomologyMatches ( ls2 );
}
int FitHitsContainer::calculateNumHomologyMatchesLS ( const FitHit& hit, const FitHit& homHit )
{
	const LinksSearch* ls = static_cast <LinksSearch*> (hit.getSfs ());
	const LinksSearch* ls2 = static_cast <LinksSearch*> (homHit.getSfs ());
	return ls->calculateNumHomologyMatches ( ls2 );
}
void FitHitsContainer::printHTMLReport ( ostream& os )
{
	os << "<hr />" << endl;
	ParameterList::printHTML ( os, "Fraction-Spot-Run ID", params.getDataSetInfo ()->getSpotRunID ( searchNumber ) );
	printNumHits ( os, "MS-Fit", numHits, params.getMaxReportedHits () );
	if ( mainHits.empty () ) {
		if ( params.getDetailedReport () ) {
			os << "MS-Fit could not fit the data to any proteins.";
			os << "<p />" << endl;
			os << "Try reducing the search specificity by:";
			os << "<ul>";
			os << "<li>decreasing Min. # peptides required to match</li>";
			os << "<li>widening the range for MW of Protein</li>";
			os << "<li>increasing the peptide-mass tolerance</li>";
			os << "</ul>";
			os << "<p />" << endl;
		}
		os << "<hr />" << endl;
	}
	else {
		printHTMLSummary ( os );

		if ( params.getDetailedReport () ) printHTMLDetail ( os );
	}
}
void FitHitsContainer::printXMLReport ( ostream& os )
{
	//printNumHits ( os, "MS-Fit", numHits, params.getMaxReportedHits () );
	printXMLDetail ( os );
}
void FitHitsContainer::printHTMLSummary ( ostream& os )
{
	int numDataSets = params.getDataSetInfo ()->getNumDataSets ();
	os << "<p>" << endl;
	ExpandableJavascriptBlock ejb ( "Results Summary", numDataSets == 1 );
	ejb.printHeader ( os );
	bool scoreFlag = params.getMowseFlag ();
	int numPeaks = peaks.size ();
	os << "<table>" << endl;
	for ( VectorPairIntIntVectorSizeType i = 0 ; i < mainHits.size () ; i++ ) {
		const FitHit& hit = hits [mainHits[i].first];
		if ( i == 0 ) hit.printHTMLHeader ( os, numPeaks, scoreFlag, false );
		hit.printHTML ( os, numPeaks, scoreFlag, false );
	}
	os << "</table>" << endl;
	printCoverageTable ( os, mainHits );

	for ( VectorPairIntIntVectorSizeType j = 0 ; j < mainHits.size () ; j++ ) {
		const IntVector& sh = mainHits [j].second;
		if ( sh.size () ) {
			const FitHit& hit = hits[mainHits [j].first];
			os << "<a name=\"HM" << hit.getHitNumber () << "\"></a>";
			os << "<div class=\"results_section_header\">Similar Matches for Protein Hit " << hit.getHitNumber () << "</div>" << endl;
			os << "<table>" << endl;
			hit.printHTMLHeader ( os, numPeaks, scoreFlag, false );
			hit.printHTML ( os, numPeaks, scoreFlag, false );
			for ( IntVectorSizeType k = 0 ; k < sh.size () ; k++ ) {
				const FitHit& hit2 = hits[sh[k]];
				if ( k == 0 ) hit2.printHTMLHeader ( os, numPeaks, scoreFlag, true );
				hit2.printHTML ( os, numPeaks, scoreFlag, true );
			}
			os << "</table>" << endl;
		}
	}
	os << "<hr />" << endl;
	ejb.printFooter ( os );
	os << "</p>" << endl;
}
void FitHitsContainer::printCoverageTable ( ostream& os, const VectorPairIntIntVector& hitList ) const
{
	os << "<table border=\"border\" cellspacing=\"1\">" << endl;
	for ( FitHitVectorSizeType i = 0 ; i < hitList.size () ; i++ ) {
		tableRowStart ( os );
			tableCell ( os, hits [hitList[i].first].getIndex (), true );
			hits [hitList[i].first].printCoverageHTML ( os );
		tableRowEnd ( os );
	}
	os << "</table>" << endl;
}
void FitHitsContainer::printHTMLDetail ( ostream& os )
{
	int numPeaks = peaks.size ();
	if ( numPeaks > 0 ) {
		os << "<p>" << endl;
		ExpandableJavascriptBlock ejb ( "Detailed Results" );
		ejb.printHeader ( os );

		for ( FitHitVectorSizeType i = 0 ; i < hits.size () ; i++ ) {
			hits [i].printHTMLDetailSummary ( os, i, numPeaks );
			hits [i].printHTMLDetail ( os, params, peaks );
		}
		ejb.printFooter ( os );
		os << "</p>" << endl;
	}
}
void FitHitsContainer::printXMLDetail ( ostream& os )
{
	int numPeaks = peaks.size ();
	if ( numPeaks > 0 ) {
		for ( FitHitVectorSizeType i = 0 ; i < hits.size () ; i++ ) {
			//hits [i].printXMLDetailSummary ( os, i, numPeaks );
			hits [i].printXMLDetail ( os, params, peaks );
		}
	}
}
