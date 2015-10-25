/******************************************************************************
*                                                                             *
*  Program    : searchCompare                                                 *
*                                                                             *
*  Filename   : sc_sres_rep.cpp                                               *
*                                                                             *
*  Created    : March 27th 2003                                               *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2003-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef VIS_C
#pragma warning( disable : 4503 )	// Disable warnings about decorated name length being too long
#endif
#ifndef VIS_C
#include <stdexcept>
#endif
#include <lgen_file.h>
#include <lgen_uncompress.h>
#include <lg_time.h>
#include <lu_aa_info.h>
#include <lu_app_gr.h>
#include <lu_blib.h>
#include <lu_html.h>
#include <lu_table.h>
#include <lu_delim.h>
#include <lu_getfil.h>
#include <lu_mat_score.h>
#include <lu_quan_ratio.h>
#include <lu_r_plot.h>
#include <lu_rep_links.h>
#include <lu_sctag_link.h>
#include <lu_scomp_form.h>
#include <lu_cgi_val.h>
#include <lu_param_list.h>
#include <lu_df_info.h>
#include <lu_proj_file.h>
#include <lu_prod_par.h>
#include <sc_anum_res.h>
#include <sc_sres_rep.h>
#include <sc_pep_xml.h>
#include <sc_mzidentml.h>
#include <sc_xlink.h>
#include <sc_quan.h>
using std::vector;
using std::count;
using std::ostream;
using std::ostringstream;
using std::string;
using std::stable_sort;
using std::cout;
using std::unique;
using std::set_intersection;
using std::set_difference;
using std::set_union;
using std::back_inserter;
using std::make_pair;
using std::endl;
using std::copy;
using std::map;
using std::runtime_error;

typedef map <int, vector <const PeptidePosition*> > MapIntVectorConstPeptidePositionPtr;
typedef MapIntVectorConstPeptidePositionPtr::const_iterator MapIntVectorConstPeptidePositionPtrConstIterator;

class SortProteinReportByTotalDiscriminantScore {
	bool taxCheck;
	vector <TaxonomyMatch*> tm;
public:
	SortProteinReportByTotalDiscriminantScore ( bool taxCheck, const vector <TaxonomyMatch*>& tm ) :
		taxCheck ( taxCheck ),
		tm ( tm ) {}
	int operator () ( const SearchResultsProteinLine* a, const SearchResultsProteinLine* b ) const
	{
		double maxA = -std::numeric_limits<double>::max();
		int indexA = -1;
		for ( ConstSearchResultsProteinInfoPtrVectorSizeType i = 0 ; i < a->searchResultsInfo.size () ; i++ ) {
			if ( a->searchResultsInfo [i]->getTotalDiscriminantScore () ) {
				maxA = genMax ( maxA, a->searchResultsInfo [i]->getTotalDiscriminantScore () );
				indexA = i;
			}
		}
		double maxB = -std::numeric_limits<double>::max();
		int indexB = -1;
		for ( ConstSearchResultsProteinInfoPtrVectorSizeType j = 0 ; j < b->searchResultsInfo.size () ; j++ ) {
			if ( b->searchResultsInfo [j]->getTotalDiscriminantScore () ) {
				maxB = genMax ( maxB, b->searchResultsInfo [j]->getTotalDiscriminantScore ());
				indexB = j;
			}
		}
		if ( indexA == -1 && indexB == -1 ) return false;
		if ( maxA == maxB ) {
			double aScore = a->searchResultsInfo [indexA]->getScore ();
			double bScore = b->searchResultsInfo [indexB]->getScore ();

			if ( aScore == bScore ) {
				if ( taxCheck ) {
					const string& sp1 = const_cast <SearchResultsProteinLine*>(a)->getSpecies ();
					const string& sp2 = const_cast <SearchResultsProteinLine*>(b)->getSpecies ();
					for ( std::vector <TaxonomyMatch*>::size_type k = 0 ; k < tm.size () ; k++ ) {
						if ( !tm [k]->equal ( sp1, sp2 ) ) {
							return tm [k]->greaterThan ( sp1, sp2 );
						}
					}
					return checkAccessionNumbers ( a, b );
				}
				else return checkAccessionNumbers ( a, b );
			}
			else
				return aScore > bScore;
		}
		else
			return maxA > maxB;
	}
	int checkAccessionNumbers ( const SearchResultsProteinLine* a, const SearchResultsProteinLine* b ) const
	{
		const string& acc1 = a->getAccessionNumber ();
		const string& acc2 = b->getAccessionNumber ();
		const string& ui1 = const_cast <SearchResultsProteinLine*>(a)->getUniprotID ();
		const string& ui2 = const_cast <SearchResultsProteinLine*>(b)->getUniprotID ();
		bool accInGName1 = ui1.find ( acc1 ) != string::npos;
		bool accInGName2 = ui2.find ( acc2 ) != string::npos;
		if ( accInGName1 && !accInGName2 ) return 0;
		if ( !accInGName1 && accInGName2 ) return 1;
		int dollarIdx1 = acc1.find ( '$' );
		int dollarIdx2 = acc2.find ( '$' );
		int idx1 = (dollarIdx1 == string::npos) ? 0 : dollarIdx1 + 1;
		int idx2 = (dollarIdx2 == string::npos) ? 0 : dollarIdx2 + 1;
		bool moreThanO1 = acc1 [idx1] >= 'O';
		bool moreThanO2 = acc2 [idx2] >= 'O';
		if ( moreThanO1 && !moreThanO2 ) return 1;
		if ( !moreThanO1 && moreThanO2 ) return 0;
		return acc1 < acc2;
	}
};
class SCMSBridgeLink {
public:
	SCMSBridgeLink () {}
	void write ( ostream& os, const string& id, const string& searchKey, const Tolerance* tol, const string& systematicError, const StringVector& accNumbers ) const;
	void printHTML ( ostream& os ) const;
};

void SCMSBridgeLink::write ( ostream& os, const string& id, const string& searchKey, const Tolerance* tol, const string& systematicError, const StringVector& accNumbers ) const
{
	ProgramLink::openLink ( os, "scMSBridgeLink", -1 );
	int end = id.find ( "-" );
	printCGIString ( os, "fraction", id.substr ( 0, end ) );
	printCGIString ( os, "spot_number", id.substr ( end+1, id.length ()-end-1 ) );
	printCGIString ( os, "search_key", searchKey );
	tol->putCGI ( os, "ms_parent_mass" );
	printCGI ( os, "ms_parent_mass_systematic_error", systematicError );
	printCGIContainer ( os, "entry_data", accNumbers );
	os << "\\\">";
	os << "MS Hits";
	ProgramLink::closeLink ( os );
}
void SCMSBridgeLink::printHTML ( ostream& os ) const
{
	os << "scMSBridgeLink" << "=\"";
	os << ProgramLink::getURLStart ( "msform" );
	os << "?";
	printCGIString ( os, "form", "msbridge" );
	printCGIString ( os, "link_search_type", "No Link" );
	printCGI ( os, "separate_proteins", true );
	os << "\";\n";
}

bool SearchResultsProteinLine::reportNumber = false;
bool SearchResultsProteinLine::reportLinks = false;
string SearchResultsProteinLine::styleID1 = "sc_stripe_1";
string SearchResultsProteinLine::styleID2 = "sc_stripe_2";
SearchResultsProteinLine::SearchResultsProteinLine ( const string& accessionNumber, const ProteinInfo& proteinInfo, int numUnique, const ConstSearchResultsProteinInfoPtrVector& searchResultsInfo ) :
	accessionNumber ( accessionNumber ),
	proteinInfo ( proteinInfo ),
	numUnique ( numUnique ),
	searchResultsInfo ( searchResultsInfo )
{
}
void SearchResultsProteinLine::printHTMLHeader ( ostream& os, const StringVector& searchNames, bool reportUniqPeps ) const
{
	tableRowStart ( os );
		if ( reportNumber ) tableHeader ( os, "Rank", "", "", false, 0, 2 );
		if ( reportUniqPeps ) tableHeader ( os, "Uniq Pep", "", "", false, 0, 2 );
		ProteinInfo::printHTMLANumHeader ( os, 2 );
		for ( ConstSearchResultsProteinInfoPtrVectorSizeType i = 0 ; i < searchResultsInfo.size () ; i++ ) {
			if ( searchResultsInfo [i]->getColspan () ) {
				tableHeaderStart ( os, ( i % 2 == 0 ) ? styleID1 : styleID2, "", false, searchResultsInfo [i]->getColspan () );
				os << ( sresMergedFlag ? "Merged" : searchNames [i] );
				tableHeaderEnd ( os );
			}
		}
		ProteinInfo::printHTMLHeader ( os, 2 );
	tableRowEnd ( os );

	tableRowStart ( os );
		for ( ConstSearchResultsProteinInfoPtrVectorSizeType j = 0 ; j < searchResultsInfo.size () ; j++ ) {
			searchResultsInfo [j]->printHeaderHTML ( os, j, ( j % 2 == 0 ) ? styleID1 : styleID2 );
		}
	tableRowEnd ( os );
}
void SearchResultsProteinLine::printHTML ( ostream& os, const SResLink& sresLink, const string& id, bool reportUniqPeps ) const
{
	tableRowStart ( os );
		if ( reportLinks ) {
			if ( reportNumber ) {
				tableDataStart ( os );
				for ( StringVectorSizeType i = 0 ; i < idStr.size () ; i++ ) {
					sresLink.write ( os, proteinInfo.getFullAccessionNumber (), id, idStr [i] );
					if ( i != idStr.size () - 1 ) os << "<br />";
					os << endl;
				}
				tableDataEnd ( os );
			}
		}
		else {
			if ( reportNumber ) tableCell ( os, getIDStrVecOutput () );
		}
		if ( reportUniqPeps ) {
			if ( gen_strcontains ( idStr [0], '-' ) ) tableCell ( os, numHomology );
			else tableCell ( os, "" );
		}
		proteinInfo.printHTMLANum ( os, false );
		for ( ConstSearchResultsProteinInfoPtrVectorSizeType i = 0 ; i < searchResultsInfo.size () ; i++ ) {
			searchResultsInfo [i]->printHTML ( os, false, ( i % 2 == 0 ) ? styleID1 : styleID2 );
		}
		proteinInfo.printHTML ( os, false );
	tableRowEnd ( os );
}
void SearchResultsProteinLine::printDelimitedHeader ( ostream& os, bool ID, const StringVector& searchNames, bool reportUniqPeps ) const
{
	delimitedRowStart ( os );
		if ( ID )								delimitedEmptyCell ( os );
		if ( reportNumber )						delimitedEmptyCell ( os );
		if ( reportUniqPeps )					delimitedEmptyCell ( os );
		if ( ProteinInfo::getReportAccession () )delimitedEmptyCell ( os );
		for ( ConstSearchResultsProteinInfoPtrVectorSizeType i = 0 ; i < searchResultsInfo.size () ; i++ ) {
			if ( searchResultsInfo [i]->getColspan () ) {
				delimitedHeader ( os, ( sresMergedFlag ? "Merged" : searchNames [i] ) );
				delimitedEmptyNCells ( os, searchResultsInfo [i]->getColspan () - 1 );
			}
		}
		delimitedEmptyNCells ( os, ProteinInfo::getColspan () );
	delimitedRowEnd ( os );

	delimitedRowStart ( os );
		if ( ID ) delimitedHeader ( os, "ID" );
		if ( reportNumber ) delimitedHeader ( os, "Rank" );
		if ( reportUniqPeps ) delimitedHeader ( os, "Uniq Pep" );
		ProteinInfo::printDelimitedANumHeader ( os );
		for ( ConstSearchResultsProteinInfoPtrVectorSizeType j = 0 ; j < searchResultsInfo.size () ; j++ ) {
			searchResultsInfo [j]->printHeaderDelimited ( os, sresMergedFlag ? -1 : j );
		}
		ProteinInfo::printDelimitedHeader ( os );
	delimitedRowEnd ( os );
}
void SearchResultsProteinLine::printDelimited ( ostream& os, const string& id, bool reportUniqPeps ) const
{
	delimitedRowStart ( os );
		if ( !id.empty () ) delimitedCell ( os, id );

		if ( reportNumber ) delimitedCell ( os, "[" + getIDStrVecOutput () + "]" );
		if ( reportUniqPeps ) {
			if ( gen_strcontains ( idStr [0], '-' ) ) delimitedCell ( os, numHomology );
			else delimitedCell ( os, "" );
		}
		proteinInfo.printDelimitedANum ( os );
		for ( ConstSearchResultsProteinInfoPtrVectorSizeType i = 0 ; i < searchResultsInfo.size () ; i++ ) {
			searchResultsInfo [i]->printDelimited ( os );
		}
		proteinInfo.printDelimited ( os );
	delimitedRowEnd ( os );
}
bool SearchResultsPeptideLine::reportNumber = false;
bool SearchResultsPeptideLine::reportLinks = false;
string SearchResultsPeptideLine::styleID1 = "sc_stripe_1";
string SearchResultsPeptideLine::styleID2 = "sc_stripe_2";

SearchResultsPeptideLine::SearchResultsPeptideLine ( const ProteinInfo& proteinInfo, const SearchResultsPeptideHitConstPtrVector& srph, const ConstSearchResultsProteinInfoPtrVector& srpi ) :
	proteinInfo ( proteinInfo ),
	srph ( srph ),
	srpi ( srpi )
{
}
SearchResultsPeptideLine::SearchResultsPeptideLine ( const SearchResultsPeptideHitConstPtrVector& srph, const ConstSearchResultsProteinInfoPtrVector& srpi ) :
	srph ( srph ),
	srpi ( srpi )
{
}
void SearchResultsPeptideLine::printHTMLTimeHeader ( ostream& os, const StringVector& searchNames ) const
{
	tableRowStart ( os );
		if ( reportNumber ) tableHeader ( os, "", "", "", false, 0, 2 );
		ProteinInfo::printHTMLANumHeader ( os, 2 );
		PeptidePosition::printHeaderHTMLPeak1 ( os, "", 2 );
		PeptidePosition::printHeaderHTMLPeak2 ( os, -1, "", 2 );

		for ( ConstPPPeptideHitInfoPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {
			int colspan = PeptidePosition::getColspan ( sresMergedFlag ? -1 : i, false, false );
			colspan += srph [i]->getPeptideHitInfo ()->getColspan ();
			colspan += srpi [i]->getColspan ();
			if ( colspan ) {
				tableHeaderStart ( os, ( i % 2 == 0 ) ? styleID1 : styleID2, "", false, colspan );
				os << ( sresMergedFlag ? "Merged" : searchNames [i] );
				tableHeaderEnd ( os );
			}
		}
		ProteinInfo::printHTMLHeader ( os, 2 );
	tableRowEnd ( os );
	tableRowStart ( os );
		for ( ConstSearchResultsProteinInfoPtrVectorSizeType j = 0 ; j < srph.size () ; j++ ) {
			string styleID = ( j % 2 == 0 ) ? styleID1 : styleID2;
			srpi [j]->printHeaderHTML ( os, sresMergedFlag ? -1 : j, styleID );
			PeptidePosition::printHeaderHTML ( os, sresMergedFlag ? -1 : j, styleID, false );
			srph [j]->getPeptideHitInfo ()->printHeaderHTML ( os, styleID );
		}
	tableRowEnd ( os );
}
void SearchResultsPeptideLine::printHTMLHeader ( ostream& os, const StringVector& searchNames ) const
{
	if ( searchNames.size () > 1 ) {
		tableRowStart ( os );
			for ( ConstPPPeptideHitInfoPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {
				int colspan = PeptidePosition::getColspan ( sresMergedFlag ? -1 : i );
				colspan += srph [i]->getPeptideHitInfo ()->getColspan ();
				if ( colspan ) {
					tableHeaderStart ( os, ( i % 2 == 0 ) ? styleID1 : styleID2, "", false, colspan );
					os << ( sresMergedFlag ? "Merged" : searchNames [i] );
					tableHeaderEnd ( os );
				}
			}
		tableRowEnd ( os );
	}
	tableRowStart ( os );
		for ( ConstSearchResultsProteinInfoPtrVectorSizeType j = 0 ; j < srph.size () ; j++ ) {
			string styleID = ( j % 2 == 0 ) ? styleID1 : styleID2;
			PeptidePosition::printHeaderHTML ( os, sresMergedFlag ? -1 : j, styleID );
			srph [j]->getPeptideHitInfo ()->printHeaderHTML ( os, styleID );
		}
	tableRowEnd ( os );
}
void SearchResultsPeptideLine::printHTMLProteinHeader ( ostream& os, const StringVector& searchNames ) const
{
	tableRowStart ( os );
		if ( searchNames.size () > 1 ) tableHeader ( os, "Search Name", styleID1 );
		srpi [0]->printHeaderHTML ( os, 0, styleID1 );
	tableRowEnd ( os );
}
void SearchResultsPeptideLine::printHTML ( ostream& os, const SCMSTagLink& smtl, const string& id ) const		// Peptide Report
{
	tableRowStart ( os );
		for ( ConstPPPeptideHitInfoPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {
			srph [i]->printHTML ( os, proteinInfo, i, ( i % 2 == 0 ) ? styleID1 : styleID2, smtl, id, sresMergedFlag );
		}
	tableRowEnd ( os );
}
void SearchResultsPeptideLine::printHTML ( ostream& os, int pline, bool empty, bool aNumEmpty, const SCMSTagLink& smtl ) const		// Time Report
{
	if ( !reportLinks ) {
		empty = false;
		aNumEmpty = false;
	}
	tableRowStart ( os );
		if ( reportNumber ) {
			if ( empty ) tableEmptyCell ( os );
			else tableCell ( os, pline );
		}
		proteinInfo.printHTMLANum ( os, aNumEmpty );
		srph [0]->printHTMLPeak1 ( os, "" );
		srph [0]->printHTMLPeak2 ( os, 0, "", smtl, true );
		for ( ConstPPPeptideHitInfoPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {
			string styleID = ( i % 2 == 0 ) ? styleID1 : styleID2;
			srpi [i]->printHTML ( os, aNumEmpty, styleID );
			srph [i]->printHTML ( os, proteinInfo, sresMergedFlag ? -1 : i, styleID, smtl, "", false, false );
		}
		proteinInfo.printHTML ( os, aNumEmpty );
	tableRowEnd ( os );
}
void SearchResultsPeptideLine::printHTMLProtein ( ostream& os, const StringVector& searchNames ) const
{
	for ( ConstPPPeptideHitInfoPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {
		bool stripe = ( i % 2 == 0 );
		string styleID = stripe ? styleID1 : styleID2;
		tableRowStart ( os );
			if ( searchNames.size () > 1 ) tableCell ( os, sresMergedFlag ? "Merged" : searchNames [i], true, true, styleID );
			srpi [i]->printHTML ( os, false, styleID );
		tableRowEnd ( os );
	}
}
void SearchResultsPeptideLine::printDelimitedHeader ( std::ostream& os, const StringVector& searchNames, bool ID, bool reportUniqPeps ) const
{
	delimitedRowStart ( os );
		if ( ID )				delimitedEmptyCell ( os );
		if ( reportNumber )		delimitedEmptyCell ( os );
		if ( reportUniqPeps )	delimitedEmptyCell ( os );
		if ( ProteinInfo::getReportAccession () )delimitedEmptyCell ( os );
		for ( ConstPPPeptideHitInfoPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {
			int colspan = PeptidePosition::getColspan ( i, true );
			colspan += srph [i]->getPeptideHitInfo ()->getColspan ();
			colspan += srpi [i]->getColspan ();
			if ( colspan ) {
				delimitedHeader ( os, searchNames [i] );
				delimitedEmptyNCells ( os, colspan - 1 );
			}
		}
		delimitedEmptyNCells ( os, ProteinInfo::getColspan () );
	delimitedRowEnd ( os );
	delimitedRowStart ( os );
		if ( ID ) delimitedHeader ( os, "ID" );
		if ( reportNumber ) delimitedEmptyCell ( os );
		if ( reportUniqPeps ) delimitedHeader ( os, "Uniq Pep" );
		ProteinInfo::printDelimitedANumHeader ( os );
		for ( ConstPPPeptideHitInfoPtrVectorSizeType j = 0 ; j < srph.size () ; j++ ) {
			srpi [j]->printHeaderDelimited ( os, sresMergedFlag ? -1 : j );
			PeptidePosition::printHeaderDelimited ( os, sresMergedFlag ? -1 : j );
			srph [j]->getPeptideHitInfo ()->printHeaderDelimited ( os );
		}
		ProteinInfo::printDelimitedHeader ( os );
	delimitedRowEnd ( os );
}
void SearchResultsPeptideLine::printDelimited ( ostream& os, const string& idStr, int numHomology, const string& id, bool reportUniqPeps ) const
{
	delimitedRowStart ( os );
		if ( !id.empty () ) delimitedCell ( os, id );

		if ( reportNumber ) delimitedCell ( os, "[" + idStr + "]" );
		if ( reportUniqPeps ) {
			if ( gen_strcontains ( idStr, '-' ) ) delimitedCell ( os, numHomology );
			else delimitedCell ( os, "" );
		}
		proteinInfo.printDelimitedANum ( os );
		for ( ConstPPPeptideHitInfoPtrVectorSizeType i = 0 ; i < srph.size () ; i++ ) {
			srpi [i]->printDelimited ( os );
			if ( srph [i] )
				srph [i]->printDelimited ( os, proteinInfo, sresMergedFlag ? -1 : i );
			else
				SearchResultsPeptideHit::printDelimitedEmpty ( os, i );
		}
		proteinInfo.printDelimited ( os );
	delimitedRowEnd ( os );
}
void SearchResultsPeptideLine::printProteinSequenceHTML ( ostream& os, bool coverage ) const
{
	proteinInfo.printHTMLLines ( os );
	if ( coverage ) {
		for ( ConstSearchResultsProteinInfoPtrVectorSizeType i = 0 ; i < srpi.size () ; i++ ) {
			htmlPrintProteinSequence ( os, proteinInfo.getProteinSequence (), IntVector ( 0 ), srpi [i]->getCoverageMap ().getAACovered (), true );
		}
	}
}
bool SearchResultsPeptideLine::outputQuanResults ( ostream& os, const StringVector& searchNames, bool area ) const
{
	bool flag = false;
	for ( ConstSearchResultsProteinInfoPtrVectorSizeType i = 0 ; i < srpi.size () ; i++ ) {
		if ( srph [i] ) {
			if ( srph [i]->getPeptidePosition ()->outputQuanResults ( os, sresMergedFlag ? "Merged" : searchNames [i], getRepeats ( i ) == 1 ? 1 : 2, area ) ) {
				flag = true;
			}
		}
	}
	return flag;
}
DoubleVectorVector SearchResultsPeptideLine::getRatios ( bool area ) const
{
	DoubleVectorVector dvv ( srpi.size () );
	for ( ConstSearchResultsProteinInfoPtrVectorSizeType i = 0 ; i < srpi.size () ; i++ ) {
		if ( srph [i] ) {
			if ( area ) dvv [i] = srph [i]->getPeptidePosition ()->getAreaRatios ();
			else		dvv [i] = srph [i]->getPeptidePosition ()->getIntensityRatios ();
		}
	}
	return dvv;
}
SearchResultsReport::SearchResultsReport ( const vector <SearchResults*>& sr, const string& reportHitsType, const string& id ) :
	searchResults ( sr ),
	reportHitsType ( reportHitsType ),
	id ( id )
{
	for ( SearchResultsPtrVectorSizeType k = 0 ; k < sr.size () ; k++ ) {
		fullSearchNames.push_back ( sr [k]->getProjectName () + '/' + sr [k]->getResultsName () );
		searchNames.push_back ( sr [k]->getResultsName () );
	}
}
SearchResultsReport::~SearchResultsReport () {}

void SearchResultsReport::printIDHTML ( ostream& os, const string& id ) const
{
	if ( id != SearchResults::getDefaultID () ) {
		if ( id.find ( "-" ) != string::npos || searchResults.size () > 1 ) {
			ParameterList::printHTML ( os, "ID", id );
		}
		else {
			ParameterList::printHTML ( os, "ID", id + " (" + searchResults [0]->getFractionName (atoi (id.c_str ())-1) + ")" );
		}
	}
}
bool SearchResultsReport::isRawForwarding () const
{
	return InfoParams::instance ().getBoolValue ( "raw_data_forwarding" ) && QuantitationRatio::getQuanReport ();
}
void SearchResultsReport::printDelimitedReport ( ostream& os, int i, bool last, const string& outputDirectory, const string& outputFilename, bool& delimHeaderPrinted ) const
{
	static PPTempFile pptf ( "", "" );
	static string fullPath;
	fullPath = outputDirectory.empty () ? pptf.getFullPath () : outputDirectory;
	string filename = outputFilename.empty () ? "report.txt" : outputFilename;
	static string actualPath = fullPath + SLASH + filename;
	static string outputPath = pptf.getURL ();

	if ( i == 0 ) {
		genCreateDirectory ( fullPath );
		GenOFStream ost ( actualPath, std::ios_base::out );
		if ( sresXLinks ) {
			printCrosslinksDelimited ( ost );
		}
		else {
			printHistogramDelimited ( ost );
			if ( !delimHeaderPrinted ) {
				if ( printDelimitedHeader ( ost ) )	delimHeaderPrinted = true;
			}
			printDelimited ( ost );
		}
	}
	else {
		GenOFStream ost ( actualPath, std::ios_base::out | std::ios_base::app );
		if ( sresXLinks ) {
			printCrosslinksDelimited ( ost );
		}
		else {
			if ( !delimHeaderPrinted ) {
				if ( printDelimitedHeader ( ost ) )	delimHeaderPrinted = true;
			}
			printDelimited ( ost );
		}
	}
	if ( last && outputDirectory.empty () ) {
		os << "<a href=\"";
		os << outputPath;
		os << "/";
		os << filename;
		os << "\">";
		if ( sresViewer )	os << "Results File";
		else				os << "Click to see report";
		os << "</a>";
		os << "<br />";
		os << endl;
	}
}
void SearchResultsReport::printHistogramDelimited ( ostream& os ) const
{
	for ( SearchResultsPtrVectorSizeType i = 0 ; i < searchResults.size () ; i++ ) {
		delimitedRowStart ( os );
			delimitedCell ( os, "Search Name:" );
			delimitedCell ( os, fullSearchNames [i] );
		delimitedRowEnd ( os );
	}
}
void SearchResultsReport::printHistogramHTML ( ostream& os ) const
{
	ParameterList::printHTMLContainer ( os, reportHitsType, fullSearchNames );
	os << "<br />" << endl;
	for ( SearchResultsPtrVectorSizeType i = 0 ; i < searchResults.size () ; i++ ) {
		searchResults [i]->drawHistogram ( os );
	}
	os << "<br />" << "<br />" << endl; 
}
bool SearchResultsProteinReport::reportUniqPeps = false;
SearchResultsProteinReport::SearchResultsProteinReport ( const vector <SearchResults*>& sr, bool remove, const MapStringToStringVector& aNumList, const string& reportHitsType, const string& reportHomologousProteins, const string& id ) :
	SearchResultsReport ( sr, reportHitsType, id ),
	reportHomologousProteins ( reportHomologousProteins )
{
	ujm->writeMessage ( cout, "Collating proteins" );
	StringVector reportedHits;
	if ( reportHitsType == "Protein Difference" || reportHitsType == "Protein Intersection" ) {
		StringVectorVector aNums;
		for ( SearchResultsPtrVectorSizeType x = 0 ; x < searchResults.size () ; x++ ) {
			aNums.push_back ( searchResults [x]->getAccessionNumbers ( id ) );
		}
		getReportedProteinHits ( reportHitsType, aNums, reportedHits );
	}
	else {
		reportedHits = SearchResults::getMergedAccNumbers ( id );
	}
	StringVector tempHits;
	if ( !aNumList.empty () ) {		// Filter by accession number
		ANumResults anr ( aNumList );
		StringVector aNumFilter = anr.getAccessionNumbers ();
		if ( remove )
			set_difference ( reportedHits.begin (), reportedHits.end (), aNumFilter.begin (), aNumFilter.end (), back_inserter ( tempHits ) );
		else
			set_intersection ( reportedHits.begin (), reportedHits.end (), aNumFilter.begin (), aNumFilter.end (), back_inserter ( tempHits ) );
		reportedHits = tempHits;
		tempHits.clear ();
	}
	for ( StringVectorSizeType i = 0 ; i < reportedHits.size () ; i++ ) {
		string aNum = reportedHits [i];
		ProteinInfo proteinInfo ( aNum );
		ConstSearchResultsProteinInfoPtrVector sri;
		const SearchResultsProteinInfo* jpi = SearchResults::getSearchResultsProteinInfo2 ( aNum, id );
		if ( sresMergedFlag )
			sri.push_back ( jpi );
		else {
			for ( SearchResultsPtrVectorSizeType j = 0 ; j < searchResults.size () ; j++ ) {
				sri.push_back ( searchResults [j]->getSearchResultsProteinInfo ( aNum, id ) );
			}
		}
		srprotl.push_back ( new SearchResultsProteinLine ( reportedHits [i], proteinInfo, jpi->getNumUnique (), sri ) );
	}
	ProteinInfo::resetANumMap ();
	ujm->writeMessage ( cout, "Sorting proteins" );
	stable_sort ( srprotl.begin (), srprotl.end (), SortProteinReportByTotalDiscriminantScore ( ProteinInfo::getTaxonomyCheck (), ProteinInfo::getTaxonomyMatch () ) );
	if ( sresMainAndSupplementary ) calculateMainAndSupplementaryHits ();
	fullAccList = getAccessionNumberListPair ();
	PeptidePosition::initialise ( searchResults );
	ujm->writeMessage ( cout, "Finished collating proteins" );
}
void SearchResultsProteinReport::getReportedProteinHits ( const string& reportedHitsType, const StringVectorVector& aNums, StringVector& reportedHits )
{
	int numSearches = aNums.size ();
	if ( reportedHitsType == "Protein Difference" ) {
		for ( int i = 0 ; i < numSearches ; i++ ) {
			StringVector tempHits1 = aNums [i];
			for ( int j = 0 ; j < numSearches ; j++ ) {
				StringVector tempHits2;
				if ( j != i ) {
					set_difference	( tempHits1.begin (), tempHits1.end (), aNums [j].begin (), aNums [j].end (), back_inserter ( tempHits2 ) );
					tempHits1 = tempHits2;
				}
			}
			copy ( tempHits1.begin (), tempHits1.end (), back_inserter ( reportedHits ) );
		}
	}
	else {						// Protein Intersection
		StringVector tempHits;
		reportedHits = aNums [0];
		for ( int a = 1 ; a < numSearches ; a++ ) {
			set_intersection( reportedHits.begin (), reportedHits.end (), aNums [a].begin (), aNums [a].end (), back_inserter ( tempHits ) );
			reportedHits = tempHits;
			tempHits.clear ();
		}
	}
}
void SearchResultsProteinReport::calculateMainAndSupplementaryHits ()
{
	vector <SearchResultsProteinLine*> mainSRL;
	set_score_matrix ( "BLOSUM62" );
	ujm->writeMessage ( cout, "Compiling main and supplementary hits" );
	for ( SearchResultsProteinLinePtrVectorSizeType i = 0 ; i < srprotl.size () ; i++ ) { // For each protein
		int num = i + 1;
		if ( num % 10 == 0 ) ujm->writeMessage ( cout, "Compiling main and supplementary hits, " + gen_itoa ( num ) + "/" + gen_itoa ( srprotl.size () ) + " proteins" );
		checkHomologyMatch ( mainSRL, srprotl [i] );
	}
	srprotl.clear ();
	makeProteinLines ( mainSRL );
}
void SearchResultsProteinReport::makeProteinLines ( const vector <SearchResultsProteinLine*>& mainSRL )
{
	for ( SearchResultsProteinLinePtrVectorSizeType i = 0 ; i < mainSRL.size () ; i++ ) {
		srprotl.push_back ( mainSRL [i] );
		if ( mainSRL [i]->supSRL.size () ) {
			reportUniqPeps = true;
			makeProteinLines ( mainSRL [i]->supSRL );
		}
	}
}
bool SearchResultsProteinReport::checkHomologyMatch ( vector <SearchResultsProteinLine*>& mainSRL, SearchResultsProteinLine* srprotl )
{
	bool newHomologyMatch = true;
	int numUnique = srprotl->getNumUnique ();
	for ( IntVectorSizeType i = 0 ; i < mainSRL.size () ; i++ ) {
		int numMatches = getNumHomologyMatches ( mainSRL [i]->getAccessionNumber (), srprotl->getAccessionNumber () );
		if ( numMatches ) {
			newHomologyMatch = false;
			if ( reportHomologousProteins == "All" || ( reportHomologousProteins == "Interesting" && numUnique != numMatches ) ) {
				srprotl->setIDStr ( mainSRL [i]->getIDStr () + "-" + gen_itoa ( mainSRL [i]->supSRL.size () + 1 ) );
				srprotl->setNumHomology ( numUnique - numMatches );
				mainSRL [i]->supSRL.push_back ( srprotl );
			}
		}
	}
	if ( newHomologyMatch ) {
		srprotl->setIDStr ( gen_itoa ( mainSRL.size () + 1 ) );
		mainSRL.push_back ( srprotl );
	}
	return newHomologyMatch;
}
int SearchResultsProteinReport::getNumHomologyMatches ( const string& aNum, const string& mainANum )
{
	int numMatches = 0;
	string lastSequence;
	PairSearchResultsPeptideHitPtrVectorConstIterator ppppvci = searchResults [0]->getSearchPeptidePositionIters ( mainANum, id );
	int num = 0;
	string sTotal = gen_itoa ( ppppvci.second-ppppvci.first );
	for ( SearchResultsPeptideHitPtrVector::const_iterator d = ppppvci.first ; d != ppppvci.second ; d++, num++ ) {
		string currentSequence = (*d)->getSequence ();
		string currentDatabaseSequence = (*d)->getDBPeptide ();
		int currentDatabaseSequenceLen = currentDatabaseSequence.length ();
		if ( num > 0 && num % 500 == 0 ) ujm->writeMessage ( cout, "Large protein. Number of peptides processed = " + gen_itoa ( num ) + "/" + sTotal );
		if ( currentSequence != lastSequence ) {
			string lastHomSequence;
			PairSearchResultsPeptideHitPtrVectorConstIterator ppppvci1 = searchResults [0]->getSearchPeptidePositionIters ( aNum, id );
			for ( SearchResultsPeptideHitPtrVector::const_iterator e = ppppvci1.first ; e != ppppvci1.second ; e++ ) {
				string homSequence = (*e)->getSequence ();
				string homDatabaseSequence = (*e)->getDBPeptide ();
				if ( homSequence != lastHomSequence ) {
					if ( currentDatabaseSequenceLen == homDatabaseSequence.length () ) {
						if ( matrixMatch ( currentDatabaseSequence, homDatabaseSequence ) ) {
							numMatches++;
							(*d)->setSubsequentOccurence ();
							break;
						}
					}
				}
				lastHomSequence = homSequence;
			}
		}
		lastSequence = currentSequence;
	}
	return numMatches;
}
void SearchResultsProteinReport::printBatchTagAccNumberList ( ostream& os, const StringVector& aNum )
{
	StringVectorVector accs;
	StringVector dbases = ProteinInfo::getDBSearchList1 ();
	BoolDeque dbFlags = ProteinInfo::getDBSearchFlagList1 ();
	accs.resize ( dbases.size () );
	for ( StringVectorSizeType i = 0 ; i < aNum.size () ; i++ ) {
		PairIntString pis = ProteinInfo::getANumPair ( aNum [i] );
		accs [pis.first].push_back ( pis.second );
	}
	int numDB = count ( dbFlags.begin (), dbFlags.end (), true );
	StringVector acc;
	StringVector up;
	StringVector dbList;
	SetPairStringString spss;
	for ( StringVectorVectorSizeType j = 0 ; j < accs.size () ; j++ ) {
		string db = dbases [j];
		if ( isSuffix ( db, "UserProtein.fasta" ) ) {
			for ( StringVectorSizeType k = 0 ; k < accs [j].size () ; k++ ) {
				ProteinInfo pi ( gen_itoa ( j+1 ) + "$" + accs [j][k] );
				PairStringString pss = make_pair ( ">" + pi.getName (),  pi.getProteinSequence () );
				PairSetPairStringStringIteratorBool flag = spss.insert ( pss );
				if ( flag.second ) {
					up.push_back ( pss.first );
					up.push_back ( pss.second );
				}
			}
		}
		else {
			if ( dbFlags [j] ) {
				dbList.push_back ( db );
				for ( StringVectorSizeType k = 0 ; k < accs [j].size () ; k++ ) {
					if ( k == 0 && numDB > 1 ) acc.push_back ( ">" + db );
					const string& ac = accs [j][k];
					if ( ac [0] != '-' ) {
						acc.push_back ( ac );
					}
				}
			}
		}
	}
	ProteinInfo::resetANumMap ();
	if ( !acc.empty () ) printHTMLFORMHiddenVector ( os, "accession_nums", acc );
	if ( !up.empty () ) {
		dbList.push_back ( "User Protein" );
		printHTMLFORMHiddenVector ( os, "user_protein_sequence", up );
	}
	if ( !dbList.empty () ) printHTMLFORMHiddenContainer ( os, "database", dbList );
}
void SearchResultsProteinReport::printHTML ( ostream& os ) const
{
	printIDHTML ( os, id );
	ParameterList::printHTML ( os, "Number of Proteins (After Homology Filtering)", fullAccList.first.size () );
	int numNegAccNums = fullAccList.second.size () - fullAccList.first.size ();
	if ( numNegAccNums ) {
		ParameterList::printHTML ( os, "Number of Decoy Proteins (After Homology Filtering)", numNegAccNums );
	}
	for ( SearchResultsPtrVectorSizeType x = 0 ; x < searchResults.size () ; x++ ) {
		os << "Search " << x+1 << " Number of Data Spectra: <b>" << searchResults [x]->getMSMSSpectrumInfoSize ( id ) << "</b><br />" << endl; 
		if ( id == SearchResults::getDefaultID () ) {
			string searchEndTime = searchResults [x]->getSearchEndTime ();
			string searchTime = searchResults [x]->getSearchTime ();
			if ( !searchEndTime.empty () )	os << "Search " << x+1 << " End Time: <b>" << searchEndTime << "</b><br />" << endl;
			if ( !searchTime.empty () )		os << "Search " << x+1 << " Search Time: <b>" << searchTime << "</b><br />" << endl; 
		}
		os << "Search " << x+1 << " Proteins (Pre Homology Filtering): " << searchResults [x]->getNumProteins ( id ) << "<br />" << endl; 
		searchResults [x]->printDatabaseInfoHTML ( os );
		PeptidePosition::printParamsHTML ( os, x );
	}
	os << "<br />" << endl;
	if ( srprotl.size () == 0 ) return;

	string idStr;
	if ( id != SearchResults::getDefaultID () ) idStr = id;

	printHTMLFORMStart ( os, "post", "msform", false, false, true );

	printHTMLFORMSubmit ( os, "Batch-Tag of Listed Accession Numbers" );
	os << "<br />" << endl;
	printHTMLFORMHidden ( os, "form", "batchtag" );
	printHTMLFORMHidden ( os, "search_key", PeptidePosition::getSearchKey ( 0 ) );
	printBatchTagAccNumberList ( os, fullAccList.first );
	printHTMLFORMEnd ( os );

	ProteinInfo::initialiseAccessionNumberLink ( os );
	SResLink sresLink ( isRawForwarding () );
	if ( SearchResultsProteinLine::getReportLinks () ) {
		startJavascript ( os );
		sresLink.printHTML ( os );
		endJavascript ( os );
	}
	tableStart ( os, true );
		if ( !srprotl.empty () ) srprotl [0]->printHTMLHeader ( os, fullSearchNames, reportUniqPeps );
		for ( SearchResultsProteinLinePtrVectorSizeType i = 0 ; i < srprotl.size () ; i++ ) {
			srprotl [i]->printHTML ( os, sresLink, idStr, reportUniqPeps );
		}
	tableEnd ( os );
}
bool SearchResultsProteinReport::printDelimitedHeader ( ostream& os ) const
{
	if ( !srprotl.empty () ) {
		srprotl [0]->printDelimitedHeader ( os, id != SearchResults::getDefaultID (), fullSearchNames, reportUniqPeps );
		return true;	// Header printed
	}
	else return false;	// Header not printed
}
void SearchResultsProteinReport::printDelimited ( ostream& os ) const
{
	string str;
	if ( id != SearchResults::getDefaultID () ) str = id;
	for ( SearchResultsProteinLinePtrVectorSizeType i = 0 ; i < srprotl.size () ; i++ ) {
		srprotl [i]->printDelimited ( os, str, reportUniqPeps );
	}
}
PairStringVectorStringVector SearchResultsProteinReport::getAccessionNumberListPair () const
{
	PairStringVectorStringVector anums;
	SetString ss;		// List without negative acc numbers
	SetString ss2;		// List with negative acc numbers
	for ( SearchResultsProteinLinePtrVectorSizeType i = 0 ; i < srprotl.size () ; i++ ) {
		string a = srprotl [i]->getAccessionNumber ();
		if ( !a.empty () ) {
			PairIntString pis = ProteinInfo::getANumPair ( a );
			int idx = pis.first;
			string an = pis.second;
			if ( an [0] != '-' && ProteinInfo::getNonDecoyFlag ( idx ) ) ss.insert ( a );
			ss2.insert ( a );
		}
	}
	for ( SetStringConstIterator j = ss.begin () ; j != ss.end () ; j++ ) {
		anums.first.push_back ( *j );
	}
	for ( SetStringConstIterator k = ss2.begin () ; k != ss2.end () ; k++ ) {
		anums.second.push_back ( *k );
	}
	return anums;
}

static inline bool srlEqualTime ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b )
{
	return a->getSpecID () == b->getSpecID ();
}
static inline bool srlEqualTime2 ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b )
{
	if ( sresKeepTimeReplicates )
		return a->getSearchIndex ( 0 ) == b->getSearchIndex ( 0 ) && a->getSpecID () == b->getSpecID () && a->getHitSequence () == b->getHitSequence ();
	else
		return a->getSearchIndex ( 0 ) == b->getSearchIndex ( 0 ) && a->getSpecID () == b->getSpecID ();
}
static inline bool srlEqualTime3 ( const SearchResultsPeptideLine* a, const SearchResultsPeptideLine* b )
{
	if ( sresKeepTimeReplicates )
		return a->getSpecID () == b->getSpecID () && a->getHitSequence () == b->getHitSequence ();
	else
		return a->getSpecID () == b->getSpecID ();
}
static inline bool vsrphEqualTime3 ( const SearchResultsPeptideHit* a, const SearchResultsPeptideHit* b )
{
	return a->getSearchIndex () == b->getSearchIndex () && a->getSpecID () == b->getSpecID () && a->getSequencePlusMods () == b->getSequencePlusMods ();
}
class SortPeptidePositionANum {
	public:
		bool operator () ( const SearchResultsPeptideHit* lhs, const string& rhs ) const
		{
			return lhs->getAccessionNumber () < rhs;
		}
		bool operator () ( const string& lhs, const SearchResultsPeptideHit* rhs ) const
		{
			return lhs < rhs->getAccessionNumber ();
		}
};
SearchResultsPeptideReport::SearchResultsPeptideReport ( const vector <SearchResults*>& sr, bool remove, const MapStringToStringVector& aNumList, const string& sortType, const string& sortType2, const string& reportHitsType, const string& reportHomologousProteins, const string& id ) :
	SearchResultsProteinReport ( sr, remove, aNumList, reportHitsType, reportHomologousProteins, id ),
	errorHistogram ( 0 ),
	errorFlag ( false )
{
	if ( sortType == "Time" || sortType == "Charge/M+H" || sortType2 == "Time" || sortType2 == "Charge/M+H" ) {
		ErrorHandler::genError ()->error ( "Sort type not currently supported for this type of report.\n" );
	}
	ujm->writeMessage ( cout, "Collating peptides" );
	formResultsLines ( sortType, sortType2 );

	if ( searchResults.size () == 1 ) {
		vector <SearchResultsPeptideLine*> pepLine = srpepl;

		stable_sort ( pepLine.begin (), pepLine.end (), sortPeptideReportByTime () );
		pepLine.erase ( unique ( pepLine.begin (), pepLine.end (), srlEqualTime ), pepLine.end () );
		doErrorHistogram ( pepLine );
		if ( PeptidePosition::getReportMModValue () ) doMModHistogram ( pepLine );
	}
	if ( sresXLinks ) {
		setQuanResultsCLink ();
		for ( int i = 0 ; i < searchResults.size () ; i++ ) {
			searchResults [i]->sortCLinkPeptideLines ( sortType );
		}
	}
	else
		setQuanResults ();
	ujm->writeMessage ( cout, "Finished collating peptides" );
}
// Time report
SearchResultsPeptideReport::SearchResultsPeptideReport ( const vector <SearchResults*>& sr, bool remove, const MapStringToStringVector& aNumList, const string& sortType, const string& sortType2, bool unmatchedSpectra, const string& reportHitsType, const string& reportHomologousProteins, const string& id, bool eraseNonUnique ) :
	SearchResultsProteinReport ( sr, remove, aNumList, reportHitsType, reportHomologousProteins, id ),
	errorHistogram ( 0 ),
	errorFlag ( false )
{
	if ( sortType == "Fraction/RT" || sortType2 == "Fraction/RT" ) {
		ErrorHandler::genError ()->error ( "Sort type not currently supported for this type of report.\n" );
	}
	ujm->writeMessage ( cout, "Collating peptides" );
	formResultsLines2 ( eraseNonUnique );
	setQuanResults ();
	if ( unmatchedSpectra && reportHitsType != "Protein Difference" && reportHitsType != "Peptide Difference" && reportHitsType != "Protein Intersection" && reportHitsType != "Peptide Intersection" && sortType != "Peptide Score" && sortType != "Discriminant Score" && sortType != "Expectation Value" && sortType != "Mass Mod" && sortType != "Error"  ) {
		addUnmatchedSpectra ();
	}
	if ( !sortType2.empty () ) sortPeptideTimesLines ( sortType2, srpepl.begin (), srpepl.end () );
	sortPeptideTimesLines ( sortType, srpepl.begin (), srpepl.end () );
	ujm->writeMessage ( cout, "Finished collating peptides" );
}
// Calibration report
SearchResultsPeptideReport::SearchResultsPeptideReport ( const vector <SearchResults*>& sr, const string& id ) :
	SearchResultsProteinReport ( sr, false, MapStringToStringVector (), "Union", "Interesting", id ),
	errorHistogram ( 0 ),
	errorFlag ( false )
{
	string lastANum;

	formResultsLines ( "", "" );

	stable_sort ( srpepl.begin (), srpepl.end (), sortPeptideReportByTime () );
	srpepl.erase ( unique ( srpepl.begin (), srpepl.end (), srlEqualTime ), srpepl.end () );
	int numFractions = searchResults [0]->getNumFractions ();
	DoubleVectorVector errors ( numFractions );
	for ( SearchResultsPeptideLinePtrVectorSizeType k = 0 ; k < srpepl.size () ; k++ ) {
		double err = srpepl [k]->getPeptideError ();
		if ( err != PeptidePosition::INVALID_ERROR ) errors [srpepl [k]->getFraction ()-1].push_back ( err );
	}
	mean.resize ( numFractions );
	sdev.resize ( numFractions );
	size.resize ( numFractions );
	for ( DoubleVectorVectorSizeType ii = 0 ; ii < errors.size () ; ii++ ) {	// For each fraction calculate mean error
		DoubleVector& e = errors [ii];
		size [ii] = e.size ();
		if ( size [ii] == 0 ) {
			mean [ii] = 0.0;
			sdev [ii] = 0.0;
		}
		else if ( size [ii] == 1 ) {	// Can't calculate stdev
			mean [ii] = e [0];
			sdev [ii] = 0.0;
		}
		else {
			double adev, svar, skew, curt;
			try {
				moment ( &e[0]-1, e.size (), &aveError, &adev, &sdevError, &svar, &skew, &curt );
			}
			catch ( lNrecMomentZeroVariance ) {}			// Not bothered
			catch ( lNrecMomentLessThanTwoDataValues ) {}	// Not bothered
			mean [ii] = aveError;
			sdev [ii] = sdevError;
		}
	}
	doErrorHistogram ( srpepl );
}
void SearchResultsPeptideReport::setQuanResultsCLink ()
{
	if ( QuantitationRatio::getQuanReport () ) {
		ujm->writeMessage ( cout, "Sorting peptides into fractions" );
		for ( SearchResultsPtrVectorSizeType i = 0 ; i < searchResults.size () ; i++ ) {		// For each compared search
			searchResults [i]->setQuanCLink ( i );
		}
		ujm->writeMessage ( cout, "Finished quantitation" );
	}
}
void SearchResultsPeptideReport::setQuanResults ()
{
	if ( QuantitationRatio::getQuanReport () ) {
		MapIntVectorConstPeptidePositionPtr ppv;
		ujm->writeMessage ( cout, "Sorting peptides into fractions" );
		if ( sresMergedFlag ) {
			for ( SearchResultsPeptideLinePtrVectorSizeType j = 0 ; j < srpepl.size () ; j++ ) {	// Collect peptide positions
				const PeptidePosition* pp = srpepl [j]->getPeptidePosition ( 0 );
				if ( pp ) ppv [pp->getSearchIndex ()].push_back ( pp );
			}
		}
		else {
			for ( SearchResultsPtrVectorSizeType i = 0 ; i < searchResults.size () ; i++ ) {		// For each compared search
				for ( SearchResultsPeptideLinePtrVectorSizeType j = 0 ; j < srpepl.size () ; j++ ) {// Collect peptide positions
					const PeptidePosition* pp = srpepl [j]->getPeptidePosition ( i );
					if ( pp ) ppv [i].push_back ( pp );
				}
			}
		}
		bool quanMSMSFlag = PeptidePositionQuan::getQuanMSMSFlag ();
		for ( MapIntVectorConstPeptidePositionPtrConstIterator ii = ppv.begin () ; ii != ppv.end () ; ii++ ) {		// For each compared search
			int searchIndex = (*ii).first;
			vector <const PeptidePosition*> pps = (*ii).second;
			stable_sort ( pps.begin (), pps.end (), SortPeptidePositionAscending3 () );			// Sort by fraction and spot
			int currentFraction = -1;
			int last = pps.size () - 1;
			int num;
			GenElapsedTime et;
			SetPairStringInt spsi;
			for ( std::vector <const PeptidePosition*>::size_type k = 0 ; k < pps.size () ; k++ ) {
				int fraction = pps [k]->getFraction ();
				if ( fraction != currentFraction ) {
					ujm->writeMessage ( cout, "Processing quantitation results for fraction " + gen_itoa ( fraction ) );
					if ( currentFraction != -1 ) PeptidePositionQuan::deleteQuan ( searchIndex );
					PeptidePositionQuan::initialiseQuan ( searchIndex, fraction );							// Init new fraction
					currentFraction = fraction;
					num = 0;
					spsi.clear ();
				}
				num++;
				if ( et.getElapsedTime () > 10 ) {
					ujm->writeMessage ( cout, "Processing quantitation results for fraction " + gen_itoa ( fraction ) + ", " + gen_itoa ( num ) + " peptides" );
					et.reset ();
				}
				PairSetPairStringIntIteratorBool pspsiib = spsi.insert ( make_pair ( pps [k]->getAccessionNumber () + QuantitationRatio::getUnmodifiedPeptide ( pps [k]->getSequence () ), pps [k]->getCharge () ) );
				if ( quanMSMSFlag || pspsiib.second ) {
					pps [k]->setQuanResults ();
				}
				if ( k == last ) PeptidePositionQuan::deleteQuan ( searchIndex );
			}
		}
		for ( SearchResultsPeptideLinePtrVectorSizeType x = 0 ; x < srpepl.size () ; x++ ) {	// Set protein stats
			if ( x == 0 || ( srpepl [x]->getFullAccessionNumber () != srpepl [x-1]->getFullAccessionNumber () ) ) {
				if ( QuantitationRatio::getAreaRatioReport () ) quanProteinStats ( x, true );
				if ( QuantitationRatio::getIntRatioReport () ) quanProteinStats ( x, false );
			}
		}
		ujm->writeMessage ( cout, "Finished quantitation" );
	}
}
void SearchResultsPeptideReport::formResultsLines ( const string& sortType, const string& sortType2 )
{
	static vector <SearchResultsPeptideHit*> srph;
	if ( srph.empty () ) {
		for ( int i = 0 ; i < searchResults.size () ; i++ ) {
			srph.push_back ( new SearchResultsPeptideHit ( 0, 0, i ) );
		}
	}
	for ( SearchResultsProteinLinePtrVectorSizeType ii = 0 ; ii < srprotl.size () ; ii++ ) {
		string an = srprotl [ii]->getAccessionNumber ();
		ProteinInfo proteinInfo ( an );
		PairSearchResultsPeptideHitPtrVectorConstIterator jphi = SearchResults::getSearchPeptidePositionIters ( an, id );
		for ( SearchResultsPeptideHitPtrVectorConstIterator jj = jphi.first ; jj != jphi.second ; jj++ ) {
			const SearchResultsProteinInfo* jpi = SearchResults::getSearchResultsProteinInfo2 ( an, id );
			vector <const SearchResultsPeptideHit*> sri;
			ConstSearchResultsProteinInfoPtrVector srpi;
			if ( sresMergedFlag ) {
				sri.push_back ( (*jj) );
				srpi.push_back ( jpi );
			}
			else {
				for ( SearchResultsPtrVectorSizeType kk = 0 ; kk < searchResults.size () ; kk++ ) {
					const SearchResultsPeptideHit* pppsrpi = searchResults [kk]->getPPPeptideHitInfoPair ( (*jj), id );

					if ( pppsrpi == 0 )	pppsrpi = srph [kk];
					sri.push_back ( pppsrpi );
					srpi.push_back ( searchResults [kk]->getSearchResultsProteinInfo ( an, id ) );
				}
			}
			srpepl.push_back ( new SearchResultsPeptideLine ( proteinInfo, sri, srpi ) );
		}
		if ( !sortType2.empty () )	sortPeptideLines ( sortType2, srpepl.end () - (jphi.second - jphi.first), srpepl.end () );
		if ( !sortType.empty () )	sortPeptideLines ( sortType, srpepl.end () - (jphi.second - jphi.first), srpepl.end () );
	}
	ProteinInfo::resetANumMap ();
}
void SearchResultsPeptideReport::formResultsLines2 ( bool eraseNonUnique )
{
	int numSearches = searchResults.size ();
	vector <const SearchResultsPeptideHit*> vsrph;
	for ( SearchResultsProteinLinePtrVectorSizeType ii = 0 ; ii < srprotl.size () ; ii++ ) {
		string an = srprotl [ii]->getAccessionNumber ();
		PairSearchResultsPeptideHitPtrVectorConstIterator jphi = SearchResults::getSearchPeptidePositionIters ( an, id );
		for ( SearchResultsPeptideHitPtrVectorConstIterator jj = jphi.first ; jj != jphi.second ; jj++ ) {
			if ( sresMergedFlag ) {
				vsrph.push_back ( (*jj) );
			}
			else {
				for ( SearchResultsPtrVectorSizeType kk = 0 ; kk < numSearches ; kk++ ) {
					const SearchResultsPeptideHit* pppsrpi = searchResults [kk]->getPPPeptideHitInfoPair ( (*jj), id );
					if ( pppsrpi != 0 )	vsrph.push_back ( pppsrpi );
				}
			}
		}
	}
	if ( sresMergedFlag ) {
		ProteinInfo proteinInfo;
		string lastANum;
		stable_sort ( vsrph.begin (), vsrph.end (), SortPeptidePositionAscendingMerged () );
		const SearchResultsProteinInfo* jpi;
		for ( int xx = 0 ; xx < vsrph.size () ; xx++ ) {
			string an = vsrph [xx]->getAccessionNumber ();
			if ( an != lastANum ) {
				proteinInfo = ProteinInfo ( an );
				lastANum = an;
				jpi = SearchResults::getSearchResultsProteinInfo2 ( an, id );
			}
			vector <const SearchResultsPeptideHit*> sri;	// Empty
			ConstSearchResultsProteinInfoPtrVector srpi;	// Empty
			sri.push_back ( vsrph [xx] );
			srpi.push_back ( jpi );
			srpepl.push_back ( new SearchResultsPeptideLine ( proteinInfo, sri, srpi ) );
		}
		if ( eraseNonUnique ) {
			if ( sresSingleProject )	srpepl.erase ( unique ( srpepl.begin (), srpepl.end (), srlEqualTime3 ), srpepl.end () );
			else						srpepl.erase ( unique ( srpepl.begin (), srpepl.end (), srlEqualTime2 ), srpepl.end () );	
		}
	}
	else {
		stable_sort ( vsrph.begin (), vsrph.end (), SortPeptidePositionAscendingSeparated () );
		if ( eraseNonUnique && sresKeepTimeReplicates ) {
			vsrph.erase ( unique ( vsrph.begin (), vsrph.end (), vsrphEqualTime3 ), vsrph.end () );
		}
		string prevSpecID;
		vector <const SearchResultsPeptideHit*> visrph;
		for ( int xx = 0 ; xx < vsrph.size () ; xx++ ) {
			bool last = ( xx == vsrph.size ()-1 );
			string specID = vsrph [xx]->getFullSpecID ();
			if ( last ) {
				if ( specID != prevSpecID && !prevSpecID.empty () ) {
					addSearchResultPeptideLines ( prevSpecID, visrph );
				}
				visrph.push_back ( vsrph [xx] );
			}
			if ( ( xx != 0 && specID != prevSpecID ) || last ) {
				if ( last ) {		// If the last acc number is different to the previous one it was being incorrectly set
					addSearchResultPeptideLines ( specID, visrph );
				}
				else
					addSearchResultPeptideLines ( prevSpecID, visrph );
			}
			if ( !last ) {
				visrph.push_back ( vsrph [xx] );
			}
			prevSpecID = specID;
		}
		if ( eraseNonUnique && !sresKeepTimeReplicates ) srpepl.erase ( unique ( srpepl.begin (), srpepl.end (), srlEqualTime3 ), srpepl.end () );
	}
	ProteinInfo::resetANumMap ();
}
void SearchResultsPeptideReport::addSearchResultPeptideLines ( const string& specID, vector <const SearchResultsPeptideHit*>& visrph )
{
	int numSearches = searchResults.size ();
	static ProteinInfo proteinInfo;
	static string lastANum;
	MSMSSpectrumInfo* m = searchResults [0]->getMSMSSpectrumInfo ( specID, id );
	SpecID* spid = new SpecID ( specID );
	for ( int i = 0 ; i < visrph.size () ; i++ ) {				// Iterate through hits
		const SearchResultsPeptideHit* srph = visrph [i];
		if ( srph ) {
			string an;
			int searchIndex = srph->getSearchIndex ();
			vector <const SearchResultsPeptideHit*> sri;
			ConstSearchResultsProteinInfoPtrVector srpi;
			for ( int j = 0 ; j < numSearches ; j++ ) {
				if ( j < searchIndex ) {
					sri.push_back ( new SearchResultsPeptideHit ( spid, m, j ) );
					srpi.push_back (searchResults [j]->getSearchResultsProteinInfo ());
				}
				else if ( j == searchIndex ) {
					sri.push_back (srph);
					srpi.push_back (searchResults [j]->getSearchResultsProteinInfo ( srph->getAccessionNumber (), id ));
					if ( an.empty () ) an = srph->getAccessionNumber ();
				}
				else {
					bool flag = false;
					for ( int k = i+1 ; k < visrph.size () ; k++ ) {
						const SearchResultsPeptideHit* srph2 = visrph [k];
						if ( srph2 && srph2->getSearchIndex () == j ) {
							if ( srph2->getSequencePlusMods () == srph->getSequencePlusMods () ) {
								sri.push_back (srph2);
								srpi.push_back (searchResults [j]->getSearchResultsProteinInfo ( srph2->getAccessionNumber (), id ));
								visrph [k] = 0;
								flag = true;
								break;
							}
						}
					}
					if ( !flag ) {
						sri.push_back ( new SearchResultsPeptideHit ( spid, m, j ) );
						srpi.push_back (searchResults [j]->getSearchResultsProteinInfo ());
					}
				}
			}
			if ( an != lastANum ) {
				proteinInfo = ProteinInfo ( an );
				lastANum = an;
			}
			srpepl.push_back ( new SearchResultsPeptideLine ( proteinInfo, sri, srpi ) );
		}
	}
	visrph.clear ();
}
void SearchResultsPeptideReport::doErrorHistogram ( vector <SearchResultsPeptideLine*>& pepLine )
{
	DoubleVector errors;
	for ( SearchResultsPeptideLinePtrVectorSizeType i = 0 ; i < pepLine.size () ; i++ ) {
		double err = pepLine [i]->getPeptideError ();
		double m = 0.0;
		if ( mean.size () ) m = mean [pepLine [i]->getFraction ()-1];
		if ( err != PeptidePosition::INVALID_ERROR ) errors.push_back ( err - m );
	}
	errorFlag = false;
	double adev, svar, skew, curt;

	if ( !errors.empty () ) {
		try {
			moment ( &errors[0]-1, errors.size (), &aveError, &adev, &sdevError, &svar, &skew, &curt );
		}
		catch ( lNrecMomentZeroVariance ) {}			// Not bothered
		catch ( lNrecMomentLessThanTwoDataValues ) {}	// Not bothered
		if ( errors.size () > 10 ) {
			errorFlag = true;
			if ( errors.size () > 50 ) {
				errorHistogram = new ErrorHistogram ();
				for ( DoubleVectorSizeType j = 0 ; j < errors.size () ; j++ ) {
					errorHistogram->add ( errors [j] );
				}
			}
		}
	}
}
void SearchResultsPeptideReport::doMModHistogram ( vector <SearchResultsPeptideLine*>& pepLine )
{
	MapDoubleToInt mModMap;
	for ( SearchResultsPeptideLinePtrVectorSizeType i = 0 ; i < pepLine.size () ; i++ ) {
		double mmod = pepLine [i]->getMModValue ( 0 );
		if ( mmod ) {
			MapDoubleToIntIterator cur = mModMap.find ( mmod );
			if ( cur == mModMap.end () )	mModMap [mmod] = 1;
			else							(*cur).second++;
		}
	}
	for ( MapDoubleToIntConstIterator j = mModMap.begin () ; j != mModMap.end () ; j++ ) {
		mModData.add ( j->first, j->second );
	}
}
void SearchResultsPeptideReport::addUnmatchedSpectra ()
{
	int num = sresSingleProject ? 1 : searchResults.size ();
	for ( int xx = 0 ; xx < num ; xx++ ) {						// If a single project is being searched the times are the same for each search
		StringVector newTimes;
		StringVector srlTimes;
		for ( SearchResultsPeptideLinePtrVectorSizeType x = 0 ; x < srpepl.size () ; x++ ) {	// Get the times that already have matches	
			srlTimes.push_back ( srpepl [x]->getFullSpecID () );
		}
		stable_sort ( srlTimes.begin (), srlTimes.end () );
		StringVector fileTimes;								// Get all the times
		for ( MapSpecIDMSMSSpectrumInfo::const_iterator y = searchResults [xx]->getMSMSSpectrumInfoBegin ( id ) ; y != searchResults [xx]->getMSMSSpectrumInfoEnd ( id ) ; y++ ) {
			fileTimes.push_back ( (*y).first );
		}
		stable_sort ( fileTimes.begin (), fileTimes.end () );
		set_difference ( fileTimes.begin (), fileTimes.end (), srlTimes.begin (), srlTimes.end (), back_inserter ( newTimes ) );
		for ( StringVectorSizeType i = 0 ; i < newTimes.size () ; i++ ) {	// The new times are the ones without initial matches
			vector <const SearchResultsPeptideHit*> sri;
			ConstSearchResultsProteinInfoPtrVector srpi;
			MSMSSpectrumInfo* m = searchResults [xx]->getMSMSSpectrumInfo ( newTimes [i], id );
			SpecID* spid = new SpecID ( newTimes [i] );
			if ( sresMergedFlag ) {
				sri.push_back ( new SearchResultsPeptideHit ( spid, m, xx ) );
				srpi.push_back ( searchResults [xx]->getSearchResultsProteinInfo () );
				srpepl.push_back ( new SearchResultsPeptideLine ( sri, srpi ) );
			}
			else {
				for ( SearchResultsPtrVectorSizeType j = 0 ; j < searchResults.size () ; j++ ) {
					const SearchResultsPeptideHit* srph = searchResults [j]->getPPPeptideHitInfoPair ();
					if ( srph == 0 ) {
						srph = new SearchResultsPeptideHit ( spid, m, j );
					}
					sri.push_back ( srph );
					srpi.push_back ( searchResults [j]->getSearchResultsProteinInfo () );
				}
				srpepl.push_back ( new SearchResultsPeptideLine ( sri, srpi ) );
			}
		}
	}
}
StringVector SearchResultsPeptideReport::getAccessionNumbers2 () const
{
	SetString ss;
	for ( SearchResultsPeptideLinePtrVectorSizeType i = 0 ; i < srpepl.size () ; i++ ) {
		ss.insert ( srpepl [i]->getAcc () );
	}
	StringVector sv;
	for ( SetStringConstIterator j = ss.begin () ; j != ss.end () ; j++ ) {
		sv.push_back ( *j );
	}
	return sv;
}
#ifdef MYSQL_DATABASE
void SearchResultsPeptideReport::printReportHeader ( ostream& os ) const
{
	ParameterList pList = *ProgramLink::getParams ();
	FormItemSCFormat fi3;
	FormItemSCReportType fi4 ( false, false, "Peptide" );
	FormItemSCCheckboxes fi5;

	printHTMLFORMStart ( os, "post", "searchCompare" );
	tableStart ( os, true );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				printHTMLFORMSubmit ( os, "Remove Selected Peptides" );
			tableHeaderEnd ( os );
		tableRowEnd ( os );

		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				fi3.printHTML ( os );
				fi4.printHTML ( os );
				fi5.printHTML ( os );
				CheckboxSettingJavascript csj ( "cb" );
				csj.print ( os );
				BoolDeque dq;
				for ( SearchResultsPeptideLinePtrVectorSizeType i = 0 ; i < srpepl.size () ; i++ ) {
					dq.push_back ( srpepl [i]->getCompositionFlag () );
				}
				CheckboxArraySettingJavascript casj ( os, "cb", "Compostion", dq );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
	tableEnd ( os );

	pList.removeName ( FormItemSCFormat::getName () );
	pList.removeName ( FormItemSCReportType::getName () );
	pList.removeName ( FormItemSCCheckboxes::getName () );
	pList.copyToHiddenFormEntry ( os );
	os << "<br />" << endl;
}
void SearchResultsPeptideReport::printReportFooter ( ostream& os ) const
{
	printHTMLFORMEnd ( os );
}
#endif
void SearchResultsPeptideReport::printHTML ( ostream& os ) const
{
	if ( sresProt && PPProteinHitQuanInfo::getQuan () ) {
		SearchResultsProteinReport::printHTML ( os );
	}
	else {
		printHTMLHeader ( os );
		if ( sresFPR ) expectationDensityPlot ( os );
		else {
			SCMSTagLink smtl;
			SResLink sresLink ( isRawForwarding (), sresXLinks );
			if ( SearchResultsProteinLine::getReportLinks () ) {
				ProteinInfo::initialiseAccessionNumberLink ( os );
				startJavascript ( os );
				smtl.printHTML ( os );
				sresLink.printHTML ( os );
				endJavascript ( os );
				if ( id != SearchResults::getDefaultID () && searchResults [0]->getMSData () && !srpepl.empty () && ProteinInfo::getNumDatabases () == 1 ) {
					os << "<p>" << endl;
					startJavascript ( os );
					SCMSBridgeLink smbl;
					smbl.printHTML ( os );
					endJavascript ( os );
					smbl.write ( os, id, PeptidePosition::getSearchKey ( 0 ), PeptidePosition::getParentTolerance ( 0 ), PeptidePosition::getSysErrorStr ( 0 ), getAccessionNumbers2 () );
					os << endl;
					os << "</p>" << endl;
				}
			}
			if ( sresTime ) printHTMLTimeTable ( os, smtl );
			else if ( sresXLinks ) printCrosslinksHTML ( os, sresLink );
			else printHTMLPeptideTables ( os, smtl, sresLink );
			os << "<br />" << endl;
		}
	}
}
void SearchResultsPeptideReport::printHTMLHeader ( ostream& os ) const
{
	printIDHTML ( os, id );
	if ( sresTime )	ParameterList::printHTML ( os, "Number of Data Spectra", searchResults [0]->getMSMSSpectrumInfoSize ( id ) ); 
	ParameterList::printHTML ( os, "Number of Proteins (After Homology Filtering)", fullAccList.first.size () );
	int numNegAccNums = fullAccList.second.size () - fullAccList.first.size ();
	if ( numNegAccNums ) {
		ParameterList::printHTML ( os, "Number of Decoy Proteins (After Homology Filtering)", numNegAccNums );
	}
	ParameterList::printHTML ( os, "Number of Peptides", srpepl.size () );
	for ( SearchResultsPtrVectorSizeType x = 0 ; x < searchResults.size () ; x++ ) {
		if ( !sresTime ) os << "Search " << x+1 << " Number of Data Spectra: <b>" << searchResults [x]->getMSMSSpectrumInfoSize ( id ) << "</b><br />" << endl; 
		if ( id == SearchResults::getDefaultID () ) {
			string searchEndTime = searchResults [x]->getSearchEndTime ();
			string searchTime = searchResults [x]->getSearchTime ();
			if ( !searchEndTime.empty () )	os << "Search " << x+1 << " End Time: <b>" << searchEndTime << "</b><br />" << endl;
			if ( !searchTime.empty () )		os << "Search " << x+1 << " Search Time: <b>" << searchTime << "</b><br />" << endl;
		}
		os << "Search " << x+1 << " Proteins (Pre Homology Filtering): <b>" << searchResults [x]->getNumProteins ( id ) << "</b> Peptides: <b>" << searchResults [x]->getNumPeptides ( id ) << "</b><br />" << endl; 
		searchResults [x]->printDatabaseInfoHTML ( os );
		PeptidePosition::printParamsHTML ( os, x );
	}
	os << "<br />" << endl;
	if ( errorHistogram ) errorHistogram->drawGraph ( os );
	if ( errorFlag ) {
		ParameterList::printHTML ( os, "Average Parent Error", aveError );
		ParameterList::printHTML ( os, "3 * StdDev", sdevError * 3.0 );
		os << "<br />" << endl;
	}
	if ( mModData.size () ) {
		os << "<p>" << endl;
		SpectrumGraph s1 ( "mmod_hist.par.txt" );
		ColouredGraphData gd ( mModData, 1 );
		s1.drawGraph ( os, gd, false );
		os << "</p>" << endl;
	}
}
void SearchResultsPeptideReport::printHTMLTimeTable ( ostream& os, const SCMSTagLink& smtl ) const
{
	int numDataLines = srpepl.size ();
	int pageNumber = 1;
	int rowsPerPage;
	if ( id != SearchResults::getDefaultID () ) {
		rowsPerPage = numDataLines;
	}
	else {
		rowsPerPage = 2000;
	}
	if ( rowsPerPage < numDataLines ) {
		const ParameterList* pList = ProgramLink::getParams ();
		pageNumber = pList->getIntValue ( "page", 1 );
		startJavascript ( os );
		ReportLinkProgram rlp ( "searchCompare", "scLink" );
		rlp.printHTML ( os, pList );
		endJavascript ( os );

		ReportLinks rLinks ( "", rowsPerPage, pageNumber, numDataLines );
		rLinks.printHTML ( os, rlp );
	}
	int startRow = (pageNumber-1) * rowsPerPage;
	int protInd = 0;
	for ( SearchResultsPeptideLinePtrVectorSizeType i = 0 ; i < srpepl.size () ; i++ ) {
		if ( i == 0 ) {
			tableStart ( os, true );
			srpepl [0]->printHTMLTimeHeader ( os, fullSearchNames );
		}
		bool empty = i != 0 && ( srpepl [i]->getSpecID () == srpepl [i-1]->getSpecID () );
		if ( !empty ) protInd++;
			
		if ( i >= startRow && i < startRow + rowsPerPage ) {
			srpepl [i]->printHTML ( os, protInd, empty, false, smtl );
		}
		if ( i == srpepl.size () - 1 ) tableEnd ( os );
	}
}
void SearchResultsPeptideReport::printCrosslinksHTML ( ostream& os, const SResLink& sresLink ) const
{
	os << "<p>" << endl;
	bool reportLinks = true;
	bool reportRepeats = PPPeptideHitInfo::getReportRepeats ();
	SearchResultsCrosslinkProteinHit::init ();
	SearchResultsCrosslinkPeptideHit::init ();
	PPPeptideHitInfo::setReportRepeats ( false );
	for ( int i = 0 ; i < searchResults.size () ; i++ ) {
		searchResults [i]->printCLinkProteinHitsLinks ( os, i );
		if ( i != searchResults.size () - 1 ) os << "<br />" << endl;
	}
	for ( int j = 0 ; j < searchResults.size () ; j++ ) {
		searchResults [j]->printCLinkProteinHitsHTML ( os, sresLink, j );
	}
	PPPeptideHitInfo::setReportRepeats ( reportRepeats );
	os << "</p>" << endl;
}
void SearchResultsPeptideReport::printCrosslinksDelimited ( ostream& os ) const
{
	bool reportRepeats = PPPeptideHitInfo::getReportRepeats ();
	SearchResultsCrosslinkProteinHit::init ();
	SearchResultsCrosslinkPeptideHit::init ();
	PPPeptideHitInfo::setReportRepeats ( false );
	for ( int i = 0 ; i < searchResults.size () ; i++ ) {
		searchResults [i]->printCLinkProteinHitsDelimited ( os, i );
	}
	PPPeptideHitInfo::setReportRepeats ( reportRepeats );
}
void SearchResultsPeptideReport::printHTMLPeptideTables ( ostream& os, const SCMSTagLink& smtl, const SResLink& sresLink ) const
{
	int protInd = 0;
	for ( SearchResultsPeptideLinePtrVectorSizeType i = 0 ; i < srpepl.size () ; i++ ) {
		if ( i == 0 || ( srpepl [i]->getFullAccessionNumber () != srpepl [i-1]->getFullAccessionNumber () ) ) {
			if ( SearchResultsProteinLine::getReportLinks () ) {
				StringVector idsv = srprotl [protInd++]->getIDStrVec ();
				for ( StringVectorSizeType j = 0 ; j < idsv.size () ; j++ ) {
					sresLink.write ( os, srpepl [i]->getFullAccessionNumber (), id, idsv [j] );
					if ( j != idsv.size () - 1 ) os << ",&nbsp;";
				}
			}
			else os << srprotl [protInd++]->getIDStrVecOutput ();
			os << " ";
			srpepl [i]->printProteinSequenceHTML ( os, fullAccList.second.size () == 1 );
			tableStart ( os, true );
				srpepl [0]->printHTMLProteinHeader ( os, fullSearchNames );
				srpepl [i]->printHTMLProtein ( os, fullSearchNames );
			tableEnd ( os );
			os << "<br />" << endl;
			if ( QuantitationRatio::getAreaRatioReport () ) quanPlot ( os, i, true );
			if ( QuantitationRatio::getIntRatioReport () ) quanPlot ( os, i, false );
			os << "<br />" << endl;
			tableStart ( os, true );
			srpepl [0]->printHTMLHeader ( os, fullSearchNames );
		}
		srpepl [i]->printHTML ( os, smtl, id );
		if ( i == srpepl.size () - 1 || ( srpepl [i]->getFullAccessionNumber () != srpepl [i+1]->getFullAccessionNumber () ) ) {
			tableEnd ( os );
			os << "<hr />" << endl;
		}
	}
}
void SearchResultsPeptideReport::printHTMLCalibration ( ostream& os ) const
{
	printIDHTML ( os, id );
	DoubleVector off;
	StringVector offStr;
	if ( errorHistogram ) errorHistogram->drawGraph ( os );
	tableStart ( os, true );
	tableRowStart ( os );
		tableHeader ( os, "Fraction" );
		tableHeader ( os, "Num Points" );
		tableHeader ( os, "Mean Error" );
		tableHeader ( os, "StDev Error" );
		tableHeader ( os, "3 StDev Error" );
	tableRowEnd ( os );
	for ( unsigned int i = 0 ; i < mean.size () ; i++ ) {
		tableRowStart ( os );
			tableCell ( os, i + 1 );
			tableCell ( os, size [i] );
			tableCellSigFig ( os, mean [i], 4 );
			tableCellSigFig ( os, sdev [i], 4 );
			tableCellSigFig ( os, sdev [i] * 3, 4 );
		tableRowEnd ( os );
	}
	tableEnd ( os );
	for ( unsigned int j = 0 ; j < mean.size () ; j++ ) {
		off.push_back ( mean [j] );
		ostringstream ostr;
		genPrintSigFig ( ostr, mean [j], 4 );
		offStr.push_back ( ostr.str () );
	}

	vector <FormItem*> formItems;
	formItems.push_back ( new FormItemText ( "", "", "form", 12, 20, "write_cal" ) );
	formItems.push_back ( new FormItemText ( "", "", "search_key", 12, 20, PeptidePosition::getParams0 ()->getStringValue ( "search_key" ) ) );
	formItems.push_back ( new FormItemText ( "", "", "msms_parent_mass_tolerance_units", 12, 20, PeptidePosition::getParams0 ()->getStringValue ( "msms_parent_mass_tolerance_units" ) ) );
	formItems.push_back ( new FormItemText ( "", "", "instrument_name", 12, 20, PeptidePosition::getParams0 ()->getStringValue ( "instrument_name" ) ) );
	formItems.push_back ( new FormItemTextArea ( "", "", "offsets", 3, 20, offStr ) );
	printHTMLFORMStart ( os, "post", "msform" );
	tableStart ( os, true, "center" );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				printHTMLFORMSubmit ( os, "Create Calibrated Project" );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
	tableEnd ( os );
	for ( int x = 0 ; x < formItems.size () ; x++ ) {
		formItems [x]->printHTMLHidden ( os );
	}
	printHTMLFORMEnd ( os );

	printHTMLFooter ( os, "2007" );
}
void SearchResultsPeptideReport::expectationDensityPlot ( ostream& os ) const
{
	if ( RPlot::getRFlag () ) {
		double infinity = std::numeric_limits<double>::infinity();
		std::vector <SearchResults*>::size_type num = sresMergedFlag ? 1 : searchResults.size ();
		RPlot rplot ( "expectationDensity.R" );
		GenOFStream ofs ( rplot.getDataFileFullPath () );
		bool ok = false;
		for ( SearchResultsPeptideLinePtrVectorSizeType i = 0 ; i < srpepl.size () ; i++ ) {
			if ( i == 0 || ( srpepl [i]->getSpecID () != srpepl [i-1]->getSpecID () ) ) {
				for ( SearchResultsPtrVectorSizeType j = 0 ; j < num ; j++ ) {
					double e = srpepl [i]->getExpectationValue ( j );
					if ( e > 0 && e != infinity ) {
						double n = -10.0 * log10 ( e );
						if ( n < 200 ) {
							genPrintSigFig ( ofs, n, 4 );
							ok = true;
						}
						else ofs << "NA";
					}
					else ofs << "NA";
					if ( searchResults [j]->isConcat () ) {
						ofs << '\t';
						if ( srpepl [i]->isDecoyHit ( j ) ) ofs << 'D';
						else ofs << 'N';
					}
					if ( j != num-1 ) ofs << '\t';
				}
				ofs << endl;
			}
		}
		ofs.close ();
		if ( ok ) {
			os << "<p>" << endl;
				rplot.printImage ( os, ".png" );
			os << "</p>" << endl;
			os << "<p>" << endl;
				rplot.printImage ( os, ".pdf" );
			os << "</p>" << endl;
		}
	}
}
void SearchResultsPeptideReport::quanPlot ( ostream& os, int n, bool area ) const
{
	if ( RPlot::getRFlag () ) {
		RPlot rplot ( "quan.R" );
		GenOFStream ofs ( rplot.getDataFileFullPath () );
		bool ok = false;
		string accNum = srpepl [n]->getFullAccessionNumber ();
		for ( SearchResultsPeptideLinePtrVectorSizeType i = n ; i < srpepl.size () ; i++ ) {
			if ( srpepl [i]->getFullAccessionNumber () != accNum ) break;
			if ( srpepl [i]->outputQuanResults ( ofs, fullSearchNames, area ) ) {
				ok = true;
			}
		}
		ofs.close ();
		if ( ok ) rplot.printImageAndLink ( os );
	}
}
void SearchResultsPeptideReport::quanProteinStats ( int n, bool area ) const
{
	string accNum = srpepl [n]->getFullAccessionNumber ();
	DoubleVectorVectorVector ratios;
	for ( SearchResultsPeptideLinePtrVectorSizeType i = n ; i < srpepl.size () ; i++ ) {
		if ( srpepl [i]->getFullAccessionNumber () != accNum ) break;
		ratios.push_back ( srpepl [i]->getRatios ( area ) );
	}
	int num = sresMergedFlag ? 1 : searchResults.size ();
	DoubleVectorVectorVector rat (num);
	static IntVector numQuanPksBySearch (num, 0);
	for ( DoubleVectorVectorVectorSizeType j = 0 ; j < ratios.size () ; j++ ) {		// Iterate through peptides
		for ( DoubleVectorVectorSizeType k = 0 ; k < num ; k++ ) {	// Iterate through searches
			DoubleVectorSizeType numQuanPeaks = ratios [j][k].size ();
			if ( numQuanPeaks > rat [k].size () ) {
				rat [k].resize ( numQuanPeaks );
				numQuanPksBySearch [k] = numQuanPeaks;
			}
			for ( DoubleVectorSizeType m = 0 ; m < numQuanPeaks ; m++ ) {	// Iterate through quan peaks
				double r = ratios [j][k][m];
				if ( r != 0.0 ) rat [k][m].push_back ( log10 ( r ) );	// k=search number, m=quan pk, j=peptides
			}
		}
	}
	for ( DoubleVectorVectorVectorSizeType x = 0 ; x < num ; x++ ) {			// Iterate through searches
		DoubleVector med ( numQuanPksBySearch [x], 0.0 );
		DoubleVector q1 ( numQuanPksBySearch [x], 0.0 );
		DoubleVector q2 ( numQuanPksBySearch [x], 0.0 );
		DoubleVector mean ( numQuanPksBySearch [x], 0.0 );
		DoubleVector stdev ( numQuanPksBySearch [x], 0.0 );
		IntVector numPks ( numQuanPksBySearch [x], 0 );
		for ( DoubleVectorVectorSizeType y = 0 ; y < rat [x].size () ; y++ ) {	// Iterate through quan peaks
			DoubleVector dv = rat [x][y];
			if ( !dv.empty () ) {
				stable_sort ( dv.begin (), dv.end () );
				med [y] = pow ( 10.0, median ( dv ) );
				if ( dv.size () > 1 ) q1 [y] = pow ( 10.0, lowerQ ( dv ) );
				if ( dv.size () > 1 ) q2 [y] = pow ( 10.0, upperQ ( dv ) );
				double av = dv [0];
				double sd = 0.0;
				double dummy;
				if ( dv.size () > 1 ) {
					try {
						moment ( &dv [0]-1, dv.size (), &av, &sd, &dummy );
					}
					catch ( lNrecMomentZeroVariance ) {}			// Not bothered
					catch ( lNrecMomentLessThanTwoDataValues ) {}	// Not bothered
				}
				mean [y] = av;
				stdev [y] = sd;
				numPks [y] = dv.size ();
			}
		}
		int idx = sresMergedFlag ? -1 : x;
		for ( SearchResultsPeptideLinePtrVectorSizeType z = n ; z < srpepl.size () ; z++ ) {
			if ( srpepl [z]->getFullAccessionNumber () != accNum ) break;
			if ( area )
				srpepl [z]->setAreaRatios ( idx, med, q1, q2, mean, stdev, numPks );
			else
				srpepl [z]->setIntensityRatios ( idx, med, q1, q2, mean, stdev, numPks );
		}
	}
}
bool SearchResultsPeptideReport::printDelimitedHeader ( ostream& os ) const
{
	if ( sresProt && PPProteinHitQuanInfo::getQuan () ) {
		return SearchResultsProteinReport::printDelimitedHeader ( os );
	}
	else {
		if ( !srpepl.empty () ) {
			srpepl [0]->printDelimitedHeader ( os, fullSearchNames, id != SearchResults::getDefaultID (), reportUniqPeps );
			return true;
		}
		else return false;	// Header not printed
	}
}
void SearchResultsPeptideReport::printMGF ( ostream& os, const string& outputDirectory ) const
{
	// Get project names, fraction names, search parameters (1 per project).
	StringVector projectNames;
	MapStringToStringVector fractionNamesMap;
	map <string, const ParameterList*> paramsMap;
	SetString projectNamesSet;
	for ( int ii = 0 ; ii < fullSearchNames.size () ; ii++ ) {
		string pName = fullSearchNames [ii].substr ( 0, fullSearchNames [ii].find ( '/' ) );
		projectNames.push_back ( pName );
		PairSetStringIteratorBool pssib = projectNamesSet.insert ( pName );
		if ( pssib.second ) {
			fractionNamesMap [pName] = searchResults [ii]->getFractionNames ();
			paramsMap [pName] = PeptidePosition::getParams ( ii );
		}
	}

	// Get a set of specID's for each project. These are automatically sorted.
	map <string, SetSpecID> specIDSets;
	if ( sresXLinks ) {
		for ( int i = 0 ; i < searchResults.size () ; i++ ) {
			const vector <SearchResultsCrosslinkProteinHit*>& cl = searchResults [i]->getCLinksHits ();
			for ( int j = 0 ; j < cl.size () ; j++ ) {
				cl [j]->addSpecIDs ( specIDSets [projectNames [i]] );
			}
		}
	}
	else {
		if ( sresMergedFlag ) {
			for ( SearchResultsPeptideLinePtrVectorSizeType i = 0 ; i < srpepl.size () ; i++ ) {	// Collect peptide positions
				const PeptidePosition* pp = srpepl [i]->getPeptidePosition ( 0 );
				if ( pp ) specIDSets [projectNames [pp->getSearchIndex ()]].insert ( pp->getSpecIDasID () );
			}
		}
		else {
			for ( SearchResultsPtrVectorSizeType i = 0 ; i < searchResults.size () ; i++ ) {		// For each compared search
				string proj = projectNames [i];
				for ( SearchResultsPeptideLinePtrVectorSizeType j = 0 ; j < srpepl.size () ; j++ ) {// Collect peptide positions
					const PeptidePosition* pp = srpepl [j]->getPeptidePosition ( i );
					if ( pp ) specIDSets [proj].insert ( pp->getSpecIDasID () );
				}
			}
		}
	}
	for ( std::map <string, SetSpecID>::const_iterator jj = specIDSets.begin () ; jj != specIDSets.end () ; jj++ ) {
		MSMSPeakListDataSetInfo* dsi = 0;
		string version;
		int fraction = -1;		// fraction not set
		StringVector actualPaths;

		PPTempFile pptf ( "", "" );
		string fullPath = outputDirectory.empty () ? pptf.getFullPath () : outputDirectory;
		genCreateDirectory ( fullPath );		// Creates a directory to dump the files
		string outputPath = pptf.getURL ();
		bool outputLink = outputDirectory.empty ();

		string projectName = (*jj).first;
		const SetSpecID& ssid = (*jj).second;
		MapStringToStringVectorConstIterator iter1 = fractionNamesMap.find ( projectName );
		const StringVector& fractionNames = (*iter1).second;
		std::map <std::string, const ParameterList*>::const_iterator iter2 = paramsMap.find ( projectName );
		const ParameterList* pl = (*iter2).second;
		string fSuffix = ".mgf";
		for ( SetSpecIDConstIterator i = ssid.begin () ; i != ssid.end () ; i++ ) {
			if ( (*i).getFraction () != fraction ) {
				if ( fraction != -1 ) delete dsi;
				fraction = (*i).getFraction ();
				string fn1 = getCentroidDataFilename ( pl, fraction );
				if ( isFileType ( fn1, "ms2" ) ) fSuffix = ".ms2";
				if ( isFileType ( fn1, "apl" ) ) fSuffix = ".apl";
				version = getProjectVersion ( pl, fraction );
				dsi = new MSMSPeakListDataSetInfo ( fn1 );
				actualPaths.push_back ( fullPath + SLASH + fractionNames [fraction-1] + fSuffix );
			}
			GenOFStream ost ( actualPaths.back (), std::ios_base::out | std::ios_base::app );
			try {
				dsi->writePeakList ( ost, *i, version );
			}
			catch ( runtime_error e ) {
				ErrorHandler::genError ()->error ( e );
			}
		}
		if ( !ssid.empty () ) {
			delete dsi;
			bool flag = gen7zaCreate ( fullPath + SLASH + projectName, fullPath + SLASH + "*" + fSuffix, "zip" );
			os << "<br /><br />" << endl;
			if ( flag && outputLink ) printArchiveLink ( os, projectName, outputPath );
			genUnlink ( actualPaths );
		}
	}
}
void SearchResultsPeptideReport::writeBiblioSpec ( ostream& os, const string& outputDirectory, const string& outputFilename, const string& id, bool norm ) const
{
	PPTempFile pptf ( "", "" );
	string fullPath = outputDirectory.empty () ? pptf.getFullPath () : outputDirectory;
	genCreateDirectory ( fullPath );
	string outputPath = pptf.getURL ();
	bool outputLink = outputDirectory.empty ();
	string outputName;
	if ( !outputFilename.empty () ) outputName = outputFilename;
	else							outputName = fullSearchNames [0].substr ( 0, fullSearchNames [0].find ( '/' ) );
	if ( id != "id" ) outputName += id;
	string actualPath = fullPath + SLASH + outputName + ".blib";
	writeBiblioSpecDB ( actualPath, outputName, norm );
	bool flag = gen7zaCreate ( fullPath + SLASH + outputName, actualPath, "zip" );
	genUnlink ( fullPath + SLASH + outputName + ".blib" );
	os << "<br /><br />" << endl;
	if ( flag && outputLink ) printArchiveLink ( os, outputName, outputPath );
}
void SearchResultsPeptideReport::writeBiblioSpecDB ( const string& actualPath, const string& outputName, bool norm ) const
{
	ujm->writeMessage ( cout, "Creating blib file" );
	BlibWrite blib ( actualPath );
	blib.create ();
	blib.insertScoreTypes ();
	MapBlibRefSpectra mbrs;
	MapStringToInt cdFilenameMap;
	StringVector cdFilenameList;
	StringVector cdVersionList;

	MapStringToPairIntDouble cdFirstRT;
	MapStringToPairIntDouble cdLastRT;
	if ( norm ) getRTRanges ( cdFirstRT, cdLastRT );
	double cfrt;
	double clrt;
	double cfrtNum = 0;
	double clrtNum = 0;
	int curFraction = -1;
	int curSearchIndex = -1;
	int sourceID = 0;
	for ( SearchResultsPeptideLinePtrVectorSizeType i = 0 ; i < srpepl.size () ; i++ ) {
		SearchResultsPeptideLine* s = srpepl [i];
		int fraction = s->getFraction ();
		int index = i+1;
		double pValue = s->getPValue ();
		double eValue = s->getExpectation ();
		string dbPeptide = s->getDBPeptide ();
		double mOverZ = s->getMOverZ ();
		string prevAA = s->getThePrevAA ();
		string nextAA = s->getTheNextAA ();
		int numPeaks = s->getNumPeaks ();
		int charge = s->getCharge ();
		double rt = s->getSpotAsNumber ();
		string msmsInfo = s->getMSMSInfo ();
		string specID = s->getSpecID ();
		string modPeptide = s->getBlibPeptide ();
		int searchIndex = s->getSearchIndex0 ();
		bool decoyFlag = s->isDecoyHit ();
		if ( searchIndex != curSearchIndex || fraction != curFraction ) {
			const ParameterList* pl = PeptidePosition::getParams ( searchIndex );
			curFraction = fraction;
			curSearchIndex = searchIndex;
			string filename = getCentroidDataFilename ( pl, fraction );
			string version = getProjectVersion ( pl, fraction );
			if ( norm ) {
				MapStringToPairIntDoubleConstIterator cur2 = cdFirstRT.find ( filename );
				cfrtNum = (*cur2).second.first;
				cfrt = (*cur2).second.second;
				MapStringToPairIntDoubleConstIterator cur3 = cdLastRT.find ( filename );
				clrtNum = (*cur3).second.first;
				clrt = (*cur3).second.second;
			}
			MapStringToIntConstIterator cur = cdFilenameMap.find ( filename );
			if ( cur != cdFilenameMap.end () ) {
				sourceID = (*cur).second;
			}
			else {
				sourceID++;
				cdFilenameMap [filename] = sourceID;
				blib.insertSpectrumSourceFiles ( filename );
				cdFilenameList.push_back ( filename );
				cdVersionList.push_back ( version );
			}
		}
		if ( !decoyFlag ) {
			if ( norm ) {
				static int minRTs = 50;
				if ( clrtNum < minRTs || cfrtNum < minRTs )
					rt = -1.0;
				else
					rt = (rt - cfrt)/(clrt - cfrt) * 100.0;
			}
			if ( rt != -1.0 ) {
				ostringstream ostr;
				ostr << sourceID;
				ostr << specID.substr ( specID.find ( '-' ) );
				specID = ostr.str ();
				PairStringInt brsek = make_pair ( modPeptide, charge );
				MapBlibRefSpectraIterator cur = mbrs.find ( brsek );
				if ( cur == mbrs.end () ) {
					VectorPairIntDouble vpid;
					s->getBlibMods ( vpid );
					mbrs [brsek] = BlibRefSpectraEntryValue ( dbPeptide, vpid, mOverZ, prevAA, nextAA, 1, numPeaks, sourceID, rt, index, s->getMSMSInfo (), specID, pValue, eValue );
				}
				else {
					BlibRefSpectraEntryValue& brsev = (*cur).second;
					brsev.incrementCopies ();
					if ( pValue < brsev.getPValue () )			// Better hit
						brsev.update ( mOverZ, prevAA, nextAA, numPeaks, sourceID, rt, index, msmsInfo, specID, pValue, eValue );
					else
						brsev.update ( sourceID, rt, index );
				}
			}
		}
	}
	typedef std::set <std::pair <SpecID, int> > SetPairSpecIDInt;
	typedef SetPairSpecIDInt::const_iterator SetPairSpecIDIntConstIterator;
	SetPairSpecIDInt specIDs;
	int idx = 0;
	blib.beginTransaction ();
	for ( MapBlibRefSpectraConstIterator j = mbrs.begin () ; j != mbrs.end () ; j++ ) {
		idx++;
		blib.insertRefSpectra ( (*j).first, (*j).second );
		blib.insertModifications ( idx, (*j).second );
		specIDs.insert ( make_pair ( (*j).second.getSpecID (), idx ) );
		blib.insertRetentionTimes ( idx, (*j).second );
	}
	blib.endTransaction ();
	int fraction = -1;
	MSMSPeakListDataSetInfo* dsi = 0;
	blib.beginTransaction ();
	string atotal = gen_itoa ( specIDs.size () );
	int num = 1;
	for ( SetPairSpecIDIntConstIterator k = specIDs.begin () ; k != specIDs.end () ; k++, num++ ) {
		if ( (*k).first.getFraction () != fraction ) {
			if ( fraction != -1 ) delete dsi;
			fraction = (*k).first.getFraction ();
			dsi = new MSMSPeakListDataSetInfo ( cdFilenameList [fraction-1] );
		}
		MSMSDataPointVector msmsDataPointList;
		try {
			dsi->getData ( msmsDataPointList, (*k).first, cdVersionList [fraction-1], ujm );
		}
		catch ( runtime_error e ) {
			ErrorHandler::genError ()->error ( e );
		}
		if ( !msmsDataPointList.empty () ) {
			blib.insertRefSpectraPeaks ( (*k).second, msmsDataPointList [0].getDataPeaks () );
		}
		if ( num % 100 == 0 ) ujm->writeMessage ( cout, "Processing spectra " + gen_itoa ( num ) + "/" + atotal );
	}
	blib.endTransaction ();
	blib.insertLibInfo ( idx, "msf.ucsf.edu", false, outputName );
	ujm->deletePreviousMessage ( cout );
}
void SearchResultsPeptideReport::getRTRanges ( MapStringToPairIntDouble& cdFirstRT, MapStringToPairIntDouble& cdLastRT ) const
{
	int curFraction = -1;
	int curSearchIndex = -1;
	MapPairIntIntToString fracName;
	PairIntInt idx;
	for ( SearchResultsPeptideLinePtrVectorSizeType i = 0 ; i < srpepl.size () ; i++ ) {
		SearchResultsPeptideLine* s = srpepl [i];
		int fraction = s->getFraction ();
		int searchIndex = s->getSearchIndex0 ();
		if ( searchIndex != curSearchIndex || fraction != curFraction ) {
			curFraction = fraction;
			curSearchIndex = searchIndex;
			const ParameterList* pl = PeptidePosition::getParams ( searchIndex );
			string filename = getCentroidDataFilename ( pl, fraction );
			idx = make_pair ( fraction, searchIndex );
			fracName [idx] = filename;
			cdFirstRT [filename].first = 0;
			cdFirstRT [filename].second = std::numeric_limits<int>::max();
			cdLastRT [filename].first = 0;
			cdLastRT [filename].second = std::numeric_limits<int>::min();
		}
	}
	curFraction = -1;
	curSearchIndex = -1;
	string curFname;
	double cfrt;
	double clrt;
	for ( SearchResultsPeptideLinePtrVectorSizeType j = 0 ; j < srpepl.size () ; j++ ) {
		SearchResultsPeptideLine* s = srpepl [j];
		double rt = s->getSpotAsNumber ();
		int fraction = s->getFraction ();
		int searchIndex = s->getSearchIndex0 ();
		bool decoyFlag = s->isDecoyHit ();
		if ( searchIndex != curSearchIndex || fraction != curFraction ) {
			curFraction = fraction;
			curSearchIndex = searchIndex;
			idx = make_pair ( fraction, searchIndex );
			MapPairIntIntToStringConstIterator cur = fracName.find ( idx );
			curFname = (*cur).second;
			MapStringToPairIntDoubleConstIterator cur2 = cdFirstRT.find ( curFname );
			cfrt = (*cur2).second.second;
			MapStringToPairIntDoubleConstIterator cur3 = cdLastRT.find ( curFname );
			clrt = (*cur3).second.second;
		}
		if ( !decoyFlag ) {
			cdFirstRT [curFname].first++;
			if ( rt < cfrt ) {
				cdFirstRT [curFname].second = rt;
				cfrt = rt;
			}
			cdLastRT [curFname].first++;
			if ( rt > clrt ) {
				cdLastRT [curFname].second = rt;
				clrt = rt;
			}
		}
	}
}
void SearchResultsPeptideReport::printMZIdentML ( ostream& os, const string& outputDirectory, const string& outputFilename ) const
{
	PPTempFile pptf ( "", "" );
	string fullPath = outputDirectory.empty () ? pptf.getFullPath () : outputDirectory;
	genCreateDirectory ( fullPath );
	string outputPath = pptf.getURL ();
	bool outputLink = outputDirectory.empty ();
	string outputName;
	if ( !outputFilename.empty () ) outputName = outputFilename;
	else							outputName = fullSearchNames [0].substr ( 0, fullSearchNames [0].find ( '/' ) );
	string actualPath = fullPath + SLASH + outputName + ".mzid";
	SCMZIdentMLReport scmzir;
	GenOFStream ost ( actualPath, std::ios_base::out );
	scmzir.printHeader ( ost );

	printMZIdentML_SequenceCollection ( ost );
	printMZIdentML_AnalysisCollection ( ost );
	printMZIdentML_AnalysisProtocolCollection ( ost );
	printMZIdentML_DataCollection ( ost );

	scmzir.printFooter ( ost );
	ost.close ();
	bool flag = gen7zaCreate ( fullPath + SLASH + outputName, actualPath, "zip" );
	genUnlink ( fullPath + SLASH + outputName + ".mzid" );
	os << "<br /><br />" << endl;
	if ( flag && outputLink ) printArchiveLink ( os, outputName, outputPath );
}
void SearchResultsPeptideReport::printMZIdentML_SequenceCollection ( ostream& ost ) const
{
	MZIdentML_SequenceCollection mzsc;
	mzsc.printOpenTag ( ost, 1 );
	for ( SearchResultsProteinLinePtrVectorSizeType i = 0 ; i < srprotl.size () ; i++ ) {
		SCMZIdentML_DBSequence dbSeq ( srprotl [i] );
		dbSeq.print ( ost, 2 );
	}
	SetString idSet;
	int idx = 1;
	int idx2 = 1;
	VectorXMLOutputItemPtr items;
	for ( SearchResultsPeptideLinePtrVectorSizeType j = 0 ; j < srpepl.size () ; j++ ) {
		SearchResultsPeptideLine* s = srpepl [j];
		const string& pepStr = s->getPrintedSequence ();
		PairSetStringIteratorBool flag = idSet.insert ( pepStr );
		string pepRef = "pep" + gen_itoa ( idx );
		string pepEvidRef = "pe"+gen_itoa ( idx2++ );
		items.push_back ( new MZIdentML_PeptideEvidence ( pepEvidRef,  s->getDatabaseMZIdentMLRef (), pepRef, s->getStartAA (), s->getEndAA (), s->getThePrevAA (), s->getTheNextAA (), s->isDecoyHit ( 0 ) ) );
		if ( flag.second ) {
			SCMZIdentML_Peptide pep ( s, pepRef );
			idx++;
			pep.print ( ost, 2 );
		}
	}
	for ( VectorXMLOutputItemPtrSizeType k = 0 ; k < items.size () ; k++ ) {
		items [k]->print ( ost, 2 );
	}
	mzsc.printCloseTag ( ost, 1 );
}
void SearchResultsPeptideReport::printMZIdentML_AnalysisCollection ( ostream& ost ) const
{
	MZIdentML_AnalysisCollection mzac;
	mzac.printOpenTag ( ost, 1 );
	mzac.printCloseTag ( ost, 1 );
}
void SearchResultsPeptideReport::printMZIdentML_AnalysisProtocolCollection ( ostream& ost ) const
{
	const ParameterList* pList = PeptidePosition::getParams0 ();
	MZIdentML_AnalysisProtocolCollection mzapc;
	mzapc.printOpenTag ( ost, 1 );
		VectorXMLOutputItemPtr subItems1;

		VectorXMLOutputItemPtr subItems1A;
		subItems1A.push_back ( new MZIdentML_CVParam ( "MS:1001083", "ms-ms search", "PSI-MS" ) );
		subItems1.push_back ( new MZIdentML_SearchType ( subItems1A ) );

		ParameterList pList2 = *(PeptidePosition::getParams0 ());
		pList2.removeName ( "msms_parent_mass_systematic_error" );
		pList2.removeName ( "msms_parent_mass_tolerance" );
		pList2.removeName ( "msms_parent_mass_tolerance_units" );
		pList2.removeName ( "fragment_masses_tolerance" );
		pList2.removeName ( "fragment_masses_tolerance_units" );
		pList2.removeName ( "msms_prot_low_mass" );
		pList2.removeName ( "msms_prot_high_mass" );
		pList2.removeName ( "msms_full_mw_range" );
		pList2.removeName ( "low_pi" );
		pList2.removeName ( "high_pi" );
		pList2.removeName ( "full_pi_range" );
		VectorXMLOutputItemPtr subItems1B = pList2.getMZIdentMLUserParameters ();
		subItems1.push_back ( new MZIdentML_AdditionalSearchParams ( subItems1B ) );

		VectorXMLOutputItemPtr subItems1C;
		subItems1.push_back ( new MZIdentML_ModificationParams ( subItems1C ) );

		VectorXMLOutputItemPtr subItems1D;
		subItems1.push_back ( new MZIdentML_Enzymes ( subItems1D ) );

		VectorXMLOutputItemPtr subItems1E;
		char aa [] = "ACDEFGHIKLMNPQRSTUVWY";
		//            123456789012345678901
		for ( int i = 0 ; i < aa [i] != '\0' ; i++ ) {
			char a = aa [i];
			double mm = a == 'U' ? 150.95363 : AAInfo::getInfo ().getMonoisotopicMass ( a );
			subItems1E.push_back ( new MZIdentML_Residue ( a, mm ) );
		}
		VectorXMLOutputItemPtr subItems1E1;
		subItems1E1.push_back ( new MZIdentML_CVParam_AlternativeSingleLetterCodes ( "D N" ) );
		subItems1E.push_back ( new MZIdentML_AmbiguousResidue ( subItems1E1, 'B' ) );

		VectorXMLOutputItemPtr subItems1E2;
		subItems1E2.push_back ( new MZIdentML_CVParam_AlternativeSingleLetterCodes ( "E Q" ) );
		subItems1E.push_back ( new MZIdentML_AmbiguousResidue ( subItems1E2, 'Z' ) );

		VectorXMLOutputItemPtr subItems1E3;
		subItems1E3.push_back ( new MZIdentML_CVParam_AlternativeSingleLetterCodes ( "A C D E F G H I K L M N O P Q R S T V W Y" ) );
		subItems1E.push_back ( new MZIdentML_AmbiguousResidue ( subItems1E3, 'X' ) );

		subItems1.push_back ( new MZIdentML_MassTable ( subItems1E, "MT", "1 2" ) );

		VectorXMLOutputItemPtr subItems1F;
		double parentTolerance = pList->getDoubleValue ( "msms_parent_mass_tolerance" );
		string parentToleranceUnits = pList->getStringValue ( "msms_parent_mass_tolerance_units" );
		double parentToleranceSysError = pList->getDoubleValue ( "msms_parent_mass_systematic_error" );
		subItems1F.push_back ( new MZIdentML_CVParam_Tolerance ( parentTolerance, parentToleranceUnits, parentToleranceSysError, "+" ) );
		subItems1F.push_back ( new MZIdentML_CVParam_Tolerance ( parentTolerance, parentToleranceUnits, parentToleranceSysError, "-" ) );
		subItems1.push_back ( new MZIdentML_FragmentTolerance ( subItems1F ) );

		VectorXMLOutputItemPtr subItems1G;
		double fragmentTolerance = pList->getDoubleValue ( "fragment_masses_tolerance" );
		string fragmentToleranceUnits = pList->getStringValue ( "fragment_masses_tolerance_units" );
		double fragmentToleranceSysError = 0.0;
		subItems1G.push_back ( new MZIdentML_CVParam_Tolerance ( fragmentTolerance, fragmentToleranceUnits, 0.0, "+" ) );
		subItems1G.push_back ( new MZIdentML_CVParam_Tolerance ( fragmentTolerance, fragmentToleranceUnits, 0.0, "-" ) );
		subItems1.push_back ( new MZIdentML_ParentTolerance ( subItems1G ) );

		VectorXMLOutputItemPtr subItems1H;
		subItems1.push_back ( new MZIdentML_Threshold ( subItems1H ) );

		VectorXMLOutputItemPtr subItems1I;
		StringVector species = pList->getStringVectorValue ( "species" );
		if ( species [0] != "All" ) {
			VectorXMLOutputItemPtr subItems1I1;
			VectorXMLOutputItemPtr subItems1I1A;
			subItems1I1A.push_back ( new MZIdentML_CVParam ( "MS:1001020", "DB filter taxonomy", "PSI-MS" ) );
			subItems1I1.push_back ( new MZIdentML_FilterType ( subItems1I1A ) );
			subItems1I.push_back ( new MZIdentML_Filter ( subItems1I1 ) );
		}
		bool fullMWRange = pList->getBoolValue ( "msms_full_mw_range" );
		if ( !fullMWRange ) {
			string msmsProtLowMass = pList->getStringValue ( "msms_prot_low_mass" );
			string msmsProtHighMass = pList->getStringValue ( "msms_prot_high_mass" );
			VectorXMLOutputItemPtr subItems1I2;
			VectorXMLOutputItemPtr subItems1I2A;
			subItems1I2A.push_back ( new MZIdentML_CVParam ( "MS:1001022", "DB MW filter", "PSI-MS" ) );
			subItems1I2.push_back ( new MZIdentML_FilterType ( subItems1I2A ) );
			VectorXMLOutputItemPtr subItems1I2B;
			subItems1I2B.push_back ( new MZIdentML_CVParam ( "MS:1001202", "DB MW filter minimum", "PSI-MS", msmsProtLowMass ) );
			subItems1I2B.push_back ( new MZIdentML_CVParam ( "MS:1001201", "DB MW filter maximum", "PSI-MS", msmsProtHighMass ) );
			subItems1I2.push_back ( new MZIdentML_Include ( subItems1I2B ) );
			subItems1I.push_back ( new MZIdentML_Filter ( subItems1I2 ) );
		}
		bool fullPIRange = pList->getBoolValue ( "full_pi_range" );
		if ( !fullPIRange ) {
			string lowPI = pList->getStringValue ( "low_pi" );
			string highPI = pList->getStringValue ( "high_pi" );
			VectorXMLOutputItemPtr subItems1I3;
			VectorXMLOutputItemPtr subItems1I3A;
			subItems1I3A.push_back ( new MZIdentML_CVParam ( "MS:1001023", "DB PI filter", "PSI-MS" ) );
			subItems1I3.push_back ( new MZIdentML_FilterType ( subItems1I3A ) );
			VectorXMLOutputItemPtr subItems1I3B;
			subItems1I3B.push_back ( new MZIdentML_CVParam ( "MS:1001204", "DB PI filter minimum", "PSI-MS", lowPI ) );
			subItems1I3B.push_back ( new MZIdentML_CVParam ( "MS:1001203", "DB PI filter maximum", "PSI-MS", highPI ) );
			subItems1I3.push_back ( new MZIdentML_Include ( subItems1I3B ) );
			subItems1I.push_back ( new MZIdentML_Filter ( subItems1I3 ) );
		}
		subItems1.push_back ( new MZIdentML_DatabaseFilters ( subItems1I ) );

		MZIdentML_SpectrumIdentificationProtocol mzsip ( subItems1, "sip1", "AS_ProteinProspector_Batch-Tag" );
		mzsip.print ( ost, 2 );

		VectorXMLOutputItemPtr subItems2;

		VectorXMLOutputItemPtr subItems2A;
		subItems2.push_back ( new MZIdentML_AnalysisParams ( subItems2A ) );

		VectorXMLOutputItemPtr subItems2B;
		subItems2.push_back ( new MZIdentML_Threshold ( subItems2B ) );

		MZIdentML_ProteinDetectionProtocol mzpdp ( subItems2, "pdp1", "AS_ProteinProspector_Search_Compare" );
		mzpdp.print ( ost, 2 );
	mzapc.printCloseTag ( ost, 1 );
}
void SearchResultsPeptideReport::printMZIdentML_DataCollection ( ostream& ost ) const
{
	MZIdentML_DataCollection mzdc;
	mzdc.printOpenTag ( ost, 1 );
		VectorXMLOutputItemPtr subItems1;

		for ( SearchResultsPtrVectorSizeType i = 0 ; i < searchResults.size () ; i++ ) {
			VectorXMLOutputItemPtr subItems1A;
			VectorXMLOutputItemPtr subItems1A1;
			subItems1A1.push_back ( new MZIdentML_CVParam ( "MS:XXXXXX", "PP Batch-Tag XML File", "PSI-MS" ) );
			subItems1A.push_back ( new MZIdentML_FileFormat ( subItems1A1 ) );
			subItems1.push_back ( new MZIdentML_SourceFile ( subItems1A, "SF_" + gen_itoa ( i+1 ), searchResults [i]->getFName () ) );
		}

		ProteinInfo::addDatabaseMZIdentMLInfo ( subItems1 );

		StringVector centroidPathList = ProjectFile::getCentroidPathList ( PeptidePosition::getParams0 () );
		for ( SearchResultsPtrVectorSizeType k = 0 ; k < searchResults.size () ; k++ ) {
			VectorXMLOutputItemPtr subItems1C;
			VectorXMLOutputItemPtr subItems1C1;
			subItems1C1.push_back ( new MZIdentML_CVParam ( "MS:1001062", "Mascot MGF file", "PSI-MS" ) );
			subItems1C.push_back ( new MZIdentML_FileFormat ( subItems1C1 ) );
			VectorXMLOutputItemPtr subItems1C2;
			subItems1C2.push_back ( new MZIdentML_CVParam ( "MS:1001528", "Mascot query number", "PSI-MS" ) );
			subItems1C.push_back ( new MZIdentML_SpectrumIDFormat ( subItems1C2 ) );
			subItems1.push_back ( new MZIdentML_SpectraData ( subItems1C, "SD_" + gen_itoa ( k+1 ), centroidPathList [k] ) );
		}

		MZIdentML_Inputs mzi ( subItems1 );
		mzi.print ( ost, 2 );

		VectorXMLOutputItemPtr subItems2;

		VectorXMLOutputItemPtr subItems2A;
		for ( SearchResultsPeptideLinePtrVectorSizeType m = 0 ; m < srpepl.size () ; m++ ) {
			SearchResultsPeptideLine* s = srpepl [m];
			VectorXMLOutputItemPtr subItems2A1;
			for ( int mm = 0 ; mm < 1 ; mm++ ) {
				VectorXMLOutputItemPtr subItems2A1A;
				for ( int mmm = 0 ; mmm < 1 ; mmm++ ) {
					subItems2A1A.push_back ( new MZIdentML_PeptideEvidenceRef ( "" ) );
				}
				subItems2A1.push_back ( new MZIdentML_SpectrumIdentificationItem ( subItems2A1A, "SII_1", s->getMOverZCalc (), s->getCharge (), s->getMOverZ (), true, s->getRank () ) );
			}
			subItems2.push_back ( new MZIdentML_SpectrumIdentificationResult ( subItems2A1, "SIR_1", "", "" ) );
		}
		subItems2.push_back ( new MZIdentML_SpectrumIdentificationList ( subItems2A, "SIL_1" ) );

		VectorXMLOutputItemPtr subItems2C;
		for ( SearchResultsProteinLinePtrVectorSizeType n = 0 ; n < srprotl.size () ; n++ ) {
		}
		subItems2.push_back ( new MZIdentML_ProteinDetectionList ( subItems2C, "id" ) );

		MZIdentML_AnalysisData mzad ( subItems2 );
		mzad.print ( ost, 2 );

	mzdc.printCloseTag ( ost, 1 );
}
void SearchResultsPeptideReport::printPepXML ( ostream& os, const string& outputDirectory, const string& outputFilename ) const
{
	PPTempFile pptf ( "", "" );
	string fullPath = outputDirectory.empty () ? pptf.getFullPath () : outputDirectory;
	genCreateDirectory ( fullPath );
	string outputPath = pptf.getURL ();
	bool outputLink = outputDirectory.empty ();

	SCPepXMLReport scpcr;

	string fractionName;
	string actualPath;
	int num = 1;
	SCPepXMLReport2* rep2 = 0;
	string specID;
	for ( SearchResultsPeptideLinePtrVectorSizeType i = 0 ; i < srpepl.size () ; i++ ) {

		bool last = ( i == srpepl.size () - 1 );
		SearchResultsPeptideLine* s = srpepl [i];
		if ( s->getFractionName () != fractionName ) {
			if ( !fractionName.empty () ) {
				GenOFStream ost ( actualPath, std::ios_base::out | std::ios_base::app );
				scpcr.printFooter ( ost );
				if ( outputLink ) printLink ( os, fractionName, outputPath, ".pep.xml" );
			}
			fractionName = s->getFractionName ();
			actualPath = fullPath + SLASH + fractionName + ".pep.xml";
			scpcr.updateFractionName ( fractionName );
			GenOFStream ost ( actualPath, std::ios_base::out );
			scpcr.printHeader ( ost );
			num = 1;
		}
		GenOFStream ost ( actualPath, std::ios_base::out | std::ios_base::app );
		if ( s->getSpecID () != specID ) {
			specID = s->getSpecID ();
			if ( rep2 != 0 ) {
				rep2->print ( ost );
				delete rep2;
			}
			rep2 = new SCPepXMLReport2 ( s, num++ );
		}
		else
			rep2->add ( s );
		if ( last ) {
			rep2->print ( ost );
			delete rep2;
		}
	}
	if ( !fractionName.empty () ) {
		GenOFStream ost ( actualPath, std::ios_base::out | std::ios_base::app );
		scpcr.printFooter ( ost );
		if ( outputLink ) printLink ( os, fractionName, outputPath, ".pep.xml" );
	}
	if ( !fractionName.empty () ) {
		os << "<br /><br />" << endl;
		string outputName;
		if ( !outputFilename.empty () ) outputName = outputFilename;
		else							outputName = fullSearchNames [0].substr ( 0, fullSearchNames [0].find ( '/' ) );
		bool flag = gen7zaCreate ( fullPath + SLASH + outputName, fullPath + SLASH + "*.pep.xml", "zip" );
		os << "<br /><br />" << endl;
		if ( flag && outputLink ) printArchiveLink ( os, outputName, outputPath );
	}
}
void SearchResultsPeptideReport::printLink ( ostream& os, const string& fraction, const string& outputPath, const string suffix ) const
{
	string filename = fraction + suffix;
	os << "<a href=\"";
	os << outputPath;
	os << "/";
	os << filename;
	os << "\">";
	os << "Fraction";
	os << " (" << filename << ")";
	os << "</a>";
	os << "<br />";
	os << endl;
}
void SearchResultsPeptideReport::printArchiveLink ( ostream& os, const string& projectName, const string& outputPath ) const
{
	string filename = projectName + ".zip";
	os << "<a href=\"";
	os << outputPath;
	os << "/";
	os << filename;
	os << "\">";
	if ( sresViewer )	os << "Peak List File";
	else				os << "Archive file (" << filename << ")";
	os << "</a>";
	os << "<br />";
	os << endl;
}
void SearchResultsPeptideReport::printDelimited ( ostream& os ) const
{
	if ( sresProt && PPProteinHitQuanInfo::getQuan () ) {
		SearchResultsProteinReport::printDelimited ( os );
	}
	else {
		string str;
		if ( id != SearchResults::getDefaultID () ) str = id;
		int protInd = 0;
		string lineStr;
		int numHomology;
		for ( SearchResultsPeptideLinePtrVectorSizeType i = 0 ; i < srpepl.size () ; i++ ) {
			if ( sresTime ) {
				lineStr = gen_itoa ( i + 1 );
			}
			else {
				if ( i == 0 || ( srpepl [i]->getFullAccessionNumber () != srpepl [i-1]->getFullAccessionNumber () ) ) {
					lineStr = srprotl [protInd]->getIDStrVecOutput ();
					numHomology = srprotl [protInd]->getNumHomology ();
					protInd++;
				}
			}
			srpepl [i]->printDelimited ( os, lineStr, numHomology, str, reportUniqPeps );
		}
	}
}
void SearchResultsPeptideReport::sortPeptideLines ( const string& sortType, const SearchResultsPeptideLinePtrVectorIterator& begin, const SearchResultsPeptideLinePtrVectorIterator& end )
{
	if		( sortType == "m/z" )				stable_sort ( begin, end, sortPeptideReportByMOverZ () );
	else if	( sortType == "M+H" )				stable_sort ( begin, end, sortPeptideReportByMPlusH () );
	else if	( sortType == "Error" )				stable_sort ( begin, end, sortPeptideReportByError () );
	else if	( sortType == "Intensity" )			stable_sort ( begin, end, sortPeptideReportByIntensity () );
	else if	( sortType == "Fraction/RT" )		stable_sort ( begin, end, sortPeptideReportBySpot () );
	else if	( sortType == "RT" )				stable_sort ( begin, end, sortPeptideReportByRT () );
	else if	( sortType == "Start Residue" )		stable_sort ( begin, end, sortPeptideReportByStartResidue () );
	else if	( sortType == "End Residue" )		stable_sort ( begin, end, sortPeptideReportByEndResidue () );
	else if	( sortType == "Peptide Score" )		stable_sort ( begin, end, sortPeptideReportByPeptideScore () );
	else if	( sortType == "Discriminant Score" )stable_sort ( begin, end, sortPeptideReportByDiscriminantScore () );
	else if	( sortType == "Expectation Value" )	stable_sort ( begin, end, sortPeptideReportByExpectationValue () );
	else if	( sortType == "Mass Mod" )			stable_sort ( begin, end, sortPeptideReportByMassMod () );
	else if	( sortType == "Crosslink AA" && sresXLinks )	{} // Do nothing but valid option
	else ErrorHandler::genError ()->error ( "Invalid sort type for peptides.\n" );
}
void SearchResultsPeptideReport::sortPeptideTimesLines ( const string& sortType, const SearchResultsPeptideLinePtrVectorIterator& begin, const SearchResultsPeptideLinePtrVectorIterator& end )
{
	if		( sortType == "Time" )				stable_sort ( begin, end, sortPeptideReportByTime () );
	else if	( sortType == "Error" )				stable_sort ( begin, end, sortPeptideTimeReportByError () );
	else if	( sortType == "RT" )				stable_sort ( begin, end, sortPeptideReportByRT () );
	else if	( sortType == "m/z" )				stable_sort ( begin, end, sortPeptideTimeReportByMOverZ () );
	else if	( sortType == "M+H" )				stable_sort ( begin, end, sortPeptideTimeReportByMPlusH () );
	else if	( sortType == "Charge/M+H" )		stable_sort ( begin, end, sortPeptideTimeReportByChargeAndMPlusH () );
	else if	( sortType == "Intensity" )			stable_sort ( begin, end, sortPeptideTimeReportByIntensity () );
	else if	( sortType == "Start Residue" )		stable_sort ( begin, end, sortPeptideReportByStartResidue () );
	else if	( sortType == "End Residue" )		stable_sort ( begin, end, sortPeptideReportByEndResidue () );
	else if	( sortType == "Peptide Score" )		stable_sort ( begin, end, sortPeptideReportByPeptideScore () );
	else if	( sortType == "Discriminant Score" )stable_sort ( begin, end, sortPeptideReportByDiscriminantScore () );
	else if	( sortType == "Expectation Value" )	stable_sort ( begin, end, sortPeptideReportByExpectationValue () );
	else if	( sortType == "Mass Mod" )			stable_sort ( begin, end, sortPeptideReportByMassMod () );
	else ErrorHandler::genError ()->error ( "Invalid sort type for peptides.\n" );
}
