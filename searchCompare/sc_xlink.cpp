/******************************************************************************
*                                                                             *
*  Program    : searchCompare                                                 *
*                                                                             *
*  Filename   : sc_xlink.cpp                                                  *
*                                                                             *
*  Created    : July 9th 2012                                                 *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2012-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef VIS_C
#pragma warning( disable : 4503 )	// Disable warnings about decorated name length being too long
#endif

#include <lg_stdlib.h>
#include <lgen_file.h>
#include <lu_cgi_val.h>
#include <lu_get_link.h>
#include <lu_delim.h>
#include <lu_getfil.h>
#include <lu_html.h>
#include <lu_param_list.h>
#include <lu_prod_par.h>
#include <lu_r_plot.h>
#include <lu_table.h>
#include <lu_sctag_link.h>
#include <lu_quan_ratio.h>
#include <sc_xlink.h>
#include <sc_quan.h>
#include <sc_sres_link.h>

using std::string;
using std::vector;
using std::ostream;
using std::ostringstream;
using std::endl;
using std::stable_sort;

StringVector SearchResultsCrosslinkPeptideHit::instrument;
vector <Tolerance*> SearchResultsCrosslinkPeptideHit::parentTolerances;
vector <Tolerance*> SearchResultsCrosslinkPeptideHit::fragmentTolerances;
BoolDeque SearchResultsCrosslinkPeptideHit::spottingPlates;
BoolDeque SearchResultsCrosslinkPeptideHit::spectrumNumber;
StringVectorVector SearchResultsCrosslinkPeptideHit::fractionNames;
bool SearchResultsCrosslinkPeptideHit::multipleFractionNames = false;
StringVectorVector SearchResultsCrosslinkPeptideHit::rawTypes;
StringVector SearchResultsCrosslinkPeptideHit::searchKey;
bool SearchResultsCrosslinkPeptideHit::reportMPlusH = false;
bool SearchResultsCrosslinkPeptideHit::reportMOverZ = false;
bool SearchResultsCrosslinkPeptideHit::reportCharge = false;
bool SearchResultsCrosslinkPeptideHit::reportIntensity = false;
bool SearchResultsCrosslinkPeptideHit::reportMSMSInfo = false;
bool SearchResultsCrosslinkPeptideHit::reportStartAA = false;
bool SearchResultsCrosslinkPeptideHit::reportEndAA = false;
bool SearchResultsCrosslinkPeptideHit::reportRepeats = false;
bool SearchResultsCrosslinkPeptideHit::reportLinks = false;
bool SearchResultsCrosslinkPeptideHit::reportError = false;
bool SearchResultsCrosslinkPeptideHit::reportTime = false;

bool SearchResultsCrosslinkPeptideHit::reportXLPeptide = false;
bool SearchResultsCrosslinkPeptideHit::reportXLScore = false;
bool SearchResultsCrosslinkPeptideHit::reportXLExpectation = false;
bool SearchResultsCrosslinkPeptideHit::reportXLMValue = false;
bool SearchResultsCrosslinkPeptideHit::reportXLRank = false;
bool SearchResultsCrosslinkPeptideHit::reportXLLowScore = false;
bool SearchResultsCrosslinkPeptideHit::reportXLLowEVal = false;
bool SearchResultsCrosslinkPeptideHit::reportXLLowMValue = false;
bool SearchResultsCrosslinkPeptideHit::reportXLAA = false;

int SearchResultsCrosslinkPeptideHit::reportPreviousAA = 0;
bool SearchResultsCrosslinkPeptideHit::reportDBPeptide = false;
int SearchResultsCrosslinkPeptideHit::reportNextAA = 0;
bool SearchResultsCrosslinkPeptideHit::extraInfo = false;

SearchResultsCrosslinkPeptideHit::SearchResultsCrosslinkPeptideHit ( const SpecID* specID, const MSMSSpectrumInfo* mmsi, const PeptideSpectralInfo* psi, double error, const HitPeptide* hitPeptide1, const HitPeptide* hitPeptide2, const XLinkDecoyInfo& xldi, int searchIndex, double xScore1, int xRank1, double xScore2, int xRank2, double xFirstScore ) :
	peptideHitInfo ( psi ),
	specID ( specID ),
	mmsi ( mmsi ),
	error ( error ),
	hitPeptide1 ( hitPeptide1 ),
	hitPeptide2 ( hitPeptide2 ),
	xldi ( xldi ),
	searchIndex ( searchIndex ),
	xScore1 ( xScore1 ),
	xRank1 ( xRank1 ),
	xScore2 ( xScore2 ),
	xRank2 ( xRank2 ),
	xFirstScore ( xFirstScore ),
	quanRatio ( 0 )
{
	if ( xRank1 != 0 ) extraInfo = true;
}
void SearchResultsCrosslinkPeptideHit::setQuanResults ( const LinkInfo* linkInfo ) const
{
	QuanPeptide qp ( hitPeptide1->getPeptide (), hitPeptide1->getNTerm (), hitPeptide1->getCTerm (), hitPeptide1->getNeutralLoss (), hitPeptide2->getPeptide (), hitPeptide2->getNTerm (), hitPeptide2->getCTerm (), linkInfo->getName () );
	quanRatio = PeptidePositionQuan::getQuanRatio ( getMOverZ (), getCharge (), qp, getSpecIDasID (), searchIndex );
}
void SearchResultsCrosslinkPeptideHit::printHeaderHTML ( ostream& os, const string& styleID, int searchIdx )
{
	tableRowStart ( os );
		if ( reportMPlusH )		tableHeader ( os, "M+H", styleID, "", false );
		if ( reportMOverZ )		tableHeader ( os, "m/z", styleID, "", false );
		if ( reportCharge )		tableHeader ( os, "z", styleID, "", false );
		if ( reportIntensity )	tableHeader ( os, "Intensity", styleID, "", false );
		if ( reportError )		tableHeader ( os, parentTolerances [searchIdx]->getUnitsString (), styleID );
		if ( reportXLPeptide )	tableHeader ( os, "Crosslinked Peptide", styleID );
		if ( reportTime ) {
			if ( fractionNames [searchIdx].size () > 1 ) tableHeader ( os, "Fraction", styleID );
			SpecID::printTableHeader ( os, spottingPlates [searchIdx], spottingPlates [searchIdx] || spectrumNumber [searchIdx], styleID );
		}
		if ( reportMSMSInfo ) SpecID::printTableMSMSInfoHeader ( os, styleID );
		PPPeptideHitInfo::printHeaderHTML ( os, styleID );
		if ( extraInfo ) {
			if ( reportXLScore )		tableHeader ( os, "MS-Tag<br />Score", styleID, "", false );
			if ( reportXLExpectation )	tableHeader ( os, "Exp", styleID, "", false );
			if ( reportXLMValue )		tableHeader ( os, "-10logP", styleID, "", false );
			if ( reportXLRank )			tableHeader ( os, "Rank", styleID, "", false );
			if ( reportXLLowScore )		tableHeader ( os, "Low<br />Score", styleID, "", false );
			if ( reportXLLowEVal )		tableHeader ( os, "Low<br />Exp", styleID, "", false );
			if ( reportXLLowMValue )	tableHeader ( os, "Low<br />-10logP", styleID, "", false );
		}
		PeptidePositionQuan::printHeaderHTML ( os, searchIdx, styleID );
		if ( reportStartAA )	tableHeader ( os, "Start", styleID );
		if ( reportEndAA )		tableHeader ( os, "End", styleID );
		if ( reportXLAA )		tableHeader ( os, "XLink AA", styleID );
		ProteinInfo::printHTMLANumHeader ( os, 0 );
		if ( reportRepeats )	tableHeader ( os, "# in DB", styleID );
	tableRowEnd ( os );
	tableEmptyRow ( os );
}
void SearchResultsCrosslinkPeptideHit::printProteinInfoHTML ( ostream& os, const SResLink& sresLink ) const
{
	xldi.printProteinInfoHTML ( os, sresLink, 0 );
}
void SearchResultsCrosslinkPeptideHit::printHeaderDelimited ( ostream& os, int searchIdx, const PPProteinHitQuanInfo& ppphqi ) const
{
	delimitedRowStart ( os );
		if ( reportMPlusH )		delimitedHeader ( os, "M+H" );
		if ( reportMOverZ )		delimitedHeader ( os, "m/z" );
		if ( reportCharge )		delimitedHeader ( os, "z" );
		if ( reportIntensity )	delimitedHeader ( os, "Intensity" );
		if ( reportError )	delimitedCell ( os, parentTolerances [searchIdx]->getUnitsString () );
		if ( reportPreviousAA )	delimitedHeader ( os, "Prev AA 1" );
		if ( reportDBPeptide )	delimitedHeader ( os, "DB Peptide 1" );
		if ( reportXLPeptide )	delimitedHeader ( os, "Peptide 1" );
		if ( reportNextAA )		delimitedHeader ( os, "Next AA 1" );
		if ( reportPreviousAA )	delimitedHeader ( os, "Prev AA 2" );
		if ( reportDBPeptide )	delimitedHeader ( os, "DB Peptide 2" );
		if ( reportXLPeptide )	delimitedHeader ( os, "Peptide 2" );
		if ( reportNextAA )		delimitedHeader ( os, "Next AA 2" );
		if ( reportTime ) {
			delimitedHeader ( os, "Fraction" );
			SpecID::printDelimitedHeader ( os, spottingPlates [searchIdx], spottingPlates [searchIdx] || spectrumNumber [searchIdx] );
		}
		if ( reportMSMSInfo )	SpecID::printDelimitedMSMSInfoHeader ( os );
		PPPeptideHitInfo::printHeaderDelimited ( os );

		if ( extraInfo ) {
			if ( reportXLScore )		delimitedHeader ( os, "Sc 1" );
			if ( reportXLExpectation )	delimitedHeader ( os, "Exp 1" );
			if ( reportXLMValue )		delimitedHeader ( os, "-10logP 1" );
			if ( reportXLRank )			delimitedHeader ( os, "Rk 1" );
			if ( reportXLScore )		delimitedHeader ( os, "Sc 2" );
			if ( reportXLExpectation )	delimitedHeader ( os, "Exp 2" );
			if ( reportXLMValue )		delimitedHeader ( os, "-10logP 2" );
			if ( reportXLRank )			delimitedHeader ( os, "Rk 2" );
			if ( reportXLLowScore )		delimitedHeader ( os, "L Sc" );
			if ( reportXLLowEVal )		delimitedHeader ( os, "L Exp" );
			if ( reportXLLowMValue )	delimitedHeader ( os, "L -10logP" );
		}
		if ( reportStartAA ) {
			delimitedHeader ( os, "Start 1" );
			delimitedHeader ( os, "Start 2" );
		}
		if ( reportEndAA ) {
			delimitedHeader ( os, "End 1" );
			delimitedHeader ( os, "End 2" );
		}
		PeptidePositionQuan::printHeaderDelimited ( os, searchIdx );
		ppphqi.printHeaderDelimited ( os );
		ProteinInfo::printDelimitedANumHeader ( os );
		if ( reportXLAA )	delimitedHeader ( os, "XLink AA 1" );
		if ( reportRepeats )XLinkDecoyInfo::printRepeatsDelimitedHeader1 ( os );
		ProteinInfo::printDelimitedANumHeader ( os );
		if ( reportXLAA )	delimitedHeader ( os, "XLink AA 2" );
		if ( reportRepeats )XLinkDecoyInfo::printRepeatsDelimitedHeader2 ( os );
		ProteinInfo::printDelimitedHeader ( os );
		ProteinInfo::printDelimitedHeader ( os );
	delimitedRowEnd ( os );
}
void SearchResultsCrosslinkPeptideHit::printHTML ( ostream& os, const string& styleID, const MSProductLink* productLink, const SCMSTagLink& smtl, const LinkInfo* linkInfo ) const
{
	tableRowStart ( os );
	MSParentLink* parentLink = 0;
	if ( reportLinks ) {
		if ( !rawTypes [searchIndex][getFraction ()-1].empty () ) parentLink = new MSParentLink;
		startJavascript ( os );
		if ( !rawTypes [searchIndex][getFraction ()-1].empty () ) parentLink->printHTML ( os );
		endJavascript ( os );
	}
	if ( reportMPlusH )		tableCell ( os, getMPlusH (), 4, false, styleID, 0, 2 );
	if ( reportMOverZ ) {
		if ( reportLinks && !rawTypes [searchIndex][getFraction ()-1].empty () ) {
			tableDataStart ( os, styleID, "", true, 0, 2 );
				string units = parentTolerances [searchIndex]->getUnitsString ();
				parentLink->write ( os, getSpecIDasID (), getMOverZ (), getCharge (), PeptidePosition::getRTIntervalStart (), PeptidePosition::getRTIntervalEnd (), PeakFitData::getSNRThreshold (), hitPeptide1->getPeptide (), hitPeptide1->getNTerm (), hitPeptide1->getCTerm (), hitPeptide1->getNeutralLoss (), hitPeptide2->getPeptide (), hitPeptide2->getNTerm (), hitPeptide2->getCTerm (), linkInfo->getName (), searchKey [searchIndex], PeptidePosition::getSysErrorStr ( searchIndex ), units );
			tableCellEnd ( os );
		}
		else
			tableCell ( os, getMOverZ (), 4, false, styleID, 0, 2 );
	}
	if ( reportCharge )		tableCell ( os, getCharge (), false, styleID, 0, 2 );
	if ( reportIntensity )	tableCellSigFig ( os, getIntensity (), 3, false, styleID, 0, 2 );
	if ( reportError )		tableCellSigFig ( os, error, 2, false, styleID, 0, 2 );
	if ( reportXLPeptide ) {
		tableHeaderStart ( os, styleID, "", true, 0, 2 );
			StringVector previousAA;
			previousAA.push_back ( getPreviousAA1 () );
			previousAA.push_back ( getPreviousAA2 () );
			StringVector nextAA;
			nextAA.push_back ( getNextAA1 () );
			nextAA.push_back ( getNextAA2 () );
			StringVector peptide;
			peptide.push_back ( hitPeptide1->getPeptide () );
			peptide.push_back ( hitPeptide2->getPeptide () );
			StringVector nt;
			nt.push_back ( hitPeptide1->getNTerm () );
			nt.push_back ( hitPeptide2->getNTerm () );
			StringVector ct;
			ct.push_back ( hitPeptide1->getCTerm () );
			ct.push_back ( hitPeptide2->getCTerm () );
			StringVector nl;
			nl.push_back ( hitPeptide1->getNeutralLoss () );
			nl.push_back ( hitPeptide2->getNeutralLoss () );
			productLink->write4 ( os, previousAA, peptide, nt, ct, nl, nextAA, !reportLinks, mmsi->getCharge (), linkInfo );
			if ( PeptidePosition::getRunMSProductFlag () ) runMSProduct ( searchIndex, peptideHitInfo.getScore (), linkInfo );
		tableHeaderEnd ( os );
	}
	if ( reportTime ) {
		if ( fractionNames [searchIndex].size () > 1 ) tableCell ( os, fractionNames [searchIndex][getFraction ()-1], false, false, styleID, 0, 2 );
		bool sPlates = spottingPlates [searchIndex];
		bool sNum = spectrumNumber [searchIndex];
		if ( reportLinks )
			specID->printTableCell2 ( os, sPlates, sPlates || sNum, styleID, smtl, searchKey [searchIndex], 0, 2 );
		else
			specID->printTableCell ( os, sPlates, sPlates || sNum, styleID, 0, 2 );
	}
	if ( reportMSMSInfo ) specID->printTableMSMSInfoCell ( os, styleID, 0, 2 );
	peptideHitInfo.printHTML ( os, styleID, 0, 2 );
	if ( extraInfo ) {
		double lowScore = peptideHitInfo.getScore ()-xFirstScore;
		if ( reportXLScore )		tableCell ( os, xScore1, 1, false, styleID );
		if ( reportXLExpectation )	tableCellSigFig ( os, peptideHitInfo.getExpectationValue ( xScore1 ), 2, false, styleID );
		if ( reportXLMValue )		tableCell ( os, peptideHitInfo.getMascotScore ( xScore1 ), 1, false, styleID );
		if ( reportXLRank )			tableCell ( os, xRank1, false, styleID );
		if ( reportXLLowScore )		tableCell ( os, lowScore, 1, false, styleID );
		if ( reportXLLowEVal )		tableCellSigFig ( os, peptideHitInfo.getExpectationValue ( lowScore ), 2, false, styleID );
		if ( reportXLLowMValue )	tableCell ( os, peptideHitInfo.getMascotScore ( lowScore ), 1, false, styleID );
	}
	if ( quanRatio == 0 )
		PeptidePositionQuan::printQuanBlankHTML ( os, searchIndex, styleID );
	else
		PeptidePositionQuan::printHTML ( os, quanRatio, styleID, false );
	if ( reportStartAA )tableCell ( os, getStartAA1 (), false, styleID );
	if ( reportEndAA )	tableCell ( os, getEndAA1 (), false, styleID );
	if ( reportXLAA )	tableCell ( os, getMModPosn1 (), false, styleID );
	xldi.printHTMLANum1 ( os, 0 );
	if ( reportRepeats )tableCell ( os, xldi.getRepeats1 ( 0 ), true, false, styleID );
	tableRowEnd ( os );
	tableRowStart ( os );
	if ( extraInfo ) {
		if ( reportXLScore )		tableCell ( os, xScore2, 1, false, styleID );
		if ( reportXLExpectation )	tableCellSigFig ( os, peptideHitInfo.getExpectationValue ( xScore2 ), 2, false, styleID );
		if ( reportXLMValue )		tableCell ( os, peptideHitInfo.getMascotScore ( xScore2 ), 1, false, styleID );
		if ( reportXLRank )			tableCell ( os, xRank2, false, styleID );
		if ( reportXLLowScore )		tableEmptyCell ( os, styleID );
		if ( reportXLLowEVal )		tableEmptyCell ( os, styleID );
		if ( reportXLLowMValue )	tableEmptyCell ( os, styleID );
	}
	PeptidePositionQuan::printQuanBlankHTML ( os, searchIndex, styleID );

	if ( reportStartAA )tableCell ( os, getStartAA2 (), false, styleID );
	if ( reportEndAA )	tableCell ( os, getEndAA2 (), false, styleID );
	if ( reportXLAA )	tableCell ( os, getMModPosn2 (), false, styleID );
	xldi.printHTMLANum2 ( os, 0 );
	if ( reportRepeats )tableCell ( os, xldi.getRepeats2 ( 0 ), true, false, styleID );
	tableRowEnd ( os );
	tableEmptyRow ( os );
}
void SearchResultsCrosslinkPeptideHit::runMSProduct ( int i, double score, const LinkInfo* linkInfo ) const
{
	static int idx = 0;
	idx++;
	string scoreStr = gen_ftoa ( score, "%.1f" );
	string outputFilename =  gen_itoa ( idx ) + "_" + stripFilenameChars ( hitPeptide1->getPeptide () + "---" + hitPeptide2->getPeptide () ) + "_" + scoreStr + "_" + gen_itoa ( getCharge () );
	PPTempFile tempFile ( "msprod", ".txt" );
	string tempFileFullPath = tempFile.getFullPath ();

	string command;
#ifdef VIS_C
	command += "\"";
#endif
	command += getSystemCall ( "mssearch.cgi" );
	command += " - ";
	command += getCommandLineNVPair ( "search_name", "msproduct" );
	command += " ";
	command += getCommandLineNVPair ( "output_filename", outputFilename );
	command += " ";
	command += getCommandLineNVPair ( "results_to_file", "1" );
	command += " ";
	command += getCommandLineNVPair ( "output_type", "XML" );
	command += " ";
	command += getCommandLineNVPair ( "report_title", "MS-Product" );
	command += " ";
	command += getCommandLineNVPair ( "version", Version::instance ().getVersion () );
	command += " ";
	command += getCommandLineNVPair ( "data_source", "List of Files" );
	command += " ";
	command += getCommandLineNVPair ( "use_instrument_ion_types", "1" );
	command += " ";
	command += getCommandLineNVPair ( "search_key", searchKey [i] );
	command += " ";
	command += getCommandLineNVPair ( "instrument_name", instrument [i] );
	command += " ";
	command += parentTolerances [i]->getCommandLineNVPair ( "msms_parent_mass" );
	command += fragmentTolerances [i]->getCommandLineNVPair ( "fragment_masses" );
	command += getCommandLineNVPair ( "parent_mass_convert", "monoisotopic" );
	command += " ";
	command += specID->getCommandLineNVPair ();
	command += getCommandLineNVPair ( "max_charge", mmsi->getCharge () );
	command += " ";
	command += getCommandLineNVPair ( "count_pos_z", "Ignore Basic AA" );
	command += " ";
	command += getCommandLineNVPair ( "link_search_type", linkInfo->getName () );
	command += " ";
	command += hitPeptide1->getCommandLineNVPair ( 1 );
	command += hitPeptide2->getCommandLineNVPair ( 2 );
	command += MSMSPeakFilterOptions::getCommandLineNVPair ( ProgramLink::getParams () );
	command += "> ";
#ifdef VIS_C
	command += "\"" + tempFileFullPath + "\"";
#else
	command += tempFileFullPath;
#endif
#ifdef VIS_C
	command +=  "\"";
#endif
	genSystem ( command, "", true );
	genUnlink ( tempFileFullPath );
}
void SearchResultsCrosslinkPeptideHit::printDelimited ( ostream& os, const PPProteinHitQuanInfo& ppphqi ) const
{
	delimitedRowStart ( os );
	if ( reportMPlusH )		delimitedCell ( os, getMPlusH (), 4 );
	if ( reportMOverZ )		delimitedCell ( os, getMOverZ (), 4 );
	if ( reportCharge )		delimitedCell ( os, getCharge () );
	if ( reportIntensity )	delimitedCellSigFig ( os, getIntensity (), 3 );
	if ( reportError )	delimitedCellSigFig ( os, error, 2 );

	if ( reportPreviousAA ) delimitedCell ( os, getPreviousAA1 () );
	if ( reportDBPeptide )	delimitedCell ( os, getDBPeptide1 () );
	if ( reportXLPeptide ) {
		ostringstream ost;
		if ( !hitPeptide1->getNTerm ().empty () ) ost << hitPeptide1->getNTerm () << '-';
		ost << hitPeptide1->getPeptide ();
		if ( !hitPeptide1->getCTerm ().empty () ) ost << '-' << hitPeptide1->getCTerm ();
		if ( !hitPeptide1->getNeutralLoss ().empty () ) ost << '+' << hitPeptide1->getNeutralLoss ();
		delimitedCell ( os, ost.str () );
	}
	if ( reportNextAA )		delimitedCell ( os, getNextAA1 () );

	if ( reportPreviousAA ) delimitedCell ( os, getPreviousAA2 () );
	if ( reportDBPeptide )	delimitedCell ( os, getDBPeptide2 () );
	if ( reportXLPeptide ) {
		ostringstream ost2;
		if ( !hitPeptide2->getNTerm ().empty () ) ost2 << hitPeptide2->getNTerm () << '-';
		ost2 << hitPeptide2->getPeptide ();
		if ( !hitPeptide2->getCTerm ().empty () ) ost2 << '-' << hitPeptide2->getCTerm ();
		if ( !hitPeptide2->getNeutralLoss ().empty () ) ost2 << '+' << hitPeptide2->getNeutralLoss ();
		delimitedCell ( os, ost2.str () );
	}
	if ( reportNextAA )		delimitedCell ( os, getNextAA2 () );

	if ( reportTime ) {
		delimitedCell ( os, fractionNames [searchIndex][getFraction ()-1] );
		specID->printDelimitedCell ( os, spottingPlates [searchIndex], spottingPlates [searchIndex] || spectrumNumber [searchIndex] );
	}
	if ( reportMSMSInfo ) specID->printDelimitedMSMSInfoCell ( os );
	peptideHitInfo.printDelimited ( os );
	if ( extraInfo ) {
		double lowScore = peptideHitInfo.getScore ()-xFirstScore;
		if ( reportXLScore )		delimitedCell ( os, xScore1, 1 );
		if ( reportXLExpectation )	delimitedCellSigFig ( os, peptideHitInfo.getExpectationValue ( xScore1 ), 2 );
		if ( reportXLMValue )		delimitedCell ( os, peptideHitInfo.getMascotScore ( xScore1 ), 1 );
		if ( reportXLRank )			delimitedCell ( os, xRank1 );

		if ( reportXLScore )		delimitedCell ( os, xScore2, 1 );
		if ( reportXLExpectation )	delimitedCellSigFig ( os, peptideHitInfo.getExpectationValue ( xScore2 ), 2 );
		if ( reportXLMValue )		delimitedCell ( os, peptideHitInfo.getMascotScore ( xScore2 ), 1 );
		if ( reportXLRank )			delimitedCell ( os, xRank2 );

		if ( reportXLLowScore )		delimitedCell ( os, lowScore, 1 );
		if ( reportXLLowEVal )		delimitedCellSigFig ( os, peptideHitInfo.getExpectationValue ( lowScore ), 2 );
		if ( reportXLLowMValue )	delimitedCell ( os, peptideHitInfo.getMascotScore ( lowScore ), 1 );
	}
	if ( reportStartAA ) {
		delimitedCell ( os, getStartAA1 () );
		delimitedCell ( os, getStartAA2 () );
	}
	if ( reportEndAA ) {
		delimitedCell ( os, getEndAA1 () );
		delimitedCell ( os, getEndAA2 () );
	}
	if ( quanRatio == 0 )
		PeptidePositionQuan::printQuanBlankDelimited ( os, searchIndex );
	else
		PeptidePositionQuan::printDelimited ( os, quanRatio, false );
	ppphqi.printDelimited ( os );
	xldi.printDelimitedANum1 ( os, 0 );
	if ( reportXLAA )	delimitedCell ( os, getMModPosn1 () );
	if ( reportRepeats )xldi.printRepeatsDelimited1 ( os, 0 );
	xldi.printDelimitedANum2 ( os, 0 );
	if ( reportXLAA )	delimitedCell ( os, getMModPosn2 () );
	if ( reportRepeats )xldi.printRepeatsDelimited2 ( os, 0 );
	xldi.printDelimited1 ( os, 0 );
	xldi.printDelimited2 ( os, 0 );
	delimitedRowEnd ( os );
}
void SearchResultsCrosslinkPeptideHit::init ()
{
	instrument				= PeptidePosition::getInstrument ();
	parentTolerances		= PeptidePosition::getParentTolerances ();
	fragmentTolerances		= PeptidePosition::getFragmentTolerances ();
	spottingPlates			= PeptidePosition::getSpottingPlates ();
	spectrumNumber			= PeptidePosition::getSpectrumNumber ();
	fractionNames			= PeptidePosition::getFractionNames ();
	multipleFractionNames	= PeptidePosition::getMultipleFractionNames ();
	rawTypes				= PeptidePosition::getRawTypes ();
	searchKey				= PeptidePosition::getSearchKey ();

	reportMPlusH	= PeptidePosition::getReportMPlusH ();
	reportMOverZ	= PeptidePosition::getReportMOverZ ();
	reportCharge	= PeptidePosition::getReportCharge ();
	reportIntensity	= PeptidePosition::getReportIntensity ();
	reportMSMSInfo	= PeptidePosition::getReportMSMSInfo ();
	reportStartAA	= PeptidePosition::getReportStartAA ();
	reportEndAA		= PeptidePosition::getReportEndAA ();
	reportError		= PeptidePosition::getReportError ();
	reportLinks		= PeptidePosition::getReportLinks ();
	reportTime		= PeptidePosition::getReportTime ();
	reportPreviousAA= PeptidePosition::getReportPreviousAA ();
	reportDBPeptide	= PeptidePosition::getReportDBPeptide ();
	reportNextAA	= PeptidePosition::getReportNextAA ();
	reportRepeats	= PPPeptideHitInfo::getReportRepeats ();

	const ParameterList* pList = ProgramLink::getParams ();
	reportXLPeptide		= pList->getBoolValue ( "report_xl_peptide" );
	reportXLScore		= pList->getBoolValue ( "report_xl_score" );
	reportXLExpectation	= pList->getBoolValue ( "report_xl_expectation" );
	reportXLMValue		= pList->getBoolValue ( "report_xl_nlog_p" );
	reportXLRank		= pList->getBoolValue ( "report_xl_rank" );
	reportXLLowScore	= pList->getBoolValue ( "report_xl_low_score" );
	reportXLLowEVal		= pList->getBoolValue ( "report_xl_low_expectation" );
	reportXLLowMValue	= pList->getBoolValue ( "report_xl_low_nlog_p" );
	reportXLAA			= pList->getBoolValue ( "report_xl_aa" );
}
bool SearchResultsCrosslinkPeptideHit::outputQuanResults ( ostream& os, bool area ) const
{
	bool flag = false;
	if ( PeptidePositionQuan::outputQuanResults ( os, quanRatio, "", 1, area ) ) {
		flag = true;
	}
	return flag;
}
DoubleVector SearchResultsCrosslinkPeptideHit::getRatios ( bool area ) const
{
	if ( area )	return PeptidePositionQuan::getAreaRatios ( quanRatio );
	else		return PeptidePositionQuan::getIntensityRatios ( quanRatio );
}

bool SearchResultsCrosslinkProteinHit::reportLinks = false;
SearchResultsCrosslinkProteinHit::SearchResultsCrosslinkProteinHit ()
{
	ppphqi.initQuan ();
}
void SearchResultsCrosslinkProteinHit::init ()
{
	reportLinks = PeptidePosition::getReportLinks ();
}
void SearchResultsCrosslinkProteinHit::printHTML ( ostream& os, const SResLink& sresLink, int searchNumber, const LinkInfo* linkInfo ) const
{
	cLinkPeptideHits [0]->printProteinInfoHTML ( os, sresLink );
	if ( PPProteinHitQuanInfo::getQuan () ) {
		tableStart ( os, true );
			tableRowStart ( os );
				ppphqi.printHeaderHTML ( os, "" );
			tableRowEnd ( os );
			tableRowStart ( os );
				ppphqi.printHTML ( os, "" );
			tableRowEnd ( os );
		tableEnd ( os );
		os << "<br />" << endl;
		if ( QuantitationRatio::getAreaRatioReport () )	quanPlot ( os, true );
		if ( QuantitationRatio::getIntRatioReport () )	quanPlot ( os, false );
		os << "<br />" << endl;
	}
	tableStart ( os, false, "", "10" );
	SearchResultsCrosslinkPeptideHit::printHeaderHTML ( os, "", searchNumber );
	for ( int i = 0 ; i < cLinkPeptideHits.size () ; i++ ) {
		SCMSTagLink smtl;
		MSProductLink* productLink = 0;
		if ( reportLinks ) {
			productLink = new MSProductLink ( PeptidePosition::getInstrument (searchNumber), PeptidePosition::getSearchKey (searchNumber), cLinkPeptideHits [i]->getSpecIDasID (), PeptidePosition::getParentTolerance (searchNumber), PeptidePosition::getFragmentTolerance (searchNumber) );
			startJavascript ( os );
			productLink->printHTML ( os );
			smtl.printHTML ( os );
			endJavascript ( os );
		}
		cLinkPeptideHits [i]->printHTML ( os, "", productLink, smtl, linkInfo );
		if ( reportLinks ) delete productLink;
	}
	tableEnd ( os );
}
void SearchResultsCrosslinkProteinHit::printHeaderDelimited ( ostream& os, int searchIdx ) const
{
	cLinkPeptideHits [0]->printHeaderDelimited ( os, searchIdx, ppphqi );
}
void SearchResultsCrosslinkProteinHit::printDelimited ( ostream& os ) const
{
	for ( int i = 0 ; i < cLinkPeptideHits.size () ; i++ ) {
		cLinkPeptideHits [i]->printDelimited ( os, ppphqi );
	}
}
void SearchResultsCrosslinkProteinHit::addSpecIDs ( SetSpecID& ssid ) const
{
	for ( int i = 0 ; i < cLinkPeptideHits.size () ; i++ ) {
		ssid.insert ( cLinkPeptideHits [i]->getSpecIDasID () );
	}
}
void SearchResultsCrosslinkProteinHit::quanPlot ( ostream& os, bool area ) const
{
	if ( RPlot::getRFlag () ) {
		RPlot rplot ( "quan.R" );
		GenOFStream ofs ( rplot.getDataFileFullPath () );
		bool ok = false;
		for ( int i = 0 ; i < cLinkPeptideHits.size () ; i++ ) {
			if ( cLinkPeptideHits [i]->outputQuanResults ( ofs, area ) ) {
				ok = true;
			}
		}
		ofs.close ();
		if ( ok ) rplot.printImageAndLink ( os );
	}
}
void SearchResultsCrosslinkProteinHit::quanProteinStats ( bool area )
{
	DoubleVectorVector ratios;
	for ( int i = 0 ; i < cLinkPeptideHits.size () ; i++ ) {
		ratios.push_back ( cLinkPeptideHits [i]->getRatios ( area ) );
	}
	DoubleVectorVector rat;
	int numQuanPeaks = 0;
	for ( DoubleVectorVectorSizeType j = 0 ; j < ratios.size () ; j++ ) {		// Iterate through peptides
		int nqp = ratios [j].size ();
		if ( nqp > rat.size () ) {
			rat.resize ( nqp );
			numQuanPeaks = nqp;
		}
		for ( DoubleVectorSizeType m = 0 ; m < nqp ; m++ ) {	// Iterate through quan peaks
			double r = ratios [j][m];
			if ( r != 0.0 ) rat [m].push_back ( log10 ( r ) );	// m=quan pk, j=peptides
		}
	}
	DoubleVector med ( numQuanPeaks, 0.0 );
	DoubleVector q1 ( numQuanPeaks, 0.0 );
	DoubleVector q2 ( numQuanPeaks, 0.0 );
	DoubleVector mean ( numQuanPeaks, 0.0 );
	DoubleVector stdev ( numQuanPeaks, 0.0 );
	IntVector numPks ( numQuanPeaks, 0 );
	for ( DoubleVectorVectorSizeType y = 0 ; y < rat.size () ; y++ ) {	// Iterate through quan peaks
		DoubleVector dv = rat [y];
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
	if ( area )		ppphqi.setAreaRatios ( med, q1, q2, mean, stdev, numPks );
	else			ppphqi.setIntensityRatios ( med, q1, q2, mean, stdev, numPks );
}
int sortXlinkPeptideHits ( const SearchResultsCrosslinkPeptideHit* a, const SearchResultsCrosslinkPeptideHit* b, int ineq )
{
	const std::string& acc1 = a->getAccessionNumbers ();
	const std::string& acc2 = b->getAccessionNumbers ();
	bool inter1 = acc1.find ( '\t' ) != std::string::npos;		
	bool inter2 = acc2.find ( '\t' ) != std::string::npos;
	if ( inter1 == inter2 ) {
		if ( acc1 == acc2 )
			return ineq;
		else
			return acc1 < acc2;
	}
	else
		return inter1 > inter2;
}

std::map <std::string, ProteinInfo*> XLinkDecoyInfo::proteinInfoMap;

void XLinkDecoyInfo::add ( const string& a, int s, const string& a2, int s2 )
{
	accNum.push_back ( a );
	startAA.push_back ( s );
	accNum2.push_back ( a2 );
	startAA2.push_back ( s2 );
	std::map <std::string, ProteinInfo*>::const_iterator cur = proteinInfoMap.find ( a );
	if ( cur == proteinInfoMap.end () ) {
		proteinInfoMap [a] = new ProteinInfo ( a );
	}
	cur = proteinInfoMap.find ( a2 );
	if ( cur == proteinInfoMap.end () ) {
		proteinInfoMap [a2] = new ProteinInfo ( a2 );
	}
}
void XLinkDecoyInfo::calc ( int num )
{
	pep1IntraCount = 1;
	pep1InterCount = 0;
	pep1DecoyCount = 0;
	pep2IntraCount = 1;
	pep2InterCount = 0;
	pep2DecoyCount = 0;
	for ( StringVectorSizeType i = 0 ; i < accNum.size () ; i++ ) {
		if ( accNum [i] == accNum [num] ) {
			if ( startAA [i] != startAA [num] ) pep1IntraCount++;
		}
		else {
			if ( ProteinInfo::isDecoy ( accNum [i] ) )	pep1DecoyCount++;
			else										pep1InterCount++;
		}
		if ( accNum2 [i] == accNum2 [num] ) {
			if ( startAA2 [i] != startAA2 [num] ) pep2IntraCount++;
		}
		else {
			if ( ProteinInfo::isDecoy ( accNum2 [i] ) )	pep2DecoyCount++;
			else										pep2InterCount++;
		}
	}
}
string XLinkDecoyInfo::getPreviousAA1 ( int num, int reportPreviousAA ) const
{
	std::map <std::string, ProteinInfo*>::const_iterator cur = proteinInfoMap.find ( accNum [num] );
	return (*cur).second->getPreviousAA ( startAA [num], reportPreviousAA );
}
string XLinkDecoyInfo::getPreviousAA2 ( int num, int reportPreviousAA ) const
{
	std::map <std::string, ProteinInfo*>::const_iterator cur = proteinInfoMap.find ( accNum2 [num] );
	return (*cur).second->getPreviousAA ( startAA2 [num], reportPreviousAA );
}
string XLinkDecoyInfo::getNextAA1 ( int num, int endAA, int reportNextAA ) const
{
	std::map <std::string, ProteinInfo*>::const_iterator cur = proteinInfoMap.find ( accNum [num] );
	return (*cur).second->getNextAA ( endAA, reportNextAA );
}
string XLinkDecoyInfo::getNextAA2 ( int num, int endAA, int reportNextAA ) const
{
	std::map <std::string, ProteinInfo*>::const_iterator cur = proteinInfoMap.find ( accNum2 [num] );
	return (*cur).second->getNextAA ( endAA, reportNextAA );
}
int XLinkDecoyInfo::getStartAA1 ( int num ) const
{
	return startAA [num];
}
int XLinkDecoyInfo::getStartAA2 ( int num ) const
{
	return startAA2 [num];
}
bool XLinkDecoyInfo::isDecoy1 ( int num ) const
{
	std::map <std::string, ProteinInfo*>::const_iterator cur = proteinInfoMap.find ( accNum [num] );
	return (*cur).second->isDecoy ();
}
bool XLinkDecoyInfo::isDecoy2 ( int num ) const
{
	std::map <std::string, ProteinInfo*>::const_iterator cur = proteinInfoMap.find ( accNum2 [num] );
	return (*cur).second->isDecoy ();
}
string XLinkDecoyInfo::getFullAccessionNumber1 ( int num ) const
{
	std::map <std::string, ProteinInfo*>::const_iterator cur = proteinInfoMap.find ( accNum [num] );
	return (*cur).second->getFullAccessionNumber ();
}
string XLinkDecoyInfo::getFullAccessionNumber2 ( int num ) const
{
	std::map <std::string, ProteinInfo*>::const_iterator cur = proteinInfoMap.find ( accNum2 [num] );
	return (*cur).second->getFullAccessionNumber ();
}
void XLinkDecoyInfo::printProteinInfoHTML ( ostream& os, const SResLink& sresLink, int num ) const
{
	std::map <std::string, ProteinInfo*>::const_iterator cur = proteinInfoMap.find ( accNum [num] );
	ProteinInfo* proteinInfo1 = (*cur).second;
	cur = proteinInfoMap.find ( accNum2 [num] );
	ProteinInfo* proteinInfo2 = (*cur).second;
	proteinInfo1->printHTMLLines ( os );
	if ( PeptidePosition::getReportLinks () ) {
		sresLink.write ( os, accNum [num], "", "Peptide Hits" );
	}
	os << "<br /><br />" << endl;
	if ( proteinInfo1->getAcc () != proteinInfo2->getAcc () || proteinInfo1->isDecoy () != proteinInfo2->isDecoy () ) {
		proteinInfo2->printHTMLLines ( os );
		if ( PeptidePosition::getReportLinks () ) {
			sresLink.write ( os, accNum2 [num], "", "Peptide Hits" );
		}
		os << "<br /><br />" << endl;
	}
}
string XLinkDecoyInfo::getRepeats1 ( int num ) const
{
	return gen_itoa ( pep1IntraCount ) + ":" + gen_itoa ( pep1InterCount ) + ":" + gen_itoa ( pep1DecoyCount );
}
string XLinkDecoyInfo::getRepeats2 ( int num ) const
{
	return gen_itoa ( pep2IntraCount ) + ":" + gen_itoa ( pep2InterCount ) + ":" + gen_itoa ( pep2DecoyCount );
}
void XLinkDecoyInfo::printRepeatsDelimitedHeader1 ( ostream& os )
{
	delimitedHeader ( os, "Intra 1" );
	delimitedHeader ( os, "Inter 1" );
	delimitedHeader ( os, "Decoy 1" );
}
void XLinkDecoyInfo::printRepeatsDelimitedHeader2 ( ostream& os )
{
	delimitedHeader ( os, "Intra 2" );
	delimitedHeader ( os, "Inter 2" );
	delimitedHeader ( os, "Decoy 2" );
}
void XLinkDecoyInfo::printRepeatsDelimited1 ( ostream& os, int num ) const
{
	delimitedCell ( os, pep1IntraCount );
	delimitedCell ( os, pep1InterCount );
	delimitedCell ( os, pep1DecoyCount );
}
void XLinkDecoyInfo::printRepeatsDelimited2 ( ostream& os, int num ) const
{
	delimitedCell ( os, pep2IntraCount );
	delimitedCell ( os, pep2InterCount );
	delimitedCell ( os, pep2DecoyCount );
}
void XLinkDecoyInfo::printHTMLANum1 ( ostream& os, int num ) const
{
	std::map <std::string, ProteinInfo*>::const_iterator cur = proteinInfoMap.find ( accNum [num] );
	(*cur).second->printHTMLANum ( os, false );
}
void XLinkDecoyInfo::printHTMLANum2 ( ostream& os, int num ) const
{
	std::map <std::string, ProteinInfo*>::const_iterator cur = proteinInfoMap.find ( accNum2 [num] );
	(*cur).second->printHTMLANum ( os, false );
}
void XLinkDecoyInfo::printDelimitedANum1 ( ostream& os, int num ) const
{
	std::map <std::string, ProteinInfo*>::const_iterator cur = proteinInfoMap.find ( accNum [num] );
	(*cur).second->printDelimitedANum ( os );
}
void XLinkDecoyInfo::printDelimitedANum2 ( ostream& os, int num ) const
{
	std::map <std::string, ProteinInfo*>::const_iterator cur = proteinInfoMap.find ( accNum2 [num] );
	(*cur).second->printDelimitedANum ( os );
}
void XLinkDecoyInfo::printDelimited1 ( ostream& os, int num ) const
{
	std::map <std::string, ProteinInfo*>::const_iterator cur = proteinInfoMap.find ( accNum [num] );
	(*cur).second->printDelimited ( os );
}
void XLinkDecoyInfo::printDelimited2 ( ostream& os, int num ) const
{
	std::map <std::string, ProteinInfo*>::const_iterator cur = proteinInfoMap.find ( accNum2 [num] );
	(*cur).second->printDelimited ( os );
}
