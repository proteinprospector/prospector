/******************************************************************************
*                                                                             *
*  Program    : searchCompare                                                 *
*                                                                             *
*  Filename   : sc_quan.cpp                                                   *
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
*  Copyright (2012-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef VIS_C
#pragma warning( disable : 4503 )	// Disable warnings about decorated name length being too long
#endif

#include <lgen_error.h>
#include <lr_main.h>
#include <lu_aa_calc.h>
#include <lu_file_type.h>
#include <lu_get_link.h>
#include <lu_inst.h>
#include <lu_spec_id.h>
#include <lu_param_list.h>
#include <lu_proj_file.h>
#include <lu_quan_multi.h>
#include <lu_quan_ratio.h>
#include <lu_tol.h>
#include <sc_quan.h>
#include <sc_search_res.h>

using std::string;
using std::vector;
using std::ostream;
using std::runtime_error;

using namespace FileTypes;

AACalculator* PeptidePositionQuan::aaCalc = 0;
bool PeptidePositionQuan::quanMSMSFlag;
bool PeptidePositionQuan::n15Flag;
bool PeptidePositionQuan::qFlag;
string PeptidePositionQuan::qType;
double PeptidePositionQuan::defaultResolution;
int PeptidePositionQuan::numIsotopePeaks = 1;
RawFile* PeptidePositionQuan::rf;
StringVector PeptidePositionQuan::instrument;
StringVectorVector PeptidePositionQuan::rawTypes;
StringVector PeptidePositionQuan::version;
VectorConstParameterListPtr PeptidePositionQuan::params;
vector <Tolerance*> PeptidePositionQuan::parentTolerances;
DoubleVectorVector PeptidePositionQuan::offsets;
bool PeptidePositionQuan::quanMultiNormalFlag = false;
double PeptidePositionQuan::rtIntervalStart = 0.0;
double PeptidePositionQuan::rtIntervalEnd = 0.0;

void PeptidePositionQuan::initialiseQuan ( int i, int fraction )
{
	resetInstrumentName ( instrument [i] );
	try {
		rf = new RawFile ( params [i], fraction, rtIntervalStart, rtIntervalEnd );
	}
	catch ( runtime_error e ) {
		ErrorHandler::genError ()->error ( e );
	}
}
void PeptidePositionQuan::deleteQuan ( int i )
{
	delete rf;
}
PeakFitData* PeptidePositionQuan::getQuanRatio ( double mOverZ, int charge, const QuanPeptide& qp, const SpecID& specID, int searchIndex )
{
	PeakFitData* quanRatio;
	vector <XYData> vXYData (1);
	double startMass = mOverZ - 6.0;
	double endMass = mOverZ + 10.0;
	int fraction = specID.getFraction ();
	if ( n15Flag ) {
		startMass -= 20.0;
		endMass += 20.0;
	}
	if ( quanMSMSFlag ) {
		startMass = QuantitationMulti::getStartMass ();
		endMass = QuantitationMulti::getEndMass ();
	}
	else if ( qFlag ) {
		QuantitationRatio::getDataRange ( qp, mOverZ, charge, startMass, endMass );
	}
	PairStringVectorString psvs;
	try {
		rf->getQuantitationXYData ( vXYData, psvs, quanMSMSFlag, specID, mOverZ, version [searchIndex], startMass, endMass );
	}
	catch ( runtime_error e ) {
		ErrorHandler::genError ()->error ( e );
	}
	if ( !quanMSMSFlag ) {
		string tolUnits = parentTolerances [searchIndex]->getUnitsString ();
		if ( offsets [searchIndex][fraction-1] != 0.0 ) {
			vXYData [0].offsetCalibration ( offsets [searchIndex][fraction-1], tolUnits );
		}
	}
	if ( quanMSMSFlag ) {
		if ( rawTypes [searchIndex][fraction-1] == WIFF || rawTypes [searchIndex][fraction-1] == RAW || rawTypes [searchIndex][fraction-1] == "" )
			quanRatio = new QuantitationMultiMassWindow ( vXYData [0] );
		else
			quanRatio = new QuantitationMultiNormal ( vXYData [0], defaultResolution );
	}
	else {
		ElementalFormula ef;
		ElementalFormula* efp;
		bool flag;
		try {
			if ( qp.isXLink () ) {
				flag = aaCalc->calculateStrippedElementalCompositionWithTerminii ( qp.getPeptide1 (), qp.getNTerm1 (), qp.getCTerm1 (), qp.getNLoss (), ef, false );
				ElementalFormula ef2;
				flag = aaCalc->calculateStrippedElementalCompositionWithTerminii ( qp.getPeptide2 (), qp.getNTerm2 (), qp.getCTerm2 (), qp.getNLoss (), ef2, false );
				ef += ef2;
				ef += qp.getBridgeFormula ();
				ef -= "H";
			}
			else {
				flag = aaCalc->calculateStrippedElementalCompositionWithTerminii ( qp.getPeptide1 (), qp.getNTerm1 (), qp.getCTerm1 (), qp.getNLoss (), ef );
			}
			efp = &ef;
		}
		catch ( AACalculatorNoElementalComposition ) {
			efp = 0;
		}
		if ( qFlag ) {
			quanRatio = new QuantitationRatio ( vXYData [0], mOverZ, charge, qp, efp, flag, defaultResolution, numIsotopePeaks );
		}
		else
			quanRatio = new ParentData ( vXYData [0], mOverZ, charge, efp, flag, defaultResolution, 1 );
	}
	return quanRatio;
}
int PeptidePositionQuan::getColspan ( int searchNumber )
{
	int colspan = 0;
	if ( quanMSMSFlag )	{
		bool qmnFlag = ( searchNumber == -1 ) ? quanMultiNormalFlag : rawTypes [searchNumber][0] != WIFF;
		colspan += QuantitationMulti::getColspan ( qmnFlag );
	}
	else if ( qFlag ) {
		for ( int i = 0 ; i < numIsotopePeaks ; i++ ) {
			colspan += QuantitationRatio::getColspan ( i );
		}
	}
	else				colspan += ParentData::getColspan ();			
	colspan += PeakFitData::getNoiseColspan ();
	return colspan;
}
void PeptidePositionQuan::printHeaderHTML ( ostream& os, int searchNumber, const string& styleID )
{
	if ( quanMSMSFlag )	{
		bool qmnFlag = ( searchNumber == -1 ) ? quanMultiNormalFlag : rawTypes [searchNumber][0] != WIFF;
		QuantitationMulti::printHTMLHeader ( os, styleID, qmnFlag );
	}
	else if ( qFlag ) {
		for ( int i = 0 ; i < numIsotopePeaks ; i++ ) QuantitationRatio::printHTMLHeader ( os, i, styleID );
	}
	else ParentData::printHTMLHeader ( os, styleID );				
	PeakFitData::printHTMLNoiseHeader ( os, styleID );
}
void PeptidePositionQuan::printQuanBlankHTML ( ostream& os, int searchNumber, const string& styleID )
{
	if ( quanMSMSFlag ) {
		bool qmnFlag = ( searchNumber == -1 ) ? quanMultiNormalFlag : rawTypes [searchNumber][0] != WIFF;
		QuantitationMulti::printHTMLBlankLine ( os, styleID, qmnFlag );
	}
	else if ( qFlag ) {
		for ( int i = 0 ; i < numIsotopePeaks ; i++ ) QuantitationRatio::printHTMLBlankLine ( os, i, styleID );
	}
	else				ParentData::printHTMLBlankLine ( os, styleID );					
	PeakFitData::printHTMLNoiseBlankLine ( os, styleID );
}
void PeptidePositionQuan::printQuanBlankDelimited ( ostream& os, int searchNumber )
{
	if ( quanMSMSFlag ) {
		bool qmnFlag = ( searchNumber == -1 ) ? quanMultiNormalFlag : rawTypes [searchNumber][0] != WIFF;
		QuantitationMulti::printDelimitedBlankLine ( os, qmnFlag );
	}
	else if ( qFlag ) {
		for ( int i = 0 ; i < numIsotopePeaks ; i++ ) QuantitationRatio::printDelimitedBlankLine ( os, i );
	}
	else				ParentData::printDelimitedBlankLine ( os );		
	PeakFitData::printDelimitedNoiseBlankLine ( os );
}
void PeptidePositionQuan::printHTML ( ostream& os, const PeakFitData* quanRatio, const string& styleID, bool joint )
{
	for ( int i = 0 ; i < ( quanMSMSFlag ? 1 : numIsotopePeaks ) ; i++ ) {
		if ( quanMSMSFlag && joint && quanMultiNormalFlag ) quanRatio->printHTMLLine ( os, -1, styleID );
		else												quanRatio->printHTMLLine ( os, i, styleID );
	}
	quanRatio->printHTMLNoiseLine ( os, styleID );
}
void PeptidePositionQuan::printHeaderDelimited ( ostream& os, int searchNumber )
{
	if ( quanMSMSFlag ) {
		bool qmnFlag = ( searchNumber == -1 ) ? quanMultiNormalFlag : rawTypes [searchNumber][0] != WIFF;
		QuantitationMulti::printDelimitedHeader ( os, qmnFlag );
	}
	else if ( qFlag ) {
		for ( int i = 0 ; i < numIsotopePeaks ; i++ ) QuantitationRatio::printDelimitedHeader ( os, i );
	}
	else				ParentData::printDelimitedHeader ( os );
	PeakFitData::printDelimitedNoiseHeader ( os );
}
void PeptidePositionQuan::printDelimited ( ostream& os, const PeakFitData* quanRatio, bool joint )
{
	for ( int i = 0 ; i < ( quanMSMSFlag ? 1 : numIsotopePeaks ) ; i++ ) {
		if ( quanMSMSFlag && joint && quanMultiNormalFlag )	quanRatio->printDelimitedLine ( os, -1 );
		else												quanRatio->printDelimitedLine ( os, i );
	}
	quanRatio->printDelimitedNoiseLine ( os );
}
bool PeptidePositionQuan::outputQuanResults ( ostream& os, const PeakFitData* quanRatio, const string& searchName, int numRepeats, bool area )
{
	bool flag = false;
	if ( quanRatio ) flag = quanRatio->outputQuanResults ( os, searchName, numRepeats, area );
	return flag;
}
DoubleVector PeptidePositionQuan::getIntensityRatios ( const PeakFitData* quanRatio )
{
	DoubleVector ret;
	if ( quanRatio ) ret = quanRatio->getIntensityRatios ();
	return ret;
}
DoubleVector PeptidePositionQuan::getAreaRatios ( const PeakFitData* quanRatio )
{
	DoubleVector ret;
	if ( quanRatio ) ret = quanRatio->getAreaRatios ();
	return ret;
}
void PeptidePositionQuan::initialiseAACalculator ( const MapStringConstModPtr& constMods, const string& rawType, const string& quanType, double resolution )
{
	qFlag = ( rawType == "Quantitation" );
	quanMSMSFlag = qFlag && isQuanMSMS ( quanType );
	n15Flag = quanType == "Label:15N";
	qType = quanType;
	defaultResolution = resolution;
	if ( quanMSMSFlag ) QuantitationMulti::init ( quanType );
	aaCalc = new AACalculator ( true, constMods );
	if ( qFlag && !quanMSMSFlag ) QuantitationRatio::setQuanResidue ( quanType );
}
void PeptidePositionQuan::initialiseParams ( const ParameterList* p )
{
	params.push_back ( p );
	instrument.push_back ( p->getStringValue ( "instrument_name" ) );
	static LinkInfo* linkInfo = new LinkInfo ( p );
	QuanPeptide::setLinkInfo ( linkInfo );
}
void PeptidePositionQuan::initialiseDataSetInfo ( double timeWindowStart, double timeWindowEnd )
{
	rtIntervalStart = timeWindowStart;
	rtIntervalEnd = timeWindowEnd;
}
int PeptidePositionQuan::getNumQuanRatioColumns ()
{
	return quanMSMSFlag ? getNumReporterIons ( qType ) - 1 : QuantitationRatio::getNumRatios ( qType );
}
void PeptidePositionQuan::initialise ( const SearchResultsPtrVector& searchResults )
{
	for ( SearchResultsPtrVectorSizeType i = 0 ; i < searchResults.size () ; i++ ) {

		ProjectFile projectFile ( params [i] );
		version.push_back ( projectFile.getProjectVersion () );
		rawTypes.push_back ( searchResults [i]->getRawTypes () );
		if ( rawTypes.back () [0] != WIFF ) quanMultiNormalFlag = true;

		parentTolerances.push_back ( searchResults [i]->getParentTolerance () );

		offsets.push_back ( searchResults [i]->getOffsets () );
	}
}
