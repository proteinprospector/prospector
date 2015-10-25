/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_inst.cpp                                                   *
*                                                                             *
*  Created    : November 3rd 1997                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1997-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#define LUCSF_INST_MAIN
#include <lg_string.h>
#include <lgen_error.h>
#include <lu_getfil.h>
#include <lu_inst.h>
#include <lu_pk_filter.h>
#include <lu_param_list.h>
using std::string;
using std::getline;

void initialiseInstrumentName ( const string& instrumentName )
{
	if ( instInf == 0 ) instInf = new InstrumentInfo ( instrumentName );
}
void resetInstrumentName ( const string& instrumentName )
{
	delete instInf;
	instInf = new InstrumentInfo ( instrumentName );
}
InstrumentInfo::InstrumentInfo ( const string& instrumentName ) :
	name					( instrumentName ),
	parentPrecision			( 4 ),
	parentErrorSigFig		( 3 ),
	parentIntensitySigFig	( 3 ),

	fragmentPrecision		( 4 ),
	fragmentErrorSigFig		( 2 ),
	fragmentIntensitySigFig	( 3 ),

	allowIncorrectCharge ( false ),
	frag ( 0 )
{
	if ( instrumentName != "" ) {		// If no instrument name specified use defaults.
		GenIFStream fromFile ( MsparamsDir::instance ().getParamPath ( "instrument.txt" ) );
		string line;
		while ( getline ( fromFile, line ) ) {
			if ( line == instrumentName ) {
				ParameterList params ( fromFile );

				params.getValue ( "parent_precision", parentPrecision );
				params.getValue ( "parent_error_significant_figures", parentErrorSigFig );
				params.getValue ( "parent_intensity_significant_figures", parentIntensitySigFig );

				params.getValue ( "fragment_precision", fragmentPrecision );
				params.getValue ( "fragment_error_significant_figures", fragmentErrorSigFig );
				params.getValue ( "fragment_intensity_significant_figures", fragmentIntensitySigFig );

				params.getValue ( "allow_incorrect_charge", allowIncorrectCharge );
				setPrecisions ();
				MSPeakFilterOptions::setDefault ( MSPeakFilterOptions ( &params ) );
				MSMSPeakFilterOptions::setDefault ( MSMSPeakFilterOptions ( &params ) );
				string fragType; 
				params.getValue ( "frag", fragType );
				frag = new MSMSFragmentation ( fragType );
				return;
			}
		}
		ErrorHandler::genError ()->error ( "Invalid or unspecified instrument name.\nFunction: InstrumentInfo.\n" );
	}
	setPrecisions ();
	frag = new MSMSFragmentation ( "" );
}
InstrumentInfo::~InstrumentInfo ()
{
	delete frag;
}
unsigned int InstrumentInfo::getAmmoniaLossMask () const
{
	return frag->getAmmoniaLossMask ();
}
unsigned int InstrumentInfo::getWaterLossMask () const
{
	return frag->getWaterLossMask ();
}
unsigned int InstrumentInfo::getPosChargeBearingMask () const
{
	return frag->getPosChargeBearingMask ();
}
unsigned int InstrumentInfo::getDIonExcludeMask () const
{
	return frag->getDIonExcludeMask ();
}
unsigned int InstrumentInfo::getWIonExcludeMask () const
{
	return frag->getWIonExcludeMask ();
}
unsigned int InstrumentInfo::getVIonExcludeMask () const
{
	return frag->getVIonExcludeMask ();
}
StringVector InstrumentInfo::getIonTypes () const
{
	return frag->getIonTypes ();
}
DoubleVector InstrumentInfo::getLossMasses () const
{
	return frag->getLossMasses ();
}
ScoreType InstrumentInfo::getUnmatchedScore () const
{
	return frag->getUnmatchedScore ();
}

ScoreType InstrumentInfo::getImmoniumScore () const
{
	return frag->getImmoniumScore ();
}
ScoreType InstrumentInfo::getRelatedIonScore () const
{
	return frag->getRelatedIonScore ();
}
ScoreType InstrumentInfo::getMScore () const
{
	return frag->getMScore ();
}

ScoreType InstrumentInfo::getAScore () const
{
	return frag->getAScore ();
}
ScoreType InstrumentInfo::getALossScore () const
{
	return frag->getALossScore ();
}
ScoreType InstrumentInfo::getAPhosLossScore () const
{
	return frag->getAPhosLossScore ();
}
ScoreType InstrumentInfo::getBScore () const
{
	return frag->getBScore ();
}
ScoreType InstrumentInfo::getBPlusH2OScore () const
{
	return frag->getBPlusH2OScore ();
}
ScoreType InstrumentInfo::getBLossScore () const
{
	return frag->getBLossScore ();
}
ScoreType InstrumentInfo::getBPhosLossScore () const
{
	return frag->getBPhosLossScore ();
}
ScoreType InstrumentInfo::getCLadderScore () const
{
	return frag->getCLadderScore ();
}
ScoreType InstrumentInfo::getCPlus2DaScore () const
{
	return frag->getCPlus2DaScore ();
}
ScoreType InstrumentInfo::getCPlus1DaScore () const
{
	return frag->getCPlus1DaScore ();
}
ScoreType InstrumentInfo::getCScore () const
{
	return frag->getCScore ();
}
ScoreType InstrumentInfo::getCMinus1DaScore () const
{
	return frag->getCMinus1DaScore ();
}
ScoreType InstrumentInfo::getDScore () const
{
	return frag->getDScore ();
}

ScoreType InstrumentInfo::getVScore () const
{
	return frag->getVScore ();
}
ScoreType InstrumentInfo::getWScore () const
{
	return frag->getWScore ();
}
ScoreType InstrumentInfo::getXScore () const
{
	return frag->getXScore ();
}
ScoreType InstrumentInfo::getNLadderScore () const
{
	return frag->getNLadderScore ();
}
ScoreType InstrumentInfo::getYScore () const
{
	return frag->getYScore ();
}
ScoreType InstrumentInfo::getYLossScore () const
{
	return frag->getYLossScore ();
}
ScoreType InstrumentInfo::getYPhosLossScore () const
{
	return frag->getYPhosLossScore ();
}
ScoreType InstrumentInfo::getBigYScore () const
{
	return frag->getBigYScore ();
}
ScoreType InstrumentInfo::getZScore () const
{
	return frag->getZScore ();
}
ScoreType InstrumentInfo::getZPlus1DaScore () const
{
	return frag->getZPlus1DaScore ();
}
ScoreType InstrumentInfo::getZPlus2DaScore () const
{
	return frag->getZPlus2DaScore ();
}
ScoreType InstrumentInfo::getZPlus3DaScore () const
{
	return frag->getZPlus3DaScore ();
}

ScoreType InstrumentInfo::getBP2Score () const
{
	return frag->getBP2Score ();
}
ScoreType InstrumentInfo::getBP2LossScore () const
{
	return frag->getBP2LossScore ();
}
ScoreType InstrumentInfo::getBP2PhosLossScore () const
{
	return frag->getBP2PhosLossScore ();
}
ScoreType InstrumentInfo::getYP2Score () const
{
	return frag->getYP2Score ();
}
ScoreType InstrumentInfo::getYP2LossScore () const
{
	return frag->getYP2LossScore ();
}
ScoreType InstrumentInfo::getYP2PhosLossScore () const
{
	return frag->getYP2PhosLossScore ();
}

ScoreType InstrumentInfo::getInternalAScore () const
{
	return frag->getInternalAScore ();
}
ScoreType InstrumentInfo::getInternalBScore () const
{
	return frag->getInternalBScore ();
}
ScoreType InstrumentInfo::getInternalLossScore () const
{
	return frag->getInternalLossScore ();
}

ScoreType InstrumentInfo::getMH3PO4Score () const
{
	return frag->getMH3PO4Score ();
}
ScoreType InstrumentInfo::getMSOCH4Score () const
{
	return frag->getMSOCH4Score ();
}

bool InstrumentInfo::getScoring () const
{
	return frag->getScoring ();
}
double InstrumentInfo::getMaximumInternalIonMass () const
{
	return frag->getMaximumInternalIonMass ();
}
bool InstrumentInfo::getChargeReducedFragmentation () const
{
	return frag->getChargeReducedFragmentation ();
}
bool InstrumentInfo::getETDTypeFragmentation () const
{
	return genToLower ( getFragName () ).find ( "etd" ) != string::npos;
}

void InstrumentInfo::setPrecisions ()
{
	parentPeakPrecision.setMassDecimalPlaces ( parentPrecision );
	parentPeakPrecision.setErrorSigFig ( parentErrorSigFig );
	parentPeakPrecision.setIntensitySigFig ( parentIntensitySigFig );
	fragmentPeakPrecision.setMassDecimalPlaces ( fragmentPrecision );
	fragmentPeakPrecision.setErrorSigFig ( fragmentErrorSigFig );
	fragmentPeakPrecision.setIntensitySigFig ( fragmentIntensitySigFig );
}
void InstrumentList::initialise ()
{
	GenCommentedIFStream ifs ( MsparamsDir::instance ().getParamPath ( "instrument.txt" ) );
	string line;
	while ( ifs.getUncommentedLine ( line ) ) {
		names.push_back ( line );
		ParameterList params ( ifs );
	}
}
InstrumentList::InstrumentList ()
{
	initialise ();
}
InstrumentList& InstrumentList::instance ()
{
	static InstrumentList instrumentList;
	return instrumentList;
}
