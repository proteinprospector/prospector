/******************************************************************************
*                                                                             *
*  Program    : msdisplay                                                     *
*                                                                             *
*  Filename   : msdisplay_params.cpp                                          *
*                                                                             *
*  Created    : December 20th 2002                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2002-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lu_get_link.h>
#include <lu_param_list.h>
#include <msdisplay_params.h>

using std::string;

MSDisplayParameters::MSDisplayParameters ( const ParameterList* params ) :
	MSProgramParameters	( params ),
	aaInitInfo			( params ),
	quanType			( params->getStringValue( "quan_type", "" ) ),
	rawType				( params->getStringValue( "raw_type", "" ) ),
	isotopePurity		( params ),
	resolution			( params->getDoubleValue( "resolution", 10000.0 ) ),
	numPeaks			( params->getIntValue	( "num_peaks", 3 ) ),
	project				( params->getStringValue( "project_name", "" ) ),
	specID				( params ),
	mOverZ				( params->getDoubleValue( "m_over_z", 0.0 ) ),
	rtIntervalStart		( params->getStringValue( "rt_int_start", "0.0" ) ),
	rtIntervalEnd		( params->getStringValue( "rt_int_end", "0.0" ) ),
	snrThreshold		( params->getDoubleValue( "snr_threshold", 0.0 ) ),
	areaThreshold		( params->getDoubleValue( "area_threshold", 0.0 ) ),
	intensityThreshold	( params->getDoubleValue( "intensity_threshold", 0.0 ) ),
	charge				( params->getIntValue	( "charge", 1 ) ),
	formula				( params->getStringValue( "formula", "" ) ),
	nterm				( params->getStringValue( "nterm", "" ) ),
	cterm				( params->getStringValue( "cterm", "" ) ),
	nloss				( params->getStringValue( "nloss", "" ) ),
	formula2			( params->getStringValue( "formula2", "" ) ),
	nterm2				( params->getStringValue( "nterm2", "" ) ),
	cterm2				( params->getStringValue( "cterm2", "" ) ),
	linkInfo			( new LinkInfo ( params ) ),
	qPeptide			( formula, nterm, cterm, nloss, formula2, nterm2, cterm2, linkInfo->getName () ),
	displayStartMass	( params->getDoubleValue ( "display_start_mass", 1000.0 ) ),
	displayEndMass		( params->getDoubleValue ( "display_end_mass", 1200.0 ) ),
	displayAllMasses	( params->getBoolValue ( "display_all_masses", true ) ),
	msType				( params->getStringValue ( "ms_type" ) ),
	msmsInfo			( params->getStringValue ( "msms_info" ) ),
	reporterIonWindow	( params->getDoubleValue ( "reporter_ion_window", 0.4 ) )
{
	if ( isQuanMSMS ( quanType ) ) {	// This is an MSMS quantitation type
		purityCorrection = new PurityCorrection ( params );
		QuantitationMulti::setPurityCorrection ( purityCorrection );
		QuantitationMultiMassWindow::setReporterIonWindow ( reporterIonWindow );
	}
}
string MSDisplayParameters::getBridgeFormula () const
{
	return linkInfo->getBridgeFormula ();
}
