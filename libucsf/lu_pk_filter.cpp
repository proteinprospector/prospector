/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_pk_filter.cpp                                              *
*                                                                             *
*  Created    : June 26th 2003                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2003-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lu_pk_filter.h>
#include <lu_param_list.h>

using std::ostream;
using std::string;

bool	MSPeakFilterOptions::peakExclusionDefault		= false;
unsigned int MSPeakFilterOptions::maxPeaksDefault		= 200;
unsigned int MSPeakFilterOptions::minPeaksDefault		= 5;
double	MSPeakFilterOptions::minIntensityDefault		= 0.0;
bool	MSPeakFilterOptions::massExclusionDefault		= false;
double	MSPeakFilterOptions::minMassDefault				= 50.0;
double	MSPeakFilterOptions::maxMassDefault				= 10000.0;
bool	MSPeakFilterOptions::matrixExclusionDefault		= false;
double	MSPeakFilterOptions::maxMatrixMassDefault		= 1300.0;

MSPeakFilterOptions::MSPeakFilterOptions ( const ParameterList* params ) :
	peakExclusion		( params->getBoolValue	( "ms_peak_exclusion",		peakExclusionDefault ) ),
	maxPeaks			( params->getUIntValue	( "ms_max_peaks",			maxPeaksDefault ) ),
	minPeaks			( params->getUIntValue	( "ms_min_peaks",			minPeaksDefault ) ),
	minIntensity		( params->getDoubleValue( "ms_min_intensity",		minIntensityDefault ) ),
	massExclusion		( params->getBoolValue	( "ms_mass_exclusion",		massExclusionDefault ) ),
	minMass				( params->getDoubleValue( "ms_min_mass",			minMassDefault ) ),
	maxMass				( params->getDoubleValue( "ms_max_mass",			maxMassDefault ) ),
	matrixExclusion		( params->getBoolValue	( "ms_matrix_exclusion",	matrixExclusionDefault ) ),
	maxMatrixMass		( params->getDoubleValue( "ms_max_matrix_mass",		maxMatrixMassDefault ) )
{
}
void MSPeakFilterOptions::copyToCGI ( ostream& os, const ParameterList* params )
{
	bool peakExclusion = params->copyToCGI ( os, "ms_peak_exclusion" );
	if ( peakExclusion ) {
		params->copyToCGI ( os, "ms_max_peaks" );
		params->copyToCGI ( os, "ms_min_peaks" );
		params->copyToCGI ( os, "ms_min_intensity" );
	}
	bool massExclusion = params->copyToCGI ( os, "ms_mass_exclusion" );
	if ( massExclusion ) {
		params->copyToCGI ( os, "ms_min_mass" );
		params->copyToCGI ( os, "ms_max_mass" );
	}
	bool matrixExclusion = params->copyToCGI ( os, "ms_matrix_exclusion" );
	if ( matrixExclusion ) {
		params->copyToCGI ( os, "ms_max_matrix_mass" );
	}
}
void MSPeakFilterOptions::copyToHiddenFormEntry ( ostream& os, const ParameterList* params )
{
	bool peakExclusion = params->copyToHiddenFormEntry ( os, "ms_peak_exclusion" );
	if ( peakExclusion ) {
		params->copyToHiddenFormEntry ( os, "ms_max_peaks" );
		params->copyToHiddenFormEntry ( os, "ms_min_peaks" );
		params->copyToHiddenFormEntry ( os, "ms_min_intensity" );
	}
	bool massExclusion = params->copyToHiddenFormEntry ( os, "ms_mass_exclusion" );
	if ( massExclusion ) {
		params->copyToHiddenFormEntry ( os, "ms_min_mass" );
		params->copyToHiddenFormEntry ( os, "ms_max_mass" );
	}
	bool matrixExclusion = params->copyToHiddenFormEntry ( os, "ms_matrix_exclusion" );
	if ( matrixExclusion ) {
		params->copyToHiddenFormEntry ( os, "ms_max_matrix_mass" );
	}
}
void MSPeakFilterOptions::setDefault ( const MSPeakFilterOptions& ms )
{
	peakExclusionDefault		= ms.peakExclusion;
	maxPeaksDefault				= ms.maxPeaks;
	minPeaksDefault				= ms.minPeaks;
	minIntensityDefault			= ms.minIntensity;
	massExclusionDefault		= ms.massExclusion;
	minMassDefault				= ms.minMass;
	maxMassDefault				= ms.maxMass;
	matrixExclusionDefault		= ms.matrixExclusion;
	maxMatrixMassDefault		= ms.maxMatrixMass;
}
string	MSMSPeakFilterOptions::peakFilterDefault		= "Max MSMS Pks";
bool	MSMSPeakFilterOptions::highResETDDeisotopeDefault	= false;
bool	MSMSPeakFilterOptions::ftPeakExclusionDefault		= false;
bool	MSMSPeakFilterOptions::ECDorETDSideChainExclusionDefault = false;
bool	MSMSPeakFilterOptions::peakExclusionDefault		= false;
unsigned int MSMSPeakFilterOptions::maxPeaksDefault		= 60;
unsigned int MSMSPeakFilterOptions::minPeaksDefault		= 5;
double	MSMSPeakFilterOptions::minIntensityDefault		= 0.0;
bool	MSMSPeakFilterOptions::massExclusionDefault		= false;
double	MSMSPeakFilterOptions::minMassDefault			= 50.0;
double	MSMSPeakFilterOptions::precursorExclusionDefault= 15.0;
double	MSMSPeakFilterOptions::minPrecursorMassDefault	= 0.0;
bool	MSMSPeakFilterOptions::matrixExclusionDefault	= false;
double	MSMSPeakFilterOptions::maxMatrixMassDefault		= 400.0;
bool	MSMSPeakFilterOptions::joinPeaksDefault			= false;
bool	MSMSPeakFilterOptions::deisotopeDefault			= false;
bool	MSMSPeakFilterOptions::deisotopeHiResDefault	= false;

MSMSPeakFilterOptions::MSMSPeakFilterOptions ( const ParameterList* params ) :
	peakFilter			( params->getStringValue( "msms_pk_filter",			peakFilterDefault ) ),
	highResETDDeisotope	( params->getBoolValue	( "msms_high_res_etd_deisotope",highResETDDeisotopeDefault ) ),
	ftPeakExclusion		( params->getBoolValue	( "msms_ft_peak_exclusion",	ftPeakExclusionDefault ) ),
	ECDorETDSideChainExclusion ( params->getBoolValue ( "msms_ecd_or_etd_side_chain_exclusion", ECDorETDSideChainExclusionDefault ) ),
	peakExclusion		( params->getBoolValue	( "msms_peak_exclusion",	peakExclusionDefault ) ),
	maxPeaks			( params->getUIntValue	( "msms_max_peaks",			maxPeaksDefault ) ),
	minPeaks			( params->getUIntValue	( "msms_min_peaks",			minPeaksDefault ) ),
	minIntensity		( params->getDoubleValue( "msms_min_intensity",		minIntensityDefault ) ),
	massExclusion		( params->getBoolValue	( "msms_mass_exclusion",	massExclusionDefault ) ),
	minMass				( params->getDoubleValue( "msms_min_mass",			minMassDefault ) ),
	precursorExclusion	( params->getDoubleValue( "msms_precursor_exclusion", precursorExclusionDefault ) ),
	minPrecursorMass	( params->getDoubleValue( "msms_min_precursor_mass",minPrecursorMassDefault ) ),
	matrixExclusion		( params->getBoolValue	( "msms_matrix_exclusion",	matrixExclusionDefault ) ),
	maxMatrixMass		( params->getDoubleValue( "msms_max_matrix_mass",	maxMatrixMassDefault ) ),
	joinPeaks			( params->getBoolValue	( "msms_join_peaks",		joinPeaksDefault ) ),
	deisotope			( params->getBoolValue	( "msms_deisotope",			deisotopeDefault ) ),
	deisotopeHiRes		( params->getBoolValue	( "msms_deisotope_hi_res",	deisotopeHiResDefault ) )
{
}
void MSMSPeakFilterOptions::copyToCGI ( ostream& os, const ParameterList* params )
{
	params->copyToCGI ( os, "msms_pk_filter" );
	params->copyToCGI ( os, "msms_ft_peak_exclusion" );
	params->copyToCGI ( os, "msms_ecd_or_etd_side_chain_exclusion" );
	bool peakExclusion = params->copyToCGI ( os, "msms_peak_exclusion" );
	params->copyToCGI ( os, "msms_max_peaks" );
	if ( peakExclusion ) {
		params->copyToCGI ( os, "msms_min_peaks" );
		params->copyToCGI ( os, "msms_min_intensity" );
	}
	bool massExclusion = params->copyToCGI ( os, "msms_mass_exclusion" );
	if ( massExclusion ) {
		params->copyToCGI ( os, "msms_min_mass" );
		params->copyToCGI ( os, "msms_precursor_exclusion" );
	}
	params->copyToCGI ( os, "msms_min_precursor_mass" );
	bool matrixExclusion = params->copyToCGI ( os, "msms_matrix_exclusion" );
	if ( matrixExclusion ) {
		params->copyToCGI ( os, "msms_max_matrix_mass" );
	}
	params->copyToCGI ( os, "msms_join_peaks" );
	params->copyToCGI ( os, "msms_deisotope" );
	params->copyToCGI ( os, "msms_deisotope_hi_res" );
}
void MSMSPeakFilterOptions::copyToHiddenFormEntry ( ostream& os, const ParameterList* params )
{
	params->copyToHiddenFormEntry ( os, "msms_pk_filter" );
	params->copyToHiddenFormEntry ( os, "msms_ft_peak_exclusion" );
	params->copyToHiddenFormEntry ( os, "msms_ecd_or_etd_side_chain_exclusion" );
	bool peakExclusion = params->copyToHiddenFormEntry ( os, "msms_peak_exclusion" );
	if ( peakExclusion ) {
		params->copyToHiddenFormEntry ( os, "msms_max_peaks" );
		params->copyToHiddenFormEntry ( os, "msms_min_peaks" );
		params->copyToHiddenFormEntry ( os, "msms_min_intensity" );
	}
	bool massExclusion = params->copyToHiddenFormEntry ( os, "msms_mass_exclusion" );
	if ( massExclusion ) {
		params->copyToHiddenFormEntry ( os, "msms_min_mass" );
		params->copyToHiddenFormEntry ( os, "msms_precursor_exclusion" );
	}
	params->copyToHiddenFormEntry ( os, "msms_min_precursor_mass" );
	bool matrixExclusion = params->copyToHiddenFormEntry ( os, "msms_matrix_exclusion" );
	if ( matrixExclusion ) {
		params->copyToHiddenFormEntry ( os, "msms_max_matrix_mass" );
	}
	params->copyToHiddenFormEntry ( os, "msms_join_peaks" );
	params->copyToHiddenFormEntry ( os, "msms_deisotope" );
	params->copyToHiddenFormEntry ( os, "msms_deisotope_hi_res" );
}
string MSMSPeakFilterOptions::getCommandLineNVPair ( const ParameterList* params )
{
	string s;
	s += params->getCommandLineNVPair ( "msms_pk_filter" );
	s += params->getCommandLineNVPair ( "msms_ft_peak_exclusion" );
	s += params->getCommandLineNVPair ( "msms_ecd_or_etd_side_chain_exclusion" );
	s += params->getCommandLineNVPair ( "msms_max_peaks" );
	string s2 = params->getCommandLineNVPair ( "msms_peak_exclusion" );
	if ( !s2.empty () ) {
		s += s2;
		s += params->getCommandLineNVPair ( "msms_min_peaks" );
		s += params->getCommandLineNVPair ( "msms_min_intensity" );
	}
	string s3 = params->getCommandLineNVPair ( "msms_mass_exclusion" );
	if ( !s3.empty () ) {
		s += s3;
		s += params->getCommandLineNVPair ( "msms_min_mass" );
		s += params->getCommandLineNVPair ( "msms_precursor_exclusion" );
	}
	s += params->getCommandLineNVPair ( "msms_min_precursor_mass" );
	string s4 = params->getCommandLineNVPair ( "msms_matrix_exclusion" );
	if ( !s4.empty () ) {
		s += s4;
		s += params->getCommandLineNVPair ( "msms_max_matrix_mass" );
	}
	s += params->getCommandLineNVPair ( "msms_join_peaks" );
	s += params->getCommandLineNVPair ( "msms_deisotope" );
	s += params->getCommandLineNVPair ( "msms_deisotope_hi_res" );

	return s;
}
void MSMSPeakFilterOptions::setDefault ( const MSMSPeakFilterOptions& msms )
{
	highResETDDeisotopeDefault	= msms.highResETDDeisotope;
	ftPeakExclusionDefault		= msms.ftPeakExclusion;
	ECDorETDSideChainExclusionDefault = msms.ECDorETDSideChainExclusion;
	peakExclusionDefault		= msms.peakExclusion;
	maxPeaksDefault				= msms.maxPeaks;
	minPeaksDefault				= msms.minPeaks;
	minIntensityDefault			= msms.minIntensity;
	massExclusionDefault		= msms.massExclusion;
	minMassDefault				= msms.minMass;
	precursorExclusionDefault	= msms.precursorExclusion;
	minPrecursorMassDefault		= msms.minPrecursorMass;
	matrixExclusionDefault		= msms.matrixExclusion;
	maxMatrixMassDefault		= msms.maxMatrixMass;
	joinPeaksDefault			= msms.joinPeaks;
	deisotopeDefault			= msms.deisotope;
	deisotopeHiResDefault		= msms.deisotopeHiRes;
}
