/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_iso_par.cpp                                                *
*                                                                             *
*  Created    : June 9th 2001                                                 *
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
#include <lg_string.h>
#include <lgen_error.h>
#include <lu_aa_calc.h>
#include <lu_iso_par.h>
#include <lu_getfil.h>
#include <lu_cgi_val.h>
#include <lu_mass_seq.h>
#include <lu_usermod.h>
#include <lu_param_list.h>
#include <lu_iso.h>
#include <lu_mass_conv.h>
#include <lu_prog.h>
using std::string;
using std::ostream;
using std::vector;

MSIsotopeParameters::MSIsotopeParameters ( const ParameterList* params ) :
	MSProgramParameters		( params ),
	aaInitInfo		( params ),
	displayGraph	( params->getBoolValue	( "display_graph" ) ),
	detailedReport	( params->getBoolValue	( "detailed_report" ) ),
	profileType		( params->getStringValue ( "profile_type", "Gaussian" ) ),
	resolution		( params->getDoubleValue ( "resolution", 10000.0 ) ),
	isotopePurity	( params )
{
	if ( resolution > 10000000 ) {
		ErrorHandler::genError ()->error ( "The maximum resolution is 10,000,000.\n" );
	}
	if ( resolution > 200000 && profileType == "Lorentzian" ) {
		ErrorHandler::genError ()->error ( "The maximum resolution for a Lorentzian distribution is 200,000.\n" );
	}
	Usermod::initialiseAllUsermodAAInfo ();
	initIsotopePeakStats ( params );
	IsotopeProfile::setDetailed ( detailedReport );
	isotopeProfile = initIsotopeProfile ();
}
void MSIsotopeParameters::initIsotopePeakStats ( const ParameterList* params )
{
	double maxIntensity = 0.0;
	double minMoverZ;
	double maxMoverZ;
	for ( int i = 0 ; ; i++ ) {
		double mOverZ;
		string num = ( i == 0 ) ? "" : gen_itoa ( i+1 );
		string dType = params->getStringValue ( "distribution_type" + num );
		if ( dType == "Off" ) continue;
		if ( dType == "" ) break;
		int z = params->getIntValue ( "parent_charge" + num, 1 );
		double intensity = params->getDoubleValue ( "parent_intensity" + num, 1.0 );
		distributionType.push_back ( dType );
		parentCharge.push_back ( z );
		intensityString.push_back ( params->getStringValue ( "parent_intensity" + num, "1.0" ) );

		if ( dType == "Peptide Sequence" ) {
			string initSeq = params->getStringValue ( "sequence" + num );
			if ( initSeq.empty () ) {
				ErrorHandler::genError ()->error ( "The peptide sequence field is empty.\n" );
			}
			string nTermName = params->getStringValue ( "nterm" + num );
			string cTermName = params->getStringValue ( "cterm" + num );
			MapStringConstModPtr constMods;
			if ( !nTermName.empty () ) constMods ["n"] = new ConstMod ( nTermName + " (N-term)" );
			if ( !cTermName.empty () ) constMods ["c"] = new ConstMod ( cTermName + " (C-term)" );
			StringVector seq = initSequence ( initSeq );
			AACalculator aaCalc ( true, constMods );
			ElementalFormula ef;
			try {
				ef = aaCalc.calculateElementalComposition ( seq );
			}
			catch ( AACalculatorNoElementalComposition ) {
				ErrorHandler::genError ()->error ( "Can't calculate the elemental composition from the peptide sequence.\n" );
			}
			mOverZ = mPlusHToMOverZ ( formula_to_monoisotopic_mass ( ef.getFormula ().c_str () ) - ELECTRON_REST_MASS, z, true );
			isotopePeakStats.push_back ( new IsotopePeakStats ( ef.getFormula (), z ) );
			string outString;
			if ( !nTermName.empty () ) outString += nTermName + '-';
			outString += initSeq;
			if ( !cTermName.empty () ) outString += '-' + cTermName;
			outputString.push_back ( outString );
		}
		else if ( dType == "Elemental Composition" ) {
			string eFormula = params->getStringValue ( "elemental_composition" + num, "" );
			if ( eFormula.empty () ) {
				ErrorHandler::genError ()->error ( "The elemental formula field is empty.\n" );
			}
			mOverZ = mPlusHToMOverZ ( formula_to_monoisotopic_mass ( eFormula.c_str () ) - ELECTRON_REST_MASS, z, true );
			isotopePeakStats.push_back ( new IsotopePeakStats ( eFormula, z ) );
			outputString.push_back ( eFormula );
		}
		else if ( dType == "Averagine" ) {
			string aMassStr = params->getStringValue ( "averagine_mass" + num );
			if ( aMassStr.empty () ) {
				ErrorHandler::genError ()->error ( "The averagine mass field is empty.\n" );
			}
			double aMass = params->getDoubleValue ( "averagine_mass" + num );
			if ( aMass < 300.0 ) {
				ErrorHandler::genError ()->error ( "The minimum averagine mass is 300 Da.\n" );
			}
			mOverZ = aMass;
			isotopePeakStats.push_back ( new IsotopePeakStats ( aMass, z ) );
			outputString.push_back ( aMassStr );
		}
		minMoverZ = ( i == 0 ) ? mOverZ : genMin ( mOverZ, minMoverZ );
		maxMoverZ = ( i == 0 ) ? mOverZ : genMax ( mOverZ, maxMoverZ );
		maxIntensity = ( i == 0 ) ? intensity : genMax ( intensity, maxIntensity );
		intensityVector.push_back ( intensity );
	}
	//if ( maxMoverZ > minMoverZ + 30.0 ) {
	//	ErrorHandler::genError ()->error ( "Monoisotopic M/Z values can't be separated by more than 30 Th.\n" );
	//}
	for ( int j = 0 ; j < intensityVector.size () ; j++ ) {
		intensityVector [j] /= maxIntensity;
	}
}
IsotopeProfile* MSIsotopeParameters::initIsotopeProfile ()
{
	IsotopeProfile* ip = 0;
	if ( profileType == "Gaussian" ) {
		ip = new GaussianIsotopeProfile ( isotopePeakStats, intensityVector, resolution );
	}
	if ( profileType == "Lorentzian" ) {
		ip = new LorentzianIsotopeProfile ( isotopePeakStats, intensityVector, resolution );
	}
	if ( profileType == "Stick" ) {
		ip = new StickIsotopeProfile ( isotopePeakStats );
	}
	return ip;
}
void MSIsotopeParameters::printHTML ( ostream& os ) const
{
	for ( int i = 0 ; i < distributionType.size () ; i++ ) {
		string num = distributionType.size () == 1 ? "" : " " + gen_itoa ( i+1 );
		string dist = distributionType [i];

		ParameterList::printHTML ( os, "Distribution Type" + num, dist );
		if ( dist == "Peptide Sequence" ) {
			ParameterList::printHTML ( os, "Sequence", outputString [i] );
		}
		else if ( dist == "Elemental Composition" ) {
			ParameterList::printHTML ( os, "Elemental Formula", outputString [i] );
		}
		else if ( dist == "Averagine" ) {
			ParameterList::printHTML ( os, "Averagine Mass", outputString [i] );
		}
		ParameterList::printHTML ( os, "Parent Charge" + num, parentCharge [i] );
		ParameterList::printHTML ( os, "Intensity" + num, intensityString [i] );
	}
	ParameterList::printHTML ( os, "Profile Type", profileType );
	if ( profileType != "Stick" ) ParameterList::printHTML ( os, "Resolution", resolution );
}
void MSIsotopeLink::putCGI ( ostream& os ) const
{
	const ParameterList* params = ProgramLink::getParams ();
	printCGIString ( os, "search_name", "msisotope" );
	printCGIString ( os, "output_type", "HTML" );
	printCGIString ( os, "report_title", "MS-Isotope" );
	printCGI ( os, "display_graph", true );
	printCGIString ( os, "version", Version::instance ().getVersion () );
	printCGIString ( os, "distribution_type", "Elemental Composition" );
	AAInitInfo::copyToCGI ( os, params );
	IsotopePurity::copyToCGI ( os, params );
}
void MSIsotopeLink::write ( ostream& os, ElementalFormula& elementalComposition, const int charge ) const
{
	ProgramLink::openLink ( os, "msisoLink", -1 );
	printCGI ( os, "parent_charge", charge );
	printCGIString ( os, "elemental_composition", elementalComposition.getFormula () );
	os << "\\\">" << "<nobr>" << elementalComposition << "</nobr>";
	ProgramLink::closeLink ( os );
}
void MSIsotopeLink::printHTML ( ostream& os ) const
{
	os << "msisoLink" << "=\"";
	os << ProgramLink::getURLStart ( "mssearch" );
	os << "?";
	putCGI ( os );
	os << "\";\n";
}
