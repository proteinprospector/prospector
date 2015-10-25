/******************************************************************************
*                                                                             *
*  Program    : msdisplay                                                     *
*                                                                             *
*  Filename   : msdisplay_main.cpp                                            *
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
#ifndef VIS_C
#include <stdexcept>
#endif
#include <lu_aa_calc.h>
#include <lu_proj_file.h>
#include <lu_param_list.h>
#include <lu_app_gr.h>
#include <lu_html.h>
#include <lu_usermod.h>
#include <lu_srch_form.h>
#include <lu_form_valid.h>
#include <lu_version.h>
#include <lu_file_type.h>
#include <lu_table.h>
#include <lu_getfil.h>
#include <lr_main.h>
#include <msdisplay_params.h>
using std::cout;
using std::string;
using std::ostream;
using std::endl;
using std::pair;
using std::vector;
using std::ostringstream;
using std::showpos;
using std::noshowpos;
using std::runtime_error;
using namespace FileTypes;

static void writeSpectrumReport ( ostream& os, const MSDisplayParameters& p, const ParameterList* params );
static void writeQuantitationReport ( ostream& os, const MSDisplayParameters& p, const ParameterList* params );
static void calibrate ( XYData& xyData, const ParameterList* params );

static pair <StringVector, string> sTimes;

int main ( int argc, char** argv )
{
	ostream& os = cout;
	initialiseProspector ();
	ParameterList paramList ( argc, argv );
	try {
		if ( paramList.empty () ) {
			ErrorHandler::genError ()->error ( "No parameters passed to Prospector program.\n" );
		}
		//string ver = paramList.getStringValue ( "version" );
		//if ( !ver.empty () && ver != Version::instance ().getVersion () ) {
		//	ErrorHandler::genError ()->error ( "Form version mismatch. Try reloading the form.\n" );
		//}
		ProgramLink::setParams ( &paramList );
		init_html ( os, "" );
		Usermod::initialiseAllUsermodAAInfo ();
		MSDisplayParameters p ( &paramList );
		if ( p.getMOverZ () != 0.0 && p.getRawType () != "MS Full Scan" ) {
			os << "<div class=\"msdisplay_header\">" << endl;
			genPrint ( os, p.getMOverZ (), 4 );
			if ( p.getCharge () != 1 ) {
				os << "<sup>" << showpos << p.getCharge () << noshowpos << "</sup>" << endl;
			}
			os << "</div>" << endl;
		}
		if ( p.getRawType () == "MS Precursor" || p.getRawType () == "Quantitation" )
			writeQuantitationReport ( os, p, &paramList );
		else
			writeSpectrumReport ( os, p, &paramList );

		FormValidatingJavascript fvj;
		FormItem* fi = new FormItemRawType ( p.getRawType () );
		FormItem* fi2;
		FormItem* fi2a;
		FormItem* fi2b;
		if ( !sTimes.first.empty () ) {
			fi2 = new FormItemSelect ( "RT", "", "spot_number", sTimes.first, sTimes.second );
			fi2a = new FormItemText ( "RT Interval (secs)", "", "rt_int_start", 5, 10, p.getRTIntervalStart (), fvj.addSignedFloatingPointValidator ( "rt_int_start", "RT Interval Start" ) );
			fi2b = new FormItemText ( "to", "", "rt_int_end", 5, 10, p.getRTIntervalEnd (), fvj.addSignedFloatingPointValidator ( "rt_int_end", "RT Interval End" ) );
		}
		FormItem* fi3 = new FormItemIsotopePurity ( "C", 13, &fvj, p.getIsotopePurity ().getPercentC13 () );
		FormItem* fi4 = new FormItemIsotopePurity ( "N", 15, &fvj, p.getIsotopePurity ().getPercentN15 () );
		FormItem* fi5 = new FormItemIsotopePurity ( "O", 18, &fvj, p.getIsotopePurity ().getPercentO18 () );
		const char* zOptions [] = { "8", "7", "6", "5", "4", "3", "2", "1", 0 };
		FormItem* fi6 = new FormItemText ( "Precursor m/z", "", "m_over_z", 10, 12, gen_ftoa ( p.getMOverZ (), "%.4f" ), fvj.addPositiveFloatingPointValidator ( "m_over_z", "Precursor m/z" ) );
		FormItem* fi7 = new FormItemSelect ( "Charge", "", "charge", zOptions, gen_itoa ( p.getCharge () ) );
		FormItem* fi8 = new FormItemText ( "Resolution", "", "resolution", 8, 10, gen_ftoa ( p.getResolution (), "%.1f" ), fvj.addPositiveFloatingPointOrExponentValidator ( "resolution", "Resolution" ) );
		fvj.print ( os );
		printHTMLFORMStart ( os, "post", "msdisplay", false, true );

		paramList.removeName ( "raw_type" );
		paramList.removeName ( "percent_C13" );
		paramList.removeName ( "percent_N15" );
		paramList.removeName ( "percent_O18" );
		paramList.removeName ( "charge" );
		paramList.removeName ( "m_over_z" );
		paramList.removeName ( "resolution" );
		tableStart ( os, true );
			tableRowStart ( os );
				tableHeaderStart ( os, "", "center", true );
					fi->printHTML ( os );
					if ( !sTimes.first.empty () ) {
							paramList.removeName ( "spot_number" );
							fi2->printHTML ( os );
							paramList.removeName ( "rt_int_start" );
							fi2a->printHTML ( os );
							paramList.removeName ( "rt_int_end" );
							fi2b->printHTML ( os );
					}
				tableHeaderEnd ( os );
			tableRowEnd ( os );
			tableRowStart ( os );
				tableHeaderStart ( os, "", "center", true );
					fi3->printHTML ( os );
					fi4->printHTML ( os );
					fi5->printHTML ( os );
					fi6->printHTML ( os );
					fi7->printHTML ( os );
					fi8->printHTML ( os );
				tableHeaderEnd ( os );
			tableRowEnd ( os );
				paramList.copyToHiddenFormEntry ( os );

			tableRowStart ( os );
				tableHeaderStart ( os, "", "center", true );
					printHTMLFORMSubmit ( os, "Draw Spectrum" );
				tableHeaderEnd ( os );
			tableRowEnd ( os );
		tableEnd ( os );

		os << "</form>" << endl << endl;

		printProgramInformationHTML ( os, "MS-Display" );
	}
	catch ( runtime_error e ) {
		paramList.writeLogError ( e.what () );
	}
	return 0;
}
static void writeSpectrumReport ( ostream& os, const MSDisplayParameters& p, const ParameterList* params )
{
	vector <XYData> vXYData ( 1 );

	double startMass = p.getDisplayAllMasses () ? 0.0 : p.getDisplayStartMass ();
	double endMass = p.getDisplayAllMasses () ? 0.0 : p.getDisplayEndMass ();
	bool msFullScan = ( p.getRawType () == "MS Full Scan" );
	SpecID specID = p.getSpecID ();
	try {
		RawFile rf ( params, specID.getFraction (), params->getDoubleValue ( "rt_int_start" ), params->getDoubleValue ( "rt_int_end" ) );
		rf.getXYData ( vXYData, sTimes, msFullScan, specID, p.getMOverZ (), startMass, endMass );
	}
	catch ( runtime_error e ) {		// Catch database login problems
		ErrorHandler::genError ()->error ( e );
	}
	if ( msFullScan ) calibrate ( vXYData [0], params );
	double xMin = vXYData [0].minX ();
	double xMax = vXYData [0].maxX ();
	for ( int i = 1 ; i < vXYData.size () ; i++ ) {
		xMin = genMin ( xMin, vXYData [i].minX () );
		xMax = genMax ( xMax, vXYData [i].maxX () );
	}
	for ( int j = 0 ; j < vXYData.size () ; j++ ) {
		if ( vXYData [j].size () ) {
			FileGraphData graphData ( vXYData [j], true );
			SpectrumGraph s ( "pr_graph.par.txt" );
			s.drawGraph ( os, graphData, true, xMin, xMax );
		}
		else
			ErrorHandler::genError ()->error ( "No data within selected mass range.\n" );
	}
}
static void writeQuantitationReport ( ostream& os, const MSDisplayParameters& p, const ParameterList* params )
{
	vector <XYData> vXYData ( 1 );
	AACalculator* aaCalc = new AACalculator ( true, p.getAAInitInfo ().getConstMods () );
	PeakFit::setGraphs ( true );
	if ( p.getRawType () == "Quantitation" && !isQuanMSMS ( p.getQuanType () ) ) QuantitationRatio::setQuanResidue ( p.getQuanType () );
	PeakFitData::setReportPeakIntensity	( true );
	PeakFitData::setReportPeakSNR		( true );
	PeakFitData::setReportPeakResolution( true );
	PeakFitData::setReportPeakFWHM		( true );
	PeakFitData::setReportPeakArea		( true );
	PeakFitData::setReportNoiseMean		( true );
	PeakFitData::setReportStdDev		( true );
	PeakFitData::setReportFormulaString	( true );
	PeakFitData::setSNRThreshold		( p.getSNRThreshold () );
	PeakFitData::setAreaThreshold		( p.getAreaThreshold () );
	PeakFitData::setIntensityThreshold	( p.getIntensityThreshold () );
	QuantitationRatio::setReportActualLightHeavyIntensityRatio		( true );
	QuantitationRatio::setReportActualLightHeavyAreaRatio			( true );
	double mOverZ = p.getMOverZ ();
	int charge = p.getCharge ();
	bool qFlag = p.getRawType () == "Quantitation";
	bool quanMSMS = qFlag && isQuanMSMS ( p.getQuanType () );
	double startMass = mOverZ - 6.0;
	double endMass = mOverZ + 10.0;
	if ( p.getQuanType () == "Label:15N" ) {
		startMass -= 20.0;
		endMass += 20.0;
	}
	if ( quanMSMS ) {
		QuantitationMulti::init ( p.getQuanType () );
		startMass = QuantitationMulti::getStartDisplayMass ();
		endMass = QuantitationMulti::getEndDisplayMass ();
	}
	else if ( qFlag ) {
		QuantitationRatio::getDataRange ( p.getQPeptide (), mOverZ, charge, startMass, endMass );
	}
	SpecID specID = p.getSpecID ();
	try {
		ProjectFile projectFile ( params );
		string version = projectFile.getProjectVersion ();
		RawFile rf ( params, specID.getFraction (), params->getDoubleValue ( "rt_int_start" ), params->getDoubleValue ( "rt_int_end" ) );
		rf.getQuantitationXYData ( vXYData, sTimes, quanMSMS, p.getSpecID (), mOverZ, version, startMass, endMass );

		if ( vXYData [0].size () ) {
			if ( !quanMSMS ) calibrate ( vXYData [0], params );
			PeakFitData* pfd;
			if ( quanMSMS ) {
				if ( rf.getRawType () == WIFF || rf.getRawType () == RAW )
					pfd = new QuantitationMultiMassWindow ( vXYData [0] );
				else
					pfd = new QuantitationMultiNormal ( vXYData [0], p.getResolution () );
			}
			else {
				ElementalFormula ef;
				ElementalFormula* efp;
				bool flag;
				try {
					if ( p.isXlink () ) {
						flag = aaCalc->calculateStrippedElementalCompositionWithTerminii ( p.getFormula (), p.getNTerm (), p.getCTerm (), p.getNLoss (), ef, false );
						ElementalFormula ef2;
						flag = aaCalc->calculateStrippedElementalCompositionWithTerminii ( p.getFormula2 (), p.getNTerm2 (), p.getCTerm2 (), p.getNLoss (), ef2, false );
						ef += ef2;
						ef += p.getBridgeFormula ();
						ef -= "H";
					}
					else {
						flag = aaCalc->calculateStrippedElementalCompositionWithTerminii ( p.getFormula (), p.getNTerm (), p.getCTerm (), p.getNLoss (), ef );
					}
					efp = &ef;
				}
				catch ( AACalculatorNoElementalComposition ) {
					efp = 0;
				}
				if ( qFlag ) {
					pfd = new QuantitationRatio ( vXYData [0], mOverZ, charge, p.getQPeptide (), efp, flag, p.getResolution (), p.getNumPeaks () );
				}
				else {
					pfd = new ParentData ( vXYData [0], mOverZ, charge, efp, flag, p.getResolution (), p.getNumPeaks () );
				}
			}
			pfd->printHTML ( os );
			delete pfd;
		}
		else
			throw runtime_error ( "No data within selected mass range." );
	}
	catch ( runtime_error e ) {
		ErrorHandler::genError ()->error ( e );
	}
}
static void calibrate ( XYData& xyData, const ParameterList* params )
{
	ProjectFile pf ( params );
	double offset = pf.getOffset ( params->getIntValue ( "fraction", 1 ) - 1 );
	string tolUnits = params->getStringValue ( "msms_parent_mass_tolerance_units" );
	if ( offset != 0.0 ) {
		xyData.offsetCalibration ( offset, tolUnits );
	}
	else {
		double sysError = params->getDoubleValue ( "msms_parent_mass_systematic_error", 1 );
		if ( sysError != 0.0 ) {
			xyData.offsetCalibration ( sysError, tolUnits );
		}
	}
}
