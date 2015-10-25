/******************************************************************************
*                                                                             *
*  Program    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_prod_srch.cpp                                              *
*                                                                             *
*  Created    : October 22nd 2001                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifndef VIS_C
#include <stdexcept>
#endif
#include <lu_aa_calc.h>
#include <lu_html.h>
#include <lu_getfil.h>
#include <lu_prod_form.h>
#include <lu_prod_srch.h>
#include <lu_mass_conv.h>
#include <lu_mass_frag.h>
#include <lu_iso_par.h>
#include <lu_table.h>
#include <lu_app_gr.h>
#include <lu_charge.h>
#include <lu_param_list.h>
#include <lu_form_valid.h>
#ifdef RAW_DATA
#include <lu_proj_file.h>
#include <lr_main.h>
#endif
using std::ostream;
using std::string;
using std::endl;
using std::set;
using std::pair;
using std::vector;
using std::runtime_error;

MSProductSearch::MSProductSearch ( MSProductParameters& params ) :
	MSProgram ( params ),
	prodParams ( params ),
	peaks ( 0 ),
	rawDataOption ( false )
{
	PeakContainer::setMultiChargeAssign ( false );
	MSMSDataSetInfo* dsi = params.getDataSetInfo ();
	if ( dsi->getNumDataSets () > 0 ) {
		MSMSDataPoint* file = dsi->getDataSet ( 0 );
		precursorMOverZ = file->getPrecursorMZ ();
		Peak parentPeak ( file->getPrecursorMZ (), file->getPrecursorTolerance (), file->getPrecursorCharge (), file->getPrecursorIntensity (), getAdductMass ( prodParams.getMonoisotopicFlag () ), prodParams.getAverageParentMonoFragments () );
		MSMSPeakFilterOptions* mmpfo = prodParams.getMSMSPeakFilterOptions ();
		peaks.push_back ( new PeakContainer ( file, mmpfo, &parentPeak, prodParams.getParentMassTolerance (), prodParams.getProductMassTolerance (), prodParams.getMonoisotopicFlag (), prodParams.getAverageParentMonoFragments () ) );
		setMSProductPrecursorCharge ( file->getPrecursorCharge () );
	}
	for ( int i = 0 ; ; i++ ) {
		StringVector s = params.getSequence ( i );
		if ( s.empty () ) break;
		AACalculator aaCalc ( prodParams.getMonoisotopicFlag (), params.getConstMods ( i ) );
		if ( params.getShow ( i ) ) pInfo.push_back ( new PepInfo ( s, aaCalc.getNTerminusWt (), aaCalc.getNTerminusString (), aaCalc.getCTerminusWt (), aaCalc.getCTerminusString (), params.getShow ( i ) ) );
	}
	int z = params.getMaxCharge ();
	IntVector zArray;
	for ( int j = 0 ; j < pInfo.size () ; j++ ) {
		zArray.push_back ( z );
	}
	biemannFragments = new BiemannFragments ( pInfo, params.getBiemannParameters (), peaks, zArray, params.getCalibrate (), params.getCalTolerance (), params.getAlternative (), params.getLinkInfo () );
#ifdef RAW_DATA
	rawDataOption = !ProjectFile::getRawTypes ( ProgramLink::getParams () ).empty ();
#endif
}
MSProductSearch::~MSProductSearch ()
{
	for ( int i = 0 ; i < peaks.size () ; i++ ) {
		delete peaks [i];
	}
	delete biemannFragments;
	for ( int j = 0 ; j < pInfo.size () ; j++ ) {
		delete pInfo [j];
	}
}
void MSProductSearch::printBodyHTML ( ostream& os )
{
	prodParams.printHTML ( os );
	if ( !peaks.empty () ) {
		if ( !peaks [0]->empty () ) {
			vector <XYData> vXYData;
			bool plotCentroids = true;
#ifdef RAW_DATA
			if ( rawDataOption && prodParams.getRawData () ) {
				vXYData.resize ( 1 );
				try {
					RawFile rf ( ProgramLink::getParams (), prodParams.getSpecID ().getFraction (), 0.0, 0.0 );
					rf.getXYData ( vXYData, prodParams.getSpecID (), precursorMOverZ );
				}
				catch ( runtime_error e ) {
					ErrorHandler::genError ()->error ( e );
				}
				plotCentroids = prodParams.getPlotCentroids ();
			}
#endif
			biemannFragments->printApplet ( os, vXYData, plotCentroids, prodParams.getDiscriminating () );
			if ( biemannFragments->getSpectrumMatch () ) printForm ( os );
		}
		else {
			os << "<p>" << endl;
			os << "<b>No fragment peaks, graph not displayed.</b>" << endl;
			os << "</p>" << endl;
		}
	}
	int i = 0;
	SetString ss;
	SetString sd;
	int precision = instInf->getFragmentPeakPrecision ().getMassDecimalPlaces ();
	string precisionStr = "%." + gen_itoa ( precision ) + "f";
	for ( ; ; i++ ) {
		if ( prodParams.getShow ( i ) ) {
			StringVector sequence = biemannFragments->getPeptide ( i );
			if ( sequence.empty () ) break;
			AACalculator monoAACalc ( true, prodParams.getConstMods ( i ) );
			bool eForm = true;
			ElementalFormula elemComp;
			bool insert = false;
			try {
				elemComp = monoAACalc.calculateElementalComposition ( sequence );
				PairSetStringIteratorBool pssib = ss.insert ( elemComp.getFormula () );
				insert = pssib.second;
			}
			catch ( AACalculatorNoElementalComposition ) {
				eForm = false;
			}
			int maxCharge = prodParams.getMaxCharge ();
			if ( eForm && insert ) printIsotopeLink ( os, elemComp, maxCharge );
			double peptideMassMi = monoAACalc.calculatePeptideMW ( sequence );
			AACalculator avAACalc ( false, prodParams.getConstMods ( i ) );
			double peptideMassAv = avAACalc.calculatePeptideMW ( sequence ); 
			modified_mass_convert ( prodParams.getMonoisotopicFlag () );
			string massStr = gen_ftoa ( peptideMassAv, precisionStr.c_str () );
			if ( sd.insert ( massStr ).second ) printPrecursor ( os, maxCharge, peptideMassMi, peptideMassAv, precision );
		}
	}
	if ( i > 0 ) {
		if ( !peaks.empty () ) biemannFragments->printSpectrumMatchHTML ( os, prodParams.getDiscriminating () );
		biemannFragments->printHTML ( os );
	}
}
void MSProductSearch::printJavascriptFunctions ( ostream& os )
{
	basicJavascriptFunctions ( os );
	startJavascript ( os );
	os << "function showLinkSearchTypeItems ( item ) {" << endl;
	os << "\t" << "var val = getSelectValue(item);" << endl;
	os << "\t" << "if ( val == 'User Defined Link' ) {" << endl;
	os << "\t\t" << "showdiv ( 'div_ls' );" << endl;
	os << "\t" << "}" << endl;
	os << "\t" << "else {" << endl;
	os << "\t\t" << "hidediv ( 'div_ls' );" << endl;
	os << "\t" << "}" << endl;
	os << "}" << endl;
	endJavascript ( os );
}
void MSProductSearch::printRawDataValidationJavascript ( ostream& os ) const
{
	startJavascript ( os );
		os << "function validateCheckRawData( form ) {" << endl;
		os << "\t" << "var rd = false;" << endl;
		os << "\t" << "if(form.data_plotted.value == \"Raw\" || form.data_plotted.value == \"Raw and Centroid\"){" << endl;
		os << "\t\t" << "rd = true;" << endl;
		os << "\t" << "}" << endl;
		os << "\t" << "if(rd){" << endl;
		os << "\t\t" << "form.action = \"" << ProgramLink::getURLStart ( "mssearchRawData" ) << "\";" << endl;
		os << "\t" << "}" << endl;
		os << "\t" << "else{" << endl;
		os << "\t\t" << "form.action = \"" << ProgramLink::getURLStart ( "mssearch" ) << "\";" << endl;
		os << "\t" << "}" << endl;
		os << "\t" << "validateForms(form);" << endl;
		os << "}" << endl;
	endJavascript ( os );
}
void MSProductSearch::printPrecursor ( ostream& os, int maxCharge, double peptideMassMi, double peptideMassAv, int precision ) const
{
	tableStart ( os, true );
		tableRowStart ( os );
			for ( int c = 1 ; c <= maxCharge ; c++ ) {
				tableHeaderStart ( os );
					os << "MH<sup>+" << c << "</sup>(av)" << endl;
				tableHeaderEnd ( os );
				tableHeaderStart ( os );
					os << "MH<sup>+" << c << "</sup>(mono)" << endl;
				tableHeaderEnd ( os );
			}
		tableRowEnd ( os );
		tableRowStart ( os );
			for ( int charge = 1 ; charge <= maxCharge ; charge++ ) {
				tableCell ( os, mPlusHToMOverZ ( peptideMassAv, charge, false ), precision );
				tableCell ( os, mPlusHToMOverZ ( peptideMassMi, charge, true ), precision );
			}
		tableRowEnd ( os );
	tableEnd ( os );
}
void MSProductSearch::printIsotopeLink ( ostream& os, ElementalFormula& elemComp, int maxCharge ) const
{
	MSIsotopeLink isotopeLink;
	startJavascript ( os );
	isotopeLink.printHTML ( os );
	endJavascript ( os );
	os << "Elemental Composition: ";
	isotopeLink.write ( os, elemComp, maxCharge );
	os << "<br />" << endl;
}
void MSProductSearch::printBodyTabDelimitedText ( ostream& os )
{
	if ( !peaks.empty () ) biemannFragments->printSpectrumMatchTabDelimitedText ( os, prodParams.getDiscriminating () );
}
void MSProductSearch::printBodyMGFSpectrum ( ostream& os )
{
	string filename = DataReader::getLastFilename ();
	string version = DataReader::getLastVersion ();
	if ( !filename.empty () ) {
		MSMSPeakListDataSetInfo dsi ( filename );
		dsi.writePeakList ( os, prodParams.getSpecID (), version );
	}
}
void MSProductSearch::setLinkSearchVisualizationFlags ( const string& val, bool& div_ls )
{
	if ( val == "User Defined Link" ) {
		div_ls = true;
	}
	else {
		div_ls = false;
	}
}
void MSProductSearch::printForm ( ostream& os )
{
	int numSequences = genMax ( prodParams.getNumSequences (), InfoParams::instance ().getIntValue ( "max_msprod_sequences", 2 ) );
	ParameterList paramList = *ProgramLink::getParams ();
	FormValidatingJavascript fvj;
	vector <FormItemText> fi;
	vector <FormItemText> fiN;
	vector <FormItemText> fiC;
	vector <FormItemText> fiL;
	vector <FormItemCheckbox> fiCB;
	int numCols = 5;
	for ( int i = 0 ; i < numSequences ; i++ ) {
		string num = ( i == 0 ) ? "" : gen_itoa ( i+1 );
		fi.push_back ( FormItemText ( "", "", "sequence" + num, 70, 1000, "" ) );
		fiN.push_back ( FormItemText ( "", "", "nterm" + num, 9, 100, "" ) );
		fiC.push_back ( FormItemText ( "", "", "cterm" + num, 9, 100, "" ) );
		fiL.push_back ( FormItemText ( "", "", "nloss" + num, 9, 100, "" ) );
		fiCB.push_back ( FormItemCheckbox ( "", "", "s" + num, true ) );
	}
	FormItemMaxCharge fi2;
	FormItemCountPosZ fi2a;
	FormItemMaxInternalLen fi2b (&fvj);
	FormItemMSMSPkFilter fi3a ( true );
	FormItemMaxMSMSPeaks fi3 ( &fvj, true );
	FormItemFragmentMassesTolerance fi4 (&fvj);
	FormItemFragmentMassesToleranceUnits fi5 ( FormItemFragmentMassesToleranceUnits::getName () );
	FormItemAlternative fi5a;
	FormItemDiscriminating fi5b;
	FormItemLinkSearchType fi5c ( "No Link", "showLinkSearchTypeItems( this.form." + FormItemLinkSearchType::getName () + " )" );
	FormItemComposition fi5d ( "Bridge", "bridge", 1, "", &fvj );
	FormItemMaxLosses fi6;
	FormItemMultiZInternal fi7;
	FormItemCalibrate fi8;
	string tolUnits = paramList.getStringValue ( FormItemFragmentMassesToleranceUnits::getName () );
	string tolVal = paramList.getStringValue ( FormItemFragmentMassesTolerance::getName () );
	FormItemCalTolerance fi9 ( tolUnits, tolVal );
	FormItemDataPlotted fi10;
	FormItemOutputType fi12 ( false, false );
	FormItemUseInstrumentIonTypes fi13 ( true );
	bool xLinkingFlag = paramList.getStringValue ( FormItemLinkSearchType::getName () ) != "No Link" && paramList.getStringValue ( FormItemLinkSearchType::getName () ) != "";
	tableStart ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center" );
				biemannFragments->printErrorChart ( os );
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center" );
				fvj.print ( os );
				if ( rawDataOption && InfoParams::instance ().getBoolValue ( "raw_data_forwarding" ) ) {
					printRawDataValidationJavascript ( os );
					os << "<form";
					os << " ";
					os << "method=\"post\"";
					os << " ";
					os << "onsubmit=\"return validateCheckRawData(this)\"";
					os << ">" << endl;
				}
				else printHTMLFORMStart ( os, "post", "mssearch", false, true );
				printJavascriptFunctions ( os );
				tableStart ( os, true );
					ExpandableJavascriptBlock* ejb = 0;
					for ( int j = 0 ; j < numSequences ; j++ ) {
						string num = ( j == 0 ) ? "" : gen_itoa ( j+1 );
						if ( j == 1 && !xLinkingFlag ) {
							tableRowStart ( os );
								tableHeaderStart ( os, "", "center", false, numCols );
									bool open = !paramList.getStrippedStringValue ( "sequence" + num ).empty ();
									ejb = new ExpandableJavascriptBlock ( "Additional Sequences", open );
									ejb->printHeader ( os );
									tableStart ( os, false );
						}
						if ( j == 0 ) {
							tableRowStart ( os );
								tableHeader ( os, "" );
								tableHeader ( os, "N Term" );
								tableHeader ( os, "Sequence" );
								tableHeader ( os, "C Term" );
								tableHeader ( os, "N Loss" );
							tableRowEnd ( os );
						}
						tableRowStart ( os );
							tableHeaderStart ( os );
								fiCB [j].setValue ( &paramList );
								fiCB [j].printHTML ( os );
								paramList.removeName ( "s" + num );
							tableHeaderEnd ( os );
							tableHeaderStart ( os );
								fiN [j].setValue ( &paramList );
								fiN [j].printHTML ( os );
								paramList.removeName ( "nterm" + num );
							tableHeaderEnd ( os );
							tableHeaderStart ( os );
								fi [j].setValue ( &paramList );
								fi [j].printHTML ( os );
								paramList.removeName ( "sequence" + num );
							tableHeaderEnd ( os );
							tableHeaderStart ( os );
								fiC [j].setValue ( &paramList );
								fiC [j].printHTML ( os );
								paramList.removeName ( "cterm" + num );
							tableHeaderEnd ( os );
							tableHeaderStart ( os );
								fiL [j].setValue ( &paramList );
								fiL [j].printHTML ( os );
								paramList.removeName ( "nloss" + num );
							tableHeaderEnd ( os );
						tableRowEnd ( os );
						if ( xLinkingFlag && genIsOdd (j) && j != numSequences - 1 ) {
							tableRowStart ( os );
								tableHeaderStart ( os, "", "center", false, numCols );
									printHTMLFORMLabel ( os, "&nbsp;", "" );
								tableHeaderEnd ( os );
							tableRowEnd ( os );
						}
						if ( j > 0 && j == numSequences - 1 ) {
										bool div_ls;
										setLinkSearchVisualizationFlags ( paramList.getStringValue ( FormItemLinkSearchType::getName () ), div_ls );
										tableRowStart ( os );
											tableHeaderStart ( os, "", "left", false, numCols );
												tableStart ( os, false );
													tableRowStart ( os );
														tableHeaderStart ( os, "", "left", true );
															fi5a.setValue ( &paramList );
															fi5a.printHTML ( os );
															paramList.removeName ( FormItemAlternative::getName () );
														tableHeaderEnd ( os );
														tableHeaderStart ( os, "", "left", true );
															fi5b.setValue ( &paramList );
															fi5b.printHTML ( os );
															paramList.removeName ( FormItemDiscriminating::getName () );
														tableHeaderEnd ( os );
														tableHeaderStart ( os, "", "left", true );
															fi5c.setValue ( &paramList );
															fi5c.printHTML ( os );
															paramList.removeName ( FormItemLinkSearchType::getName () );
														tableHeaderEnd ( os );
														tableHeaderStart ( os, "", "left", true );
															divStart ( os, "div_ls", div_ls );
																fi5d.setValue ( &paramList );
																fi5d.printHTML ( os );
																paramList.removeName ( FormItemComposition::getName ("bridge") );
															divEnd ( os );
														tableHeaderEnd ( os );
													tableRowEnd ( os );
												tableEnd ( os );
											tableHeaderEnd ( os );
										tableRowEnd ( os );
							if ( ejb ) {
											tableEnd ( os );
											ejb->printFooter ( os );
											delete ejb;
									tableHeaderEnd ( os );
								tableRowEnd ( os );
							}
						}
					}
					tableRowStart ( os );
						tableHeaderStart ( os, "", "center", true, numCols );
							fi2.setValue ( &paramList );
							fi2.printHTML ( os );
							paramList.removeName ( FormItemMaxCharge::getName () );
							fi2a.setValue ( &paramList );
							fi2a.printHTML ( os );
							paramList.removeName ( FormItemCountPosZ::getName () );
							fi6.setValue ( &paramList );
							fi6.printHTML ( os );
							paramList.removeName ( FormItemMaxLosses::getName () );
							fi7.setValue ( &paramList );
							fi7.printHTML ( os );
							paramList.removeName ( FormItemMultiZInternal::getName () );
							fi2b.setValue ( &paramList );
							fi2b.printHTML ( os );
							paramList.removeName ( FormItemMaxInternalLen::getName () );
						tableHeaderEnd ( os );
					tableRowEnd ( os );
					tableRowStart ( os );
						tableHeaderStart ( os, "", "center", false, numCols );
							fi3a.setValue ( &paramList );
							fi3a.printHTML ( os );
							paramList.removeName ( FormItemMSMSPkFilter::getName () );
							fi3.setValue ( &paramList );
							fi3.printHTML ( os );
							paramList.removeName ( FormItemMaxMSMSPeaks::getName () );
							fi4.setValue ( &paramList );
							fi4.printHTML ( os );
							paramList.removeName ( FormItemFragmentMassesTolerance::getName () );
							fi5.setValue ( &paramList );
							fi5.printHTML ( os );
							paramList.removeName ( FormItemFragmentMassesToleranceUnits::getName () );
							fi8.setValue ( &paramList );
							fi8.printHTML ( os );
							paramList.removeName ( FormItemCalibrate::getName () );
							fi9.setValue ( &paramList );
							fi9.printHTML ( os );
							paramList.removeName ( FormItemCalTolerance::getName () );
						tableHeaderEnd ( os );
					tableRowEnd ( os );
					tableRowStart ( os );
						tableHeaderStart ( os, "", "center", false, numCols );
							tableStart ( os );
								tableRowStart ( os );
									if ( rawDataOption ) {
										tableHeaderStart ( os );
											fi10.setValue ( &paramList );
											fi10.printHTML ( os );
											paramList.removeName ( FormItemDataPlotted::getName () );
										tableHeaderEnd ( os );
									}
									tableHeaderStart ( os );
										printHTMLFORMSubmit ( os, "MS-Product" );
									tableHeaderEnd ( os );
									tableHeaderStart ( os );
										fi12.setValue ( &paramList );
										fi12.printHTML ( os );
										paramList.removeName ( FormItemOutputType::getName () );
									tableHeaderEnd ( os );
								tableRowEnd ( os );
							tableEnd ( os );
						tableHeaderEnd ( os );
					tableRowEnd ( os );
					tableRowStart ( os );
						tableHeaderStart ( os, "", "center", false, numCols );
							bool open2 = !paramList.getBoolValue ( FormItemUseInstrumentIonTypes::getName () );
							ExpandableJavascriptBlock ejb2 ( "Ion Types", open2 );
							ejb2.printHeader ( os );
							fi13.setValue ( &paramList );
							fi13.printHTML ( os );
							paramList.removeName ( FormItemUseInstrumentIonTypes::getName () );
							VectorConstParameterListPtr p1;
							p1.push_back ( &paramList );
							IonTypeForm itf ( p1 );
							itf.printHTML ( os );
							itf.removeName ( &paramList );
							ejb2.printFooter ( os );
						tableHeaderEnd ( os );
					tableRowEnd ( os );
				tableEnd ( os );
				paramList.copyToHiddenFormEntry ( os );
				printHTMLFORMEnd ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
	tableEnd ( os );
}
void MSProductSearch::printBodyXML ( ostream& os )
{
	biemannFragments->printXML ( os, prodParams.getDiscriminating () );
}
