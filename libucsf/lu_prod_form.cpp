/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_prod_form.cpp                                              *
*                                                                             *
*  Created    : January 10th 2005                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2005-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lu_getfil.h>
#include <lu_html.h>
#include <lu_table.h>
#include <lu_prod_form.h>
#include <lu_form_valid.h>
#include <lu_param_list.h>
#include <lu_usermod.h>
using std::string;
using std::vector;
using std::ostream;
using std::endl;

namespace {
void printInstructions ( ostream& os, bool massesAllowed )
{
	os << "Enter&nbsp;Sequence&nbsp;in&nbsp;<b>Capital&nbsp;letters</b>&nbsp;(B,&nbsp;J,&nbsp;O,&nbsp;X,&nbsp;Z&nbsp;not&nbsp;allowed)&nbsp;except:" << endl;
	os << "<br />" << endl;
	os << "|&nbsp;m&nbsp;-&nbsp;Met-ox&nbsp;|&nbsp;h&nbsp;-&nbsp;Homoserine&nbsp;lactone&nbsp;|&nbsp;U&nbsp;-&nbsp;Selenocysteine&nbsp;|" << endl;
	os << "<br />" << endl;
	os << "|&nbsp;s,&nbsp;t,&nbsp;y&nbsp;-&nbsp;Phosphorylated&nbsp;S,&nbsp;T,&nbsp;Y&nbsp;|&nbsp;u,&nbsp;v,&nbsp;w,&nbsp;x&nbsp;-&nbsp;user&nbsp;specified&nbsp;amino&nbsp;acids&nbsp;|>" << endl;
	os << "<br />" << endl;
	os << "Modified&nbsp;amino&nbsp;acids&nbsp;may&nbsp;be&nbsp;entered&nbsp;using&nbsp;PSI&nbsp;notation&nbsp;-&nbsp;eg.&nbsp;M(Oxidation),&nbsp;S(Phospho)" << endl;
	os << "<br />" << endl;
	ExpandableJavascriptBlock ejb ( "Click + to see list of available PSI modifications (enter exactly as shown)" );
	ejb.printHeader ( os );
	StringVector sv = Usermod::getMSProdNames ();
	for ( StringSizeType i = 0 ; i < sv.size () ; i++ ) {
		os << sv [i] << "<br />" << endl;
	}
	ejb.printFooter ( os );
	os << "<br />" << endl;
	if ( massesAllowed ) {
		os << "An&nbsp;amino&nbsp;acid&nbsp;can&nbsp;also&nbsp;be&nbsp;followed&nbsp;by&nbsp;an&nbsp;exact&nbsp;mass&nbsp;-&nbsp;eg.&nbsp;P(-27.9949)&nbsp;or&nbsp;N(0.9840)" << endl;
		os << "<br />" << endl;
	}
	os << "Modifications&nbsp;may&nbsp;be&nbsp;added&nbsp;together&nbsp;-&nbsp;eg.&nbsp;X(MappingN+Deamidated)" << endl;
	os << "<br />" << endl;
	os << "Modified&nbsp;N&nbsp;and&nbsp;C&nbsp;termini&nbsp;must&nbsp;be&nbsp;selected&nbsp;from&nbsp;the&nbsp;menus" << endl;
	os << "<br />" << endl;
	os << "<br />" << endl;
}
}

class FormItemSequence : public FormItemText {
public:
	FormItemSequence ( int num = 1 ) :
	  FormItemText ( "", "", getName ( num ), 60, 1000, num == 1 ? "SAMPLER" : "" ) {}
	static string getName ( int num = 1 )
	{
		string n = ( num == 1 ) ? "" : gen_itoa ( num );
		return "sequence" + n;
	}
};

class FormItemElementalComposition : public FormItemText {
public:
	FormItemElementalComposition ( FormValidatingJavascript* fvj, int num ) :
		FormItemText ( "", "", getName ( num ), 26, 60, num == 1 ? "C60 H86 N13 O13 S2" : "", fvj->addElementalFormulaAllowBlankValidator (getName ( num ), "Elemental Composition") ) {}
	static string getName ( int num )
	{
		std::string n = ( num == 1 ) ? "" : gen_itoa ( num );
		return "elemental_composition" + n;
	}
};

class FormItemAveragineMass : public FormItemText {
public:
	FormItemAveragineMass ( FormValidatingJavascript* fvj, int num ) :
		FormItemText ( "", "", getName ( num ), 10, 12, num == 1 ? "1000.0" : "", fvj->addPositiveFloatingPointAllowBlankValidator (getName ( num ), "Averagine Mass") ) {}
	static string getName ( int num )
	{
		std::string n = ( num == 1 ) ? "" : gen_itoa ( num );
		return "averagine_mass" + n;
	}
};

class FormItemDistributionType : public FormItemSelect {
	static const char* options [];
public:
	FormItemDistributionType ( int num ) :
		FormItemSelect ( "", "", getName ( num ), options, num == 1 ? "Peptide Sequence" : "Off", 1, "showDistributionTypeItems( this.form." + getName ( num ) + ", " + gen_itoa ( num ) + ")" ) {}
	static string getName ( int num )
	{
		string n = ( num == 1 ) ? "" : gen_itoa ( num );
		return "distribution_type" + n;
	}
};
const char* FormItemDistributionType::options [] = { "Peptide Sequence", "Elemental Composition", "Averagine", "Off", 0 };

class FormItemProfileType : public FormItemSelect {
	static const char* options [];
public:
	FormItemProfileType () :
		FormItemSelect ( "Profile Type", "isoman.htm#profile_type", getName (), options, "Gaussian" ) {}
	static string getName () { return "profile_type"; }
};
const char* FormItemProfileType::options [] = { "Stick", "Gaussian", "Lorentzian", 0 };

class FormItemParentCharge : public FormItemText {
public:
	FormItemParentCharge ( FormValidatingJavascript* fvj, int num = 1, bool showLabel = true ) :
		FormItemText ( showLabel ? "Charge" : "", "", getName ( num ), 2, 2, "1", fvj->addPositiveNonZeroIntegerValidator (getName ( num ), "Charge") ) {}
	static string getName ( int num = 1 )
	{
		string n = ( num == 1 ) ? "" : gen_itoa ( num );
		return "parent_charge" + n;
	}
};

class FormItemIntensity : public FormItemText {
public:
	FormItemIntensity ( FormValidatingJavascript* fvj, int num ) :
		FormItemText ( "", "", getName ( num ), 6, 10, "1", fvj->addPositiveFloatingPointValidator (getName ( num ), "Intensity") ) {}
	static string getName ( int num = 1 )
	{
		string n = ( num == 1 ) ? "" : gen_itoa ( num );
		return "parent_intensity" + n;
	}
};

class FormItemResolution : public FormItemText {
public:
	FormItemResolution ( FormValidatingJavascript* fvj ) :
		FormItemText ( "Resolution", "isoman.htm#resolution", getName (), 8, 20, "10000.0", fvj->addPositiveFloatingPointValidator (getName (), "Resolution") ) {}
	static string getName () { return "resolution"; }
};

class FormItemIsoDetailedReport : public FormItemCheckbox {
public:
	FormItemIsoDetailedReport () :
		FormItemCheckbox ( "Detailed Report", "isoman.htm#detailed_report", getName (), false ) {}
	static string getName () { return "detailed_report"; }
};

MSIsotopeForm::MSIsotopeForm ( const VectorConstParameterListPtr& params ) :
	userAAForm ( params, &fvj ),
	saveResultsForm ( params, &fvj, false, true ),
	numSequences ( InfoParams::instance ().getIntValue ( "max_msprod_sequences", 2 ) )
{
	create ( params );
}
void MSIsotopeForm::printJavascriptFunctions ( ostream& os )
{
	basicJavascriptFunctions ( os );
	startJavascript ( os );
	os << "function showDistributionTypeItems ( item, num ) {" << endl;
	os << "\t" << "var val = getSelectValue(item);" << endl;
	os << "\t" << "if ( val == 'Peptide Sequence' ) {" << endl;
	os << "\t\t" << "showdiv ( 'ps' + num );" << endl;
	os << "\t\t" << "hidediv ( 'ec' + num );" << endl;
	os << "\t\t" << "hidediv ( 'am' + num );" << endl;
	os << "\t" << "}" << endl;
	os << "\t" << "else if ( val == 'Elemental Composition' ) {" << endl;
	os << "\t\t" << "hidediv ( 'ps' + num );" << endl;
	os << "\t\t" << "showdiv ( 'ec' + num );" << endl;
	os << "\t\t" << "hidediv ( 'am' + num );" << endl;
	os << "\t" << "}" << endl;
	os << "\t" << "else if ( val == 'Averagine' ) {" << endl;
	os << "\t\t" << "hidediv ( 'ps' + num );" << endl;
	os << "\t\t" << "hidediv ( 'ec' + num );" << endl;
	os << "\t\t" << "showdiv ( 'am' + num );" << endl;
	os << "\t" << "}" << endl;
	os << "\t" << "else {" << endl;				// Off
	os << "\t\t" << "hidediv ( 'ps' + num );" << endl;
	os << "\t\t" << "hidediv ( 'ec' + num );" << endl;
	os << "\t\t" << "hidediv ( 'am' + num );" << endl;
	os << "\t" << "}" << endl;
	os << "}" << endl;
	endJavascript ( os );
}
void MSIsotopeForm::createItems ()
{
	formItemMap [FormItemSearchName::getName ()] = new FormItemSearchName ( "msisotope" );
	formItemMap [FormItemReportTitle::getName ()]	= new FormItemReportTitle ( "MS-Isotope" );
	formItemMap [FormItemVersion::getName ()]	= new FormItemVersion ();

	for ( int i = 1 ; i <= numSequences ; i++ ) {
		formItemMap [FormItemDistributionType::getName (i)]		= new FormItemDistributionType ( i );
		formItemMap [FormItemParentCharge::getName (i)]			= new FormItemParentCharge ( &fvj, i, false );
		formItemMap [FormItemIntensity::getName (i)]			= new FormItemIntensity ( &fvj, i );
		formItemMap [FormItemNTerm::getName (i)]				= new FormItemNTerm ( i, false );
		formItemMap [FormItemSequence::getName (i)]				= new FormItemSequence ( i );
		formItemMap [FormItemCTerm::getName (i)]				= new FormItemCTerm ( i, false );
		formItemMap [FormItemElementalComposition::getName (i)]	= new FormItemElementalComposition ( &fvj, i );
		formItemMap [FormItemAveragineMass::getName (i)]		= new FormItemAveragineMass ( &fvj, i );
	}
	formItemMap [FormItemProfileType::getName ()]		= new FormItemProfileType ();
	formItemMap [FormItemResolution::getName ()]		= new FormItemResolution ( &fvj );
	formItemMap [FormItemDisplayGraph::getName ()]		= new FormItemDisplayGraph ( true );
	formItemMap [FormItemIsoDetailedReport::getName ()]	= new FormItemIsoDetailedReport ();

	formItemMap [FormItemIsotopePurity::getName ( "C", 13 )] = new FormItemIsotopePurity ( "C", 13, &fvj );
	formItemMap [FormItemIsotopePurity::getName ( "N", 15 )] = new FormItemIsotopePurity ( "N", 15, &fvj );
	formItemMap [FormItemIsotopePurity::getName ( "O", 18 )] = new FormItemIsotopePurity ( "O", 18, &fvj );
}
void MSIsotopeForm::printHTML ( ostream& os )
{
	init_html ( os, "MS-Isotope", "isoman.htm" );
	fvj.print ( os );
	printHTMLFORMStart ( os, "post", "mssearch", false, true );
	int numCols = 4;
	printJavascriptFunctions ( os );
	tableStart ( os, true );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true, 2 );
				tableStart ( os, true );
					tableRowStart ( os );
						tableHeaderStart ( os, "", "center", false, numCols );
							os << "<div class=\"form_section_header\">Instructions For Entering Peptide Sequences</div>" << endl;
							printInstructions ( os, false );
						tableHeaderEnd ( os );
					tableRowEnd ( os );
					ExpandableJavascriptBlock* ejb = 0;
					for ( int i = 1 ; i <= numSequences ; i++ ) { 
						if ( i == 2 ) {
							tableRowStart ( os );
								tableHeaderStart ( os, "", "left", false, numCols );
									ejb = new ExpandableJavascriptBlock ( "Additional Distributions" );
									ejb->printHeader ( os );
									tableStart ( os, true );
						}
						if ( i == 1 || i == 2 ) {
							tableRowStart ( os );
								tableHeaderStart ( os );
									printHTMLFORMLabel ( os, "Charge", "isoman.htm#charge" );
								tableHeaderEnd ( os );
								tableHeaderStart ( os );
									printHTMLFORMLabel ( os, "Intensity", "isoman.htm#intensity" );
								tableHeaderEnd ( os );
								tableHeaderStart ( os, "", "center" );
									printHTMLFORMLabel ( os, "Distribution Type", "isoman.htm#distribution_type" );
								tableHeaderEnd ( os );
								tableHeaderStart ( os, "", "center" );
									printHTMLFORMLabel ( os, "Value", "" );
								tableHeaderEnd ( os );
							tableRowEnd ( os );
						}
						tableRowStart ( os );
							tableHeaderStart ( os );
								formItemMap [FormItemParentCharge::getName (i)]->printHTML ( os );
							tableHeaderEnd ( os );
							tableHeaderStart ( os );
								formItemMap [FormItemIntensity::getName (i)]->printHTML ( os );
							tableHeaderEnd ( os );
							tableHeaderStart ( os, "", "center", true );
								formItemMap [FormItemDistributionType::getName (i)]->printHTML ( os );
							tableHeaderEnd ( os );
							tableHeaderStart ( os, "", "left", true );
								if ( i == 1 )	os << "<div id=\"ps" << i << "\" style=\"display:block\">" << endl;
								else			os << "<div id=\"ps" << i << "\" style=\"display:none\">" << endl;
									formItemMap [FormItemNTerm::getName (i)]->printHTML ( os );
									formItemMap [FormItemSequence::getName (i)]->printHTML ( os );
									formItemMap [FormItemCTerm::getName (i)]->printHTML ( os );
								os << "</div>" << endl;
								os << "<div id=\"ec" << i << "\" style=\"display:none\">" << endl;
									formItemMap [FormItemElementalComposition::getName (i)]->printHTML ( os );
								os << "</div>" << endl;
								os << "<div id=\"am" << i << "\" style=\"display:none\">" << endl;
									formItemMap [FormItemAveragineMass::getName (i)]->printHTML ( os );
								os << "</div>" << endl;
							tableHeaderEnd ( os );
						tableRowEnd ( os );
						if ( ejb && i == numSequences ) {
									tableEnd ( os );
									ejb->printFooter ( os );
									delete ejb;
								tableHeaderEnd ( os );
							tableRowEnd ( os );
						}
					}
				tableEnd ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center" );
				userAAForm.printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );

		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true, 2 );
				tableStart ( os );
					tableRowStart ( os );
						tableHeaderStart ( os, "", "center", true );
							printHTMLFORMSubmit ( os, "Calculate Isotope Distribution" );
						tableHeaderEnd ( os );
					tableRowEnd ( os );
				tableEnd ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );

		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true, 2 );
				formItemMap [FormItemProfileType::getName ()]->printHTML ( os );
				formItemMap [FormItemResolution::getName ()]->printHTML ( os );
				formItemMap [FormItemDisplayGraph::getName ()]->printHTML ( os );
				formItemMap [FormItemIsoDetailedReport::getName ()]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );

		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true, 2 );
				saveResultsForm.printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );

		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true, 2 );
				formItemMap [FormItemIsotopePurity::getName ( "C", 13 )]->printHTML ( os );
				formItemMap [FormItemIsotopePurity::getName ( "N", 15 )]->printHTML ( os );
				formItemMap [FormItemIsotopePurity::getName ( "O", 18 )]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
	tableEnd ( os );
	showHiddenItems ( os );
	printHTMLFORMEnd ( os );
	printHTMLFooter ( os, "1998" );
}
void MSIsotopeForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];

		for ( int i = 1 ; i <= numSequences ; i++ ) { 
			formItemMap [FormItemDistributionType::getName (i)]->setValue ( p );
			formItemMap [FormItemParentCharge::getName (i)]->setValue ( p );
			formItemMap [FormItemIntensity::getName (i)]->setValue ( p );
			formItemMap [FormItemNTerm::getName (i)]->setValue ( p );
			formItemMap [FormItemSequence::getName (i)]->setValue ( p );
			formItemMap [FormItemCTerm::getName (i)]->setValue ( p );
			formItemMap [FormItemElementalComposition::getName (i)]->setValue ( p );
			formItemMap [FormItemAveragineMass::getName (i)]->setValue ( p );
		}
		formItemMap [FormItemProfileType::getName ()]->setValue ( p );
		formItemMap [FormItemDisplayGraph::getName ()]->setValue ( p );
		formItemMap [FormItemResolution::getName ()]->setValue ( p );
		formItemMap [FormItemIsoDetailedReport::getName ()]->setValue ( p );
		formItemMap [FormItemIsotopePurity::getName ( "C", 13 )]->setValue ( p );
		formItemMap [FormItemIsotopePurity::getName ( "N", 15 )]->setValue ( p );
		formItemMap [FormItemIsotopePurity::getName ( "O", 18 )]->setValue ( p );
	}
}

const char* FormItemMaxCharge::options [] = {
	"No Limit", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
	"13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", 0 };

const char* FormItemCountPosZ::options [] = { "Count Basic AA", "Ignore Basic AA", 0 };

FormItemMaxInternalLen::FormItemMaxInternalLen ( FormValidatingJavascript* fvj ) :
	FormItemText ( "Max Internal Len", "prodman.htm#Max_internal_len", getName (), 4, 4, "200", fvj->addPositiveNonZeroIntegerValidator ( getName (), "Max Internal Len" ) )
{
}

MSProductForm::MSProductForm ( const VectorConstParameterListPtr& params, ProspectorForm* dataForm ) :
	ionTypeForm ( params ),
	userAAForm ( params, &fvj ),
	saveResultsForm ( params, &fvj ),
	dataForm ( dataForm ),
	numSequences ( InfoParams::instance ().getIntValue ( "max_msprod_sequences", 2 ) )
{
	create ( params );
}
void MSProductForm::createItems ()
{
	formItemMap [FormItemSearchName::getName ()]	= new FormItemSearchName ( "msproduct" );
	formItemMap [FormItemReportTitle::getName ()]	= new FormItemReportTitle ( "MS-Product" );
	formItemMap [FormItemVersion::getName ()]	= new FormItemVersion ();

	for ( int i = 1 ; i <= numSequences ; i++ ) { 
		formItemMap [FormItemS::getName (i)]		= new FormItemS ( i );
		formItemMap [FormItemNTerm::getName (i)]	= new FormItemNTerm ( i, false );
		formItemMap [FormItemSequence::getName (i)]	= new FormItemSequence ( i );
		formItemMap [FormItemCTerm::getName (i)]	= new FormItemCTerm ( i, false );
	}
	formItemMap [FormItemAlternative::getName ()]	= new FormItemAlternative ();
	formItemMap [FormItemDiscriminating::getName ()]= new FormItemDiscriminating ();
	formItemMap [FormItemComposition::getName("bridge")]= new FormItemComposition ( "Bridge", "bridge", 1, "", &fvj );

	formItemMap [FormItemUseInstrumentIonTypes::getName ()]	= new FormItemUseInstrumentIonTypes ( false );

	formItemMap [FormItemMaxCharge::getName ()]	= new FormItemMaxCharge ();
	formItemMap [FormItemCountPosZ::getName ()] = new FormItemCountPosZ ();
	formItemMap [FormItemMaxLosses::getName ()] = new FormItemMaxLosses ();
	formItemMap [FormItemMaxInternalLen::getName ()] = new FormItemMaxInternalLen (&fvj);

	formItemMap [FormItemDisplayGraph::getName()] = new FormItemDisplayGraph ( true );

	formItemMap [FormItemParentMassConvert::getName ()]			= new FormItemParentMassConvert ( "MS" );
	formItemMap [FormItemFragmentMassesTolerance::getName ()]	= new FormItemFragmentMassesTolerance (&fvj);
	formItemMap [FormItemFragmentMassesToleranceUnits::getName ()]= new FormItemFragmentMassesToleranceUnits ();

	formItemMap [FormItemMSMSPkFilter::getName ()]			= new FormItemMSMSPkFilter ( true );
	formItemMap [FormItemMaxMSMSPeaks::getName ()]			= new FormItemMaxMSMSPeaks ( &fvj, true );
}
void MSProductForm::printHTML ( ostream& os )
{
	init_html ( os, "MS-Product", "prodman.htm" );
	fvj.print ( os );
	printHTMLFORMStart ( os, "post", "mssearch", false, true );
	int numCols = 4;
	tableStart ( os, true );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				tableStart ( os, true );
					tableRowStart ( os );
						tableHeaderStart ( os, "", "center", false, numCols );
							os << "<div class=\"form_section_header\">Peptide Sequence</div>" << endl;
						tableHeaderEnd ( os );
					tableRowEnd ( os );
					tableRowStart ( os );
						tableHeaderStart ( os, "", "center", false, numCols );
							printInstructions ( os, true );
						tableHeaderEnd ( os );
					tableRowEnd ( os );
					ExpandableJavascriptBlock* ejb = 0;
					for ( int i = 1 ; i <= numSequences ; i++ ) { 
						if ( i == 2 ) {
							tableRowStart ( os );
								tableHeaderStart ( os, "", "left", false, numCols );
									ejb = new ExpandableJavascriptBlock ( "Additional Sequences" );
									ejb->printHeader ( os );
									tableStart ( os, true );
						}
						if ( i == 1 ) {
							tableRowStart ( os );
								tableHeader ( os, "" );
								tableHeaderStart ( os );
									printHTMLFORMLabel ( os, "N Term", "allman.htm#term_groups" );
								tableHeaderEnd ( os );
								tableHeader ( os, "Sequence" );
								tableHeaderStart ( os );
									printHTMLFORMLabel ( os, "C Term", "allman.htm#term_groups" );
								tableHeaderEnd ( os );
							tableRowEnd ( os );
						}
						tableRowStart ( os );
							tableHeaderStart ( os );
								formItemMap [FormItemS::getName (i)]->printHTML ( os );
							tableHeaderEnd ( os );
							tableHeaderStart ( os );
								formItemMap [FormItemNTerm::getName (i)]->printHTML ( os );
							tableHeaderEnd ( os );
							tableHeaderStart ( os );
								formItemMap [FormItemSequence::getName (i)]->printHTML ( os );
							tableHeaderEnd ( os );
							tableHeaderStart ( os );
								formItemMap [FormItemCTerm::getName (i)]->printHTML ( os );
							tableHeaderEnd ( os );
						tableRowEnd ( os );
						if ( ejb && i == numSequences ) {
										tableRowStart ( os );
											tableHeaderStart ( os, "", "left", false, numCols );
												formItemMap [FormItemAlternative::getName ()]->printHTML ( os );
												formItemMap [FormItemDiscriminating::getName ()]->printHTML ( os );
												formItemMap [FormItemComposition::getName("bridge")]->printHTML ( os );
											tableHeaderEnd ( os );
										tableRowEnd ( os );
									tableEnd ( os );
									ejb->printFooter ( os );
									delete ejb;
								tableHeaderEnd ( os );
							tableRowEnd ( os );
						}
					}
					tableRowStart ( os );
						tableHeaderStart ( os, "", "center", false, numCols );
							userAAForm.printHTML ( os );
						tableHeaderEnd ( os );
					tableRowEnd ( os );
				tableEnd ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", false, numCols );
				formItemMap [FormItemUseInstrumentIonTypes::getName ()]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", false, numCols );
				ionTypeForm.printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				printHTMLFORMSubmit ( os, "Induce Fragmentation" );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				formItemMap [FormItemDisplayGraph::getName()]->printHTML ( os );
				formItemMap [FormItemMaxCharge::getName ()]->printHTML ( os );
				formItemMap [FormItemCountPosZ::getName ()]->printHTML ( os );
				formItemMap [FormItemMaxLosses::getName ()]->printHTML ( os );
				formItemMap [FormItemMaxInternalLen::getName ()]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				saveResultsForm.printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				formItemMap [FormItemParentMassConvert::getName ()]->printHTML ( os );
				formItemMap [FormItemFragmentMassesTolerance::getName ()]->printHTML ( os );
				formItemMap [FormItemFragmentMassesToleranceUnits::getName ()]->printHTML ( os );
				formItemMap [FormItemMSMSPkFilter::getName ()]->printHTML ( os );
				formItemMap [FormItemMaxMSMSPeaks::getName ()]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		dataForm->printHTML ( os );
	tableEnd ( os );
	showHiddenItems ( os );
	printHTMLFORMEnd ( os );
	printHTMLFooter ( os, "1995" );
}
void MSProductForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];
		for ( int i = 1 ; i <= numSequences ; i++ ) { 
			formItemMap [FormItemNTerm::getName (i)]->setValue ( p );
			formItemMap [FormItemSequence::getName (i)]->setValue ( p );
			formItemMap [FormItemCTerm::getName (i)]->setValue ( p );
		}
		formItemMap [FormItemAlternative::getName ()]->setValue ( p );
		formItemMap [FormItemDiscriminating::getName ()]->setValue ( p );
		formItemMap [FormItemComposition::getName("bridge")]->setValue ( p );
		formItemMap [FormItemUseInstrumentIonTypes::getName ()]->setValue ( p );
		//ionTypeForm.printHTML ( os );
		formItemMap [FormItemDisplayGraph::getName()]->setValue ( p );
		formItemMap [FormItemMaxCharge::getName ()]->setValue ( p );
		formItemMap [FormItemCountPosZ::getName ()]->setValue ( p );
		formItemMap [FormItemMaxInternalLen::getName ()]->setValue ( p );
		formItemMap [FormItemMaxLosses::getName ()]->setValue ( p );
		formItemMap [FormItemParentMassConvert::getName ()]->setValue ( p );
		formItemMap [FormItemFragmentMassesTolerance::getName ()]->setValue ( p );
		formItemMap [FormItemFragmentMassesToleranceUnits::getName ()]->setValue ( p );
		formItemMap [FormItemMSMSPkFilter::getName ()]->setValue ( p );
		formItemMap [FormItemMaxMSMSPeaks::getName ()]->setValue ( p );
	}
}
MSProductPasteForm::MSProductPasteForm ( const VectorConstParameterListPtr& params ) :
	MSProductForm ( params, new PasteAreaDataForm ( params, "msproduct" ) )
{
	create ( params );
}
MSProductUploadForm::MSProductUploadForm ( const VectorConstParameterListPtr& params ) :
	MSProductForm ( params, new UploadDataForm ( params ) )
{
	create ( params );
}
IonTypeForm::IonTypeForm ( const VectorConstParameterListPtr& params )
{
	create ( params );
}
void IonTypeForm::createItems ()
{
	formItemMap [FormItemIonType::getName ()+string("i")]	= new FormItemIonType ( "", "i", true );
	formItemMap [FormItemIonType::getName ()+string("m")]	= new FormItemIonType ( "", "m", false );

	formItemMap [FormItemIonType::getName ()+string("a")]	= new FormItemIonType ( "", "a", true );
	formItemMap [FormItemIonType::getName ()+string("b")]	= new FormItemIonType ( "", "b", true );
	formItemMap [FormItemIonType::getName ()+string("c_plus_2")] = new FormItemIonType ( "", "c+2", false );
	formItemMap [FormItemIonType::getName ()+string("c_plus_1")] = new FormItemIonType ( "", "c+1", false );
	formItemMap [FormItemIonType::getName ()+string("c")]	= new FormItemIonType ( "", "c", false );
	formItemMap [FormItemIonType::getName ()+string("c_minus_1")] = new FormItemIonType ( "", "c-1", false );

	formItemMap [FormItemIonType::getName ()+string("x")]	= new FormItemIonType ( "", "x", false );
	formItemMap [FormItemIonType::getName ()+string("y")]	= new FormItemIonType ( "", "y", true );
	formItemMap [FormItemIonType::getName ()+string("Y")]	= new FormItemIonType ( "", "Y", false );
	formItemMap [FormItemIonType::getName ()+string("z")]	= new FormItemIonType ( "", "z", false );
	formItemMap [FormItemIonType::getName ()+string("z_plus_1")] = new FormItemIonType ( "", "z+1", false );
	formItemMap [FormItemIonType::getName ()+string("z_plus_2")] = new FormItemIonType ( "", "z+2", false );
	formItemMap [FormItemIonType::getName ()+string("z_plus_3")] = new FormItemIonType ( "", "z+3", false );

	formItemMap [FormItemIonType::getName ()+string("I")]	= new FormItemIonType ( "", "I", false );

	formItemMap [FormItemIonType::getName ()+string("N")]	= new FormItemIonType ( "", "N", false );
	formItemMap [FormItemIonType::getName ()+string("C")]	= new FormItemIonType ( "", "C", false );

	formItemMap [FormItemIonType::getName ()+string("d")]	= new FormItemIonType ( "", "d", false );
	formItemMap [FormItemIonType::getName ()+string("v")]	= new FormItemIonType ( "", "v", false );
	formItemMap [FormItemIonType::getName ()+string("w")]	= new FormItemIonType ( "", "w", false );
	formItemMap [FormItemIonType::getName ()+string("h")]	= new FormItemIonType ( "", "h", true );
	formItemMap [FormItemIonType::getName ()+string("n")]	= new FormItemIonType ( "", "n", true );
	formItemMap [FormItemIonType::getName ()+string("P")]	= new FormItemIonType ( "", "P", false );
	formItemMap [FormItemIonType::getName ()+string("S")]	= new FormItemIonType ( "", "S", false );
	formItemMap [FormItemIonType::getName ()+string("B")]	= new FormItemIonType ( "", "B", true );
}
void IonTypeForm::printHTML ( ostream& os )
{
	tableStart ( os, true, "", "1" );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true, 2 );
				os << "AA<br />Composition" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 6 );
				os << "N-term<br />Sequence" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 7 );
				os << "C-term<br />Sequence" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				os << "Internal<br />Fragment" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 2 );
				os << "Ladder<br />Sequencing" << endl;
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("i")]->printHTML ( os );
				os << "<br />i" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("m")]->printHTML ( os );
				os << "<br />m" << endl;
			tableHeaderEnd ( os );

			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("a")]->printHTML ( os );
				os << "<br />a" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("b")]->printHTML ( os );
				os << "<br />b" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("c_plus_2")]->printHTML ( os );
				os << "<br />c+2" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("c_plus_1")]->printHTML ( os );
				os << "<br />c+1" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("c")]->printHTML ( os );
				os << "<br />c" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("c_minus_1")]->printHTML ( os );
				os << "<br />c-1" << endl;
			tableHeaderEnd ( os );

			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("x")]->printHTML ( os );
				os << "<br />x" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("y")]->printHTML ( os );
				os << "<br />y" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("Y")]->printHTML ( os );
				os << "<br />Y" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("z")]->printHTML ( os );
				os << "<br />z" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("z_plus_1")]->printHTML ( os );
				os << "<br />z+1" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("z_plus_2")]->printHTML ( os );
				os << "<br />z+2" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("z_plus_3")]->printHTML ( os );
				os << "<br />z+3" << endl;
			tableHeaderEnd ( os );

			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("I")]->printHTML ( os );
				os << "<br />Internal" << endl;
			tableHeaderEnd ( os );

			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("N")]->printHTML ( os );
				os << "<br />N-term" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("C")]->printHTML ( os );
				os << "<br />C-term" << endl;
			tableHeaderEnd ( os );
		tableRowEnd ( os );
	tableEnd ( os );
	tableStart ( os, true, "", "1" );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true, 3 );
				os << "Satellite Sequence<br />(side-chain loss)" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 4 );
				os << "Neutral-loss<br />Sequence" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				os << "Peeling<br />Sequence" << endl;
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("d")]->printHTML ( os );
				os << "<br />d" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("v")]->printHTML ( os );
				os << "<br />v" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("w")]->printHTML ( os );
				os << "<br />w" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("h")]->printHTML ( os );
				os << "<br />-H<sub>2</sub>O" << endl;
				os << "<br />S, T, E, D" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("n")]->printHTML ( os );
				os << "<br />-NH<sub>3</sub>" << endl;
				os << "<br />R, K, Q, N" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("P")]->printHTML ( os );
				os << "<br />-H<sub>3</sub>PO<sub>4</sub>" << endl;
				os << "<br />(S, T, Y + PO<sub>4</sub>)" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("S")]->printHTML ( os );
				os << "<br />-SOCH<sub>4</sub>" << endl;
				os << "<br />(M + Ox)" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("B")]->printHTML ( os );
				os << "<br />b+H<sub>2</sub>O" << endl;
				os << "<br />R, H, K" << endl;
			tableHeaderEnd ( os );
		tableRowEnd ( os );
	tableEnd ( os );
}
void IonTypeForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];

		formItemMap [FormItemIonType::getName()+string("i")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("m")]->setValue ( p );

		formItemMap [FormItemIonType::getName()+string("a")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("b")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("c_plus_2")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("c_plus_1")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("c")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("c_minus_1")]->setValue ( p );

		formItemMap [FormItemIonType::getName()+string("x")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("y")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("Y")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("z")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("z_plus_1")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("z_plus_2")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("z_plus_3")]->setValue ( p );

		formItemMap [FormItemIonType::getName()+string("I")]->setValue ( p );

		formItemMap [FormItemIonType::getName()+string("N")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("C")]->setValue ( p );

		formItemMap [FormItemIonType::getName()+string("d")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("v")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("w")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("h")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("n")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("P")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("S")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("B")]->setValue ( p );
	}
}
void IonTypeForm::removeName ( ParameterList* p )
{
	p->removeName ( FormItemIonType::getName () );
}
