/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_filter_form.cpp                                            *
*                                                                             *
*  Created    : September 3rd 2012                                            *
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
#include <lu_data_form.h>
#include <lu_html.h>
#include <lu_table.h>
#include <lu_filter_form.h>

using std::ostream;
using std::endl;
using std::string;

class FormItemAllCharges : public FormItemCheckbox {
public:
	FormItemAllCharges () :
		FormItemCheckbox ( "All", "", getName (), true ) {}
	static string getName () { return "all_charges"; }
};

class FormItemChargeFilter : public FormItemSelectMultiple {
	static const char* options [];
public:
	FormItemChargeFilter () :
		FormItemSelectMultiple ( "Charge Filter", "filterman.htm#precursor_charge", getName (), options, StringVector (), 4 ) {}
	static string getName () { return "charge_filter"; }
};
const char* FormItemChargeFilter::options [] = {
	"1", "2", "3", "4", "5", "6", "7", "8", "9", "10 and above", 0 };

class FormItemLossComposition : public FormItemText {
public:
	FormItemLossComposition ( FormValidatingJavascript* fvj ) :
		FormItemText ( "Loss Composition", "filterman.htm#loss_composition", getName (), 16, 60, "", fvj->addElementalFormulaAllowBlankValidator (getName (), "Loss Composition") ) {}
	static string getName () { return "loss_composition"; }
};

class FormItemFragmentMZ : public FormItemText {
public:
	FormItemFragmentMZ ( FormValidatingJavascript* fvj ) :
		FormItemText ( "Fragment M/Z", "filterman.htm#fragment_mz", getName (), 6, 10, "", fvj->addPositiveFloatingPointValidator (getName (), "Fragment M/Z") ) {}
	static string getName () { return "fragment_mz"; }
};

class FormItemFragmentMZs : public FormItemTextArea {
public:
	FormItemFragmentMZs () :
		FormItemTextArea ( "Fragment M/Zs", "filterman.htm#fragment_mz", getName (), 6, 10, StringVector () ) {}
	static string getName () { return "fragment_mzs"; }
};

class FormItemMinFragMatches : public FormItemText {
public:
	FormItemMinFragMatches ( FormValidatingJavascript* fvj ) :
		FormItemText ( "Min Matches", "filterman.htm#min_matches", getName (), 2, 2, "0", fvj->addPositiveIntegerValidator ( getName (), "Min Matches" ) ) {}
	static std::string getName () { return "min_frag_matches"; }
};

class FormItemMSMSPrecursorExclusion : public FormItemText {
public:
	FormItemMSMSPrecursorExclusion ( FormValidatingJavascript* fvj, bool allowBlank );
	static string getName () { return "msms_precursor_exclusion"; }
};

FormItemMSMSPrecursorExclusion::FormItemMSMSPrecursorExclusion ( FormValidatingJavascript* fvj, bool allowBlank ) :
	FormItemText ( "", "", getName (), 4, 6, "0", allowBlank ? fvj->addPositiveNonZeroIntegerAllowBlankValidator ( getName (), "MSMS Precursor Exclusion" ) : fvj->addPositiveNonZeroIntegerValidator ( getName (), "MSMS Precursor Exclusion" ) )
{
}

class FormItemKeepOrRemove : public FormItemSelect {
	static const char* options [];
public:
	FormItemKeepOrRemove ();
	static string getName () { return "keep_or_remove"; }
};

const char* FormItemKeepOrRemove::options [] = {
	"Keep Spectra Matching Criteria",
	"Remove Spectra Matching Criteria",
	"Create Both Sets of Spectra", 0
};

FormItemKeepOrRemove::FormItemKeepOrRemove () :
	FormItemSelect ( "", "", getName (), options, options [0] )
{
}

MSFilterForm::MSFilterForm ( const VectorConstParameterListPtr& params ) :
	fvj ( fvj ),
	mPlusHForm ( &fvj, params )
{
	create ( params );
}
void MSFilterForm::createItems ()
{
	formItemMap [FormItemSearchName::getName ()]	= new FormItemSearchName ( "msfilter" );
	formItemMap [FormItemReportTitle::getName ()]	= new FormItemReportTitle ( "MS-Filter" );
	formItemMap [FormItemVersion::getName ()]		= new FormItemVersion ();

	formItemMap [FormItemUploadPeakList::getName ()]= new FormItemUploadPeakList ();

	formItemMap [FormItemChargeFilter::getName ()]	= new FormItemChargeFilter ();
	formItemMap [FormItemAllCharges::getName ()]	= new FormItemAllCharges ();

	formItemMap [FormItemLossComposition::getName ()] = new FormItemLossComposition (&fvj);
	formItemMap [FormItemFragmentMZs::getName ()] = new FormItemFragmentMZs ();

	formItemMap [FormItemMinFragMatches::getName ()] = new FormItemMinFragMatches (&fvj);

	formItemMap [FormItemParentMassConvert::getName ()]			= new FormItemParentMassConvert ( "MS" );
	//formItemMap [FormItemMSMSParentMassTolerance::getName ()]		= new FormItemMSMSParentMassTolerance ( &fvj );
	formItemMap [FormItemMSMSParentMassToleranceUnits::getName ()]	= new FormItemMSMSParentMassToleranceUnits ();
	formItemMap [FormItemMassSystematicError::getName ("msms_")]	= new FormItemMassSystematicError (&fvj, "msms_");
	formItemMap [FormItemFragmentMassesTolerance::getName ()]	= new FormItemFragmentMassesTolerance (&fvj);
	formItemMap [FormItemFragmentMassesToleranceUnits::getName ()]= new FormItemFragmentMassesToleranceUnits ();
	formItemMap [FormItemInstrumentName::getName ()]= new FormItemInstrumentName ();

	formItemMap [FormItemMSMSPkFilter::getName ()]	= new FormItemMSMSPkFilter ( true );
	formItemMap [FormItemMaxMSMSPeaks::getName ()]	= new FormItemMaxMSMSPeaks ( &fvj, true );
	formItemMap [FormItemMSMSMinPrecursorMass::getName ()]		= new FormItemMSMSMinPrecursorMass ( &fvj, true );
	formItemMap [FormItemMSMSPrecursorExclusion::getName ()]	= new FormItemMSMSPrecursorExclusion ( &fvj, true );

	formItemMap [FormItemKeepOrRemove::getName ()]	= new FormItemKeepOrRemove ();
}
void MSFilterForm::printHTML ( ostream& os )
{
	init_html ( os, "MS-Filter", "filterman.htm" );
	printHTMLFORMStart ( os, "post", "mssearch", true, true );
	tableStart ( os, true );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				mPlusHForm.printHTML ( os );
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true );
				formItemMap [FormItemChargeFilter::getName ()]->printHTML ( os );
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true );
				formItemMap [FormItemAllCharges::getName ()]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true );
				formItemMap [FormItemFragmentMZs::getName ()]->printHTML ( os );
				os << "<br /><br />" << endl;
				formItemMap [FormItemMinFragMatches::getName ()]->printHTML ( os );
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "left", true, 2 );
				formItemMap [FormItemLossComposition::getName ()]->printHTML ( os );
				os << "<br /><br />" << endl;
				formItemMap [FormItemInstrumentName::getName ()]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true, 3 );
				formItemMap [FormItemParentMassConvert::getName ()]->printHTML ( os );
				//formItemMap [FormItemMSMSParentMassTolerance::getName ()]->printHTML ( os );
				formItemMap [FormItemFragmentMassesTolerance::getName ()]->printHTML ( os );
				formItemMap [FormItemFragmentMassesToleranceUnits::getName ()]->printHTML ( os );
				formItemMap [FormItemMassSystematicError::getName ("msms_")]->printHTML ( os );
				formItemMap [FormItemMSMSParentMassToleranceUnits::getName ()]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true, 3 );
				formItemMap [FormItemMSMSPkFilter::getName ()]->printHTML ( os );
				formItemMap [FormItemMaxMSMSPeaks::getName ()]->printHTML ( os );
				formItemMap [FormItemKeepOrRemove::getName ()]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true, 3 );
				printHTMLFORMSubmit ( os, "Upload File" );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true, 3 );
				os << "<div class=\"form_large_label\">Peak List File</div>" << endl;
				formItemMap [FormItemUploadPeakList::getName ()]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
	tableEnd ( os );
	showHiddenItems ( os );
	printHTMLFORMEnd ( os );
	printHTMLFooter ( os, "2012" );
}
void MSFilterForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];

		formItemMap [FormItemParentMassConvert::getName ()]->setValue ( p );
		//formItemMap [FormItemMSMSParentMassTolerance::getName ()]->setValue ( p );
		formItemMap [FormItemMSMSParentMassToleranceUnits::getName ()]->setValue ( p );
		formItemMap [FormItemFragmentMassesTolerance::getName ()]->setValue ( p );
		formItemMap [FormItemFragmentMassesToleranceUnits::getName ()]->setValue ( p );
		formItemMap [FormItemMassSystematicError::getName ("msms_")]->setValue ( p );

		formItemMap [FormItemChargeFilter::getName ()]->setValue ( p );
		formItemMap [FormItemAllCharges::getName ()]->setValue ( p );

		formItemMap [FormItemLossComposition::getName ()]->setValue ( p );
		formItemMap [FormItemFragmentMZs::getName ()]->setValue ( p );
		formItemMap [FormItemMinFragMatches::getName ()]->setValue ( p );
		formItemMap [FormItemMSMSMinPrecursorMass::getName ()]->setValue ( p );
		formItemMap [FormItemKeepOrRemove::getName ()]->setValue ( p );
	}
}
