/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_comp_form.cpp                                              *
*                                                                             *
*  Created    : January 11th 2005                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2005-2010) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lu_html.h>
#include <lu_table.h>
#include <lu_comp_form.h>
#include <lu_data_form.h>

using std::vector;
using std::ostream;
using std::string;
using std::endl;

class FormItemMSMSParentMass : public FormItemText {
public:
	FormItemMSMSParentMass () :
		FormItemText ( "m/z (Da)", "", getName (), 9, 12, "1260.5992" ) {}
	static string getName () { return "msms_parent_mass"; }
};

class FormItemMSMSParentCharge : public FormItemSelect {
	static const char* options [];
public:
	FormItemMSMSParentCharge () :
		FormItemSelect ( "Charge", "", getName (), options, "1" ) {}
	static string getName () { return "msms_parent_charge"; }
};
const char* FormItemMSMSParentCharge::options [] = { "6", "5", "4", "3", "2", "1", "0", "-1", "-2", "-3", "-4", "-5", "-6", 0 };

class FormItemIT : public FormItemSelectMultiple {
	static const char* options [];
public:
	FormItemIT () :
		FormItemSelectMultiple ( "Ion Types", "", getName (), options, StringVector (), 6 )
	{
		select.push_back ( "MH+" );
	}
	static string getName () { return "it"; }
};
const char* FormItemIT::options [] = {
	"MH+", "y", "b", "a", "c", "b+H2O",
	"y-NH3", "b-NH3", "a-NH3", "y-H2O", "b-H2O", "a-H2O",
	"y-H3PO4", "b-H3PO4", "a-H3PO4", "y-SOCH4", "b-SOCH4",
	0
};

class FormItemMSMSMaxReportedHits : public FormItemText {
public:
	FormItemMSMSMaxReportedHits ( FormValidatingJavascript* fvj ) :
		FormItemText ( "Maximum Reported Compositions", "", getName (), 6, 10, "1000", fvj->addPositiveNonZeroIntegerValidator ( getName (), "Maximum Reported Compositions" ) ) {}
	static string getName () { return "msms_max_reported_hits"; }
};

class FormItemCombinationType : public FormItemSelect {
	static const char* options [];
public:
	FormItemCombinationType () :
		FormItemSelect ( "Combination Type", "compman.htm#combination_types", getName (), options, "Amino Acid" ) {}
	static string getName () { return "combination_type"; }
};
const char* FormItemCombinationType::options [] = { "Amino Acid", "Peptide Elemental", "Elemental", "Any Elemental", 0 };

class FormItemCompositionSearch : public FormItemCheckbox {
public:
	FormItemCompositionSearch () :
		FormItemCheckbox ( "", "", getName (), true ) {}
	static string getName () { return "composition_search"; }
};

MSCompForm::MSCompForm ( const VectorConstParameterListPtr& params ) :
	userAAForm ( params, &fvj ),
	immoniumIonForm ( params ),
	extraIonForm ( params ),
	absentIonForm ( params ),
	modifiedIonForm ( params ),
	saveResultsForm ( params, &fvj )
{
	create ( params );
}
void MSCompForm::createItems ()
{
	formItemMap [FormItemSearchName::getName ()]	= new FormItemSearchName ( "mscomp" );
	formItemMap [FormItemReportTitle::getName ()]	= new FormItemReportTitle ( "MS-Comp" );
	formItemMap [FormItemVersion::getName ()]	= new FormItemVersion ();

	formItemMap [FormItemMSMSParentMass::getName ()]	= new FormItemMSMSParentMass ();
	formItemMap [FormItemMSMSParentCharge::getName ()]	= new FormItemMSMSParentCharge ();
	formItemMap [FormItemParentMassConvert::getName ()]	= new FormItemParentMassConvert ( "MS" );

	formItemMap [FormItemMSMSParentMassTolerance::getName ()]	= new FormItemMSMSParentMassTolerance ( &fvj, "10" );
	formItemMap [FormItemMSMSParentMassToleranceUnits::getName ()]= new FormItemMSMSParentMassToleranceUnits ();

	formItemMap [FormItemInstrumentName::getName ()]= new FormItemInstrumentName ();
	formItemMap [FormItemIT::getName ()]			= new FormItemIT ();

	formItemMap [FormItemConstMod::getName ()]	= new FormItemConstMod ( true, true );

	formItemMap [FormItemMSMSMaxReportedHits::getName ()]	= new FormItemMSMSMaxReportedHits (&fvj);
	formItemMap [FormItemCombinationType::getName ()]		= new FormItemCombinationType ();

	formItemMap [FormItemCompositionSearch::getName ()] = new FormItemCompositionSearch ();
}
void MSCompForm::printHTML ( ostream& os )
{
	init_html ( os, "MS-Comp", "compman.htm" );
	fvj.print ( os );
	printHTMLFORMStart ( os, "post", "mssearch", false, true );
	tableStart ( os, true );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true );
				formItemMap [FormItemMSMSParentMass::getName()]->printHTML ( os );
				formItemMap [FormItemMSMSParentCharge::getName()]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [FormItemMSMSParentMassTolerance::getName ()]->printHTML ( os );
				formItemMap [FormItemMSMSParentMassToleranceUnits::getName ()]->printHTML ( os );
				formItemMap [FormItemParentMassConvert::getName ()]->printHTML ( os );
				os << "<br />" << endl;
				saveResultsForm.printHTML ( os );
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 1, 2 );
				formItemMap [FormItemInstrumentName::getName ()]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [FormItemIT::getName ()]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true );
				formItemMap [FormItemConstMod::getName ()]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true );
				formItemMap [FormItemMSMSMaxReportedHits::getName ()]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [FormItemCombinationType::getName ()]->printHTML ( os );
				os << "<br />" << endl;
				os << "Note: If the mass tolerance is large and few ions below are<br />"; 
				os << "checked then the program will take a long time to run.";
				os << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "left", true );
				userAAForm.printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true, 2 );
				printHTMLFORMSubmit ( os, "Calculate Compositions" );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true, 2 );
				immoniumIonForm.printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true, 2 );
				extraIonForm.printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true, 2 );
				absentIonForm.printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true, 2 );
				modifiedIonForm.printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
	tableEnd ( os );
	showHiddenItems ( os );
	printHTMLFORMEnd ( os );
	printHTMLFooter ( os, "1996" );
}
void MSCompForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];

		formItemMap [FormItemMSMSParentMass::getName()]->setValue ( p );
		formItemMap [FormItemMSMSParentCharge::getName()]->setValue ( p );
		formItemMap [FormItemMSMSParentMassTolerance::getName ()]->setValue ( p );
		formItemMap [FormItemMSMSParentMassToleranceUnits::getName ()]->setValue ( p );
		formItemMap [FormItemParentMassConvert::getName ()]->setValue ( p );
		formItemMap [FormItemInstrumentName::getName ()]->setValue ( p );
		formItemMap [FormItemIT::getName ()]->setValue ( p );
		formItemMap [FormItemConstMod::getName ()]->setValue ( p );
		formItemMap [FormItemMSMSMaxReportedHits::getName ()]->setValue ( p );
		formItemMap [FormItemCombinationType::getName ()]->setValue ( p );
	}
}
