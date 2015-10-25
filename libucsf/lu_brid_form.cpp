/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_brid_form.cpp                                              *
*                                                                             *
*  Created    : January 7th 2005                                              *
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
#include <lu_html.h>
#include <lu_table.h>
#include <lu_srch_form.h>
#include <lu_brid_form.h>
#include <lu_form_valid.h>

using std::string;
using std::vector;
using std::ostream;
using std::endl;

class FormItemAccessMethod : public FormItemSelect {
	static const char* options [];
public:
	FormItemAccessMethod () :
		FormItemSelect ( "Retrieve Entry by", "allman.htm#retrieval_method", getName (), options, "Accession Number" ) {}
	static string getName () { return "access_method"; }
};
const char* FormItemAccessMethod::options [] = { "Index Number", "Accession Number", 0 };

class FormItemEntryData : public FormItemTextArea {
public:
	FormItemEntryData () :
		FormItemTextArea ( "", "", getName (), 4, 20, StringVector () )
		{
			value.push_back ( "P15497" );
		}
	static string getName () { return "entry_data"; }
};


class FormItemEndTerminus : public FormItemCheckbox {
public:
	FormItemEndTerminus () :
		FormItemCheckbox ( "End Terminus", "allman.htm#end_terminus", getName (), false ) {}
	static string getName () { return "end_terminus"; }
};

class FormItemStrippingTerminal : public FormItemSelect {
	static const char* options [];
public:
	FormItemStrippingTerminal () :
		FormItemSelect ( "Stripping Terminal", "allman.htm#end_terminus", getName (), options, "N" ) {}
	static string getName () { return "stripping_terminal"; }
};
const char* FormItemStrippingTerminal::options [] = { "N", "C", 0 };

class FormItemStartStrip : public FormItemText {
public:
	FormItemStartStrip ( FormValidatingJavascript* fvj ) :
		FormItemText ( "Stripping Range", "allman.htm#end_terminus", getName (), 2, 2, "0", fvj->addPositiveIntegerValidator ( getName (), "Start Stripping Range" ) ) {}
	static string getName () { return "start_strip"; }
};

class FormItemEndStrip : public FormItemText {
public:
	FormItemEndStrip ( FormValidatingJavascript* fvj ) :
		FormItemText ( "to", "", getName (), 2, 2, "4", fvj->addPositiveNonZeroIntegerValidator ( getName (), "End Stripping Range" ) ) {}
	static string getName () { return "end_strip"; }
};

class FormItemSeparateProteins : public FormItemCheckbox {
public:
	FormItemSeparateProteins () :
		FormItemCheckbox ( "Separate Proteins", "allman.htm#separate_proteins", getName (), false ) {}
	static string getName () { return "separate_proteins"; }
};


class FormItemHideHTMLLinks : public FormItemCheckbox {
public:
	FormItemHideHTMLLinks () :
		FormItemCheckbox ( "Hide HTML Links", "allman.htm#hide_html_links", getName (), false ) {}
	static string getName () { return "hide_html_links"; }
};

string MSSingleForm::defaultUserSequence = "\
MKAVVLTLAVLFLTGSQARHFWQQDDPQSSWDRVKDFATVYVEAIKDSGRDYVAQFEASALGKQLNLKLLDNWDTLASTL\n\
SKVREQLGPVTQEFWDNLEKETASLRQEMHKDLEEVKQKVQPYLDEFQKKWHEEVEIYRQKVAPLGEEFREGARQKVQEL\n\
QDKLSPLAQELRDRARAHVETLRQQLAPYSDDLRQRLTARLEALKEGGGSLAEYHAKASEQLKALGEKAKPVLEDLRQGL\n\
LPVLESLKVSILAAIDEASKKLNAQ";

const string MSSingleForm::getDefaultUserSequence ()
{
	return defaultUserSequence;
}
MSSingleForm::MSSingleForm ( const VectorConstParameterListPtr& params, bool enzymeOption, bool tabDelimitedTextOption, FormValidatingJavascript* fvj ) :
	fvj ( fvj ),
	enzymeOption ( enzymeOption ),
	saveResultsForm ( params, fvj, false, tabDelimitedTextOption ),
	compIonForm ( 0 )
{
	if ( enzymeOption ) compIonForm = new PresentCompIonForm ( params, "ACDEFGHIKLMNPQRSTVWY" );
	create ( params );
}
MSSingleForm::~MSSingleForm ()
{
	if ( enzymeOption ) delete compIonForm;
}
void MSSingleForm::createItems ()
{
	formItemMap [FormItemDatabase::getName ()]	= new FormItemDatabase ( true );

	formItemMap [FormItemAccessMethod::getName ()]= new FormItemAccessMethod ();

	formItemMap [FormItemNTerminusAALimit::getName ()] = new FormItemNTerminusAALimit ( fvj );

	formItemMap [FormItemEntryData::getName ()]= new FormItemEntryData ();

	if ( enzymeOption ) {
#ifdef CHEM_SCORE
		formItemMap [FormItemChemScore::getName ()]			= new FormItemChemScore ();
		formItemMap [FormItemMetOxFactor::getName ()]		= new FormItemMetOxFactor (fvj);
#endif
		formItemMap [FormItemEnzyme::getName ()]			= new FormItemEnzyme ( "Trypsin", true );
		formItemMap [FormItemEndTerminus::getName ()]		= new FormItemEndTerminus ();
		formItemMap [FormItemStrippingTerminal::getName ()]	= new FormItemStrippingTerminal ();
		formItemMap [FormItemStartStrip::getName ()]		= new FormItemStartStrip (fvj);
		formItemMap [FormItemEndStrip::getName ()]			= new FormItemEndStrip (fvj);
		formItemMap [FormItemMissedCleavages::getName ()]	= new FormItemMissedCleavages (fvj);
		formItemMap [FormItemCompMaskType::getName ()]		= new FormItemCompMaskType ();
	}
	formItemMap [FormItemConstMod::getName ()] = new FormItemConstMod ( true, true );
}
void MSSingleForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];
		formItemMap [FormItemDatabase::getName ()]->setValue ( p );
		formItemMap [FormItemAccessMethod::getName ()]->setValue ( p );
		formItemMap [FormItemNTerminusAALimit::getName ()]->setValue ( p );
		formItemMap [FormItemEntryData::getName ()]->setValue ( p );
		if ( enzymeOption ) {
			formItemMap [FormItemEnzyme::getName ()]->setValue ( p );
			formItemMap [FormItemEndTerminus::getName ()]->setValue ( p );
			formItemMap [FormItemStrippingTerminal::getName ()]->setValue ( p );
			formItemMap [FormItemStartStrip::getName ()]->setValue ( p );
			formItemMap [FormItemEndStrip::getName ()]->setValue ( p );
			formItemMap [FormItemMissedCleavages::getName ()]->setValue ( p );
			formItemMap [FormItemCompMaskType::getName ()]->setValue ( p );
		}
		formItemMap [FormItemConstMod::getName ()]->setValue ( p );
	}
}
void MSSingleForm::printHTML ( ostream& os )
{
	tableRowStart ( os );
		tableHeaderStart ( os, "", "left" );
			formItemMap [FormItemDatabase::getName ()]->printHTML ( os );
			os << "<br />" << endl;
			saveResultsForm.printHTML ( os );
		tableHeaderEnd ( os );
		tableHeaderStart ( os, "", "left", true );
			formItemMap [FormItemAccessMethod::getName ()]->printHTML ( os );
			os << "<br />" << endl;
			formItemMap [FormItemNTerminusAALimit::getName ()]->printHTML ( os );
		tableHeaderEnd ( os );
	tableRowEnd ( os );
	tableRowStart ( os );
		tableHeaderStart ( os, "", "left" );
			if ( enzymeOption ) {
				formItemMap [FormItemEnzyme::getName ()]->printHTML ( os );
				formItemMap [FormItemMissedCleavages::getName ()]->printHTML ( os );
				os << "<br />" << endl;
#ifdef CHEM_SCORE
				formItemMap [FormItemChemScore::getName ()]->printHTML ( os );
				formItemMap [FormItemMetOxFactor::getName ()]->printHTML ( os );
				os << "<br />" << endl;
#endif
				ExpandableJavascriptBlock ejb ( "End Terminus Parameters" );
				ejb.printHeader ( os );
					formItemMap [FormItemEndTerminus::getName ()]->printHTML ( os );
					os << "<br />" << endl;
					formItemMap [FormItemStrippingTerminal::getName ()]->printHTML ( os );
					formItemMap [FormItemStartStrip::getName ()]->printHTML ( os );
					formItemMap [FormItemEndStrip::getName ()]->printHTML ( os );
				ejb.printFooter ( os );
			}
			formItemMap [FormItemConstMod::getName ()]->printHTML ( os );
			if ( enzymeOption ) {
				ExpandableJavascriptBlock ejb2 ( "Present Amino Acids<a href=\"../html/instruct/allman.htm#present_amino_acids\">(see instructions)</a>" );
				ejb2.printHeader ( os );
					formItemMap [FormItemCompMaskType::getName ()]->printHTML ( os );
					compIonForm->printHTML ( os );
				ejb2.printFooter ( os );
			}
		tableHeaderEnd ( os );
		tableHeaderStart ( os, "", "center" );
			os << "List of Entries<a href=\"../html/instruct/allman.htm#list_of_entries\">(see instructions)</a>" << endl;
			os << "<br />" << endl;
			formItemMap [FormItemEntryData::getName ()]->printHTML ( os );
		tableHeaderEnd ( os );
	tableRowEnd ( os );
	tableRowStart ( os );
		tableHeaderStart ( os, "", "center", true, 2 );
			printHTMLFORMSubmit ( os, "Perform Digest" );
		tableHeaderEnd ( os );
	tableRowEnd ( os );
}
MSBridgeForm::MSBridgeForm ( const VectorConstParameterListPtr& params ) :
	fvj (),
	singleForm ( params, true, false, &fvj ),
	userAAForm ( params, &fvj ),
	variableModsForm ( params, "", &fvj ),
	msToleranceForm ( &fvj, params ),
	crosslinkingForm ( params, &fvj )
{
	create ( params );
}
void MSBridgeForm::createItems ()
{
	formItemMap [FormItemSearchName::getName ()] = new FormItemSearchName ( "msbridge" );
	formItemMap [FormItemReportTitle::getName ()]	= new FormItemReportTitle ( "MS-Bridge" );
	formItemMap [FormItemVersion::getName ()]	= new FormItemVersion ();
	
	formItemMap [FormItemParentContaminantMasses::getName()] = new FormItemParentContaminantMasses ();

	formItemMap [FormItemParentMassConvert::getName ()]			= new FormItemParentMassConvert ( "MS" );

	formItemMap [FormItemHideProteinSequence::getName()]= new FormItemHideProteinSequence ( false );
	formItemMap [FormItemHideHTMLLinks::getName()]		= new FormItemHideHTMLLinks ();
	formItemMap [FormItemDisplayGraph::getName()]		= new FormItemDisplayGraph ();
	formItemMap [FormItemSeparateProteins::getName ()]	= new FormItemSeparateProteins ();

	static char defaultSeq [] = "RVCMGKSQHHSFPCISDRLCSNECVKEEGGWTAGYCHLRYCRCQKAC";
	formItemMap [FormItemUserProteinSequence::getName ()]	= new FormItemUserProteinSequence ( false, defaultSeq );
}
void MSBridgeForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];

		formItemMap [FormItemParentContaminantMasses::getName ()]->setValue ( p );
		formItemMap [FormItemParentMassConvert::getName ()]->setValue ( p );
		formItemMap [FormItemHideProteinSequence::getName()]->setValue ( p );
		formItemMap [FormItemHideHTMLLinks::getName()]->setValue ( p );
		formItemMap [FormItemDisplayGraph::getName()]->setValue ( p );
		formItemMap [FormItemSeparateProteins::getName ()]->setValue ( p );
		formItemMap [FormItemUserProteinSequence::getName ()]->setValue ( p );
	}
}
MSBridgeFileForm::MSBridgeFileForm ( const VectorConstParameterListPtr& params ) :
	MSBridgeForm ( params ),
	dataForm ( params, "msbridge" )
{
}
void MSBridgeFileForm::printHTML ( ostream& os )
{
	init_html ( os, "MS-Bridge", "bridgeman.htm" );
	fvj.print ( os );
	printHTMLFORMStart ( os, "post", "mssearch", false, true );
	tableStart ( os, true );
		MSBridgeForm::printHTML ( os );
		dataForm.printHTML ( os );
	tableEnd ( os );
	showHiddenItems ( os );
	printHTMLFORMEnd ( os );
	printHTMLFooter ( os, "1999" );
}
MSBridgeStandardForm::MSBridgeStandardForm ( const VectorConstParameterListPtr& params ) :
	MSBridgeForm ( params ),
	dataForm ( params, "msbridge" )
{
}
void MSBridgeStandardForm::printHTML ( ostream& os )
{
	init_html ( os, "MS-Bridge", "bridgeman.htm" );
	fvj.print ( os );
	printHTMLFORMStart ( os, "post", "mssearch", false, true );
	tableStart ( os, true );
		MSBridgeForm::printHTML ( os );
		dataForm.printHTML ( os );
	tableEnd ( os );
	showHiddenItems ( os );
	printHTMLFORMEnd ( os );
	printHTMLFooter ( os, "1999" );
}
void MSBridgeForm::printHTML ( ostream& os )
{
	basicJavascriptFunctions ( os );
	msBridgeJavascriptFunctions ( os );
	singleForm.printHTML ( os );
	tableRowStart ( os );
		tableHeaderStart ( os, "", "left", true );
			variableModsForm.printHTML ( os );
			formItemMap [FormItemParentContaminantMasses::getName ()]->printHTML ( os );
		tableHeaderEnd ( os );
		tableHeaderStart ( os, "", "left", true );
			crosslinkingForm.printHTML ( os );
		tableHeaderEnd ( os );
	tableRowEnd ( os );
	tableRowStart ( os );
		tableHeaderStart ( os, "", "left", true );
			formItemMap [FormItemParentMassConvert::getName ()]->printHTML ( os );
			os << "<br />" << endl;
			msToleranceForm.printHTML ( os );
		tableHeaderEnd ( os );
		tableHeaderStart ( os, "", "left", true );
			formItemMap [FormItemHideProteinSequence::getName()]->printHTML ( os );
			os << "<br />" << endl;
			formItemMap [FormItemHideHTMLLinks::getName()]->printHTML ( os );
			os << "<br />" << endl;
			formItemMap [FormItemDisplayGraph::getName()]->printHTML ( os );
			formItemMap [FormItemSeparateProteins::getName ()]->printHTML ( os );
		tableHeaderEnd ( os );
	tableRowEnd ( os );
	tableRowStart ( os );
		tableHeaderStart ( os, "", "left", true, 2 );
			os << "<div class=\"form_medium_label\">For digestion of a user supplied sequence select <b>User Protein</b> above.</div>" << endl;
		tableHeaderEnd ( os );
	tableRowEnd ( os );
	tableRowStart ( os );
		tableHeaderStart ( os, "", "left", true, 2 );
			formItemMap [FormItemUserProteinSequence::getName ()]->printHTML ( os );
		tableHeaderEnd ( os );
	tableRowEnd ( os );
	tableRowStart ( os );
		tableHeaderStart ( os, "", "left", true, 2 );
			userAAForm.printHTML ( os );
		tableHeaderEnd ( os );
	tableRowEnd ( os );
}
class FormItemMinDigestFragmentMass : public FormItemText {
public:
	FormItemMinDigestFragmentMass ( FormValidatingJavascript* fvj ) :
		FormItemText ( "Peptide Mass", "digestman.htm#min_fragment_mass", getName (), 6, 10, "800.0", fvj->addPositiveFloatingPointValidator (getName (), "Min Fragment Mass") ) {}
	static string getName () { return "min_digest_fragment_mass"; }
};
class FormItemMaxDigestFragmentMass : public FormItemText {
public:
	FormItemMaxDigestFragmentMass ( FormValidatingJavascript* fvj ) :
		FormItemText ( "to", "digestman.htm#max_fragment_mass", getName (), 6, 10, "4000.0", fvj->addPositiveFloatingPointValidator (getName (), "Max Fragment Mass") ) {}
	static string getName () { return "max_digest_fragment_mass"; }
};
class FormItemMinFragmentLength : public FormItemText {
public:
	FormItemMinFragmentLength ( FormValidatingJavascript* fvj ) :
		FormItemText ( "Min Peptide Length", "digestman.htm#min_fragment_length", getName (), 2, 2, "5", fvj->addPositiveNonZeroIntegerValidator (getName (), "Min Fragment Length") ) {}
	static string getName () { return "min_digest_fragment_length"; }
};
class FormItemReportMultCharge : public FormItemCheckbox {
public:
	FormItemReportMultCharge () :
		FormItemCheckbox ( "Report Multiple Charges", "digestman.htm#report_mult_charge", getName (), false ) {}
	static string getName () { return "report_mult_charge"; }
};
class FormItemBullBreese : public FormItemCheckbox {
public:
	FormItemBullBreese () :
		FormItemCheckbox ( "Bull Breese Indicies", "digestman.htm#bull_breese", getName (), false ) {}
	static string getName () { return "bull_breese"; }
};
class FormItemHPLCIndex : public FormItemCheckbox {
public:
	FormItemHPLCIndex () :
		FormItemCheckbox ( "HPLC Indicies", "digestman.htm#hplc_indicies", getName (), false ) {}
	static string getName () { return "hplc_index"; }
};
MSDigestForm::MSDigestForm ( const VectorConstParameterListPtr& params ) :
	fvj (),
	singleForm ( params, true, true, &fvj ),
	userAAForm ( params, &fvj ),
	variableModsForm ( params, "", &fvj )
{
	create ( params );
}
void MSDigestForm::createItems ()
{
	formItemMap [FormItemSearchName::getName ()]	= new FormItemSearchName ( "msdigest" );
	formItemMap [FormItemReportTitle::getName ()]	= new FormItemReportTitle ( "MS-Digest" );
	formItemMap [FormItemVersion::getName ()]	= new FormItemVersion ();

	formItemMap [FormItemMinDigestFragmentMass::getName ()]	= new FormItemMinDigestFragmentMass (&fvj);
	formItemMap [FormItemMaxDigestFragmentMass::getName ()]	= new FormItemMaxDigestFragmentMass (&fvj);
	formItemMap [FormItemMinFragmentLength::getName ()]		= new FormItemMinFragmentLength (&fvj);

	formItemMap [FormItemHideProteinSequence::getName()]= new FormItemHideProteinSequence ( false );
	formItemMap [FormItemHideHTMLLinks::getName()]		= new FormItemHideHTMLLinks ();
	formItemMap [FormItemReportMultCharge::getName()]	= new FormItemReportMultCharge ();

	formItemMap [FormItemBullBreese::getName ()]		= new FormItemBullBreese ();
	formItemMap [FormItemHPLCIndex::getName ()]			= new FormItemHPLCIndex ();
	formItemMap [FormItemSeparateProteins::getName ()]	= new FormItemSeparateProteins ();

	formItemMap [FormItemUserProteinSequence::getName ()]	= new FormItemUserProteinSequence ( false, MSSingleForm::getDefaultUserSequence () );

	formItemMap [FormItemInstrumentName::getName ()]= new FormItemInstrumentName ();
	formItemMap [FormItemMaxHits::getName ()]	= new FormItemMaxHits ("50000");
}
void MSDigestForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];

		formItemMap [FormItemSearchName::getName ()]->setValue ( p );
		formItemMap [FormItemReportTitle::getName ()]->setValue ( p );
		formItemMap [FormItemVersion::getName ()]->setValue ( p );

		formItemMap [FormItemMinDigestFragmentMass::getName ()]->setValue ( p );
		formItemMap [FormItemMaxDigestFragmentMass::getName ()]->setValue ( p );
		formItemMap [FormItemMinFragmentLength::getName ()]->setValue ( p );

		formItemMap [FormItemHideProteinSequence::getName()]->setValue ( p );
		formItemMap [FormItemHideHTMLLinks::getName()]->setValue ( p );
		formItemMap [FormItemReportMultCharge::getName()]->setValue ( p );

		formItemMap [FormItemBullBreese::getName ()]->setValue ( p );
		formItemMap [FormItemHPLCIndex::getName ()]->setValue ( p );
		formItemMap [FormItemSeparateProteins::getName ()]->setValue ( p );

		formItemMap [FormItemUserProteinSequence::getName ()]->setValue ( p );

		formItemMap [FormItemInstrumentName::getName ()]->setValue ( p );
		formItemMap [FormItemMaxHits::getName ()]->setValue ( p );
	}
}
void MSDigestForm::printHTML ( ostream& os )
{
	init_html ( os, "MS-Digest", "digestman.htm" );
	fvj.print ( os );
	printHTMLFORMStart ( os, "post", "mssearch", false, true );
	basicJavascriptFunctions ( os );
	tableStart ( os, true );
		singleForm.printHTML ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true );
				variableModsForm.printHTML ( os );
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "left", true );
				formItemMap [FormItemMinDigestFragmentMass::getName ()]->printHTML ( os );
				formItemMap [FormItemMaxDigestFragmentMass::getName ()]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [FormItemMinFragmentLength::getName ()]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [FormItemHideProteinSequence::getName()]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [FormItemHideHTMLLinks::getName()]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [FormItemReportMultCharge::getName()]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [FormItemBullBreese::getName ()]->printHTML ( os );
				formItemMap [FormItemHPLCIndex::getName ()]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [FormItemSeparateProteins::getName ()]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true, 2 );
			os << "<div class=\"form_medium_label\">For digestion of a user supplied sequence select <b>User Protein</b> above.</div>" << endl;
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true, 2 );
				formItemMap [FormItemUserProteinSequence::getName ()]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true, 2 );
				userAAForm.printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true, 2 );
				formItemMap [FormItemInstrumentName::getName ()]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
	tableEnd ( os );
	showHiddenItems ( os );
	printHTMLFORMEnd ( os );
	printHTMLFooter ( os, "1995" );
}
MSNonSpecificForm::MSNonSpecificForm ( const VectorConstParameterListPtr& params ) :
	fvj (),
	singleForm ( params, false, false, &fvj ),
	userAAForm ( params, &fvj ),
	dataForm ( params, "msnonspecific" ),
	msToleranceForm ( &fvj, params )
{
	create ( params );
}
void MSNonSpecificForm::createItems ()
{
	formItemMap [FormItemSearchName::getName ()]	= new FormItemSearchName ( "msnonspecific" );
	formItemMap [FormItemReportTitle::getName ()]	= new FormItemReportTitle ( "MS-NonSpecific" );
	formItemMap [FormItemVersion::getName ()]	= new FormItemVersion ();

	formItemMap [FormItemHideProteinSequence::getName()]= new FormItemHideProteinSequence ( false );
	formItemMap [FormItemHideHTMLLinks::getName()]		= new FormItemHideHTMLLinks ();

	formItemMap [FormItemDisplayGraph::getName()]		= new FormItemDisplayGraph ( false );
	formItemMap [FormItemSeparateProteins::getName ()]	= new FormItemSeparateProteins ();

	formItemMap [FormItemParentMassConvert::getName ()]			= new FormItemParentMassConvert ( "MS" );

	formItemMap [FormItemParentContaminantMasses::getName()] = new FormItemParentContaminantMasses ();

	formItemMap [FormItemUserProteinSequence::getName ()]	= new FormItemUserProteinSequence ( false, MSSingleForm::getDefaultUserSequence () );
	formItemMap [FormItemMaxHits::getName ()]	= new FormItemMaxHits ("50000");
}
void MSNonSpecificForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];

		formItemMap [FormItemSearchName::getName ()]->setValue ( p );
		formItemMap [FormItemReportTitle::getName ()]->setValue ( p );
		formItemMap [FormItemVersion::getName ()]->setValue ( p );

		formItemMap [FormItemHideProteinSequence::getName()]->setValue ( p );
		formItemMap [FormItemHideHTMLLinks::getName()]->setValue ( p );

		formItemMap [FormItemDisplayGraph::getName()]->setValue ( p );
		formItemMap [FormItemSeparateProteins::getName ()]->setValue ( p );

		formItemMap [FormItemParentMassConvert::getName ()]->setValue ( p );

		formItemMap [FormItemParentContaminantMasses::getName()]->setValue ( p );

		formItemMap [FormItemUserProteinSequence::getName ()]->setValue ( p );
		formItemMap [FormItemMaxHits::getName ()]->setValue ( p );
	}
}
void MSNonSpecificForm::printHTML ( ostream& os )
{
	init_html ( os, "MS-NonSpecific", "nonspecificman.htm" );
	fvj.print ( os );
	printHTMLFORMStart ( os, "post", "mssearch", false, true );
	tableStart ( os, true );
		singleForm.printHTML ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true );
				formItemMap [FormItemHideProteinSequence::getName()]->printHTML ( os );
				formItemMap [FormItemHideHTMLLinks::getName()]->printHTML ( os );
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "left", true );
				formItemMap [FormItemDisplayGraph::getName()]->printHTML ( os );
				formItemMap [FormItemSeparateProteins::getName ()]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true );
				formItemMap [FormItemParentMassConvert::getName ()]->printHTML ( os );
				os << "<br />" << endl;
				msToleranceForm.printHTML ( os );
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "left", true );
				formItemMap [FormItemParentContaminantMasses::getName ()]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true, 2 );
				os << "<div class=\"form_medium_label\">For digestion of a user supplied sequence select <b>User Protein</b> above.</div>" << endl;
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true, 2 );
				formItemMap [FormItemUserProteinSequence::getName ()]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true, 2 );
				userAAForm.printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		dataForm.printHTML ( os );
	tableEnd ( os );
	showHiddenItems ( os );
	printHTMLFORMEnd ( os );
	printHTMLFooter ( os, "2002" );
}
