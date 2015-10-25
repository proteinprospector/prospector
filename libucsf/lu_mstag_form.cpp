/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_mstag_form.cpp                                             *
*                                                                             *
*  Created    : December 7th 2004                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2004-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lu_html.h>
#include <lu_table.h>
#include <lu_mstag_form.h>
using std::ostream;
using std::string;
using std::vector;
using std::endl;

static string SCORE_HISTOGRAM_ONLY = "score_histogram_only";

class FormItemMSSeqSearchType : public FormItemSelect {
	static const char* options [];
public:
	FormItemMSSeqSearchType () :
		FormItemSelect ( "Search Type", "", getName (), options, "no errors" ) {}
	static string getName () { return "msms_search_type"; }
};
const char* FormItemMSSeqSearchType::options [] = {
	"no errors", "high mass error", "low mass error", "middle masses error", "parent mass", 0
};

class FormItemMSSeqIonType : public FormItemSelect {
	static const char* options [];
public:
	FormItemMSSeqIonType () :
		FormItemSelect ( "Fragment ion series", "", getName (), options, "b" ) {}
	static string getName () { return "ion_type"; }
};
const char* FormItemMSSeqIonType::options [] = { "a", "b", "c", "y", 0 };

class FormItemCompositionExclude : public FormItemText {
public:
	FormItemCompositionExclude () :
		FormItemText ( "Composition Exclude", "", getName (), 40, 40, "" ) {}
	static string getName () { return "composition_exclude"; }
};

class FormItemRegularExpression : public FormItemText {
public:
	FormItemRegularExpression () :
		FormItemText ( "Regular Expression", "", getName (), 40, 100, "" ) {}
	static string getName () { return "regular_expression"; }
};

MSSeqForm::MSSeqForm ( const VectorConstParameterListPtr& params ) :
	fvj (),
	headerForm ( params, "msseq", "MS-Seq" ),
	searchForm ( params, "msseq", &fvj ),
	msmsToleranceForm ( &fvj, params, false, false ),
	immoniumIonForm ( params ),
	dataForm ( params, "msseq" )
{
	create ( params );
}
void MSSeqForm::createItems ()
{
	formItemMap [FormItemComment::getName()]				= new FormItemComment ();
	formItemMap [FormItemMaxReportedHits::getName("msms_")]	= new FormItemMaxReportedHits ( "msms_", "25", &fvj, false );

	formItemMap [FormItemMSSeqSearchType::getName()]	= new FormItemMSSeqSearchType ();
	formItemMap [FormItemMSSeqIonType::getName()]		= new FormItemMSSeqIonType ();
	formItemMap [FormItemCompositionExclude::getName()]	= new FormItemCompositionExclude ();
	formItemMap [FormItemRegularExpression::getName()]	= new FormItemRegularExpression ();

	formItemMap [FormItemMSMSPkFilter::getName ()]		= new FormItemMSMSPkFilter ( true );
}
void MSSeqForm::printHTML ( ostream& os )
{
	init_html ( os, "MS-Seq", "seqman.htm" );
	fvj.print ( os );
	printHTMLFORMStart ( os, "post", "mssearch", false, true );
	headerForm.printHTML ( os );
	tableStart ( os, true );
		searchForm.printHTML ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				printHTMLFORMSubmit ( os, "Start Search" );
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "left", true );
				formItemMap [FormItemComment::getName()]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [FormItemMaxReportedHits::getName("msms_")]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true );
				formItemMap [FormItemMSSeqSearchType::getName()]->printHTML ( os );
				formItemMap [FormItemMSSeqIonType::getName()]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [FormItemCompositionExclude::getName()]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [FormItemRegularExpression::getName()]->printHTML ( os );
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "left", true );
				msmsToleranceForm.printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true, 2 );
				immoniumIonForm.printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		dataForm.printHTML ( os );
	tableEnd ( os );
	showHiddenItems ( os );
	printHTMLFORMEnd ( os );
	printHTMLFooter ( os, "1995" );
}
void MSSeqForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];

		formItemMap [FormItemComment::getName()]->setValue ( p );
		formItemMap [FormItemMaxReportedHits::getName("msms_")]->setValue ( p );
		formItemMap [FormItemMSSeqSearchType::getName()]->setValue ( p );
		formItemMap [FormItemMSSeqIonType::getName()]->setValue ( p );
		formItemMap [FormItemCompositionExclude::getName()]->setValue ( p );
		formItemMap [FormItemRegularExpression::getName()]->setValue ( p );
		formItemMap [FormItemMSMSPkFilter::getName ()]->setValue ( p );
	}
}
MSTagForm::MSTagForm ( const VectorConstParameterListPtr& params, bool filterPeaks ) :
	fvj (),
	headerForm ( params, "mstag", "MS-Tag" ),
	searchForm ( params, "mstag", &fvj ),
	msmsToleranceForm ( &fvj, params, true, false ),
	variableModsForm ( params, "msms_", &fvj ),
	extraUserModsForm ( params, &fvj, true ),
	massModificationForm ( params, &fvj ),
	crosslinkingForm ( params, &fvj, true ),
	matrixModeForm ( params, "msms" ),
	filterPeaks ( filterPeaks )
{
	create ( params );
}
void MSTagForm::createItems ()
{
	formItemMap [FormItemComment::getName()]				= new FormItemComment ();
	formItemMap [FormItemMaxReportedHits::getName("msms_")]	= new FormItemMaxReportedHits ( "msms_", "5",  &fvj, false  );
	formItemMap [FormItemExpectationCalculationMethod::getName ()]	= new FormItemExpectationCalculationMethod ( "None" );
	formItemMap [FormItemMaxHits::getName ()]				= new FormItemMaxHits ();
	formItemMap [FormItemDisplayGraph::getName()]			= new FormItemDisplayGraph ();

	formItemMap [FormItemMSMSMaxModifications::getName ()]	= new FormItemMSMSMaxModifications ( &fvj );
	formItemMap [FormItemMSMSMaxPeptidePermutations::getName ()]= new FormItemMSMSMaxPeptidePermutations ( &fvj );

	formItemMap [FormItemUseInstrumentIonTypes::getName ()]	= new FormItemUseInstrumentIonTypes ();
	formItemMap [FormItemMSMSPkFilter::getName ()]			= new FormItemMSMSPkFilter ( !filterPeaks );
	formItemMap [FormItemMaxMSMSPeaks::getName ()]			= new FormItemMaxMSMSPeaks ( &fvj, true );
	formItemMap [FormItemMSMSMinPrecursorMass::getName ()]	= new FormItemMSMSMinPrecursorMass ( &fvj, true );

	formItemMap [SCORE_HISTOGRAM_ONLY] = new FormItemCheckbox ( "", "", SCORE_HISTOGRAM_ONLY, false );
}
void MSTagForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];
		formItemMap [FormItemComment::getName()]->setValue ( p );
		formItemMap [FormItemMaxReportedHits::getName("msms_")]->setValue ( p );
		formItemMap [FormItemExpectationCalculationMethod::getName ()]->setValue ( p );
		formItemMap [FormItemMSMSPkFilter::getName ()]->setValue ( p );
		formItemMap [FormItemMaxMSMSPeaks::getName ()]->setValue ( p );
		formItemMap [FormItemMSMSMinPrecursorMass::getName ()]->setValue ( p );
		formItemMap [FormItemMaxHits::getName ()]->setValue ( p );
		formItemMap [FormItemDisplayGraph::getName()]->setValue ( p );

		formItemMap [FormItemMSMSMaxModifications::getName ()]->setValue ( p );
		formItemMap [FormItemMSMSMaxPeptidePermutations::getName ()]->setValue ( p );

		formItemMap [FormItemUseInstrumentIonTypes::getName ()]->setValue ( p );

		formItemMap [SCORE_HISTOGRAM_ONLY]->setValue ( p );
	}
}
void MSTagForm::printHTMLFirstTable ( ostream& os )
{
	tableRowStart ( os );
		tableHeaderStart ( os, "", "center", true );
			printHTMLFORMSubmit ( os, "Start Search" );
		tableHeaderEnd ( os );
		tableHeaderStart ( os, "", "left", true );
			formItemMap [FormItemComment::getName()]->printHTML ( os );
			os << "<br />" << endl;
			formItemMap [FormItemMaxReportedHits::getName("msms_")]->printHTML ( os );
			formItemMap [FormItemMSMSPkFilter::getName ()]->printHTML ( os );
			formItemMap [FormItemMaxMSMSPeaks::getName ()]->printHTML ( os );
			os << "<br />" << endl;
			formItemMap [FormItemExpectationCalculationMethod::getName ()]->printHTML ( os );
		tableHeaderEnd ( os );
	tableRowEnd ( os );
}
void MSTagForm::printHTMLSecondTable ( ostream& os )
{
	basicJavascriptFunctions ( os );
	massModCrosslinkingJavascriptFunctions ( os );
	tableRowStart ( os );
		tableHeaderStart ( os, "", "left", true );
			tableStart ( os );
				tableRowStart ( os );
					tableHeaderStart ( os, "", "left", true );
						variableModsForm.printHTML ( os );
					tableHeaderEnd ( os );
				tableRowEnd ( os );
				tableRowStart ( os );
					tableHeaderStart ( os, "", "left", true );
						extraUserModsForm.printHTML ( os );
					tableHeaderEnd ( os );
				tableRowEnd ( os );
				tableRowStart ( os );
					tableHeaderStart ( os, "", "left", true );
						formItemMap [FormItemMSMSMaxModifications::getName ()]->printHTML ( os );
						formItemMap [FormItemMSMSMaxPeptidePermutations::getName ()]->printHTML ( os );
					tableHeaderEnd ( os );
				tableRowEnd ( os );
			tableEnd ( os );
			massModificationForm.printHTML ( os );
			os << "<br />" << endl;
			divStart ( os, "div_cl", massModificationForm.getDivCL () );
				crosslinkingForm.printHTML ( os );
				os << "<br />" << endl;
			divEnd ( os );
			matrixModeForm.printHTML ( os );
		tableHeaderEnd ( os );
		tableHeaderStart ( os, "", "left", true );
			msmsToleranceForm.printHTML ( os );
		tableHeaderEnd ( os );
	tableRowEnd ( os );
}
MSTagFileForm::MSTagFileForm ( const VectorConstParameterListPtr& params ) :
	MSTagForm ( params, true ),
	dataForm ( params, "mstag" )
{
	create ( params );
}
void MSTagFileForm::printHTML ( ostream& os )
{
	init_html ( os, "MS-Tag from File", "batchtagman.htm#mstag_from_file" );
	fvj.print ( os );
	printHTMLFORMStart ( os, "post", "mssearch", false, true );
	headerForm.printHTML ( os );
	tableStart ( os, true );
		searchForm.printHTML ( os );
		printHTMLFirstTable ( os );
		printHTMLSecondTable ( os );
		dataForm.showHiddenItems ( os );
	tableEnd ( os );
	showHiddenItems ( os );
	printHTMLFORMEnd ( os );
	printHTMLFooter ( os, "1995" );
}
MSTagFileFormFromViewer::MSTagFileFormFromViewer ( const VectorConstParameterListPtr& params ) :
	MSTagForm ( params, true ),
	dataForm ( params )
{
	create ( params );
}
void MSTagFileFormFromViewer::printHTML ( ostream& os )
{
	init_html ( os, "MS-Tag from File", "batchtagman.htm#mstag_from_file" );
	fvj.print ( os );
	printHTMLFORMStart ( os, "post", "mssearch", false, true );
	headerForm.printHTML ( os );
	tableStart ( os, true );
		searchForm.printHTML ( os );
		printHTMLFirstTable ( os );
		printHTMLSecondTable ( os );
		dataForm.showHiddenItems ( os );
	tableEnd ( os );
	showHiddenItems ( os );
	printHTMLFORMEnd ( os );
	printHTMLFooter ( os, "1995" );
}
MSTagStandardForm::MSTagStandardForm ( const VectorConstParameterListPtr& params ) :
	MSTagForm ( params, false ),
	immoniumIonForm ( params ),
	ionTypeForm ( params ),
	dataForm ( params, "mstag" )
{
	create ( params );
}
void MSTagStandardForm::printHTML ( ostream& os )
{
	init_html ( os, "MS-Tag", "tagman.htm" );
	fvj.print ( os );
	printHTMLFORMStart ( os, "post", "mssearch", false, true );
	headerForm.printHTML ( os );
	tableStart ( os, true );
		searchForm.printHTML ( os );
		printHTMLFirstTable ( os );
		printHTMLSecondTable ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true, 2 );
				ExpandableJavascriptBlock ejb ( "Instrument Ion Types" );
				ejb.printHeader ( os );
				os << "<br />" << endl;
				formItemMap [FormItemUseInstrumentIonTypes::getName ()]->printHTML ( os );
				os << "<br />" << endl;
				os << "<br />" << endl;
				ionTypeForm.printHTML ( os );
				ejb.printFooter ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true, 2 );
				ExpandableJavascriptBlock ejb2 ( "Composition Ions" );
				ejb2.printHeader ( os );
				os << "<br />" << endl;
				immoniumIonForm.printHTML ( os );
				ejb2.printFooter ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		dataForm.printHTML ( os );
	tableEnd ( os );
	showHiddenItems ( os );
	printHTMLFORMEnd ( os );
	printHTMLFooter ( os, "1995" );
}
MSTagIonTypeForm::MSTagIonTypeForm ( const VectorConstParameterListPtr& params )
{
	create ( params );
}
void MSTagIonTypeForm::createItems ()
{
	formItemMap [FormItemIonType::getName ()+string("i")]	= new FormItemIonType ( "", "i", false );
	formItemMap [FormItemIonType::getName ()+string("a")]	= new FormItemIonType ( "", "a", false );
	formItemMap [FormItemIonType::getName ()+string("b")]	= new FormItemIonType ( "", "b", false );
	formItemMap [FormItemIonType::getName ()+string("B")]	= new FormItemIonType ( "", "B", false );
	formItemMap [FormItemIonType::getName ()+string("c_plus_2")] = new FormItemIonType ( "", "c+2", false );
	formItemMap [FormItemIonType::getName ()+string("c_plus_1")] = new FormItemIonType ( "", "c+1", false );
	formItemMap [FormItemIonType::getName ()+string("c")]	= new FormItemIonType ( "", "c", false );
	formItemMap [FormItemIonType::getName ()+string("c_minus_1")] = new FormItemIonType ( "", "c-1", false );
	formItemMap [FormItemIonType::getName ()+string("x")]	= new FormItemIonType ( "", "x", false );
	formItemMap [FormItemIonType::getName ()+string("y")]	= new FormItemIonType ( "", "y", false );
	formItemMap [FormItemIonType::getName ()+string("Y")]	= new FormItemIonType ( "", "Y", false );
	formItemMap [FormItemIonType::getName ()+string("z")]	= new FormItemIonType ( "", "z", false );
	formItemMap [FormItemIonType::getName ()+string("z_plus_1")] = new FormItemIonType ( "", "z+1", false );
	formItemMap [FormItemIonType::getName ()+string("z_plus_2")] = new FormItemIonType ( "", "z+2", false );
	formItemMap [FormItemIonType::getName ()+string("z_plus_3")] = new FormItemIonType ( "", "z+3", false );
	formItemMap [FormItemIonType::getName ()+string("d")]	= new FormItemIonType ( "", "d", false );
	formItemMap [FormItemIonType::getName ()+string("v")]	= new FormItemIonType ( "", "v", false );
	formItemMap [FormItemIonType::getName ()+string("w")]	= new FormItemIonType ( "", "w", false );
	formItemMap [FormItemIonType::getName ()+string("n")]	= new FormItemIonType ( "", "n", false );
	formItemMap [FormItemIonType::getName ()+string("h")]	= new FormItemIonType ( "", "h", false );
	formItemMap [FormItemIonType::getName ()+string("P")]	= new FormItemIonType ( "", "P", false );
	formItemMap [FormItemIonType::getName ()+string("S")]	= new FormItemIonType ( "", "S", false );
	formItemMap [FormItemIonType::getName ()+string("I")]	= new FormItemIonType ( "", "I", false );
	formItemMap [FormItemIonType::getName ()+string("N")]	= new FormItemIonType ( "", "N", false );
	formItemMap [FormItemIonType::getName ()+string("C")]	= new FormItemIonType ( "", "C", false );

}
void MSTagIonTypeForm::printHTML ( ostream& os )
{
	tableStart ( os, true, "", "1" );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				os << "AA<br />Composition" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 7 );
				os << "N-term" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 7 );
				os << "C-term" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 3 );
				os << "Satellite Sequence<br />(side-chain loss)" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 4 );
				os << "Neutral-loss<br />(from ions also selected)" << endl;
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
				formItemMap [FormItemIonType::getName()+string("a")]->printHTML ( os );
				os << "<br />a" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("b")]->printHTML ( os );
				os << "<br />b" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("B")]->printHTML ( os );
				os << "<br />b+H<sub>2</sub>O" << endl;
				os << "<br />RHK" << endl;
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
				formItemMap [FormItemIonType::getName()+string("n")]->printHTML ( os );
				os << "<br />";
				os << "-";
				os << getElementalFormulaHTML ( "NH3" ) << endl;
				os << "<br />";
				os << "RKQN" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("h")]->printHTML ( os );
				os << "<br />";
				os << "-";
				os << getElementalFormulaHTML ( "H2O" ) << endl;
				os << "<br />";
				os << "STED" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("P")]->printHTML ( os );
				os << "<br />";
				os << "-";
				os << getElementalFormulaHTML ( "H3PO4" ) << endl;
				os << "<br />";
				os << "STY";
				os << "+";
				os << getElementalFormulaHTML ( "PO4" ) << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				formItemMap [FormItemIonType::getName()+string("S")]->printHTML ( os );
				os << "<br />";
				os << "-";
				os << getElementalFormulaHTML ( "SOCH4" ) << endl;
				os << "<br />";
				os << "M+Ox" << endl;
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
}
void MSTagIonTypeForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];

		formItemMap [FormItemIonType::getName()+string("i")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("a")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("b")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("B")]->setValue ( p );
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
		formItemMap [FormItemIonType::getName()+string("d")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("v")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("w")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("n")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("h")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("P")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("S")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("I")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("N")]->setValue ( p );
		formItemMap [FormItemIonType::getName()+string("C")]->setValue ( p );
	}
}
