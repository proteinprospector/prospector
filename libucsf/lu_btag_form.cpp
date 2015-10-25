/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_btag_form.cpp                                              *
*                                                                             *
*  Created    : December 13th 2004                                            *
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
#ifdef MYSQL_DATABASE
#include <lu_cookie.h>
#include <lu_html.h>
#include <lu_table.h>
#include <lu_btag_form.h>
#include <lu_param_list.h>

using std::ostream;
using std::string;
using std::vector;
using std::endl;

BatchTagForm::BatchTagForm ( const VectorConstParameterListPtr& params, const string& formType ) :
	fvj (),
	headerForm ( params, "batchtag", "BatchTag" ),
	searchForm ( params, "batchtag", &fvj ),
	variableModsForm ( params, "msms_", &fvj ),
	extraUserModsForm ( params, &fvj, true ),
	massModificationForm ( params, &fvj ),
	crosslinkingForm ( params, &fvj, true ),
	matrixModeForm ( params, "msms" ),
	msmsToleranceForm ( &fvj, params, false, true ),
	formType ( formType )
{
}
void BatchTagForm::createItems ()
{
	formItemMap [FormItemComment::getName()]	= new FormItemComment ();

	formItemMap [FormItemMaxReportedHits::getName("msms_")]		= new FormItemMaxReportedHits ("msms_", "5", &fvj, false);
	formItemMap [FormItemMaxHits::getName ()]					= new FormItemMaxHits ();
	formItemMap [FormItemMSMSMaxModifications::getName ()]		= new FormItemMSMSMaxModifications ( &fvj );
	formItemMap [FormItemMSMSMaxPeptidePermutations::getName ()]= new FormItemMSMSMaxPeptidePermutations ( &fvj );

	formItemMap [FormItemSaveParameters::getName ()]	= new FormItemSaveParameters ();
	formItemMap [FormItemExpectationCalculationMethod::getName ()]	= new FormItemExpectationCalculationMethod ( "Linear Tail Fit" );

	formItemMap [FormItemUseInstrumentIonTypes::getName ()] = new FormItemUseInstrumentIonTypes ();

	formItemMap [FormItemCreateParams::getName ()]	= new FormItemCreateParams ();
	formItemMap [FormItemScriptFilename::getName ()]= new FormItemScriptFilename ();

	formItemMap [FormItemMSMSPkFilter::getName ()]			= new FormItemMSMSPkFilter ( false );
	formItemMap [FormItemMaxMSMSPeaks::getName ()]			= new FormItemMaxMSMSPeaks ( &fvj, true );
}
void BatchTagForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];
		formItemMap [FormItemComment::getName()]->setValue ( p );

		formItemMap [FormItemMaxReportedHits::getName("msms_")]->setValue ( p );
		formItemMap [FormItemMaxHits::getName ()]->setValue ( p );
		formItemMap [FormItemMSMSMaxModifications::getName ()]->setValue ( p );
		formItemMap [FormItemMSMSMaxPeptidePermutations::getName ()]->setValue ( p );

		formItemMap [FormItemExpectationCalculationMethod::getName ()]->setValue ( p );

		formItemMap [FormItemUseInstrumentIonTypes::getName ()]->setValue ( p );

		formItemMap [FormItemMSMSPkFilter::getName ()]->setValue ( p );
		formItemMap [FormItemMaxMSMSPeaks::getName ()]->setValue ( p );
	}
}
void BatchTagForm::printHTML ( ostream& os )
{
	printHTMLFormInit ( os );
	headerForm.printHTML ( os );
	basicJavascriptFunctions ( os );
	massModCrosslinkingJavascriptFunctions ( os );
	tableStart ( os, true );
		searchForm.printHTML ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				tableStart ( os );
					tableRowStart ( os );
						tableHeaderStart ( os, "", "center", true );
							printHTMLFORMSubmit ( os, "Start Search" );
						tableHeaderEnd ( os );
						tableHeaderStart ( os, "", "center", true );
							formItemMap [FormItemSaveParameters::getName ()]->printHTML ( os );
						tableHeaderEnd ( os );
					tableRowEnd ( os );
				tableEnd ( os );
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "left", true );
				formItemMap [FormItemExpectationCalculationMethod::getName ()]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [FormItemMSMSPkFilter::getName ()]->printHTML ( os );
				formItemMap [FormItemMaxMSMSPeaks::getName ()]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true );
				msmsToleranceForm.printHTML ( os );
			tableHeaderEnd ( os );
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
		tableRowEnd ( os );
		printDataForm ( os );
	tableEnd ( os );
	showHiddenItems ( os );
	printHTMLFORMEnd ( os );
	printHTMLFooter ( os, "2002" );
}
BatchTagFullForm::BatchTagFullForm ( const VectorConstParameterListPtr& params, const string& formType ) :
	BatchTagForm ( params, formType ),
	dataForm ( params, "batchtag" )
{
	create ( params );
	if ( formType == "Save Parameters" ) {
		formItemMap [FormItemCreateParams::getName()]->setValue ( true );
	}
}
void BatchTagFullForm::printHTMLFormInit ( ostream& os ) const
{
	init_html ( os, "Batch-Tag", "batchtagman.htm#batch_tag" );
	fvj.print ( os );
	printHTMLFORMStart ( os, "post", "msbatch", false, true );
}
void BatchTagFullForm::printDataForm ( ostream& os )
{
	dataForm.showHiddenItems ( os );
}
BatchTagWebForm::BatchTagWebForm ( const VectorConstParameterListPtr& params, bool setCookie ) :
	BatchTagForm ( params, "Search" ),
	uploadDataForm ( params ),
	user ( params [0]->getStringValue ( "user" ) ),
	setCookie ( setCookie )
{
	create ( params );
}
void BatchTagWebForm::printHTMLFormInit ( ostream& os ) const
{
	init_html ( os, "Batch-Tag Web", "batchtagman.htm#batch_tag" );
	if ( setCookie ) setUserCookie ( os, user );
	fvj.print ( os );
	printHTMLFORMStart ( os, "post", "msbatch", true, true );
}
void BatchTagWebForm::printDataForm ( ostream& os )
{
	uploadDataForm.printHTML ( os );
}
#endif
