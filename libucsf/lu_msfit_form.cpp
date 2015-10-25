/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_msfit_form.cpp                                             *
*                                                                             *
*  Created    : December 9th 2004                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2004-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lu_html.h>
#include <lu_table.h>
#include <lu_msfit_form.h>
using std::ostream;
using std::string;
using std::vector;
using std::endl;

class FormItemReportHomologousProteins : public FormItemSelect {
	static const char* options [];
public:
	FormItemReportHomologousProteins () :
		FormItemSelect ( "Report Homologous Proteins", "fitman.htm#report_homologous_proteins", getName (), options, "Interesting" ) {}
	static string getName () { return "ms_report_homologous_proteins"; }
};
const char* FormItemReportHomologousProteins::options [] = { "All", "Interesting", "None", 0 };

class FormItemSortType : public FormItemSelect {
	static const char* options [];
public:
	FormItemSortType () :
		FormItemSelect ( "Sort By", "", getName (), options, "Score Sort" ) {}
	static string getName () { return "sort_type"; }
};
const char* FormItemSortType::options [] = { "Score Sort", "MW Sort", "pI Sort", 0 };

class FormItemMSMaxModifications : public FormItemText {
public:
	FormItemMSMaxModifications ( FormValidatingJavascript* fvj ) :
		FormItemText ( "Max Mods", "", getName (), 2, 2, "1", fvj->addPositiveIntegerValidator ( getName (), "Max Mods" ) ) {}
	static string getName () { return "ms_max_modifications"; }
};

class FormItemMinParentIonMatches : public FormItemText {
public:
	FormItemMinParentIonMatches ( FormValidatingJavascript* fvj ) :
		FormItemText ( "Min. # match with NO AA subs", "", getName (), 2, 2, "1", fvj->addPositiveIntegerValidator ( getName (), "Min. # match with NO AA subs" ) ) {}
	static string getName () { return "min_parent_ion_matches"; }
};


MSFitForm::MSFitForm ( const VectorConstParameterListPtr& params ) :
	fvj (),
	headerForm ( params, "msfit", "MS-Fit" ),
	searchForm ( params, "msfit", &fvj ),
	matrixModeForm ( params, "ms" ),
	msToleranceForm ( &fvj, params )
{
	create ( params );
}
void MSFitForm::createItems ()
{
	formItemMap [FormItemComment::getName()]		= new FormItemComment ();
	formItemMap [FormItemDisplayGraph::getName()]	= new FormItemDisplayGraph ( false );

	formItemMap [FormItemMaxReportedHits::getName ("ms_")]	= new FormItemMaxReportedHits ("ms_", "5", &fvj, false);
	formItemMap [FormItemSortType::getName ()]				= new FormItemSortType ();
	formItemMap [FormItemReportHomologousProteins::getName ()]= new FormItemReportHomologousProteins ();
	formItemMap [FormItemMinMatches::getName ()]			= new FormItemMinMatches (&fvj);

	formItemMap [FormItemMowseOn::getName ()]		= new FormItemMowseOn ();
	formItemMap [FormItemMowsePfactor::getName ()]	= new FormItemMowsePfactor (&fvj);
	formItemMap [FormItemParentMassConvert::getName ()]			= new FormItemParentMassConvert ( "MS" );
	formItemMap [FormItemParentContaminantMasses::getName()]	= new FormItemParentContaminantMasses ();

	formItemMap [FormItemModAA::getName ()]			= new FormItemModAA ( "Fit" );
	formItemMap [FormItemUserName::getName ("1")]	= new FormItemUserName ("1");
	formItemMap [FormItemUserName::getName ("2")]	= new FormItemUserName ("2");
	formItemMap [FormItemUserName::getName ("3")]	= new FormItemUserName ("3");
	formItemMap [FormItemUserName::getName ("4")]	= new FormItemUserName ("4");
	formItemMap [FormItemMSMaxModifications::getName ()]	= new FormItemMSMaxModifications (&fvj);
	formItemMap [FormItemMinParentIonMatches::getName ()]	= new FormItemMinParentIonMatches (&fvj);

	formItemMap [FormItemDetailedReport::getName ()]= new FormItemDetailedReport ();
}
void MSFitForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];
		formItemMap [FormItemComment::getName()]->setValue ( p );
		formItemMap [FormItemDisplayGraph::getName()]->setValue ( p );

		formItemMap [FormItemMaxReportedHits::getName ("ms_")]->setValue ( p );
		formItemMap [FormItemSortType::getName ()]->setValue ( p );
		formItemMap [FormItemReportHomologousProteins::getName ()]->setValue ( p );
		formItemMap [FormItemMinMatches::getName ()]->setValue ( p );
		formItemMap [FormItemMowseOn::getName ()]->setValue ( p );
		formItemMap [FormItemMowsePfactor::getName ()]->setValue ( p );
		formItemMap [FormItemParentMassConvert::getName ()]->setValue ( p );
		formItemMap [FormItemParentContaminantMasses::getName()]->setValue ( p );

		formItemMap [FormItemModAA::getName ()]->setValue ( p );
		formItemMap [FormItemUserName::getName ("1")]->setValue ( p );
		formItemMap [FormItemUserName::getName ("2")]->setValue ( p );
		formItemMap [FormItemUserName::getName ("3")]->setValue ( p );
		formItemMap [FormItemUserName::getName ("4")]->setValue ( p );
		formItemMap [FormItemMSMaxModifications::getName ()]->setValue ( p );
		formItemMap [FormItemMinParentIonMatches::getName ()]->setValue ( p );
	}
}
void MSFitForm::printHTMLMSFitItems ( ostream& os )
{
	tableRowStart ( os );
		tableHeaderStart ( os, "", "center", true );
			printHTMLFORMSubmit ( os, "Start Search" );
		tableHeaderEnd ( os );
		tableHeaderStart ( os, "", "left", true );
			formItemMap [FormItemComment::getName()]->printHTML ( os );
			os << "<br />" << endl;
			formItemMap [FormItemDisplayGraph::getName()]->printHTML ( os );
			headerForm.printHTML ( os );
		tableHeaderEnd ( os );
	tableRowEnd ( os );
	tableRowStart ( os );
		tableHeaderStart ( os, "", "left", true );
			formItemMap [FormItemMaxReportedHits::getName ("ms_")]->printHTML ( os );
			os << "<br />" << endl;
			formItemMap [FormItemSortType::getName ()]->printHTML ( os );
			os << "<br />" << endl;
			formItemMap [FormItemReportHomologousProteins::getName ()]->printHTML ( os );
			os << "<br />" << endl;
			formItemMap [FormItemMinMatches::getName ()]->printHTML ( os );
			os << "<br />" << endl;
			formItemMap [FormItemMowseOn::getName ()]->printHTML ( os );
			formItemMap [FormItemMowsePfactor::getName ()]->printHTML ( os );
			os << "<br />" << endl;
			formItemMap [FormItemParentMassConvert::getName ()]->printHTML ( os );
			os << "<br />" << endl;
			msToleranceForm.printHTML ( os );
			formItemMap [FormItemParentContaminantMasses::getName()]->printHTML ( os );
		tableHeaderEnd ( os );
		tableHeaderStart ( os, "", "left", true );
			tableStart ( os, true );
				tableRowStart ( os );
					tableHeaderStart ( os, "", "left", true );
						formItemMap [FormItemModAA::getName ()]->printHTML ( os );
						formItemMap [FormItemUserName::getName ("1")]->printHTML ( os );
						os << "<br />" << endl;
						formItemMap [FormItemUserName::getName ("2")]->printHTML ( os );
						os << "<br />" << endl;
						formItemMap [FormItemUserName::getName ("3")]->printHTML ( os );
						os << "<br />" << endl;
						formItemMap [FormItemUserName::getName ("4")]->printHTML ( os );
					tableHeaderEnd ( os );
				tableRowEnd ( os );
				tableRowStart ( os );
					tableHeaderStart ( os, "", "center", true );
						os << "OR<br />" << endl;
					tableHeaderEnd ( os );
				tableRowEnd ( os );
				tableRowStart ( os );
					tableHeaderStart ( os, "", "left", true );
						matrixModeForm.printHTML ( os );
						os << "<br />" << endl;
						formItemMap [FormItemMSMaxModifications::getName ()]->printHTML ( os );
						formItemMap [FormItemMinParentIonMatches::getName ()]->printHTML ( os );
					tableHeaderEnd ( os );
				tableRowEnd ( os );
			tableEnd ( os );
		tableHeaderEnd ( os );
	tableRowEnd ( os );
}
MSFitStandardForm::MSFitStandardForm ( const VectorConstParameterListPtr& params ) :
	MSFitForm ( params ),
	dataForm ( params, "msfit" )
{
	create ( params );
}
void MSFitStandardForm::printHTML ( ostream& os )
{
	init_html ( os, "MS-Fit", "fitman.htm" );
	fvj.print ( os );
	printHTMLFORMStart ( os, "post", "mssearch", false, true );
	tableStart ( os, true );
		searchForm.printHTML ( os );
		printHTMLMSFitItems ( os );
		dataForm.printHTML ( os );
	tableEnd ( os );
	showHiddenItems ( os );
	printHTMLFORMEnd ( os );
	printHTMLFooter ( os, "1995" );
}
MSFitUploadForm::MSFitUploadForm ( const VectorConstParameterListPtr& params ) :
	MSFitForm ( params ),
	dataForm ( params, "msfit" )
{
	create ( params );
}
void MSFitUploadForm::printHTML ( ostream& os )
{
	init_html ( os, "MS-Fit", "fitman.htm" );
	fvj.print ( os );
	printHTMLFORMStart ( os, "post", "mssearch", true, true );
	tableStart ( os, true );
		searchForm.printHTML ( os );
		printHTMLMSFitItems ( os );
		dataForm.printHTML ( os );
	tableEnd ( os );
	showHiddenItems ( os );
	printHTMLFORMEnd ( os );
	printHTMLFooter ( os, "1995" );
}
