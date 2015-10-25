/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_patt_form.cpp                                              *
*                                                                             *
*  Created    : January 4th 2005                                              *
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
#include <lu_patt_form.h>
using std::string;
using std::vector;
using std::ostream;
using std::endl;

class FormItemRegularExpression : public FormItemText {
public:
	FormItemRegularExpression ( FormValidatingJavascript* fvj ) :
		FormItemText ( "Reg Expression (Use CAPS)", "patternman.htm#reg_expression", getName (), 40, 100, "[AK][AK]PV[IL]ED[IL]R", fvj->addMSPatternValidator ( getName (), "Reg Expression" ) ) {}
	static string getName () { return "regular_expression"; }
};

class FormItemMaxAASubstitutions : public FormItemText {
public:
	FormItemMaxAASubstitutions ( FormValidatingJavascript* fvj ) :
		FormItemText ( "Max. # of Mismatched AAs", "patternman.htm#mismatched_AA", getName (), 2, 2, "2", fvj->addPositiveIntegerValidator ( getName (), "Max. # of Mismatched AAs" ) ) {}
	static string getName () { return "max_aa_substitutions"; }
};

class FormItemPreSearchOnly : public FormItemCheckbox {
public:
	FormItemPreSearchOnly () :
		FormItemCheckbox ( "Pre Search Only", "patternman.htm#pre_search_only", getName (), false ) {}
	static string getName () { return "pre_search_only"; }
};

MSPatternForm::MSPatternForm ( const VectorConstParameterListPtr& params ) :
	fvj (),
	headerForm ( params, "mspattern", "MS-Pattern" ),
	searchForm ( params, "mspattern", &fvj )
{
	create ( params );
}
void MSPatternForm::createItems ()
{
	formItemMap [FormItemComment::getName()]			= new FormItemComment ();
	formItemMap [FormItemRegularExpression::getName ()]	= new FormItemRegularExpression (&fvj);
	formItemMap [FormItemMaxAASubstitutions::getName ()]= new FormItemMaxAASubstitutions (&fvj);
	formItemMap [FormItemMaxReportedHits::getName ()]	= new FormItemMaxReportedHits ("", "200", &fvj, false);
	formItemMap [FormItemPreSearchOnly::getName ()]		= new FormItemPreSearchOnly ();
	formItemMap [FormItemMaxHits::getName ()]	= new FormItemMaxHits ("100000");
}
void MSPatternForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];
		formItemMap [FormItemComment::getName()]->setValue ( p );
		formItemMap [FormItemRegularExpression::getName ()]->setValue ( p );
		formItemMap [FormItemMaxAASubstitutions::getName ()]->setValue ( p );
		formItemMap [FormItemMaxReportedHits::getName ()]->setValue ( p );
		formItemMap [FormItemPreSearchOnly::getName ()]->setValue ( p );
		formItemMap [FormItemMaxHits::getName ()]->setValue ( p );
	}
}
void MSPatternForm::printHTML ( ostream& os )
{
	init_html ( os, "MS-Pattern", "patternman.htm" );
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
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				formItemMap [FormItemRegularExpression::getName ()]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [FormItemMaxAASubstitutions::getName ()]->printHTML ( os );
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "left", true );
				formItemMap [FormItemMaxReportedHits::getName ()]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [FormItemPreSearchOnly::getName ()]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
	tableEnd ( os );
	showHiddenItems ( os );
	printHTMLFORMEnd ( os );
	printHTMLFooter ( os, "1996" );
}
