/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_hom_form.cpp                                               *
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
*  Copyright (2005-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lu_html.h>
#include <lu_table.h>
#include <lu_hom_form.h>
#include <lu_mat_score.h>
using std::string;
using std::vector;
using std::ostream;
using std::endl;

class FormItemPossibleSequences : public FormItemTextAreaSimple {
public:
	FormItemPossibleSequences () :
		FormItemTextAreaSimple ( "Possible Sequences (Use CAPITALS)", "homologyman.htm#possible_sequences", getName (), 10, 50, StringVector () )
		{
			value.push_back ("1 [213]ENFAGVGV[I|L]DFES 6");
		}
	static string getName () { return "possible_sequences"; }
};

class FormItemHomMinMatches : public FormItemText {
public:
	FormItemHomMinMatches ( FormValidatingJavascript* fvj ) :
		FormItemText ( "Min. # peptides required to match", "homologyman.htm#min_peptides", getName (), 2, 2, "1", fvj->addPositiveNonZeroIntegerValidator ( getName (), "Min. # peptides required to match" ) ) {}
	static string getName () { return "min_matches"; }
};

class FormItemScoreMatrix : public FormItemSelect {
public:
	FormItemScoreMatrix () :
		FormItemSelect ( "Score Matrix", "homologyman.htm#score_matrix", getName (), getScoreMatrixNameList (), "BLOSUM62" ) {}
	static string getName () { return "score_matrix"; }
};

MSHomologyForm::MSHomologyForm ( const VectorConstParameterListPtr& params ) :
	fvj (),
	headerForm ( params, "mshomology", "MS-Homology" ),
	searchForm ( params, "mshomology", &fvj )
{
	create ( params );
}
void MSHomologyForm::createItems ()
{
	formItemMap [FormItemComment::getName()]			= new FormItemComment ();
	formItemMap [FormItemPossibleSequences::getName ()]	= new FormItemPossibleSequences ();
	formItemMap [FormItemMaxReportedHits::getName ()]	= new FormItemMaxReportedHits ("", "20", &fvj, false);
	formItemMap [FormItemHomMinMatches::getName ()]		= new FormItemHomMinMatches (&fvj);
	formItemMap [FormItemScoreMatrix::getName ()]		= new FormItemScoreMatrix ();
	formItemMap [FormItemFragmentMassesTolerance::getName ()]		= new FormItemFragmentMassesTolerance (&fvj,"0.5");
	formItemMap [FormItemFragmentMassesToleranceUnits::getName ()]	= new FormItemFragmentMassesToleranceUnits ("Da");
	formItemMap [FormItemMaxHits::getName ()]	= new FormItemMaxHits ("2000000");
	formItemMap [FormItemPreviousAA::getName()]			= new FormItemPreviousAA ( &fvj );
	formItemMap [FormItemNextAA::getName()]				= new FormItemNextAA ( &fvj );
}
void MSHomologyForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];
		formItemMap [FormItemComment::getName()]->setValue ( p );
		formItemMap [FormItemPossibleSequences::getName ()]->setValue ( p );
		formItemMap [FormItemMaxReportedHits::getName ()]->setValue ( p );
		formItemMap [FormItemHomMinMatches::getName ()]->setValue ( p );
		formItemMap [FormItemScoreMatrix::getName ()]->setValue ( p );
		formItemMap [FormItemFragmentMassesTolerance::getName ()]->setValue ( p );
		formItemMap [FormItemFragmentMassesToleranceUnits::getName ()]->setValue ( p );
		formItemMap [FormItemMaxHits::getName ()]->setValue ( p );
		formItemMap [FormItemPreviousAA::getName ()]->setValue ( p );
		formItemMap [FormItemNextAA::getName ()]->setValue ( p );
	}
}
void MSHomologyForm::printHTML ( ostream& os )
{
	init_html ( os, "MS-Homology", "homologyman.htm" );
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
				formItemMap [FormItemPreviousAA::getName()]->printHTML ( os );
				formItemMap [FormItemNextAA::getName()]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				formItemMap [FormItemPossibleSequences::getName ()]->printHTML ( os );
				os << "<br />" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "left", true );
				formItemMap [FormItemMaxReportedHits::getName ()]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [FormItemHomMinMatches::getName ()]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [FormItemScoreMatrix::getName ()]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [FormItemFragmentMassesTolerance::getName ()]->printHTML ( os );
				formItemMap [FormItemFragmentMassesToleranceUnits::getName ()]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
	tableEnd ( os );
	showHiddenItems ( os );
	printHTMLFORMEnd ( os );
	printHTMLFooter ( os, "1999" );
}
