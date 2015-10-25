/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_dbst_form.cpp                                              *
*                                                                             *
*  Created    : January 3rd 2005                                              *
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
#include <lu_dbst_form.h>
using std::string;
using std::vector;
using std::ostream;
using std::endl;

class FormItemShowAAStatistics : public FormItemCheckbox {
public:
	FormItemShowAAStatistics () :
		FormItemCheckbox ( "Amino Acid Statistics", "dbstatman.htm#AA_Statistics", getName (), false ) {}
	static string getName () { return "show_aa_statistics"; }
};

class FormItemMinHistogramMass : public FormItemText {
public:
	FormItemMinHistogramMass ( FormValidatingJavascript* fvj ) :
		FormItemText ( "Histogram Range", "", getName (), 8, 8, "600.0", fvj->addPositiveFloatingPointValidator ( getName (), "Histogram Start Mass" ) ) {}
	static std::string getName () { return "min_histogram_mass"; }
};

class FormItemMaxHistogramMass : public FormItemText {
public:
	FormItemMaxHistogramMass ( FormValidatingJavascript* fvj ) :
		FormItemText ( "to", "", getName (), 8, 8, "15000.0", fvj->addPositiveFloatingPointValidator ( getName (), "Histogram End Mass" ) ) {}
	static std::string getName () { return "max_histogram_mass"; }
};

class FormItemBandwidth : public FormItemText {
public:
	FormItemBandwidth ( FormValidatingJavascript* fvj ) :
		FormItemText ( "Bandwidth", "", getName (), 8, 8, "1.0", fvj->addPositiveFloatingPointValidator ( getName (), "Bandwidth" ) ) {}
	static std::string getName () { return "density_bandwidth"; }
};

DBStatForm::DBStatForm ( const VectorConstParameterListPtr& params ) :
	fvj (),
	headerForm ( params, "dbstat", "DB-Stat" ),
	searchForm ( params, "dbstat", &fvj )
{
	create ( params );
}
void DBStatForm::createItems ()
{
	formItemMap [FormItemShowAAStatistics::getName ()]	= new FormItemShowAAStatistics ();
	formItemMap [FormItemMaxHistogramMass::getName ()]	= new FormItemMaxHistogramMass (&fvj);
	formItemMap [FormItemMinHistogramMass::getName ()]	= new FormItemMinHistogramMass (&fvj);
	formItemMap [FormItemBandwidth::getName ()]			= new FormItemBandwidth (&fvj);
}
void DBStatForm::printHTML ( ostream& os )
{
	init_html ( os, "DB-Stat", "dbstatman.htm" );
	fvj.print ( os );
	printHTMLFORMStart ( os, "post", "mssearch", false, true );
	headerForm.printHTML ( os );
	tableStart ( os, true );
		searchForm.printHTML ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				printHTMLFORMSubmit ( os, "Start" );
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "left", true );
				formItemMap [FormItemShowAAStatistics::getName ()]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [FormItemMinHistogramMass::getName ()]->printHTML ( os );
				formItemMap [FormItemMaxHistogramMass::getName ()]->printHTML ( os );
				os << "<br />" << endl;
				formItemMap [FormItemBandwidth::getName ()]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
	tableEnd ( os );
	showHiddenItems ( os );
	printHTMLFORMEnd ( os );
	printHTMLFooter ( os, "1996" );
}
void DBStatForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];

		formItemMap [FormItemShowAAStatistics::getName ()]->setValue ( p );
		formItemMap [FormItemMinHistogramMass::getName ()]->setValue ( p );
		formItemMap [FormItemMaxHistogramMass::getName ()]->setValue ( p );
		formItemMap [FormItemBandwidth::getName ()]->setValue ( p );
	}
}
