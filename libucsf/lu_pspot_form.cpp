/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_pspot_form.cpp                                             *
*                                                                             *
*  Created    : December 15th 2004                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2004-2010) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lu_table.h>
#include <lu_pspot_form.h>
using std::ostream;
using std::string;
using std::vector;

static string SERVER_NAME		= "server_name";
static string USERNAME			= "username";
static string PASSWORD			= "password";

static string SPOT_SET_NAMES	= "spot_set_names";

static string RUN_NUMBER		= "run_number";
static string ALL_RUNS			= "all_runs";
static string WRITE_RAW_MS		= "write_raw_ms";
static string WRITE_RAW_MSMS	= "write_raw_msms";

static string INTENSITY_TYPE	= "intensity_type";
static string MINIMUM_AREA		= "minimum_area";
static string RETAIN_ISOTOPES	= "retain_isotopes";

static string CENTROID_DIR	= "centroid_dir";
static string RAW_DIR		= "raw_dir";

static string DIAGNOSTICS	= "diagnostics";

PeakSpotterForm::PeakSpotterForm ( const VectorConstParameterListPtr& params, const StringVector& spotSetNamesOptions ) :
	spot_set_names_options ( spotSetNamesOptions )
{
	create ( params );
}
void PeakSpotterForm::createItems ()
{
	formItemMap [SERVER_NAME]	= new FormItemText ( "Server Name", "", SERVER_NAME, 12, 20, "" );
	formItemMap [USERNAME]		= new FormItemText ( "Username", "", USERNAME, 12, 20, "" );
	formItemMap [PASSWORD]		= new FormItemText ( "Password", "", PASSWORD, 12, 20, "" );

	int spot_set_names_size = genMin ( 16, (int)spot_set_names_options.size () );
	formItemMap [SPOT_SET_NAMES]= new FormItemSelectMultiple ( "Spot Set<br />Names", "", SPOT_SET_NAMES, spot_set_names_options, StringVector (), spot_set_names_size );

	formItemMap [RUN_NUMBER]	= new FormItemText ( "Run Number", "", RUN_NUMBER, 3, 3, "" );
	formItemMap [ALL_RUNS]		= new FormItemCheckbox ( "All Runs", "", ALL_RUNS, true );
	formItemMap [WRITE_RAW_MS]	= new FormItemCheckbox ( "Write Raw MS", "", WRITE_RAW_MS, false );
	formItemMap [WRITE_RAW_MSMS]= new FormItemCheckbox ( "Write Raw MSMS", "", WRITE_RAW_MSMS, false );

	formItemMap [INTENSITY_TYPE]= new FormItemSelect ( "Intensity Type", "", INTENSITY_TYPE, intensity_type_options, "Height" );
	formItemMap [MINIMUM_AREA]	= new FormItemText ( "Minimum Peak Area", "", MINIMUM_AREA, 6, 6, "100.0" );
	formItemMap [RETAIN_ISOTOPES]= new FormItemCheckbox ( "Retain Isotopes", "", RETAIN_ISOTOPES, false );

	formItemMap [CENTROID_DIR]	= new FormItemText ( "Centroid Dir", "", CENTROID_DIR, 60, 100, "" );
	formItemMap [RAW_DIR]		= new FormItemText ( "Raw Dir", "", RAW_DIR, 60, 100, "" );
	formItemMap [DIAGNOSTICS]	= new FormItemCheckbox ( "Diagnostics", "", DIAGNOSTICS, false );
}
void PeakSpotterForm::setOptions ()
{
	intensity_type_options.push_back ( "Height" );
	intensity_type_options.push_back ( "Area" );
	intensity_type_options.push_back ( "Cluster Area" );
	intensity_type_options.push_back ( "Signal to Noise" );
}
void PeakSpotterForm::setValues ( const VectorConstParameterListPtr& params ) 
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];
		formItemMap [SERVER_NAME]->setValue ( p );
		formItemMap [USERNAME]->setValue ( p );
		formItemMap [PASSWORD]->setValue ( p );

		formItemMap [CENTROID_DIR]->setValue ( p );
		formItemMap [RAW_DIR]->setValue ( p );
	}
}
void PeakSpotterForm::printHTML ( ostream& os )
{
	printHTMLFORMStart ( os, "post", "peakSpotter" );
	tableStart ( os, true );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true );
				formItemMap [SPOT_SET_NAMES]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true );
				formItemMap [RUN_NUMBER]->printHTML ( os );
				formItemMap [ALL_RUNS]->printHTML ( os );
				formItemMap [WRITE_RAW_MS]->printHTML ( os );
				formItemMap [WRITE_RAW_MSMS]->printHTML ( os );
				formItemMap [DIAGNOSTICS]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true );
				formItemMap [INTENSITY_TYPE]->printHTML ( os );
				formItemMap [MINIMUM_AREA]->printHTML ( os );
				formItemMap [RETAIN_ISOTOPES]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				printHTMLFORMSubmit ( os, "Create Files" );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
	tableEnd ( os );
	showHiddenItems ( os );
	printHTMLFORMEnd ( os );
	printHTMLFooter ( os, "2002" );
}
