/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_data_form.cpp                                              *
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
*  Copyright (2004-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lg_string.h>
#include <lu_table.h>
#include <lu_data_form.h>
#include <lu_inst.h>
#include <lu_param_list.h>
using std::ostream;
using std::string;
using std::vector;
using std::endl;

class FormItemUser : public FormItemText {
public:
	FormItemUser ( const string& defaultVal );
	static string getName () { return "user"; }
};

class FormItemSearchKey : public FormItemText {
public:
	FormItemSearchKey ();
	static string getName () { return "search_key"; }
};

class FormItemDataSource : public FormItemSelect {
	static const char* options [];
public:
	FormItemDataSource ( int i = 0 ) :
		FormItemSelect ( "", "allman.htm#data_sets", getName (), options, options [i] ) {}
	static string getName () { return "data_source"; }
};
const char* FormItemDataSource::options [] = { "List of Files", "Upload Data From File", "Data Paste Area", 0 };

class FormItemProjectName : public FormItemText {
public:
	FormItemProjectName ( const string& name = "" ) :
		FormItemText ( "Project Name", "", getName (), 40, 256, name ) {}
	static string getName () { return "project_name"; }
};

class FormItemFraction : public FormItemText {
public:
	FormItemFraction () :
		FormItemText ( "Fraction", "", getName (), 3, 10, "1" ) {}
	static string getName () { return "fraction"; }
};

class FormItemSpotNumber : public FormItemText {
public:
	FormItemSpotNumber () :
		FormItemText ( "Spot", "", getName (), 8, 20, "1" ) {}
	static string getName () { return "spot_number"; }
};

class FormItemRun : public FormItemText {
public:
	FormItemRun () :
		FormItemText ( "Run", "", getName (), 3, 10, "1" ) {}
	static string getName () { return "run"; }
};

class FormItemSpectrumNumber : public FormItemText {
public:
	FormItemSpectrumNumber () :
		FormItemText ( "Spec No", "", getName (), 4, 10, "1" ) {}
	static string getName () { return "spectrum_number"; }
};
FormItemInstrumentName::FormItemInstrumentName ( int num ) :
	FormItemSelect ( "Instrument", "allman.htm#instrument", getName ( num ), InstrumentList::instance ().getNames (), "ESI-Q-TOF" )
{
}
FormItemUser::FormItemUser ( const string& defaultVal ) :
	FormItemText ( "User", "", getName (), 20, 60, defaultVal )
{
}
FormItemSearchKey::FormItemSearchKey () :
	FormItemText ( "Search Key", "", getName (), 20, 60, "" )
{
}
DataForm::DataForm ( const VectorConstParameterListPtr& params, const string& searchName ) :
	searchName ( searchName )
{
	if ( searchName == "mstag" || searchName == "msfit" || searchName == "msbridge" || searchName == "msproduct" ) filterValue = true;
	else filterValue = false;
	create ( params );
}
void DataForm::createItems ()
{
	formItemMap [FormItemDataSource::getName ()]= new FormItemDataSource ();
	formItemMap [FormItemInstrumentName::getName ()]= new FormItemInstrumentName ();
	if ( searchName == "mstag" || searchName == "msbridge" ) {
		formItemMap [FormItemSearchKey::getName ()]	= new FormItemSearchKey ();
	}
	else {
		formItemMap [FormItemProjectName::getName ()]	= new FormItemProjectName ();
		formItemMap [FormItemUser::getName ()]			= new FormItemUser ( "" );
	}
	if ( filterValue ) {
		formItemMap [FormItemFraction::getName ()]		= new FormItemFraction ();
		formItemMap [FormItemSpotNumber::getName ()]	= new FormItemSpotNumber ();
		formItemMap [FormItemRun::getName ()]			= new FormItemRun ();
		formItemMap [FormItemSpectrumNumber::getName ()]= new FormItemSpectrumNumber ();
	}
}
void DataForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];
		//formItemMap [FormItemDataSource::getName ()]->setValue ( p );	// An uploaded file should be a list of files when run subsequently
		formItemMap [FormItemInstrumentName::getName ()]->setValue ( p );
		if ( searchName == "mstag" || searchName == "msbridge" ) {
			formItemMap [FormItemSearchKey::getName ()]->setValue ( p );
		}
		else {
			formItemMap [FormItemProjectName::getName ()]->setValue ( p );
			formItemMap [FormItemUser::getName ()]->setValue ( p );
		}
		if ( filterValue ) {
			formItemMap [FormItemFraction::getName ()]->setValue ( p );
			formItemMap [FormItemSpotNumber::getName ()]->setValue ( p );
			formItemMap [FormItemRun::getName ()]->setValue ( p );
			formItemMap [FormItemSpectrumNumber::getName ()]->setValue ( p );
		}
	}
}
void DataForm::printHTML ( ostream& os )
{
	tableRowStart ( os );
		tableHeaderStart ( os, "", "left", true, 2 );
			formItemMap [FormItemInstrumentName::getName ()]->printHTML ( os );
			if ( searchName != "batchtag" && searchName != "mstag" && searchName != "msbridge" ) {
				formItemMap [FormItemProjectName::getName ()]->printHTML ( os );
			}
			os << "<br />" << endl;
			if ( searchName == "msfit" || searchName == "msbridge" || searchName == "mstag" ) {
				formItemMap [FormItemFraction::getName ()]->printHTML ( os );
				formItemMap [FormItemSpotNumber::getName ()]->printHTML ( os );
				formItemMap [FormItemRun::getName ()]->printHTML ( os );
				if ( searchName == "mstag" ) {
					formItemMap [FormItemSpectrumNumber::getName ()]->printHTML ( os );
				}
			}
			showHiddenItems ( os );
		tableHeaderEnd ( os );
	tableRowEnd ( os );
}
DataFormFromViewer::DataFormFromViewer ( const VectorConstParameterListPtr& params )
{
	create ( params );
}
void DataFormFromViewer::createItems ()
{
	formItemMap ["data_source"]	= new FormItemText ( "", "", "data_source", 12, 20, "Data From File" );
	formItemMap ["data_filename"]	= new FormItemText ( "", "", "data_filename", 12, 20, "" );
	formItemMap ["title"]		= new FormItemText ( "", "", "title", 12, 20, "" );
	formItemMap ["index"]		= new FormItemText ( "", "", "index", 12, 20, "" );
	formItemMap ["m_over_z"]	= new FormItemText ( "", "", "m_over_z", 12, 20, "" );
	formItemMap ["scan_number"]	= new FormItemText ( "", "", "scan_number", 12, 20, "" );
	formItemMap [FormItemInstrumentName::getName ()]= new FormItemInstrumentName ();
	formItemMap [FormItemFraction::getName ()]		= new FormItemFraction ();
	formItemMap [FormItemSpotNumber::getName ()]	= new FormItemSpotNumber ();
	formItemMap [FormItemRun::getName ()]			= new FormItemRun ();
	formItemMap [FormItemSpectrumNumber::getName ()]= new FormItemSpectrumNumber ();
}
void DataFormFromViewer::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];
		formItemMap ["data_filename"]->setValue ( p );
		formItemMap [FormItemInstrumentName::getName ()]->setValue ( p );
		formItemMap [FormItemFraction::getName ()]->setValue ( p );
		formItemMap [FormItemSpotNumber::getName ()]->setValue ( p );
		formItemMap ["title"]->setValue ( p );
		formItemMap ["index"]->setValue ( p );
		formItemMap ["m_over_z"]->setValue ( p );
		formItemMap ["scan_number"]->setValue ( p );
		formItemMap [FormItemRun::getName ()]->setValue ( p );
		formItemMap [FormItemSpectrumNumber::getName ()]->setValue ( p );
	}
}
void DataFormFromViewer::printHTML ( ostream& os )
{
}
class FormItemUploadData : public FormItemFile {
public:
	FormItemUploadData ( int i = 0 ) :
		FormItemFile ( "", "", getName ( i ), 100, 256, "" ) {}
	static string getName ( int i = 0 ) { return i == 0 ? "upload_data" : "upload_data" + gen_itoa ( i ); }
};
class FormItemMultiUploadData : public FormItemMultiFile {
public:
	FormItemMultiUploadData ( int i = 0 ) :
		FormItemMultiFile ( getName ( i ) ) {}
	static string getName ( int i = 0 ) { return i == 0 ? "upload_data" : "upload_data" + gen_itoa ( i ); }
};
UploadDataForm::UploadDataForm ( const VectorConstParameterListPtr& params )
{
	create ( params );
}
void UploadDataForm::createItems ()
{
	formItemMap [FormItemDataSource::getName ()]	= new FormItemDataSource ( 1 );
	formItemMap [FormItemInstrumentName::getName ()]= new FormItemInstrumentName ();
	formItemMap [FormItemUser::getName ()]			= new FormItemUser ("");
	formItemMap [FormItemUploadData::getName ()]	= new FormItemUploadData ();
}
void UploadDataForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];
		formItemMap [FormItemInstrumentName::getName ()]->setValue ( p );
		formItemMap [FormItemUser::getName ()]->setValue ( p );
	}
}
void UploadDataForm::printHTML ( ostream& os )
{
	tableRowStart ( os );
		tableHeaderStart ( os, "", "left", true, 2 );
			formItemMap [FormItemInstrumentName::getName ()]->printHTML ( os );
			os << "<br />" << endl;
			os << "<div class=\"form_large_label\">Upload Data From File</div>" << endl;
			formItemMap [FormItemUploadData::getName ()]->printHTML ( os );
			showHiddenItems ( os );
		tableHeaderEnd ( os );
	tableRowEnd ( os );
}

class FormItemDataFormat : public FormItemSelect {
	static const char* options [];
	static const char* optionsMSMS [];
public:
	FormItemDataFormat ( const string& s, int i = 0 ) :
		FormItemSelect ( "Data Format", "allman.htm#data_format", getName (), s == "mstag" || s == "msproduct" ? optionsMSMS : options, s == "mstag" ? optionsMSMS [i] : options [i] ) {}
	static string getName () { return "data_format"; }
};
const char* FormItemDataFormat::options [] = { "PP M/Z Charge", "PP M/Z Intensity Charge", 0 };
const char* FormItemDataFormat::optionsMSMS [] = { "PP M/Z Charge", "PP M/Z Intensity Charge", "mgf", "pkl", "dta", "ms2", "apl", 0 };

MSUploadDataForm::MSUploadDataForm ( const VectorConstParameterListPtr& params, const string& searchName ) :
	searchName ( searchName )
{
	create ( params );
}
void MSUploadDataForm::createItems ()
{
	formItemMap [FormItemDataFormat::getName ()]	= new FormItemDataFormat ( searchName, 0 );
	formItemMap [FormItemDataSource::getName ()]	= new FormItemDataSource ( 1 );
	formItemMap [FormItemInstrumentName::getName ()]= new FormItemInstrumentName ();
	formItemMap [FormItemUploadData::getName ()]	= new FormItemUploadData ();
}
void MSUploadDataForm::setValues ( const VectorConstParameterListPtr& params )
{
}
void MSUploadDataForm::printHTML ( ostream& os )
{
	tableRowStart ( os );
		tableHeaderStart ( os, "", "left", true, 2 );
			formItemMap [FormItemInstrumentName::getName ()]->printHTML ( os );
			formItemMap [FormItemDataFormat::getName ()]->printHTML ( os );
			os << "<br />" << endl;
			os << "<div class=\"form_large_label\">Upload Data From File</div>" << endl;
			formItemMap [FormItemUploadData::getName ()]->printHTML ( os );
			showHiddenItems ( os );
		tableHeaderEnd ( os );
	tableRowEnd ( os );
}

class FormItemData : public FormItemTextArea {
	static const char* msfitDefaultData [];
	static const char* msBridgeDefaultData [];
	static const char* msTagDefaultData [];
	static const char* msSeqDefaultData [];
	static const char* defaultData [];
	static const char** getData ( const std::string& searchName )
	{
		if ( searchName == "msfit" )		return msfitDefaultData;
		if ( searchName == "msnonspecific" )return msfitDefaultData;
		if ( searchName == "msbridge" )		return msBridgeDefaultData;
		if ( searchName == "mstag" )		return msTagDefaultData;
		if ( searchName == "msseq" )		return msSeqDefaultData;
		else return defaultData;
	}
public:
	FormItemData ( const string& searchNameValue ) :
	  FormItemTextArea ( "", "", getName (), 8, searchNameValue == "msbridge" ? 100 : 40, getData ( searchNameValue ) ) {}
	static string getName () { return "data"; }
};
const char* FormItemData::msfitDefaultData [] = {
"842.5100","856.5220","864.4733","870.5317","940.4754","943.4885",
"959.4934","970.4308","975.4785","1045.5580","1048.5716","1063.5712",
"1064.5892","1098.6185","1147.5876","1163.5996","1178.6280","1179.6014",
"1187.6316","1193.5461","1211.6607","1248.5664","1280.5561","1289.7670",
"1314.7019","1328.6521","1332.7121","1360.6820","1406.6617","1447.7010",
"1459.7311","1475.7471","1508.8107","1576.7986","1624.7649","1699.9255",
"1721.9134","1767.9147","1776.8961","1783.9077","1794.8293","1799.9017",
"1816.9798","1859.8805","2088.9872","2211.1046","2240.1851","2256.2412",
"2284.2079","2299.2019","2808.4450","3156.6352",0
};
const char* FormItemData::msBridgeDefaultData [] = {
"1347.5309","2222.8939","2435.2433","3344.3744","3573.6221","3620.6058",
"3683.7694","3698.5582","3725.6378","3848.7644",0
};
const char* FormItemData::msTagDefaultData [] = {
"1784.8375","60.0260","70.0468","84.0637","86.0800","87.0435","102.0332",
"112.0623","120.0685","129.0865","136.0690","157.0753","158.0879","159.0808",
"175.1034","185.0819","202.0876","217.1499","233.1270","234.1454","240.1426",
"245.1380","250.1414","258.1417","261.1298","262.1239","272.1949","278.1367",
"300.1896","320.1692","358.2232","365.1832","391.1774","392.1980","417.1967",
"434.2274","462.2164","479.2096","486.3239","543.3594","544.3416","545.3406",
"598.3135","599.4026","609.2665","626.2820","658.4249","859.5782","900.5853",
"902.5632","1013.6333","1030.6641","1084.7102","1159.7124","1306.7697",
"1403.8115",0
};
const char* FormItemData::msSeqDefaultData [] = {
"1181.6","279.1","407.2",0
};
const char* FormItemData::defaultData [] = { 0 };

class FormItemMSPeakExclusion: public FormItemText {
public:
	FormItemMSPeakExclusion () : FormItemText ( "", "", getName (), 1, 1, "0" ) {}
	static std::string getName () { return "ms_peak_exclusion"; }
};
class FormItemMSMassExclusion: public FormItemText {
public:
	FormItemMSMassExclusion () : FormItemText ( "", "", getName (), 1, 1, "0" ) {}
	static std::string getName () { return "ms_mass_exclusion"; }
};
class FormItemMSMatrixExclusion: public FormItemText {
public:
	FormItemMSMatrixExclusion () : FormItemText ( "", "", getName (), 1, 1, "0" ) {}
	static std::string getName () { return "ms_matrix_exclusion"; }
};
class FormItemMSMSPeakExclusion: public FormItemText {
public:
	FormItemMSMSPeakExclusion () : FormItemText ( "", "", getName (), 1, 1, "0" ) {}
	static std::string getName () { return "msms_peak_exclusion"; }
};
class FormItemMSMSMassExclusion: public FormItemText {
public:
	FormItemMSMSMassExclusion () : FormItemText ( "", "", getName (), 1, 1, "0" ) {}
	static std::string getName () { return "msms_mass_exclusion"; }
};
class FormItemMSMSMatrixExclusion: public FormItemText {
public:
	FormItemMSMSMatrixExclusion () : FormItemText ( "", "", getName (), 1, 1, "0" ) {}
	static std::string getName () { return "msms_matrix_exclusion"; }
};
class FormItemMSMSJoinPeaks: public FormItemText {
public:
	FormItemMSMSJoinPeaks () : FormItemText ( "", "", getName (), 1, 1, "0" ) {}
	static std::string getName () { return "msms_join_peaks"; }
};
class FormItemMSMSDeisotope: public FormItemText {
public:
	FormItemMSMSDeisotope () : FormItemText ( "", "", getName (), 1, 1, "0" ) {}
	static std::string getName () { return "msms_deisotope"; }
};
PasteAreaDataForm::PasteAreaDataForm ( const VectorConstParameterListPtr& params, const string& searchNameValue ) :
	searchNameValue ( searchNameValue )
{
	create ( params );
}
void PasteAreaDataForm::createItems ()
{
	formItemMap [FormItemDataFormat::getName ()]	= new FormItemDataFormat ( searchNameValue, 0 );
	formItemMap [FormItemDataSource::getName ()]	= new FormItemDataSource ( 2 );
	formItemMap [FormItemInstrumentName::getName ()]= new FormItemInstrumentName ();
	formItemMap [FormItemData::getName ()]			= new FormItemData ( searchNameValue );

	formItemMap [FormItemMSPeakExclusion::getName ()]		= new FormItemMSPeakExclusion ();
	formItemMap [FormItemMSMassExclusion::getName ()]		= new FormItemMSMassExclusion ();
	formItemMap [FormItemMSMatrixExclusion::getName ()]		= new FormItemMSMatrixExclusion ();
}
void PasteAreaDataForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];
		formItemMap [FormItemDataFormat::getName ()]->setValue ( p );
		formItemMap [FormItemDataSource::getName ()]->setValue ( p );
		formItemMap [FormItemInstrumentName::getName ()]->setValue ( p );
		formItemMap [FormItemData::getName ()]->setValue ( p );

		formItemMap [FormItemMSPeakExclusion::getName ()]->setValue ( p );
		formItemMap [FormItemMSMassExclusion::getName ()]->setValue ( p );
		formItemMap [FormItemMSMatrixExclusion::getName ()]->setValue ( p );
	}
}
void PasteAreaDataForm::printHTML ( ostream& os )
{
	tableRowStart ( os );
		tableHeaderStart ( os, "", "left", true, 2 );
			formItemMap [FormItemInstrumentName::getName ()]->printHTML ( os );
			formItemMap [FormItemDataFormat::getName ()]->printHTML ( os );
			os << "<br />" << endl;
			os << "<div class=\"form_large_label\">Data Paste Area</div>" << endl;
			if ( searchNameValue == "msbridge" ) printHTMLFORMLabel ( os, "(Information on Additional Filters)", "bridgeman.htm#add_filter" );
			formItemMap [FormItemData::getName ()]->printHTML ( os );
			showHiddenItems ( os );
		tableHeaderEnd ( os );
	tableRowEnd ( os );
}
