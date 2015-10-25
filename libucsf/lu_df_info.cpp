/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_df_info.cpp                                                *
*                                                                             *
*  Created    : August 8th 2005                                               *
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
#include <lg_string.h>
#include <lgen_error.h>
#include <lgen_file.h>
#include <lu_df_info.h>
#include <lu_param_list.h>
#include <lu_proj_file.h>
#include <lu_file_type.h>
#include <lu_cgi_val.h>
#include <lu_spec_id.h>
#include <lu_tol.h>
#include <lu_getfil.h>
#include <lu_version.h>
using std::string;
using std::stable_sort;
using std::ostream;
using std::vector;
using std::istringstream;
using namespace FileTypes;

class SortDataByParentMass {
public:
	bool operator () ( const MSMSDataPoint& a, const MSMSDataPoint& b ) const
	{
		return a.getParentMPlusH () < b.getParentMPlusH ();
	}
};

DataSetInfo::DataSetInfo ( const string& dataFormat, const string& dataSource, const string& searchKey ) :
	dataFormat ( dataFormat ),
	dataSource ( dataSource ),
	searchKey ( searchKey )
{
}
DataSetInfo::DataSetInfo ( const ParameterList* params ) :
	dataFormat ( params->getStringValue	( "data_format" ) ),
	dataSource ( params->getStringValue	( "data_source", "Data Paste Area" ) ),
	searchKey ( params->getStringValue ( "search_key" ) )
{
}
DataSetInfo::DataSetInfo () :
	dr ( 0 )
{
}
DataSetInfo::~DataSetInfo () {}
MSDataSetInfo::MSDataSetInfo ( const ParameterList* params ) :
	DataSetInfo ( params )
{
	parseData ( params, false, SpecID ( params ), false );
}
MSDataSetInfo::~MSDataSetInfo () {}
MSMSPeakListDataSetInfo::MSMSPeakListDataSetInfo ( const string& filepath ) :
	DataSetInfo ()
{
	initDataReader ( filepath, true );
}
MSMSPeakListDataSetInfo::~MSMSPeakListDataSetInfo ()
{
	if ( dr ) delete dr;
}
bool MSMSPeakListDataSetInfo::writePeakList ( ostream& os, const SpecID& specID, const string& version )
{
	if ( !dr->writePeakList ( os, specID, version ) ) {		// Failed to find scan
		return false;
	}
	return true;		// Scan found
}
void MSMSPeakListDataSetInfo::getData ( MSMSDataPointVector& msmsDataPointList, const SpecID& specID, const string& version, UpdatingJavascriptMessage* ujm )
{
	dr->getData ( msmsDataPointList, specID, 0, version, ujm );
}
MSMSPeakListDataFilterInfo::MSMSPeakListDataFilterInfo ( const string& filepath ) :
	DataSetInfo ()
{
	initDataReader ( filepath, true );
}
MSMSPeakListDataFilterInfo::~MSMSPeakListDataFilterInfo ()
{
	if ( dr ) delete dr;
}
void MSMSPeakListDataFilterInfo::readPeakList ( MSMSDataPointVector& msmsDataPointList, int index )
{
	dr->getDataFromIndex ( msmsDataPointList, -1, SpectrumRange (), index, 0 );
}
void MSMSPeakListDataFilterInfo::writePeakList ( ostream& os )
{
	dr->writePeakList ( os );
}
void MSMSPeakListDataFilterInfo::rewind ()
{
	dr->rewind ();
}
bool MSMSPeakListDataFilterInfo::isEOF ()
{
	return dr->isEOF ();
}
MSMSQuanDataSetInfo::MSMSQuanDataSetInfo ( const string& filepath ) :
	DataSetInfo ()
{
	initDataReader ( filepath, true );
}
MSMSQuanDataSetInfo::~MSMSQuanDataSetInfo ()
{
	if ( dr ) delete dr;
}
bool MSMSQuanDataSetInfo::getMSMSQuanData ( XYData& xyData, const SpecID& specID, double startMass, double endMass, const string& version )
{
	if ( !dr->getMSMSQuanData ( xyData, specID, startMass, endMass, version ) ) {		// Failed to find scan
		return false;
	}
	return true;		// Scan found
}
MSMSDataSetInfo::MSMSDataSetInfo ( const ParameterList* params, bool noFilesOK ) :
	DataSetInfo ( params ),
	spectrumRange ( params ),
	charge ( setCharge ( params ) ),
	chargeRange ( setChargeRange ( params ) ),
	title ( params->getStringValue ( "title", "" ) ),
	index ( params->getStringValue ( "index", "-1" ) ),
	mOverZ ( params->getStringValue ( "m_over_z", "0.0" ) ),
	scanNumber ( params->getStringValue ( "scan_number", "-1" ) )
{
	int ind = params->getIntValue ( "index", -1 );
	double mz = params->getDoubleValue ( "m_over_z", 0.0 );
	int offset = params->getIntValue ( "offset", 0 );
	if ( !title.empty () )
		parseData ( params, title, offset );
	else if ( ind != -1 )
		parseData ( params, ind, offset );
	else if ( mz != 0.0 )
		parseData ( params, mz, offset );
	else if ( scanNumber != "-1" )
		parseData ( params, scanNumber, false, offset );
	else
		parseData ( params, noFilesOK, SpecID ( params ), true, offset );
	calibrate ( params );
	stable_sort ( dataPointList.begin (), dataPointList.end (), SortDataByParentMass () );	// sort causes occasional problems with g++ compiler
}
MSMSDataSetInfo::~MSMSDataSetInfo () {}

int MSMSDataSetInfo::setCharge ( const ParameterList* params ) const
{
	string zString = params->getStringValue ( "msms_precursor_charge", "Automatic" );
	if ( zString == "Automatic" ) {	// This checks if the charge info is carried in the spectrum number
		int s = params->getIntValue ( "spectrum_number" );
		return s / 10000;
	}
	else return atoi ( zString.c_str () );
}
IntVector MSMSDataSetInfo::setChargeRange ( const ParameterList* params )
{
	IntVector iv;
	string zRange = params->getStringValue ( "msms_precursor_charge_range" );
	istringstream istr ( zRange );
	string s;
	while ( istr >> s ) {
		int idx = s.find ( '-' );
		if ( idx == string::npos )
			iv.push_back ( atoi ( s.c_str () ) );	// Single integer
		else {
			int r1 = atoi ( s.substr ( 0, idx ).c_str () );
			int r2 = atoi ( s.substr ( idx+1 ).c_str () );
			for ( int ch = r1 ; ch <= r2 ; ch++ ) {
				iv.push_back ( ch );
			}
		}
	}
	return iv;
}
void DataSetInfo::parseData ( const ParameterList* params, bool noFilesOK, const SpecID& specID, bool msmsData, int off )
{
	if ( dataSource == "List of Files" ) {
#ifdef BATCHTAG
		ProjectFile projectFile ( params );
		string version = projectFile.getProjectVersion ();
		for ( StringVectorSizeType i = 0 ; i < projectFile.getNumFiles () ; i++ ) {
			int fraction = i+1;
			int f = specID.getFraction ();
			if ( f == -1 || f == fraction ) {
				getData ( projectFile.getCentroidPath ( i ), true, fraction, specID, version );
			}
		}
#endif
	}
	else if ( dataSource == "Data From File" ) {
		dataFilename = params->getStringValue ( "data_filename", "" );
		if ( !genFileExists ( dataFilename ) ) {
			ErrorHandler::genError ()->error ( "The data no longer exists.\n" );
		}
		string version;
		if ( dataFilename.find ( "msviewer" ) != string::npos ) {	// This is an msviewer file read the version number
			string pFilename = genDirectoryFromPath ( genDirectoryFromPath ( dataFilename ) ) + SLASH + "params.xml";
			version = getVersionFromPPXMLFile ( pFilename );
		}
		getData ( dataFilename, true, 1, specID, version, off );
	}
	else if ( dataSource == "Upload Data From File" ) {
		string f = params->getStringValue ( "upload_data_filepath" );
		getData ( f, true, 1, specID, "" );
		genUnlink ( f );
	}
	else {	// dataSource == "Data Paste Area"
		string data = params->getFileStringValue ( "data", "" );
		if ( gen_strstrip ( data ) == "" ) {
			if ( !noFilesOK ) ErrorHandler::genError ()->error ( "No data entered.\n" );
		}
		else {
			if ( msmsData ) {	// Currently the data is only written out for MS-Product
				PPTempFile ppTempFile ( "data", ".txt" );
				dataFilename = ppTempFile.getAdjustedPath ();
				GenOFStream ost ( ppTempFile.getAdjustedPath () );
				ost << data;
			}
			getData ( data, false, 1, specID , "" );
		}
	}
	if ( specID.isFilterSpectrum () && getNumDataSets () == 0 ) {
		ErrorHandler::genError ()->error ( "Spectrum not found in data file.\n" );
	}
}
void DataSetInfo::parseData ( const ParameterList* params, const string& title, int off )
{
	dataFilename = params->getStringValue ( "data_filename", "" );
	if ( !genFileExists ( dataFilename ) ) {
		ErrorHandler::genError ()->error ( "The data no longer exists.\n" );
	}
	getDataFromTitle ( dataFilename, true, 1, title, off );
}
void DataSetInfo::parseData ( const ParameterList* params, int index, int off )
{
	dataFilename = params->getStringValue ( "data_filename", "" );
	if ( !genFileExists ( dataFilename ) ) {
		ErrorHandler::genError ()->error ( "The data no longer exists.\n" );
	}
	getDataFromIndex ( dataFilename, true, 1, index, off );
}
void DataSetInfo::parseData ( const ParameterList* params, double mOverZ, int off )
{
	dataFilename = params->getStringValue ( "data_filename", "" );
	if ( !genFileExists ( dataFilename ) ) {
		ErrorHandler::genError ()->error ( "The data no longer exists.\n" );
	}
	getDataFromMOverZ ( dataFilename, true, 1, mOverZ, off );
}
void DataSetInfo::parseData ( const ParameterList* params, const string& scanNumber, bool dummy, int off )
{
	dataFilename = params->getStringValue ( "data_filename", "" );
	if ( !genFileExists ( dataFilename ) ) {
		ErrorHandler::genError ()->error ( "The data no longer exists.\n" );
	}
	getDataFromScanNumber ( dataFilename, true, 1, scanNumber, off );
}
void DataSetInfo::initDataReader ( const string& s, bool fileFlag )
{
	bool illegal = false;
	if ( dataFormat == "PP M/Z Intensity Charge" || dataFormat == "PP M/Z Charge" ) {
		if ( fileFlag ) {
			if ( !isTextFile ( s, 256 ) ) {
				ErrorHandler::genError ()->error ( "The data is in an illegal format.\n" );
			}
			if ( isMGFFile ( s ) ) {
				ErrorHandler::genError ()->error ( "The selected data format is '" + dataFormat + "'.\nThis appears to be an mgf file.\n" );
			}
			if ( isMS2File ( s ) ) {
				ErrorHandler::genError ()->error ( "The selected data format is '" + dataFormat + "'.\nThis appears to be an ms2 file.\n" );
			}
			if ( isAPLFile ( s ) ) {
				ErrorHandler::genError ()->error ( "The selected data format is '" + dataFormat + "'.\nThis appears to be an apl file.\n" );
			}
		}
		dr = new PPDataReader2 ( s, fileFlag, dataFormat );
	}
	else if ( dataFormat == MGF )
		dr = new MGFDataReader ( s, fileFlag );
	else if ( dataFormat == PKL )
		dr = new PKLDataReader ( s, fileFlag );
	else if ( dataFormat == DTA )
		dr = new DTADataReader ( s, fileFlag );
	else if ( dataFormat == MS2 )
		dr = new MS2DataReader ( s, fileFlag );
	else if ( dataFormat == APL )
		dr = new APLDataReader ( s, fileFlag );
	else if ( isFileType ( s, MGF ) )
		dr = new MGFDataReader ( s, fileFlag );
	else if ( isFileType ( s, MS2 ) )
		dr = new MS2DataReader ( s, fileFlag );
	else if ( isFileType ( s, APL ) )
		dr = new APLDataReader ( s, fileFlag );
#ifdef BATCHTAG
	else if ( isFileType ( s, MZXML ) )
		dr = new XMLDataReader ( s, fileFlag, MZXML );
	else if ( isFileType ( s, MZDATA ) )
		dr = new XMLDataReader ( s, fileFlag, MZDATA );
	else if ( isFileType ( s, XML ) ) {
		if ( isMZXMLFile ( s ) )
			dr = new XMLDataReader ( s, fileFlag, MZXML );
		else if ( isMZDataFile ( s ) )
			dr = new XMLDataReader ( s, fileFlag, MZDATA );
		else
			illegal = true;
	}
#endif
	else if ( isTextFile ( s, 256 ) ) {		// Check first 256 bytes to see if it is a text file
		if ( isPPSingleFile ( s ) )
			dr = new PPDataReader ( s, fileFlag );
		else if ( isMGFFile ( s ) )
			dr = new MGFDataReader ( s, fileFlag );
		else if ( isMS2File ( s ) )
			dr = new MS2DataReader ( s, fileFlag );
		else if ( isAPLFile ( s ) )
			dr = new APLDataReader ( s, fileFlag );
		else
			illegal = true;
	}
	else
		illegal = true;
	if ( illegal ) {
		ErrorHandler::genError ()->error ( "Illegal Data Format.\n" );
	}
}
void MSDataSetInfo::getData ( const string& s, bool fileFlag, int fraction, const SpecID& specID, const string& version, int off )
{
	initDataReader ( s, fileFlag );
	dr->getData ( dataPointList, fraction, specID );
	delete dr;
}
void MSMSDataSetInfo::getData ( const string& s, bool fileFlag, int fraction, const SpecID& specID, const string& version, int off )
{
	if ( spectrumRange.readFile ( fraction ) ) {
		initDataReader ( s, fileFlag );
		dr->setChargeRange ( chargeRange );
		dr->getData ( dataPointList, fraction, spectrumRange, specID, charge, version, off );
		delete dr;
	}
}
void MSMSDataSetInfo::getDataFromTitle ( const string& s, bool fileFlag, int fraction, const string& title, int off )
{
	if ( spectrumRange.readFile ( fraction ) ) {
		initDataReader ( s, fileFlag );
		dr->setChargeRange ( chargeRange );
		dr->getDataFromTitle ( dataPointList, fraction, spectrumRange, title, charge, off );
		delete dr;
	}
}
void MSMSDataSetInfo::getDataFromIndex ( const string& s, bool fileFlag, int fraction, int index, int off )
{
	if ( spectrumRange.readFile ( fraction ) ) {
		initDataReader ( s, fileFlag );
		dr->setChargeRange ( chargeRange );
		dr->getDataFromIndex ( dataPointList, fraction, spectrumRange, index, charge, off );
		delete dr;
	}
}
void MSMSDataSetInfo::getDataFromMOverZ ( const string& s, bool fileFlag, int fraction, double mOverZ, int off )
{
	if ( spectrumRange.readFile ( fraction ) ) {
		initDataReader ( s, fileFlag );
		dr->setChargeRange ( chargeRange );
		dr->getDataFromMOverZ ( dataPointList, fraction, spectrumRange, mOverZ, charge, off );
		delete dr;
	}
}
void MSMSDataSetInfo::getDataFromScanNumber ( const string& s, bool fileFlag, int fraction, const string& scanNumber, int off )
{
	if ( spectrumRange.readFile ( fraction ) ) {
		initDataReader ( s, fileFlag );
		dr->setChargeRange ( chargeRange );
		dr->getDataFromScanNumber ( dataPointList, fraction, spectrumRange, scanNumber, charge, off );
		delete dr;
	}
}
void MSMSDataSetInfo::calibrate ( const ParameterList* params )
{
	double defaultOffset = params->getDoubleValue ( "msms_parent_mass_systematic_error", 0.0 );
	double defaultTolValue = params->getDoubleValue ( "msms_parent_mass_tolerance" );
	string defaultTolUnits = params->getStringValue ( "msms_parent_mass_tolerance_units" );
	ToleranceInfo formTol ( "msms_parent_mass", params );
	vector <ToleranceInfo*> ti;
	DoubleVector offsets;
	if ( searchKey.empty () ) {
		ti.push_back ( new ToleranceInfo ( "msms_parent_mass", params ) );	// Default tolerance
		offsets.push_back ( defaultOffset );
	}
#ifdef BATCHTAG
	else {
		ProjectFile projectFile ( params );
		string tolUnits = projectFile.getToleranceUnits ();
		if ( tolUnits.empty () ) tolUnits = defaultTolUnits;
		for ( int i = 0 ; i < projectFile.getNumFiles () ; i++ ) {
			double fileOffset = projectFile.getOffset ( i );
			if ( fileOffset == 0.0 ) {				// No offset in the file, use the form value.
				ti.push_back ( new ToleranceInfo ( "msms_parent_mass", params ) );
				offsets.push_back ( defaultOffset );
			}
			else {									// Offset is read from file
				ti.push_back ( new ToleranceInfo ( "msms_parent_mass", fileOffset, tolUnits ) );
				offsets.push_back ( fileOffset );
			}
		}
	}
#endif
	for ( MSMSDataPointVectorSizeType j = 0 ; j < dataPointList.size () ; j++ ) {
		int index = dataPointList [j].getFraction ()-1;
		if ( index >= ti.size () ) index = 0;
		dataPointList [j].calibrate ( ti [index]->getTolerance (), offsets [index] );
		dataPointList [j].setTolerance ( formTol.getTolerance () );
	}
	for ( std::vector <ToleranceInfo*>::size_type k = 0 ; k < ti.size () ; k++ ) delete ti [k];
}
void MSMSDataSetInfo::putProductCGI ( ostream& os, int i ) const
{
	if ( !dataFormat.empty () ) printCGIString ( os, "data_format", dataFormat );
	if ( dataSource == "List of Files" ) {
		printCGIString ( os, "data_source", "List of Files" );
		printCGIString ( os, "search_key", searchKey );
		dataPointList [i].putCGI ( os );
	}
	else if ( dataSource == "Data From File" ) {		// Viewer
		printCGIString ( os, "data_source", "Data From File" );
		printCGIString ( os, "data_filename", dataFilename );
		if ( !title.empty () )
			printCGIString ( os, "title", title );
		else if ( index != "-1" )
			printCGIString ( os, "index", index );
		else if ( mOverZ != "0.0" )
			printCGIString ( os, "m_over_z", mOverZ );
		else if ( scanNumber != "-1" )
			printCGIString ( os, "scan_number", scanNumber );
		else
			dataPointList [i].putCGI ( os );
	}
	else {
		printCGIString ( os, "data_source", "Data From File" );
		printCGIString ( os, "data_filename", dataFilename );
		dataPointList [i].putCGI ( os );
	}
}
