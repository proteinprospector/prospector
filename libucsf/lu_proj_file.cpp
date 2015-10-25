/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_proj_file.cpp                                              *
*                                                                             *
*  Created    : September 27th 2004                                           *
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
#ifdef BATCHTAG
#ifndef VIS_C
#include <stdexcept>
#endif
#include <lg_string.h>
#include <lgen_error.h>
#include <lgen_file.h>
#include <lu_getfil.h>
#include <lu_proj_file.h>
#include <lu_file_type.h>
#include <lu_xml.h>
#include <lu_xml_data.h>
#include <lu_param_list.h>
#include <lu_version.h>
#ifdef MYSQL_DATABASE
#include <ld_init.h>
#endif
using std::string;
using std::endl;
using std::runtime_error;
using namespace FileTypes;

string getProjectDir ( const ParameterList* params )
{
#ifdef MYSQL_DATABASE
	string searchKey = params->getStringValue ( "search_key" );
	if ( !searchKey.empty () ) {
		BatchJobItem* bji = MySQLPPSDDBase::instance ().getBatchJobByKey ( searchKey );
		return bji->getProjectPath ();
	}
	else {
#endif
		return params->getStringValue ( "project_filepath" );
#ifdef MYSQL_DATABASE
	}
#endif
}
static string getProjectFilePath ( const ParameterList* params )
{
#ifdef MYSQL_DATABASE
	string searchKey = params->getStringValue ( "search_key" );
	if ( !searchKey.empty () ) {
		static string sk;
		static BatchJobItem* bji;
		if ( searchKey != sk ) {
			bji = MySQLPPSDDBase::instance ().getBatchJobByKey ( searchKey );
			sk = searchKey;
		}
		return bji->getProjectFullPath ();
	}
	else {
#endif
		string projectPath = params->getStringValue ( "project_filepath" );
		string projectName = params->getStringValue ( "project_filename" );
		if ( !projectPath.empty () && !projectName.empty () )
			return projectPath + SLASH + projectName + ".xml";
		else
			return "";
#ifdef MYSQL_DATABASE
	}
#endif
}
string ProjectFile::centroidBaseDir = adjustPPOutputPath ( InfoParams::instance ().getCentroidDir () );
string ProjectFile::rawBaseDir = adjustPPOutputPath ( InfoParams::instance ().getRawDir () );
string ProjectFile::userBaseDir = InfoParams::instance ().getUserRepository ();

ProjectFile::ProjectFile ( const string& path )
{
	init ( path );
}
ProjectFile::ProjectFile ( const ParameterList* params )
{
	init ( getProjectFilePath ( params ) );
}
void ProjectFile::init ( const string& projectFilePath )
{
	if ( !genFileExists ( projectFilePath ) ) {
		ErrorHandler::genError ()->error ( "The project file " + projectFilePath + " does not exist.\n" );
	}
	else {
		projectVersion = getVersionFromPPXMLFile ( projectFilePath );
		char* info = getFileInfo ( projectFilePath, '\n', 1, false );
		toleranceUnits = XMLParser::getStringValue ( info, "msms_parent_mass_tolerance_units", "" );
		StringVector files = XMLParser::getStringVectorValue ( info, "file" );
		for ( StringVectorSizeType i = 0 ; i < files.size () ; i++ ) {
			centroidFiles.push_back ( XMLParser::getStringValue ( files [i], "centroid" ) );
			string r = XMLParser::getStringValue ( files [i], "raw", "" );
			if ( !r.empty () ) rawFiles.push_back ( r );
			else rawFiles.push_back ( "" );
			string nStr = XMLParser::getStringValue ( files [i], "num_msms_spectra", "" );
			int n;
			if ( !nStr.empty () ) n = atoi ( nStr.c_str () );
			else {
				ErrorHandler::genError ()->error ( "The project file does not contain the number of msms spectra.\n" );
			}
			numMSMSSpectra.push_back ( n );
			offset.push_back ( XMLParser::getDoubleValue ( files [i], "offset", "0.0" ) );
			tolerance.push_back ( XMLParser::getDoubleValue ( files [i], "tolerance", "0.0" ) );
		}
		delete [] info;
	}
}
bool ProjectFile::isUploadProject () const
{
	return centroidFiles [0][0] == '#';
}
string ProjectFile::getUploadDirectory () const
{
	return userBaseDir + SLASH + genDirectoryFromPath ( centroidFiles [0].substr ( 1 ) );
}
string ProjectFile::getCentroidPath ( StringVectorSizeType i ) const
{
	if ( i < centroidFiles.size () ) {
		if ( genIsFullPath ( centroidFiles [i] ) ) return genTranslateSlash ( centroidFiles [i] );
		else if ( centroidFiles [i][0] == '$' ) return genTranslateSlash ( centroidBaseDir + centroidFiles [i].substr ( 1 ) );
		else if ( centroidFiles [i][0] == '#' ) return genTranslateSlash ( userBaseDir + SLASH + centroidFiles [i].substr ( 1 ) );
		else return "";
	}
	else throw runtime_error ( "Centroid data not present." );
}
string ProjectFile::getAbsoluteRawPath ( StringVectorSizeType i ) const
{
	if ( i < rawFiles.size () ) {
		if ( genIsFullPath ( rawFiles [i] ) ) return rawFiles [i];
		else if ( rawFiles [i][0] == '$' ) return rawBaseDir + rawFiles [i].substr ( 1 );
		else if ( rawFiles [i][0] == '#' ) return userBaseDir + SLASH + rawFiles [i].substr ( 1 );
		else return "";
	}
	else return "";
}
string ProjectFile::updateProjectFile ( const ParameterList* params )
{
	string projectName;
	string projectPath;
	string projectFullPath;
	int calIndex;
#ifdef MYSQL_DATABASE		// Needs refining
	string searchKey = params->getStringValue ( "search_key" );
	BatchJobItem* bji = 0;
	if ( !searchKey.empty () ) {
		bji = MySQLPPSDDBase::instance ().getBatchJobByKey ( searchKey );
		projectFullPath = bji->getProjectFullPath ();
		projectName = bji->getProjectName ();
		projectPath = bji->getProjectPath ();
		calIndex = bji->getCalibrationIndex () + 1;
	}
	else {
#endif
		projectName = params->getStringValue ( "project_filename" );
		projectPath = params->getStringValue ( "project_filepath" );
		projectFullPath = projectPath + SLASH + projectName + ".xml";
		FileList fileList ( projectPath, projectName + ".cal.", ".xml", true );
		StringVector nList = fileList.getNameList ();
		calIndex = 1;
		int maxCal = 0;
		size_t start = projectName.length () + 5;	// 5 is length of .cal.
		for ( StringVectorSizeType j = 0 ; j < nList.size () ; j++ ) {
			string n = nList [j];
			maxCal = genMax ( maxCal, atoi ( n.substr ( start ).c_str () ) );
		}
		maxCal += 1;
#ifdef MYSQL_DATABASE
	}
#endif
	string calProjectName = projectName + ".cal." + gen_itoa ( calIndex );
	string calProjectFilename = calProjectName + ".xml";
	char* info = getFileInfo ( projectFullPath, '\n', 1, false );
	StringVector files = XMLParser::getStringVectorValue ( info, "file" );
	delete [] info;
	StringVector offset = params->getPQStringVectorValue ( "offsets" );
	GenOFStream os ( projectPath + SLASH + calProjectFilename );
	printXMLHeader ( os );
	printXMLVersion ( os );
	os << "<project>" << endl;
		ParameterList::printXML ( os, "project_name", calProjectName );
		ParameterList::printXML ( os, "msms_parent_mass_tolerance_units", params->getStringValue ( "msms_parent_mass_tolerance_units" ) );
		for ( StringVectorSizeType i = 0 ; i < files.size () ; i++ ) {
			os << "<file>";
				os << files [i];
				if ( !offset.empty () ) ParameterList::printXML ( os, "offset", offset [i] );
			os << "</file>" << endl;
		}
	os << "</project>" << endl;
#ifdef MYSQL_DATABASE
	if ( bji ) {
		MySQLPPSDDBase::instance ().updateCalIndex ( bji->getProjectID (), calIndex );
		MySQLPPSDDBase::instance ().submitProject ( bji->getUserID (), calProjectName, calProjectFilename, bji->getProjectRelativePath (), params->getStringValue ( "instrument_name" ) );
	}
#endif
	return calProjectName;
}
StringVector ProjectFile::getFractionNameList ( const ParameterList* params )
{
	string projectFilePath = getProjectFilePath ( params );
	PPExpatStringList ppesx ( "centroid_name" );
	try {
		ppesx.parseXMLFromFile ( projectFilePath );
	}
	catch ( runtime_error e ) {
		ErrorHandler::genError ()->error ( e );
	}
	StringVector names = ppesx.getNames ();
	if ( names.empty () ) {
		PPExpatStringList ppes2 ( "centroid" );
		try {
			ppes2.parseXMLFromFile ( projectFilePath );
		}
		catch ( runtime_error e ) {
			ErrorHandler::genError ()->error ( e );
		}
		StringVector names2 = ppes2.getNames ();
		for ( StringVectorSizeType i = 0 ; i < names2.size () ; i++ ) {
			names2 [i] = genFilenameFromPath ( names2 [i] );
			names2 [i] = names2 [i].substr ( 0, names2 [i].find_first_of ( "." ) );
		}
		return names2;
	}
	else
		return names;
}
StringVector ProjectFile::getCentroidPathList ( const ParameterList* params )
{
	string projectFilePath = getProjectFilePath ( params );
	PPExpatStringList ppesx ( "centroid" );
	try {
		ppesx.parseXMLFromFile ( projectFilePath );
	}
	catch ( runtime_error e ) {
		ErrorHandler::genError ()->error ( e );
	}
	StringVector names = ppesx.getNames ();
	for ( StringVectorSizeType i = 0 ; i < names.size () ; i++ ) {
		if ( names [i][0] == '$' )		names [i] = centroidBaseDir + names [i].substr ( 1 );
		else if ( names [i][0] == '#' ) names [i] = userBaseDir + SLASH + names [i].substr ( 1 );
	}
	return names;
}
IntVector ProjectFile::getNumMSSpectra ( const ParameterList* params )
{
	PPExpatStringList ppesx ( "num_ms_spectra" );
	try {
		ppesx.parseXMLFromFile ( getProjectFilePath ( params ) );
	}
	catch ( runtime_error e ) {
		ErrorHandler::genError ()->error ( e );
	}
	return ppesx.getNamesAsInt ();
}
StringVector ProjectFile::getRawTypes ( const ParameterList* params )
{
	string filePath = getProjectFilePath ( params );
	StringVector files;
	if ( filePath.empty () ) return files;	// No project file - return empty list
	try {
		char* info = getFileInfo ( filePath, '\n', 1, false );
		files = XMLParser::getStringVectorValue ( info, "file" );
		delete [] info;
	}
	catch ( runtime_error e ) {
		ErrorHandler::genError ()->error ( e );
	}
	StringVector sv;
	for ( StringVectorSizeType i = 0 ; i < files.size () ; i++ ) {
		string r = XMLParser::getStringValue ( files [i], "raw", "" );
		if ( !r.empty () ) {
			if ( isFileType ( r, WIFF ) )		sv.push_back ( WIFF );
			else if ( isFileType ( r, RAW ) )	sv.push_back ( RAW );
			else								sv.push_back ( T2D );
		}
		else sv.push_back ( "" );
	}
	return sv;
}
bool ProjectFile::getSpottingPlate ( const ParameterList* params )
{
	PPExpatStringList ppesx ( "spotting_plate" );
	try {
		ppesx.parseXMLFromFile ( getProjectFilePath ( params ) );
	}
	catch ( runtime_error e ) {
		ErrorHandler::genError ()->error ( e );
	}
	return !ppesx.getNames ().empty ();
}
DoubleVector ProjectFile::getOffsets ( const ParameterList* params, double defaultOffset )
{
	char* info = getFileInfo ( getProjectFilePath ( params ), '\n', 1, false );
	StringVector files = XMLParser::getStringVectorValue ( info, "file" );
	DoubleVector dv;
	for ( StringVectorSizeType i = 0 ; i < files.size () ; i++ ) {
		dv.push_back ( XMLParser::getDoubleValue ( files [i], "offset", defaultOffset ) );
	}
	delete [] info;
	return dv;
}
string getCentroidDataFilename ( const ParameterList* params, int fraction )
{
	ProjectFile pf ( params );
	return pf.getCentroidPath (fraction-1);
}
string getRawDataFilename ( const ParameterList* params, int fraction )
{
	ProjectFile pf ( params );
	return pf.getAbsoluteRawPath (fraction-1);
}
string getProjectVersion ( const ParameterList* params, int fraction )
{
	ProjectFile pf ( params );
	return pf.getProjectVersion ();
}
bool isCalibratedProject ( const string& project )
{
	string::size_type index = project.rfind ( ".cal." );
	return index != string::npos;
}
string getUncalibratedProject ( const string& project )
{
	string::size_type index = project.rfind ( ".cal." );
	if ( index == string::npos ) return project;
	else return project.substr ( 0, index );
}
#endif
