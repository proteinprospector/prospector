/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_viewer_par.cpp                                             *
*                                                                             *
*  Created    : March 8th 2011                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2011-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lg_stdlib.h>
#include <lg_string.h>
#include <lgen_error.h>
#include <lgen_file.h>
#include <lgen_uncompress.h>
#include <lu_file_type.h>
#include <lu_param_list.h>
#include <lu_viewer_form.h>
#include <lu_viewer_par.h>
#include <lu_getfil.h>
#include <lu_prog.h>
#include <lu_html.h>
#include <lu_repository.h>
#include <lu_usermod.h>

using std::ostringstream;
using std::ostream;
using std::string;
using std::getline;
using std::istream;
using std::endl;
using std::cout;
using namespace FileTypes;

MSViewerFilter::MSViewerFilter ()
{
}
void MSViewerFilter::initParameters ( const ParameterList* params )
{
	filterColumn.clear ();
	filterType.clear ();
	filterValue.clear ();
	for ( int i = 1 ; i <= MSViewerForm::getMaxFilterLevels () ; i++ ) {
		string num = gen_itoa (i);
		int colNumber = getColumnNumber ( params->getStringValue ( "column_num_filter_" + num ) );
		if ( colNumber == 0 ) break;
		filterColumn.push_back ( colNumber );
		filterType.push_back ( params->getStringValue ( "filter_type_" + num, "Equals" ) );
		filterValue.push_back ( params->getPQStringVectorValue ( "filter_value_" + num ) );
	}
}
int MSViewerFilter::getColumnNumber ( const string& name )
{
	if ( name == "Undefined" )
		return 0;
	else
		return atoi ( name.c_str () );
}
void MSViewerFilter::printHTML ( ostream& os ) const
{
	for ( int i = 0 ; i < filterColumn.size () ; i++ ) {
		int num = i+1;
		os << "column_num_filter_" << num << "=" << filterColumn [i] << '&';
		os << "filter_type_" << num << "=" << filterType [i] << '&';

		os << "filter_value_" << num << "=";
		for ( int j = 0 ; j < filterValue [i].size () ; j++ ) {
			os << filterValue [i][j];
			os << "%0D%0A";
		}
		os << "&";
	
	}
}
MSViewerSortOrder::MSViewerSortOrder ()
{
}
void MSViewerSortOrder::initParameters ( const ParameterList* params )
{
	sortLevel.clear ();
	sortOrderType.clear ();
	sortOrderDirection.clear ();
	for ( int i = 1 ; i <= MSViewerForm::getMaxSortLevels () ; i++ ) {
		string num = gen_itoa (i);
		int colNumber = getColumnNumber ( params->getStringValue ( "column_num_sort_level_" + num ) );
		if ( colNumber == 0 ) break;
		sortLevel.push_back ( colNumber );
		sortOrderType.push_back ( params->getStringValue ( "sort_order_type_" + num, "Alphabetic" ) );
		sortOrderDirection.push_back ( params->getStringValue ( "sort_order_direction_" + num, "Ascending" ) );
	}
}
int MSViewerSortOrder::getColumnNumber ( const string& name )
{
	if ( name == "Undefined" )
		return 0;
	else
		return atoi ( name.c_str () );
}
void MSViewerSortOrder::printHTML ( ostream& os ) const
{
	for ( int i = 0 ; i < sortLevel.size () ; i++ ) {
		int num = i+1;
		os << "column_num_sort_level_" << num << "=" << sortLevel [i] << '&';
		os << "sort_order_type_" << num << "=" << sortOrderType [i] << '&';
		os << "sort_order_direction_" << num << "=" << sortOrderDirection [i] << '&';
	}
}
MSViewerRemoveReplicate::MSViewerRemoveReplicate ()
{
}
void MSViewerRemoveReplicate::initParameters ( const ParameterList* params )
{
	replicateTest.clear ();
	for ( int i = 1 ; i <= MSViewerForm::getMaxRemoveReplicateLevels () ; i++ ) {
		string num = gen_itoa (i);
		int colNumber = getColumnNumber ( params->getStringValue ( "column_num_replicate_test_" + num ) );
		if ( colNumber == 0 ) break;
		replicateTest.push_back ( colNumber );
	}
}
int MSViewerRemoveReplicate::getColumnNumber ( const string& name )
{
	if ( name == "Undefined" )
		return 0;
	else
		return atoi ( name.c_str () );
}
void MSViewerRemoveReplicate::printHTML ( ostream& os ) const
{
	for ( int i = 0 ; i < replicateTest.size () ; i++ ) {
		int num = i+1;
		os << "column_num_replicate_test_" << num << "=" << replicateTest [i] << '&';
	}
}
MSViewerParameters::MSViewerParameters ( const ParameterList* params ) :
	MSProgramParameters	( params ),
	saveSettings ( params->getBoolValue ( "save_params" ) ),
	searchKey ( getParamsSearchKey ( params ) ),
	deleteFlag ( params->getBoolValue ( "delete" ) ),
	viewerOutputType ( params->getStringValue ( "viewer_output_type" ) ),
	page ( params->getIntValue ( "page", 1 ) ),
	rowsPerPage ( params->getStringValue ( "rows_per_page", "20" ) ),
	commandLine ( params->getStringValue ( "cl_results_filepath" ) != "" ),
	numPepXML ( 0 ),
	numMSF ( 0 )
{
	init_html ( cout, "MS-Viewer Report" );
	if ( !deleteFlag ) {
		if ( !searchKey.empty () ) {	// A saved dataset
			string f1 = getViewerRepositoryPath ( searchKey ) + SLASH + "params.xml";
			if ( !genFileExists ( f1 ) ) {
				ErrorHandler::genError ()->error ( "No data set was found for the search key '" + searchKey + "'." );
			}
			pList = new ParameterList ( f1, false, false, false, false );
			const_cast <ParameterList*>(pList)->appendParameters ( params );
			pListFlag = false;
			ProgramLink::setParams ( pList );
			init ( pList );
		}
		else {
			pList = params;
			pListFlag = true;
			init ( params );
		}
	}
}
MSViewerParameters::~MSViewerParameters ()
{
	if ( !deleteFlag ) {
		delete parentTolerance;
		delete fragmentTolerance;
		delete aaInitInfo;
	}
}
string MSViewerParameters::getParamsSearchKey ( const ParameterList* params ) const
{
	string key = params->getStringValue ( "search_key" );
	key = gen_strstrip ( key );
	key = gen_strstripchar ( key, '.' );	// Strip out all the . characters
	return key;
}
string MSViewerParameters::getKey () const
{
	if ( resultsFpath.empty () ) return "";
	else {
		string::size_type idx2 = resultsFpath.find_last_of ( "\\/" );
		string::size_type idx1 = resultsFpath.find_last_of ( "\\/", idx2-1 );
		return resultsFpath.substr ( idx1+1, idx2-idx1-1 );
	}
}
string MSViewerParameters::getViewerRepositoryContainerPath ()
{
	string path = InfoParams::instance ().getStringValue ( "viewer_repository" );
	if ( path.empty () ) {
		string bDir = genCurrentWorkingDirectory ();
		if ( genFilenameFromPath ( bDir ) == "cgi-bin" ) {
			bDir = genDirectoryFromPath ( bDir );
		}
		path = bDir + SLASH + string ( "results" ) + SLASH + string ( "msviewer" );
	}
	if ( path.find_last_of ( "\\/" ) != path.length () - 1 ) path += SLASH;
	return path;
}
string MSViewerParameters::getViewerRepositoryContainerPath ( const string& searchKey )
{
	string path = getViewerRepositoryContainerPath ();
	path += searchKey [0];
	path += SLASH;
	path += searchKey [1];
	return path;
}
string MSViewerParameters::getViewerRepositoryPath ( const string& searchKey )
{
	return getViewerRepositoryContainerPath ( searchKey ) + SLASH + searchKey;
}
void MSViewerParameters::init ( const ParameterList* params )
{
	UpdatingJavascriptMessage ujm;
	resultsFileFormat = params->getStringValue ( "results_file_format" );
	instrumentFilter = params->getStringValue ( "instrument_filter" );
	probabilityLimit = params->getDoubleValue ( "probability_limit", 0.05 );
	initResults ( params, ujm );
	initPeakList ( params, ujm );
	initParams ( params );
	ujm.deletePreviousMessage ( cout );
}
void MSViewerParameters::initResults ( const ParameterList* params, UpdatingJavascriptMessage& ujm )
{
	resultsFpath = params->getStringValue ( "results_filepath" );
	if ( resultsFpath.empty () ) {
		ujm.writeMessage ( cout, "Starting process of results file." );
		bool multi = false;
		string uploadResultsFname;
		string uploadResultsFpath;
		uploadResultsFpath = params->getStringValue ( "cl_results_filepath" );
		string resultsFname;
		if ( !uploadResultsFpath.empty () ) {
			uploadResultsFname = genFilenameFromPath ( uploadResultsFpath );
			resultsFname = uploadResultsFname;
		}
		else {
			uploadResultsFname = params->getStringValue ( "upload_temp_results_filename" );
			uploadResultsFpath = params->getStringValue ( "upload_temp_results_filepath" );
			resultsFname = genFilenameFromPath ( uploadResultsFname );		// IE gives the full path whereas Mozilla give the filename (what we want)
		}
		if ( uploadResultsFpath.empty () ) {
			ErrorHandler::genError ()->error ( "The results file must be specified." );
		}
		string path = genPreprocessFile ( uploadResultsFpath );
		string origPath;
		if ( genIsDirectory ( path ) ) {	// File extracted from an archive
			origPath = path;
			FileList fList ( path, "", "", false );
			StringVector n1 = fList.getNameList ();
			if ( n1.size () != 3 ) {
				for ( StringVectorSizeType i = 0 ; i < n1.size () ; i++ ) {
					string n = n1 [i];
					if ( n == "." || n == ".." ) continue;
					if ( isFileType ( n, MSF ) ) {
						numMSF++;
					}
					else if ( isPepXMLFile ( path + SLASH + n ) ) {
						numPepXML++;
					}
					else {
						ErrorHandler::genError ()->error ( "Only a single results file or a set of msf or pepXML files may be uploaded." );
					}
				}
				if ( numPepXML && numMSF ) {
					ErrorHandler::genError ()->error ( "A mixture of msf and pepXML files cannot be processed in a single upload." );
				}
				multi = true;
				resName = genShortFilenameFromPath ( resultsFname );
			}
			else {
				string f;
				for ( StringVectorSizeType i = 0 ; i < n1.size () ; i++ ) {
					string n = n1 [i];
					if ( n != "." && n != ".." ) {
						f = n;
					}
				}
				path += SLASH + f;
				resultsFname = f;
			}
		}
		else {
			string suffix = genSuffixFromPath ( resultsFname );
			const char* csuffix = suffix.c_str ();
			if ( !genStrcasecmp ( csuffix, "zip" ) || !genStrcasecmp ( csuffix, "7z" ) || !genStrcasecmp ( csuffix, "rar" ) || !genStrcasecmp ( csuffix, "gz" ) || !genStrcasecmp ( csuffix, "z" ) || !genStrcasecmp ( csuffix, "bz2" ) || !genStrcasecmp ( csuffix, "cmn" ) ) {
				resultsFname = resultsFname.substr ( 0, resultsFname.length () - strlen ( csuffix ) - 1 );
			}
		}
		PPTempFile pptf ( "", "" );					// Move the directory to the temporary directory
		tempFileFullPath = pptf.getFullPathDir () + SLASH + genToLower ( genRandomString ( 10 ) );
		if ( multi ) {
			resultsFpath = tempFileFullPath;
		}
		else {
			genCreateNewDirectory ( tempFileFullPath );
			resultsFpath = tempFileFullPath + SLASH + resultsFname;
		}

		genRename ( path, resultsFpath );
		if ( !origPath.empty () ) {
			genUnlinkDirectory ( origPath );
		}
		scriptConversion ( ujm );
		ujm.writeMessage ( cout, "Ending process of results file." );
	}
	else {
		if ( genIsFullPath ( resultsFpath ) ) {
			tempFileFullPath = resultsFpath.substr ( 0, resultsFpath.find_last_of ( "\\/" ) );
		}
	}
}
void MSViewerParameters::scriptConversion ( UpdatingJavascriptMessage& ujm )
{
	bool prospector = getProspector ();
	bool prospectorXL = getProspectorXL ();
	bool other = getOther ();
	bool automatic = getAuto ();
	if ( !prospector && !prospectorXL && !other && !automatic ) {
		string script = getScript ();
		if ( script != "N/A" ) {
			ujm.writeMessage ( cout, "Starting script conversion of results file." );
			bool maxQuant = getMaxQuant ();
			string command = getSystemCall ( script );
			string resultsFpath2 = resultsFpath + "_1";
			command += " ";
			command += "\"" + resultsFpath +"\"";
			command += " ";
			command += "\"" + resultsFpath2 +"\"";
			if ( maxQuant ) {
				command += " ";
				command += "\"" + gen_ftoa ( probabilityLimit, "%.4f" ) +"\"";
				if ( !instrumentFilter.empty () ) {
					command += " ";
					command += "\"" + instrumentFilter +"\"";
				}
			}
			int ret = genSystem ( command, "", true );
			if ( ret != 0 ) {
				ErrorHandler::genError ()->error ( "Problems converting the results file." );
			}
			ujm.writeMessage ( cout, "Ending script conversion of results file." );
			genUnlink ( resultsFpath );
			genRename ( resultsFpath2, resultsFpath );
		}
	}
}
void MSViewerParameters::initPeakList ( const ParameterList* params, UpdatingJavascriptMessage& ujm )
{
	peakListFpath = params->getStringValue ( "peak_list_filepath" );
	if ( peakListFpath.empty () ) {
		string uploadPeakListFpath;
		string uploadPeakListFname;
		uploadPeakListFpath = params->getStringValue ( "cl_peak_list_filepath" );
		if ( !uploadPeakListFpath.empty () ) {
			uploadPeakListFname = genFilenameFromPath ( uploadPeakListFpath );
		}
		else {
			uploadPeakListFname = params->getStringValue ( "upload_temp_peak_list_filename" );
			uploadPeakListFpath = params->getStringValue ( "upload_temp_peak_list_filepath" );
			uploadPeakListFname = genFilenameFromPath ( uploadPeakListFname );	// IE gives the full path whereas Mozilla give the filename (what we want)
		}
		if ( uploadPeakListFname.empty () ) {
			return;
		}
		ujm.writeMessage ( cout, "Starting preprocessing of results file." );
		string uploadName = genPreprocessFile ( uploadPeakListFpath );
		ujm.writeMessage ( cout, "Ending preprocessing of results file." );
		if ( uploadName == uploadPeakListFpath || isCompressedUpload ( uploadName ) ) {		// If the file hasn't been preprocessed or it has just been uncompressed it must be a single file
			string shortName = genShortFilenameFromPath ( uploadName );
			string newDir = genDirectoryFromPath ( uploadName ) + SLASH + shortName;
			if ( newDir == uploadName ) {
				newDir += "_1";
			}
			bool ret = genCreateDirectory ( newDir );
			string newUploadName;
			if ( isCompressedUpload ( uploadName ) )
				newUploadName = newDir + SLASH + genShortFilenameFromPath ( uploadPeakListFname );
			else
				newUploadName = newDir + SLASH + uploadPeakListFname;
			genRename ( uploadName, newUploadName );
			uploadName = newDir;
		}
		peakListFpath = tempFileFullPath + SLASH + genShortFilenameFromPath ( uploadPeakListFname );
		genRename ( uploadName, peakListFpath );
		processAPLFiles ( ujm );		// This can be left out but best to process these before PPProject
	}
}
void MSViewerParameters::initParams ( const ParameterList* params )
{
	parentTolerance		= new ToleranceInfo ( "msms_parent_mass", params );
	fragmentTolerance	= new ToleranceInfo ( "fragment_masses", params );

	instrumentName = params->getStringValue ( "instrument_name" );

	fractionColumnNumber	= getColumnNumber ( params->getStringValue ( "column_num_fraction" ) );
	spectrumIdentifier		= params->getStringValue ( "spectrum_identifier" );
	scanIDColumnNumber		= getColumnNumber ( params->getStringValue ( "column_num_scan_id" ) );

	peptideColumnNumber		= getColumnNumber ( params->getStringValue ( "column_num_peptide" ) );
	zColumnNumber			= getColumnNumber ( params->getStringValue ( "column_num_z" ) );
	modificationReporting	= params->getStringValue ( "modifications" );
	constantModColumnNumber	= getColumnNumber ( params->getStringValue ( "column_num_constant_mod" ) );
	variableModColumnNumber	= getColumnNumber ( params->getStringValue ( "column_num_variable_mod" ) );
	allModColumnNumber		= getColumnNumber ( params->getStringValue ( "column_num_all_mod" ) );

	separator		= getColumnSeparator ( params->getStringValue ( "column_separator" ) );
	numTitleLines	= params->getIntValue ( "num_title_lines" );
	numHeaderLines	= params->getIntValue ( "num_header_lines" );
	aaInitInfo = new AAInitInfo ( params );

	removeColumn = params->getStringVectorValue ( "remove_column" );

	msvSortOrder.initParameters ( params );
	msvFilter.initParameters ( params );
	msvRemoveReplicate.initParameters ( params );

	if ( genIsFullPath ( resultsFpath ) ) {
		const_cast <ParameterList*> (params)->removeName ( "upload_temp_peak_list_filename" );
		const_cast <ParameterList*> (params)->removeName ( "upload_temp_peak_list_filepath" );
		const_cast <ParameterList*> (params)->removeName ( "upload_temp_results_filename" );
		const_cast <ParameterList*> (params)->removeName ( "upload_temp_results_filepath" );
		const_cast <ParameterList*> (params)->addOrReplaceName ( "results_filepath", resultsFpath );
		if ( !peakListFpath.empty () ) const_cast <ParameterList*> (params)->addOrReplaceName ( "peak_list_filepath", peakListFpath );
		params->XMLParameterFile ( tempFileFullPath + SLASH + "params.xml" );
	}
	linkSearchType = params->getStringValue ( "link_search_type" );
	bridgeFormula = params->getStringValue ( "bridge_composition" );
	ExtraUserMods::instance ().addUserMods2 ( params, ExtraUserModsForm::getNumUserMods () );
}
int MSViewerParameters::getColumnNumber ( const string& name )
{
	if ( name == "Undefined" )
		return 0;
	else
		return atoi ( name.c_str () );
}
string MSViewerParameters::getColumnSeparator ( const string& name )
{
	if ( name == "Tab Delimited" )		return "\t";
	else if ( name == "CSV" )			return ",";
	else								return name;
}
void MSViewerParameters::deleteDataSet () const
{
	if ( searchKey.length () == 10 ) {
		string viewerDir = getViewerRepositoryPath ( searchKey );
		genUnlinkDirectory ( viewerDir );
	}
}
void MSViewerParameters::saveDataSet () const
{
	string searchKey = getKey ();
	if ( !searchKey.empty () ) {
		string originalDir = genDirectoryFromPath ( resultsFpath );
		string resFile = genFilenameFromPath ( resultsFpath );
		string pkListFile = genFilenameFromPath ( peakListFpath );
		const_cast <ParameterList*> (pList)->addOrReplaceName ( "results_filepath", resFile );
		const_cast <ParameterList*> (pList)->addOrReplaceName ( "peak_list_filepath", pkListFile );
		const_cast <ParameterList*> (pList)->removeName ( "save_params" );
		if ( commandLine ) {
			const_cast <ParameterList*> (pList)->removeName ( "cl_results_filepath" );
			const_cast <ParameterList*> (pList)->removeName ( "cl_peak_list_filepath" );
		}
		pList->XMLParameterFile ( originalDir + SLASH + "params.xml" );

		string outputPath = getViewerRepositoryContainerPath ( searchKey );
#ifndef VIS_C
		outputPath += SLASH;
		outputPath += searchKey;
#endif
		genCreateDirectoryPath ( outputPath );
		genRename ( originalDir, outputPath, true );
	}
}
string MSViewerParameters::getScript () const
{
	return ViewerResultsFileConverter::instance ().getScript ( resultsFileFormat );
}
void MSViewerParameters::getScriptParameters ( int& numTitleLines, int& numHeaderLines, string& delimiter ) const
{
	ViewerResultsFileConverter::instance ().getScriptParameters ( resultsFileFormat, numTitleLines, numHeaderLines, delimiter );
}
void MSViewerParameters::getScriptParameters2 ( string& spectrumIdentifier, string& fraction, string& scanID, string& peptide, string& charge, string& modifications ) const
{
	ViewerResultsFileConverter::instance ().getScriptParameters2 ( resultsFileFormat, spectrumIdentifier, fraction, scanID, peptide, charge, modifications );
}
void MSViewerParameters::setPeakListFpath ( const string& f )
{
	const_cast <ParameterList*> (pList)->addOrReplaceName ( "peak_list_filepath", f );
	peakListFpath = f;
}
void MSViewerParameters::setResultsFpath ( const string& f )
{
	const_cast <ParameterList*> (pList)->addOrReplaceName ( "results_filepath", f );
	resultsFpath = f;
}
void MSViewerParameters::setResultsFileFormat ( const string& f )
{
	const_cast <ParameterList*> (pList)->addOrReplaceName ( "results_file_format", f );
	const_cast <ParameterList*> (pList)->removeName ( "instrument_filter" );
	const_cast <ParameterList*> (pList)->removeName ( "probability_limit" );
	resultsFileFormat = f;
}
void MSViewerParameters::setSpectrumIdentifier ( const string& f )
{
	const_cast <ParameterList*> (pList)->addOrReplaceName ( "spectrum_identifier", f );
	spectrumIdentifier = f;
}
BoolDeque MSViewerParameters::getRemoveColumn ( int numCols ) const
{
	BoolDeque bd;
	for ( int i = 0 ; i < numCols ; i++ ) {
		bd.push_back ( false );
	}
	for ( StringVectorSizeType j = 0 ; j < removeColumn.size () ; j++ ) {
		int num = atoi ( removeColumn [j].c_str () );
		bd [num-1] = true;
	}
	return bd;
}
void MSViewerParameters::setScanIDColumnNumber ( int n )
{
	const_cast <ParameterList*> (pList)->addOrReplaceName ( "column_num_scan_id", n );
	scanIDColumnNumber = n;
}
void MSViewerParameters::setPeptideColumnNumber ( int n )
{
	const_cast <ParameterList*> (pList)->addOrReplaceName ( "column_num_peptide", n );
	peptideColumnNumber = n;
}
void MSViewerParameters::setVariableModColumnNumber ( int n )
{
	const_cast <ParameterList*> (pList)->addOrReplaceName ( "column_num_variable_mod", n );
	variableModColumnNumber = n;
}
void MSViewerParameters::setZColumnNumber ( int n )
{
	const_cast <ParameterList*> (pList)->addOrReplaceName ( "column_num_z", n );
	zColumnNumber = n;
}
void MSViewerParameters::setFractionColumnNumber ( int n )
{
	const_cast <ParameterList*> (pList)->addOrReplaceName ( "column_num_fraction", n );
	fractionColumnNumber = n;
}
void MSViewerParameters::setSeparator ( const string& s )
{
	const_cast <ParameterList*> (pList)->addOrReplaceName ( "column_separator", s );
	separator = s;
}
void MSViewerParameters::setNumHeaderLines ( int n )
{
	const_cast <ParameterList*> (pList)->addOrReplaceName ( "num_header_lines", n );
	numHeaderLines = n;
}
void MSViewerParameters::setPage ( int n )
{
	const_cast <ParameterList*> (pList)->addOrReplaceName ( "page", n );
	page = n;
}
void MSViewerParameters::setConstMod ( const StringVector& s )
{
	const_cast <ParameterList*> (pList)->addOrReplaceName ( "const_mod", s );
}
void MSViewerParameters::processAPLFiles ( UpdatingJavascriptMessage& ujm )
{
	FileList fList ( peakListFpath, "", "", false );
	StringVector sv = fList.getNameList ();
	for ( StringVectorSizeType i = 0 ; i < sv.size () ; i++ ) {
		string df2 = sv [i];
		string path2 =  peakListFpath + SLASH + df2;
		if ( !genIsDirectory ( path2 ) && isFileType ( path2, APL ) ) {
			ujm.writeMessage ( cout, "Processing file " + df2 + "." );
			PPProject::parseAPLFile ( peakListFpath, df2 );
		}
	}
}
