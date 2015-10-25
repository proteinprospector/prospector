/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_expec_par.cpp                                              *
*                                                                             *
*  Created    : July 20st 2006                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2006-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef BATCHTAG
#include <lgen_file.h>
#ifdef MYSQL_DATABASE
#include <ld_init.h>
#endif
#include <lu_getfil.h>
#include <lu_proj_file.h>
#include <lu_xml.h>
#include <lu_param_list.h>
#include <lu_version.h>
using std::string;

namespace {
bool versionCheck ( const string& verEFile, const string& version )
{		// Previously this was done at 5.4.2 and 5.6.0
	if ( Version::isOlderVersion ( verEFile, "5.14.3" ) && !Version::isOlderVersion ( version, "5.14.3" ) ) {
		return false;		// A new evalue calculation method introduced
	}
	return true;
}
}
ParameterList* getExpectationParams ( const ParameterList* searchParamList, string& outputName, bool force )
{
	char* info = getFileInfo ( MsparamsDir::instance ().getParamPath ( "expectation.xml" ), '\n', 1, false );
	double parentTolSP = searchParamList->getDoubleValue ( "msms_parent_mass_tolerance" );				// precursor tol on form
	string parentTolUnitsSP = searchParamList->getStringValue ( "msms_parent_mass_tolerance_units" );
	int mCleavesSP = searchParamList->getIntValue ( "missed_cleavages" );

// Get parameters for a .exp search
	string xparams = XMLParser::getStringValue ( info, "parameters" );

// Get parameters to copy from the normal search
	StringVector cparams = XMLParser::getStringVectorValue ( info, "copy_parameter" );
	delete [] info;

// Make a parameter list from these 2 sources
	ParameterList* expParamList = new ParameterList ( xparams, false, false, false );
	double parentTolEP = expParamList->getDoubleValue ( "msms_parent_mass_tolerance" );
	string parentTolUnitsEP = expParamList->getStringValue ( "msms_parent_mass_tolerance_units" );
	if ( parentTolUnitsEP == parentTolUnitsSP && parentTolSP > parentTolEP ) {
		cparams.push_back ( "msms_parent_mass_tolerance" );
		cparams.push_back ( "msms_parent_mass_tolerance_units" );
	}
	double mCleavesEP = expParamList->getIntValue ( "missed_cleavages" );
	if ( mCleavesSP > mCleavesEP ) {
		cparams.push_back ( "missed_cleavages" );
	}
	expParamList->copyName ( searchParamList, "data_source" );
	for ( StringVectorSizeType ii = 0 ; ii < cparams.size () ; ii++ ) {
		expParamList->copyName ( searchParamList, cparams [ii] );
	}
	if ( expParamList->getStringValue ( "data_source" ) != "List of Files" ) {
		expParamList->copyFileContents ( searchParamList, "upload_data" );
	}
// This is the output directory
	string outputDir;
	string projectName;
#ifdef MYSQL_DATABASE
	string searchKey = searchParamList->getStringValue ( "search_key" );
	if ( !searchKey.empty () ) {
		BatchJobItem* bji = MySQLPPSDDBase::instance ().getBatchJobByKey ( searchKey );
		outputDir = bji->getProjectPath ();
		projectName = bji->getProjectName ();
	}
	else {
#endif
		outputDir = searchParamList->getStringValue ( "project_filepath" );
		projectName = searchParamList->getStringValue ( "project_filename" );
		expParamList->copyName ( searchParamList, "project_filename" );
		expParamList->copyName ( searchParamList, "project_filepath" );
		expParamList->copyName ( searchParamList, "output_filepath" );
#ifdef MYSQL_DATABASE
	}
#endif
	if ( isCalibratedProject ( projectName ) ) {
		projectName = getUncalibratedProject ( projectName );
	}
// This is the output filename
	string outputFilename = projectName + ".exp.";
// Is this file there
	FileList fList ( outputDir, outputFilename, "", false );
	StringVector names = fList.getNameList ();
	bool found = false;
	string version = Version::instance ().getVersion ();
	for ( StringVectorSizeType i = 0 ; i < names.size () ; i++ ) {
		string n = names [i];
		int len = n.length ();
		int chop = 4;
		if ( n [len-2] == '_' ) {
			if ( n [len-1] == '0' ) chop = 6;
			else continue;

		}
		string filepath = outputDir + SLASH + n;
		ParameterList fParams ( filepath, false, false, false, false );
		string verExpFile = getVersionFromPPXMLFile ( filepath );
		if ( versionCheck ( verExpFile, version ) ) {	// The versions are compatible
			StringVectorSizeType j = 0;
			for ( ; j < cparams.size () ; j++ ) {
				if ( !fParams.isEqual ( expParamList, cparams [j] ) ) break;	// These parameters are different
			}
			if ( j == cparams.size () ) {	// The exp search has already been done
				outputName = names [i].substr ( 0, n.length () - chop );
				found = true;
				break;
			}
		}
	}
	if ( !found ) {
		outputName = genNextFreeFileName ( outputDir, outputFilename, "" );
	}
	expParamList->setValue ( "output_filename", outputName );
	expParamList->copyName ( searchParamList, "search_key" );

	if ( found && !force ) {
		delete expParamList;
		expParamList = 0;
	}
	return expParamList;
}
#endif
