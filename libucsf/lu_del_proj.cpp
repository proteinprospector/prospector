/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_del_proj.cpp                                               *
*                                                                             *
*  Created    : September 25th 2007                                           *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2007-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef MYSQL_DATABASE
#include <lg_string.h>
#include <lgen_file.h>
#include <ld_init.h>
#include <lu_getfil.h>
#include <lu_proj_file.h>
#include <lgen_uncompress.h>
using std::string;

void deleteResults ( const string& user, const string& project, const string& results )
{
	string resultsFullPath = MySQLPPSDDBase::instance ().getResultsFullPath ( user, project, results );
	string searchKey = genShortFilenameFromPath ( resultsFullPath );
	MySQLPPSDDBase::instance ().deleteSearchKey ( searchKey );		// Delete the search database entry.
	string d = genDirectoryFromPath ( resultsFullPath );
	string f = genFilenameFromPath (resultsFullPath);
	StringVector sv;
	sv.push_back ( searchKey + "*" );
	FileList f1 ( d, sv );
	f1.unlink ();													// Delete the results file(s).
	genRmdir ( d );													// Try to delete the results date directory. This will fail if the directory isn't empty.
}
bool deleteProject ( const string& user, const string& project )
{
	bool calProject = project.rfind ( ".cal." ) != string::npos;

	if ( MySQLPPSDDBase::instance ().projectHasQueuedOrRunningJobs ( user, project ) ) {
		return false;
	}
	StringVector resultsList = MySQLPPSDDBase::instance ().getResultsList ( user, project );	// Delete all the results
	for ( StringVectorSizeType i = 0 ; i < resultsList.size () ; i++ ) {
		deleteResults ( user, project, resultsList [i] );
	}

	ProjectInfo* pi = MySQLPPSDDBase::instance ().getProjectInfo ( user, project );

	string projectPath = pi->getProjectFullPath ();
	if ( genFileExists ( projectPath + ".7z" ) ) {		// Check that the project isn't compressed
		gen7zaUncompress ( projectPath + ".7z" );
	}
	if ( !calProject ) {
		ProjectFile pf ( projectPath );
		if ( pf.isUploadProject () ) {
			string dataPath = pf.getUploadDirectory ();
			if ( genFileExists ( dataPath ) ) genUnlinkDirectory ( dataPath );				// Delete the data
			genRmdir ( genDirectoryFromPath ( dataPath ) );	// Try to delete the data date directory. This will fail if the directory isn't empty.
		}
	}
	genUnlink ( projectPath );													// Delete the project file.

	if ( !calProject ) {
		FileList f1 ( pi->getProjectDirectory (), pi->getProjectName () + ".exp.", "", false );	// Delete the .exp files.
		f1.unlink ();
	}

	MySQLPPSDDBase::instance ().deleteProjectID ( pi->getProjectID () );						// Delete the database entry.

	if ( !calProject ) {																		// Delete the calibrated projects.
		StringVector sv = MySQLPPSDDBase::instance ().getProjectList ( user, true );
		for ( StringVectorSizeType j = 0 ; j < sv.size () ; j++ ) {
			if ( isPrefix ( sv [j], project + ".cal." ) ) {
				deleteProject ( user, sv [j] );
			}
		}
	}
	genRmdir ( genDirectoryFromPath ( projectPath ) );	// Try to delete the project date directory. This will fail if the directory isn't empty.
	delete pi;
	return true;
}
#endif
