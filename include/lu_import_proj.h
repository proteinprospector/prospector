/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_import_proj.h                                              *
*                                                                             *
*  Created    : June 28th 2012                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2012-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_import_proj_h
#define __lu_import_proj_h

#ifdef MYSQL_DATABASE
#include <string>
#include <lgen_define.h>

class ParameterList;
class Repository;
class ImportProjectFile;
class DatabaseProjectImport;

class ImportProject {
	std::string newProjectName;
	Repository* reposit;
	std::string userID;
	std::string uploadName;

	StringVector dataFiles;

	std::string projName;
	std::string projFile;
	StringVector expFiles;

	StringVector resFiles;
	StringVector discFiles;

	DatabaseProjectImport* dpi;
	ImportProjectFile* ipf;

	void init ( const ParameterList& paramList );
	bool getDataFiles ();
	void getProjectFiles ();
	void getResultsFiles ();
	void parseResultsFiles ();
	static bool importDataStreamEditor ( const std::string& file, const std::string& newFile, const PairStringString& edit );
	static bool importDataStreamEditor ( const std::string& file, const std::string& newFile, const std::string& searchKey, const std::string& newProject );
	void writeDataFiles ( const std::string& key ) const;
	void writeResultsFiles () const;
	void writeProjectFiles ( const std::string& key ) const;
public:
	ImportProject ( const ParameterList& paramList );
	~ImportProject ();
	void write ( std::ostream& os ) const;
};
#endif

#endif /* ! __lu_import_proj_h */
