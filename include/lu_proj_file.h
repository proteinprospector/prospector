/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_proj_file.h                                                *
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
*  Copyright (2004-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_proj_file_h
#define __lu_proj_file_h

#include <string>
#include <lgen_define.h>

class ParameterList;

class ProjectFile {
	std::string projectVersion;
	std::string toleranceUnits;
	StringVector centroidFiles;
	StringVector rawFiles;
	IntVector numMSSpectra;
	IntVector numMSMSSpectra;
	DoubleVector tolerance;
	DoubleVector offset;
	static std::string centroidBaseDir;
	static std::string rawBaseDir;
	static std::string userBaseDir;
	void init ( const std::string& projectFilePath );
public:
	ProjectFile ( const std::string& path );
	ProjectFile ( const ParameterList* params );
	bool isUploadProject () const;
	std::string getUploadDirectory () const;
	std::string getCentroidPath ( StringVectorSizeType i ) const;
	std::string getAbsoluteRawPath ( StringVectorSizeType i ) const;
	std::string getToleranceUnits () const { return toleranceUnits; }
	StringVector getRawFileList () const { return rawFiles; }
	StringVectorSizeType getNumFiles () const { return centroidFiles.size (); }
	std::string getProjectVersion () const { return projectVersion; }
	double getTolerance ( int i ) const { return tolerance [i]; }
	double getOffset ( int i ) const { return offset [i]; }
	DoubleVector getTolerance () const { return tolerance; }
	DoubleVector getOffset () const { return offset; }
	IntVector getNumMSMSSpectra () const { return numMSMSSpectra; }
	int getNumMSMSSpectra ( int i ) const { return numMSMSSpectra [i]; }
	static std::string updateProjectFile ( const ParameterList* params );
	static StringVector getFractionNameList ( const ParameterList* params );
	static StringVector getCentroidPathList ( const ParameterList* params );
	static DoubleVector getOffsets ( const ParameterList* params, double defaultOffset );
	static IntVector getNumMSSpectra ( const ParameterList* params );
	static StringVector getRawTypes ( const ParameterList* params );
	static bool getSpottingPlate ( const ParameterList* params );
};
std::string getProjectDir ( const ParameterList* params );
std::string getCentroidDataFilename ( const ParameterList* params, int fraction );
std::string getRawDataFilename ( const ParameterList* params, int fraction );
std::string getProjectVersion ( const ParameterList* params, int fraction );
bool isCalibratedProject ( const std::string& project );
std::string getUncalibratedProject ( const std::string& project );

#endif /* ! __lu_proj_file_h */
