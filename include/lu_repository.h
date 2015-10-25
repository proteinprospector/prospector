/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_repository.h                                               *
*                                                                             *
*  Created    : March 5th 2007                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2007-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_repository_h
#define __lu_repository_h

#include <string>

class Repository {
	std::string baseDir;
	std::string searchName;
protected:
	std::string dataPath;
	std::string projectPath;
	std::string resultsPath;
public:
	Repository ( const std::string& baseDir, const std::string& program );
	std::string getBaseDir () const { return baseDir; }
	std::string getSearchName () const { return searchName; }
	std::string getDataPath () const { return dataPath; }
	std::string getProjectFileDataPath ( bool upload ) const;
	std::string getProjectPath () const { return projectPath; }
	std::string getResultsPath () const { return resultsPath; }
	std::string getFullDataPath () const { return baseDir + SLASH + dataPath; }
	std::string getFullProjectPath () const { return baseDir + SLASH + projectPath; }
	std::string getFullResultsPath () const { return baseDir + SLASH + resultsPath; }
};

class UserRepository : public Repository {
	static std::string uploadRepository;
public:
	UserRepository ( const std::string& program, const std::string& userDir );
};

class FractionFile;
class MSMSDataPoint;
class PPExpat;
class SpecID;
typedef std::vector <MSMSDataPoint> MSMSDataPointVector;
typedef std::vector <const FractionFile*> VectorFractionFilePtr;
typedef VectorFractionFilePtr::size_type VectorFractionFilePtrSizeType;

class PPProject {
	bool upload;
	StringVector t2ds;
	StringVector wiffs;
	StringVector wiffScans;
	StringVector mgfs;
	StringVector ms2s;
	StringVector apls;
	StringVector ppsfs;
	StringVector dtas;
	StringVector pkls;
	StringVector finnRaws;
	StringVector xmls;

	int nFiles;

	StringVector fractionNames;
	StringVector centroidFiles;
	StringVector rawFiles;

	std::string name;
	VectorFractionFilePtr fractionFiles;
	unsigned int totalMSSpectra;
	unsigned int totalMSMSSpectra;
	unsigned int maxMSMSSpectra;
	bool rawDataAllowed;
	bool duplicateScans;
	IntVector chargeRange;
	bool spottingPlate;
	bool init;
	std::string errMessage;
	bool deleteFlag;
	bool parseFile ( const std::string& fullPath, const std::string& f );
	bool parseCentroidFile ( const std::string& fullPath, const std::string& f );
	bool getDTAPKLFractionMap ( MapStringToStringVector& fNameMap, const StringVector& fNames );
	bool parseDTAs ( const std::string& uploadName );
	bool parseDTAs2 ( const std::string& uploadName );
	bool parsePKLs ( const std::string& uploadName );
	bool parsePKLs2 ( const std::string& uploadName );
	void processDTAPKLHeader ( std::ostream& osf, const std::string& titleInfo, double mOZ, int charge, double intensity );
	bool processDTAPKLFragmentIons ( std::ostream& ofs, std::istream& ifs, bool flag );
	bool parseAPLs ( const std::string& uploadName );
	static void parseAPLHeaderLines ( std::istream& istr, double& mz, int& z, std::string& title, std::string& file );
	bool parseXMLs ( const std::string& uploadName, StringVector& xmlFiles );
	static std::string getXMLFractionName ( const std::string& f );
	PPExpat* getPPExpatPtr ( const std::string& fullPath, MSMSDataPointVector& msmsDataPointList, int fraction, const SpecID& specID ) const;
	bool parseWiffFiles ( const std::string& uploadName );
	bool parseThermoFiles ();
	bool parseT2DFiles ();
	bool parseCentroidFiles ();
	bool parseCentroidFiles ( const StringVector& files );
	bool createFractionFiles ( const std::string& uploadName, const std::string& searchKey, const Repository* reposit );
	bool createFractionFiles ();
	bool createFractionFiles ( const std::string& dir );
	void printTotalSpectraHTML ( std::ostream& os ) const;
	void removeMGFZeroIntensities ( const std::string& name );
	void removeMGFZeroIntensities ( const std::string& name1, const std::string& name2 );
public:
	PPProject ( const std::string& name, const StringVector& fileList );
	PPProject ( const std::string& name, const std::string& dir );
	PPProject ( unsigned int maxMSMSSpectra, bool rawDataAllowed, const std::string& name, const StringVector& files );
	PPProject ( unsigned int maxMSMSSpectra, bool rawDataAllowed, const std::string& name, const std::string& uploadName, const std::string& searchKey, const Repository* reposit, const IntVector& chargeRange );
	PPProject ( const std::string& uploadName );			// MS-Viewer
	void createProjectFile ( const std::string& filename, bool fullPath = false );
	void addFractionFile ( const FractionFile* fFile ) { fractionFiles.push_back ( fFile ); }
	bool initialised () const { return init; }
	std::string getErrMessage () const { return errMessage; }
	bool getDeleteFlag () const { return deleteFlag; }
	StringVector getFractionNames () const { return fractionNames; }
	StringVector getCentroidFiles () const { return centroidFiles; }
	static void writeMGFScan ( std::ostream& os, MSMSDataPoint& msms, int charge );
	static void parseAPLFile ( const std::string& uploadName, const std::string& f );
};

#endif /* ! __lu_repository_h */
