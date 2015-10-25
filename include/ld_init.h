/******************************************************************************
*                                                                             *
*  Library    : libdbase                                                      *
*                                                                             *
*  Filename   : ld_init.h                                                     *
*                                                                             *
*  Created    : February 27th 2007                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2007-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __ld_init_h
#define __ld_init_h

#include <ostream>
#include <string>
#include <lgen_define.h>

struct st_mysql;						// forward definitions so mysql header files not required
typedef struct st_mysql MYSQL;
typedef char **MYSQL_ROW;		/* return data as array of strings */

class UserInfo {
	std::string userID;
	std::string userName;
	std::string password;
	std::string email;
	std::string directoryName;
	std::string isInactive;
	unsigned int maxMSMSSpectra;
	unsigned int maxSearchTime;
	std::string allowRaw;
public:
	UserInfo ( const MYSQL_ROW& row );
	std::string getDirectoryName () const { return directoryName; }
	std::string getUserID () const { return userID; }
	bool getIsInactive () const { return isInactive == "1"; }
	bool getIsGuest () const { return isInactive == "2"; }
	unsigned int getMaxMSMSSpectra () const { return maxMSMSSpectra; }
	unsigned int getMaxSearchTime () const { return maxSearchTime; }
	bool getAllowRaw () const { return allowRaw == "1"; }
	std::string getPassword () const { return password; }
};

class ProjectInfo {
	std::string projectID;
	std::string projectName;
	std::string projectFile;
	std::string projectPath;
	int calibrationIndex;
	std::string projectDirectory;
	std::string projectFullPath;
public:
	ProjectInfo ( const MYSQL_ROW& row );
	std::string getProjectID () const { return projectID; }
	std::string getProjectName () const { return projectName; }
	std::string getProjectFile () const { return projectFile; }
	std::string getProjectPath () const { return projectPath; }
	int getCalibrationIndex () const { return calibrationIndex; }
	std::string getProjectDirectory () const { return projectDirectory; }
	std::string getProjectFullPath () const { return projectFullPath; }
};

class DaemonJobItem {
	std::string searchJobID;
	std::string searchJobKey;
	int priority;
	std::string resultsName;
	std::string resultsPath;
	std::string resultsFile;
	std::string nodeName;
	int nodePID;
	int jobStatus;
	int jobSignal;

	std::string projectPath;
	std::string projectFile;
	std::string projectName;
	std::string projectID;

	std::string email;
public:
	DaemonJobItem ( const MYSQL_ROW& row );
	std::string getSearchJobID () const { return searchJobID; }
	std::string getSearchJobKey () const { return searchJobKey; }
	std::string getResultsName () const { return resultsName; }
	std::string getResultsPath () const { return resultsPath; }
	std::string getResultsFile () const { return resultsFile; }
	std::string getNodeName () const { return nodeName; }
	int getPID () const { return nodePID; }
	bool getJobSignalAbort () const;
	bool isSubmitted () const;
	bool isRunningOrStarted () const;
	bool isRunning () const;

	std::string getProjectPath ()	const { return projectPath; }
	std::string getProjectFile ()	const { return projectFile; }
	std::string getProjectName ()	const { return projectName; }
	std::string getProjectID ()		const { return projectID; }

	std::string getJobStatus () const;
	std::string getJobSignal () const;

	std::string getEmail () const { return email; }

	std::string getProjDir () const { return resultsPath.substr ( 4, 10 ); }
};

typedef std::vector <DaemonJobItem> VectorDaemonJobItem;
typedef VectorDaemonJobItem::size_type VectorDaemonJobItemSizeType;

class DaemonJobQueue {
	VectorDaemonJobItem jobItems;
	VectorDaemonJobItem cleanUpJobItems;
	SetString activeJobKeys;
	std::string nextSearchJobKey;
	int numRunning;
public:
	DaemonJobQueue ();
	void addItem ( const MYSQL_ROW& row );
	void setActions ( const std::string& host, int maxJobsPerUser );
	VectorDaemonJobItem getCleanUpJobItems () const { return cleanUpJobItems; }
	SetString getActiveJobKeys () const { return activeJobKeys; }
	std::string getNextSearchJobKey () const { return nextSearchJobKey; }
	int getNumRunning () const { return numRunning; }
};

class BatchJobItem {
	std::string ppUserID;
	std::string projectID;
	std::string projectName;
	std::string projectPath;
	std::string projectRelativePath;
	std::string projectFullPath;
	std::string resultsName;
	std::string resultsPath;
	std::string resultsFullPath;
	int calibrationIndex;
public:
	BatchJobItem ( const MYSQL_ROW& row );
	std::string getUserID () const { return ppUserID; }
	std::string getProjectID () const { return projectID; }
	std::string getProjectName () const { return projectName; }
	std::string getProjectPath () const { return projectPath; }
	std::string getProjectRelativePath () const { return projectRelativePath; }
	std::string getProjectFullPath () const { return projectFullPath; }
	std::string getResultsName () const { return resultsName; }
	std::string getResultsPath () const { return resultsPath; }
	std::string getResultsFullPath () const { return resultsFullPath; }
	int getCalibrationIndex () const { return calibrationIndex; }
};

class JobItem {
	std::string searchJobID;
	std::string searchJobKey;
	std::string priority;
	int numSerial;
	int searchStage;
	int searchNumber;
	int percentComplete;
	std::string resultsName;
	int jobStatus;
	std::string searchSubmitted;
	std::string searchTime;
public:
	JobItem ( const MYSQL_ROW& row );
	bool isSubmittedOrStart () const;
	bool isRunning () const;
	bool jobDone () const;
	void printJobProgress ( std::ostream& os ) const;
	std::string getSearchJobID () const { return searchJobID; }
	std::string getSearchJobKey () const { return searchJobKey; }
	std::string getJobStatus () const;
	std::string getSearchSubmitted () const { return searchSubmitted; }
	int getNumSerial () const { return numSerial; }
	int getSearchStage () const { return searchStage; }
	int getSearchNumber () const { return searchNumber; }
	int getStartSerial () const { return searchNumber == 0 ? 1 : searchNumber; }
};

class JobItemForSearchTable {
	std::string searchJobID;
	std::string searchJobKey;
	int numSerial;
	int searchStage;
	int searchNumber;
	int percentComplete;
	std::string nodeName;
	std::string resultsName;
	int jobStatus;
	std::string searchTime;
	std::string projectName;
	std::string userName;
	std::string email;
	std::string directoryName;
public:
	JobItemForSearchTable ( const MYSQL_ROW& row );
	static void printHTMLHeader ( std::ostream& os, const std::string& user );
	void printHTML ( std::ostream& os, const std::string& user ) const;
};

class JobQueue {
	std::vector <JobItem> jobItems;
public:
	JobQueue () {}
	void addItem ( const MYSQL_ROW& row );
	int numRunningJobs () const;
	int numJobsBeforeJob ( const std::string& searchJobKey ) const;
	void printJobsBeforeJob ( std::ostream& os, const std::string& searchJobKey ) const;
};

class JobQueueForSearchTable {
	std::vector <JobItemForSearchTable> jobItems;
public:
	JobQueueForSearchTable () {}
	void addItem ( const MYSQL_ROW& row );
	void printHTML ( std::ostream& os, const std::string& user ) const;
};

class ProjectDateFilter {
	std::string dateFilter;
	std::string startYear;
	std::string endYear;
	std::string startMonth;
	std::string endMonth;
public:
	ProjectDateFilter ( const std::string& dateFilter, const std::string& startYear, const std::string& endYear, const std::string& startMonth, const std::string& endMonth );
	std::string getDateFilter ()const { return dateFilter; }
	std::string getStartYear ()	const { return startYear; }
	std::string getEndYear ()	const { return endYear; }
	std::string getStartMonth ()const { return startMonth; }
	std::string getEndMonth ()	const { return endMonth; }
};

struct DatabaseImportSearchJobInfo;

class MySQLPPSDDBase {
	static const std::string QUOTE;
	static const std::string EQUOTE;
	static const std::string searchProgram;
	static unsigned int port;
	static char* sock;
	static unsigned int flags;
	MYSQL* conn;
	MySQLPPSDDBase ();
	MySQLPPSDDBase ( unsigned int conn );
	unsigned int submitSQL ( const std::string& sql ) const;
	unsigned int deleteQuery ( const std::string& table, const std::string& field, const std::string& value ) const;
	unsigned int updateQuery ( const std::string& table, const PairStringString& setClause, const PairStringString& whereClause ) const;
	unsigned int insertQuery ( const std::string& table, const StringVector& names, const StringVector& values ) const;
	bool checkKeyPresence ( const std::string& sql ) const;
	static JobItem* getJobItem ( MYSQL* conn, const std::string& sql );
	static std::string getSingleString ( MYSQL* conn, const std::string& sql );
	static std::string getFullPath ( MYSQL* conn, const std::string& sql );
	std::pair <StringVector, StringVectorVector> getSQLResults ( const std::string& sql ) const;
	void printSQLResultsHTML ( std::ostream& os, const std::pair <StringVector, StringVectorVector>& p ) const;
	void printSQLResultsTabDelimitedText ( std::ostream& os, const std::pair <StringVector, StringVectorVector>& p ) const;
	void printSQLResultsXML ( std::ostream& os, const std::pair <StringVector, StringVectorVector>& p ) const;
	std::pair <StringVector, StringVectorVector> getUserInfoTable () const;
public:
	~MySQLPPSDDBase ();
	static MySQLPPSDDBase& instance ( bool reset = false, bool close = false );
	static std::string getJobStatus ( int jobStatus );
	static std::string getJobSignal ( int jobSignal );
	static std::string getCellValue ( const MYSQL_ROW& r, int i );
	bool checkSearchKey ( const std::string& searchJobKey ) const;
	bool checkProject ( const std::string& userID, const std::string& projectName ) const;
	bool checkResults ( const std::string& userID, const std::string& projectName, const std::string& resultsName ) const;
	std::string getProjectID ( const std::string& userID, const std::string& projectName ) const;
	UserInfo* getUserInfo ( const std::string& userName ) const;
	ProjectInfo* getProjectInfo ( const std::string& user, const std::string& project ) const;
	std::string getUserID ( const std::string& userName ) const;
	std::string getUserName ( const std::string& key ) const;
	bool submitUser ( const std::string& userName, const std::string& secretPhrase, const std::string& email, const std::string& directoryName, const std::string& firstName, const std::string& lastName ) const;
	bool submitProject ( const std::string& userID, const std::string& projectName, const std::string& projectFile, const std::string& projectPath, const std::string& instrument ) const;
	bool submitCopiedProject ( const std::string& userID, const std::string& projectName, const std::string& projectFile, const std::string& projectPath, const std::string& instrument, const std::string& calibrationIndex, const std::string& recordCreated ) const;
	void submitSearch ( const std::string& searchJobKey, const std::string& projectID, const std::string& resultsName, const std::string& resultsFile, const std::string& resultsPath ) const;
	void submitCopiedSearch ( const std::string& searchJobKey, const std::string& projectID, const std::string& resultsFile, const std::string& resultsPath, const DatabaseImportSearchJobInfo* disji ) const;
	bool submitSession ( const std::string& sessionKey, const std::string& ppUserID ) const;
	bool submitAbortedJob ( const std::string& searchJobKey, const std::string& message ) const;
	JobItem* getSearchJobByKey ( const std::string& searchJobKey ) const;
	std::string getSearchIDByKey ( const std::string& searchJobKey ) const;
	std::string getUserNameByKey ( const std::string& searchJobKey ) const;
	std::string getSearchTimeByKey ( const std::string& searchJobKey ) const;
	std::string getSearchEndTimeByKey ( const std::string& searchJobKey ) const;
	BatchJobItem* getBatchJobByKey ( const std::string& searchJobKey ) const;
	void setNumSerial ( const std::string& searchJobID, int numSerial ) const;
	void updateSearchStage ( const std::string& searchJobKey, int searchStage ) const;
	void updateSearchNumber ( const std::string& searchJobID, int searchNumber ) const;
	void updatePercentComplete ( const std::string& searchJobID, int percentComplete ) const;
	void setAbortSignal ( const std::string& searchJobID ) const;
	void setPassword ( const std::string& user, const std::string& password ) const;
	void setJobSubmitted ( const std::string& searchKey ) const;
	bool setJobStart ( const std::string& searchKey ) const;
	void setJobSearching ( const std::string& searchKey, const std::string& host, int pid ) const;
	void setJobDone ( const std::string& searchJobID ) const;
	void setJobAborted ( const std::string& searchJobID ) const;
	bool setJobAbortedUnknown ( const std::string& searchJobID ) const;
	bool processAbortSignal ( const std::string& searchJobID ) const;
	JobQueue getJobQueue () const;
	JobQueueForSearchTable getSearchTableJobQueue () const;
	DaemonJobQueue getDaemonJobQueue () const;
	DaemonJobItem* getDaemonJobItemByKey ( const std::string& searchJobKey ) const;
	StringVector getProjectList ( const std::string& user, bool noResults, const ProjectDateFilter* pdf = 0, bool compressInfo = false ) const;
	bool projectHasQueuedOrRunningJobs ( const std::string& user, const std::string& project ) const;
	StringVector getResultsList ( const std::string& projectID ) const;
	StringVector getResultsList ( const std::string& user, const std::string& project ) const;
	std::string getResultsFullPath ( const std::string& user, const std::string& project, const std::string& results ) const;
	std::string getResultsFullPath ( const std::string& searchKey ) const;
	std::string getProjectFullPath ( const std::string& user, const std::string& project ) const;
	std::string getInstrument ( const std::string& user, const std::string& project ) const;
	std::string getErrorMessage ( const std::string& searchJobKey ) const;
	void updateProjectRecordUpdated ( const std::string& projectID );
	void updateCalIndex ( const std::string& projectID, int calIndex );
// delete functions
	unsigned int deleteSearchKey ( const std::string& searchKey ) const;
	unsigned int deleteSearchID ( const std::string& searchID ) const;
	unsigned int deleteProjectID ( const std::string& projectID ) const;
	unsigned int deleteUserID ( const std::string& userID ) const;
	unsigned int deleteUserName ( const std::string& userName ) const;
	unsigned int deleteSession ( const std::string& sessionKey ) const;
	unsigned int deleteOldSessions ( int maxDays ) const;
	unsigned int deleteOldAbortedJobs ( int maxDays ) const;
	void printSQLResultsHTML ( std::ostream& os, const std::string& sql ) const;
	void printSQLResultsTabDelimitedText ( std::ostream& os, const std::string& sql ) const;
	void printSQLResultsXML ( std::ostream& os, const std::string& sql ) const;
	void printUserInfoHTML ( std::ostream& os ) const;
	void printUserInfoTabDelimitedText ( std::ostream& os ) const;
};

class Repository;

typedef std::map <std::string, DatabaseImportSearchJobInfo*> MapStringToDatabaseImportSearchJobInfoPtr;

class DatabaseProjectImport {
	StringVector instrument;
	StringVector calibrationIndex;
	StringVector recordCreated;

	MapStringToDatabaseImportSearchJobInfoPtr msdisji;

	void readProjectExportFile ( const std::string& file );
	void readSearchJobsExportFile ( const std::string& file );
public:
	DatabaseProjectImport ( const std::string& dir );
	~DatabaseProjectImport ();
	void submit ( const std::string& userID, const std::string& projName, const std::string& projFile, const Repository* reposit, const StringVector& resFiles ) const;
	void updateSearchKey ( const std::string& oldSearchKey, const std::string& newSearchKey );
};

#endif /* ! __ld_init_h */
