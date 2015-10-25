/******************************************************************************
*                                                                             *
*  Library    : liboracle                                                     *
*                                                                             *
*  Filename   : lo_job_run.h                                                  *
*                                                                             *
*  Created    :                                                               *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2002-2009) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifndef __lo_job_run_h
#define __lo_job_run_h

#include <lgen_define.h>
#include <string>
#include <map>
#include <oci.h>
#include <lo_init.h>

class JobRunItemInfo {
	unsigned int jobRunID;
	bool msFlag;
	double msmsMass;
	std::string spotLabel;
public:
	JobRunItemInfo () {};
	JobRunItemInfo ( unsigned int jobRunID, bool msFlag, double msmsMass, std::string spotLabel ) :
		jobRunID ( jobRunID ),
		msFlag ( msFlag ),
		msmsMass ( msmsMass ),
		spotLabel ( spotLabel ) {}
	unsigned int getJobRunID () const { return jobRunID; }
	bool getMSFlag () const { return msFlag; }
	double getMsmsMass () const { return msmsMass; }
	std::string getSpotLabel () const { return spotLabel; }
};

class JobRunItems : public OracleStatement {
	typedef std::map <unsigned int, JobRunItemInfo> MapUIntToJobRunItemInfo;
	typedef MapUIntToJobRunItemInfo::const_iterator MapUIntToJobRunItemInfoConstIterator;

	MapUIntToJobRunItemInfo items;
	MapUIntToUInt runMap;
	static bool verbose;
	static OCILobLocator* blob;
	static void printHTMLHeader ( std::ostream& os );
public:
	JobRunItems ( OracleConnection& oc, unsigned int nSpotSetID, bool rawFlag, const std::string& rawDir, bool writeMS, bool writeMSMS );
	virtual ~JobRunItems ();
	int getRunListSize () { return runMap.size (); }
	unsigned int getRunNumber ( unsigned int runID ) const
	{
		MapUIntToUIntConstIterator cur = runMap.find ( runID );
		if ( cur != runMap.end () ) {
			return (*cur).second;
		}
		else return 0;
	}
	const bool getJobRunItemInfo ( JobRunItemInfo& jrii, unsigned int runItemID ) const
	{
		MapUIntToJobRunItemInfoConstIterator cur = items.find ( runItemID );
		if ( cur != items.end () ) {
			jrii = (*cur).second;
			return true;
		}
		else return false;	// bug
	}
	static void setVerbose ( bool flag ) { verbose = flag; }
	static OCILobLocator* getBlob () { return blob; };
};

#endif /* ! __lo_job_run_h */
