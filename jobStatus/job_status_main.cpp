/******************************************************************************
*                                                                             *
*  Program    : jobStatus                                                     *
*                                                                             *
*  Filename   : job_status_main.cpp                                           *
*                                                                             *
*  Created    : March 12th 2007                                               *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2007-2009) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifndef VIS_C
#include <stdexcept>
#endif
#include <ld_init.h>
#include <lu_param_list.h>
#include <lu_getfil.h>
#include <lu_prog_par.h>
#include <lu_prog.h>
#include <lu_html.h>
#include <lu_html_form.h>
using std::ostream;
using std::cout;
using std::string;
using std::endl;
using std::runtime_error;

int main ( int argc, char** argv )
{
	initialiseProspector ();

	ParameterList paramList ( argc, argv );
	try {
		ostream& os = cout;
		ProgramLink::setParams ( &paramList );
		int jobStatusRefresh = InfoParams::instance ().getIntValue ( "job_status_refresh_time", 5 );
		string searchJobKey = paramList.getStringValue ( "search_key" );
		if ( searchJobKey.empty () ) {
			init_html ( os, "Search Table" );
			try {
				JobQueueForSearchTable jq = MySQLPPSDDBase::instance ().getSearchTableJobQueue ();
				string user;
				if ( !ParameterList::getCookieValue ( "key" ).empty () ) {
					user = ParameterList::getCookieValue ( "user" );
				}
				jq.printHTML ( os, user );
			}
			catch ( runtime_error e ) {		// Catch database login problems
				ErrorHandler::genError ()->error ( e );
			}
			refreshJavascript ( os, jobStatusRefresh*1000, ProgramLink::getURLStart ( "jobStatus" ), false, true );
		}
		else {
			init_html ( os, "Job Status" );
			JobItem* jobItem;
			try {
				jobItem = MySQLPPSDDBase::instance ().getSearchJobByKey ( searchJobKey );
			}
			catch ( runtime_error e ) {		// Catch database login problems
				ErrorHandler::genError ()->error ( e );
			}
			if ( jobItem ) {
				bool abort = paramList.getBoolValue ( "abort" );
				if ( jobItem->jobDone () ) {		// Done
					ParameterList pList ( "msform", false, false, false, false, false );
					pList.addName ( "form", "search_compare" );
					pList.addName ( "search_key", searchJobKey );
					refreshJavascript ( os, 500, pList.getURL (), true );
				}
				else if ( jobItem->isRunning () || jobItem->isSubmittedOrStart () ) {
					if ( abort ) {
						MySQLPPSDDBase::instance ().setAbortSignal ( jobItem->getSearchJobID () );
						back2Javascript ( os );
					}
					else {
						printAbortButton ( os, searchJobKey );
						if ( jobItem->isRunning () ) {	// Running
							jobItem->printJobProgress ( os );
						}
						else {							// Submitted
							JobQueue jq = MySQLPPSDDBase::instance ().getJobQueue ();
							jq.printJobsBeforeJob ( os, searchJobKey );
						}
						os << "<p>This page will refresh in " << jobStatusRefresh << " seconds.</p>" << endl;
						//os << "<p>When the search completes, an email will be sent to you.</p>" << endl;
						os << "<p>Please do not resubmit the search if your connection to this page is lost.</p>" << endl;
						refreshJavascript ( os, jobStatusRefresh*1000, ProgramLink::getURLStart ( "jobStatus" ) + "?search_key=" + searchJobKey, false, true );
					}
				}
				else {
					os << "<p>Search job " << searchJobKey << " is not running and not done.</p>" << endl;
					os << "<p>The search daemon reports the status as " << jobItem->getJobStatus () << ".</p>" << endl;
				}
			}
			else {
				string message = MySQLPPSDDBase::instance ().getErrorMessage ( searchJobKey );
				if ( !message.empty () ) {
					os << "<p>Search job " << searchJobKey << " aborted.</p>" << endl;
					ErrorHandler::genError ()->error ( message );
				}
				else
					os << "<p>Search job " << searchJobKey << " was not found.</p>" << endl;
			}
		}
		printHTMLFooter ( os, "2007" );
	}
	catch ( runtime_error e ) {
		paramList.writeLogError ( e.what () );
	}
	MySQLPPSDDBase::instance ( false, true );
	return 0;
}
