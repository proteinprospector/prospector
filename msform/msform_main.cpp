/******************************************************************************
*                                                                             *
*  Program    : msform                                                        *
*                                                                             *
*  Filename   : msform_main.cpp                                               *
*                                                                             *
*  Created    : December 1st 2004                                             *
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
#ifndef VIS_C
#include <stdexcept>
#endif
#include <lgen_file.h>
#include <lgen_process.h>
#include <lg_email.h>
#include <lu_cookie.h>
#include <lu_prog.h>
#include <lu_prog_par.h>
#include <lu_scsel_form.h>
#include <lu_btsel_form.h>
#include <lu_scomp_form.h>
#include <lu_btag_form.h>
#include <lu_dbst_form.h>
#include <lu_patt_form.h>
#include <lu_file_type.h>
#include <lu_hom_form.h>
#include <lu_brid_form.h>
#include <lu_msfit_form.h>
#include <lu_mstag_form.h>
#include <lu_filter_form.h>
#include <lu_viewer_form.h>
#include <lu_prod_form.h>
#include <lu_comp_form.h>
#include <lu_faind_form.h>
#include <lu_proj_file.h>
#include <lu_proj_form.h>
#include <lu_getfil.h>
#include <lu_param_list.h>
#include <lu_html.h>
#include <lu_login_form.h>
#include <lu_repository.h>
#include <lu_add_user.h>
#include <lu_dlsel_form.h>
#include <lu_del_proj.h>
#include <lu_xml.h>
#ifdef MYSQL_DATABASE
#include <lu_repos_info.h>
#include <ld_init.h>
#include <lu_import_proj.h>
#include <lu_export_proj.h>
#endif
using namespace FileTypes;
using std::ostream;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::copy;
using std::back_inserter;
using std::runtime_error;

namespace {
#ifdef MYSQL_DATABASE
	void importProject ( ostream& os, const ParameterList& paramList );
	void addUserToDatabase ( ostream& os, const ParameterList& paramList );
	string getUserFromCookieKey ();
	bool checkGuest ( const string& form, const UserInfo* userInfo );
	bool checkUser ( ParameterList& paramList, bool& setCookie );
	void deleteData ( const ParameterList& paramList );
#endif
	void makeProject ( ostream& os, const ParameterList& paramList );
string getDefaultResultsPath ( const string& form, const string& defaultFile = "default.xml" )
{
	return adjustPPOutputPath ( string ( "params" ) + SLASH + form ) + defaultFile;
}
ParameterList* getDefaultParameters ( const string& form )
{
#ifdef MYSQL_DATABASE
	if ( form == "batchtag" ) {
		string params = ParameterList::getCookieValue ( form + "_params" );
		if ( !params.empty () ) {
			return new ParameterList ( params, false, false );
		}
	}
#endif
	return new ParameterList ( getDefaultResultsPath ( form ), false, false, false, false );
}
}

int main ( int argc, char** argv )
{
	initialiseProspector ();
	ParameterList paramList ( argc, argv );
#ifdef MYSQL_DATABASE
	bool mySQLLoginFault = false;
#endif
	try {
		ProgramLink::setParams ( &paramList );
		ostream& os = cout;

		bool selectProject = paramList.getBoolValue ( "select_project" );
		bool selectResults = paramList.getBoolValue ( "select_results" );
		bool addUser = paramList.getBoolValue ( "add_user" );
		string form = paramList.getStringValue ( "form" );

		VectorConstParameterListPtr pv;
		ProspectorForm* f = 0;
		if ( selectProject || ( selectResults && paramList.getStringValue ( "project_name" ).empty () ) || paramList.getBoolValue ( "filter_project_list" ) ) {
#ifdef MYSQL_DATABASE
			try {			// Check the database connection
				MySQLPPSDDBase::instance ();
			}
			catch ( runtime_error e ) {		// Catch database login problems
				mySQLLoginFault = true;
				ErrorHandler::genError ()->error ( e );
			}
			bool setCookie = false;
			if ( checkUser ( paramList, setCookie ) ) {
				if ( form == "batchtag" )
					f = new BTSelProjectForm ( &paramList, setCookie );
				else if ( form == "search_compare" )
					f = new SCSelProjectForm ( &paramList, setCookie );
				else
					f = new DeleteSelForm ( &paramList, setCookie );
			}
			else
				f = new LoginForm ( &paramList );
#endif
		}
#ifdef MYSQL_DATABASE
		else if ( selectResults ) {
			string project = paramList.getStringValue ( "project_name" );
			if ( project == "Make Project...." ) {
				f = new MakeProjectForm ( &paramList );
			}
			else if ( project == "Import Project...." ) {
				f = new ImportProjectForm ( &paramList );
			}
			else {
				bool flag = uncompressProjects2 ( paramList );
				if ( flag ) {
					ParameterList pL ( "msform", false, false, false, false, false );
					pL.addName ( "form", "search_compare" );
					pL.addName ( "select_results", "1" );
					pL.copyName ( &paramList, "user" );
					pL.copyName ( &paramList, "project_name" );
					refreshJavascript ( os, 500, pL.getURL (), true );
					printHTMLFooter ( os, "2012" );
				}
				else {
					if ( paramList.getStringValue ( "form" ) == "batchtag" )
						f = new BTSelResultsForm ( &paramList );
					else
						f = new SCSelResultsForm ( &paramList );
				}
			}
		}
#endif
		else {
			if ( form == "batchtag" ) {
#ifdef MYSQL_DATABASE
				string searchKey = paramList.getStringValue ( "search_key" );
				string resultsFullPath;
				string instrument;
				string cookieParams;
				if ( !searchKey.empty () ) {
					resultsFullPath = MySQLPPSDDBase::instance ().getResultsFullPath ( searchKey );
					if ( !genFileExists ( resultsFullPath ) ) resultsFullPath += "_0";
					string u = MySQLPPSDDBase::instance ().getUserNameByKey ( searchKey );
					BatchJobItem* bti = MySQLPPSDDBase::instance ().getBatchJobByKey ( searchKey );
					paramList.addName ( "user", u );
					paramList.addName ( "project_name", bti->getProjectName () );
				}
				else {
					string resultsFile = paramList.getStringValue ( "results_file", "Default...." );
					string u = paramList.getStringValue ( "user" );
					string p = paramList.getStringValue ( "project_name" );
					if ( resultsFile == "Default...." ) {
						cookieParams = ParameterList::getCookieValue ( "batchtag_params" );
						resultsFullPath = getDefaultResultsPath ( "batchtag" );
						instrument = MySQLPPSDDBase::instance ().getInstrument ( u, p );
					}
					else {
						resultsFullPath = MySQLPPSDDBase::instance ().getResultsFullPath ( u, p, resultsFile );
						if ( !genFileExists ( resultsFullPath ) ) resultsFullPath += "_0";
					}
				}
				if ( !cookieParams.empty () )
					pv.push_back ( new ParameterList ( cookieParams, false, false ) );
				else
					pv.push_back ( new ParameterList ( resultsFullPath, false, false, false, false ) );

				ParameterList pl = paramList;
				pl.removeName ( "form" );
				pl.removeName ( "results_file" );
				pl.removeName ( "search_key" );
				const_cast <ParameterList*>(pv [0])->appendParameters ( pl );
				if ( !instrument.empty () ) {
					const_cast <ParameterList*>(pv [0])->addOrReplaceName ( "instrument_name", instrument );
				}
				f = new BatchTagFullForm ( pv, "Search" );
#endif
			}
			else if ( form == "search_table" ) {
#ifdef MYSQL_DATABASE
				bool setCookie = false;
				if ( checkUser ( paramList, setCookie ) ) {
					if ( setCookie ) setUserCookie ( os, paramList.getStringValue ( "user" ) );
					refreshJavascript ( os, 0, ProgramLink::getURLStart ( "jobStatus" ), true );
					return 0;
				}
				else
					f = new LoginForm ( &paramList );
#endif
			}
			else if ( form == "mstag" ) {	// MS-Tag from File
#ifdef MYSQL_DATABASE
				string resultsFullPath = MySQLPPSDDBase::instance ().getResultsFullPath ( paramList.getStringValue ( "search_key" ) );
				if ( !genFileExists ( resultsFullPath ) ) resultsFullPath += "_0";
				pv.push_back ( new ParameterList ( resultsFullPath, false, false, false, false ) );
				ParameterList pl = paramList;
				pl.removeName ( "form" );
				const_cast <ParameterList*>(pv [0])->appendParameters ( pl );
				f = new MSTagFileForm ( pv );
#endif
			}
			else if ( form == "mstagfromviewer" ) {	// MS-Tag from File
				pv.push_back ( getDefaultParameters ( "mstag" ) );
				ParameterList pl = paramList;
				pl.removeName ( "form" );
				const_cast <ParameterList*>(pv [0])->removeName ( "const_mod" );
				const_cast <ParameterList*>(pv [0])->appendParameters ( pl );
				f = new MSTagFileFormFromViewer ( pv );
			}
#ifdef MYSQL_DATABASE
			else if ( form == "msbridge" ) {	// MS-Bridge
				string resultsFullPath = MySQLPPSDDBase::instance ().getResultsFullPath ( paramList.getStringValue ( "search_key" ) );
				if ( !genFileExists ( resultsFullPath ) ) resultsFullPath += "_0";
				pv.push_back ( new ParameterList ( resultsFullPath, false, false, false, false ) );
				ParameterList pl = paramList;
				pl.removeName ( "form" );
				const_cast <ParameterList*>(pv [0])->appendParameters ( pl );
				f = new MSBridgeFileForm ( pv );
			}
			else if ( form == "search_compare" ) {
				string searchKey = paramList.getStringValue ( "search_key" );
				if ( !searchKey.empty () ) {
					BatchJobItem* bji = MySQLPPSDDBase::instance ().getBatchJobByKey ( searchKey );
					string resultsFullPath = bji->getResultsFullPath ();
					if ( !genFileExists ( resultsFullPath ) ) resultsFullPath += "_0";
					pv.push_back ( new ParameterList ( resultsFullPath, false, false, false, false ) );
				}
				else {
					StringVector files;
					files.push_back ( paramList.getStringValue ( "results_file" ) );
					StringVector compareFiles = paramList.getStringVectorValue ( "compare_files" );
					copy ( compareFiles.begin (), compareFiles.end (), back_inserter ( files ) );
					for ( StringVectorSizeType i = 0 ; i < files.size () ; i++ ) {
						StringVector sv = genGetSubstrings ( files [i], '/' );
						string resultsFullPath;
						if ( sv.size () == 3 )
							resultsFullPath = MySQLPPSDDBase::instance ().getResultsFullPath ( sv [0], sv [1], sv [2] );
						else
							resultsFullPath = MySQLPPSDDBase::instance ().getResultsFullPath ( paramList.getStringValue ( "user" ), sv [0], sv [1] );
						if ( !genFileExists ( resultsFullPath ) ) resultsFullPath += "_0";
						pv.push_back ( new ParameterList ( resultsFullPath, false, false, false, false ) );
					}
				}
				string cookieParams = ParameterList::getCookieValue ( "search_compare_params" );
				if ( !cookieParams.empty () ) {
					ParameterList p2 ( cookieParams, false, false );
					const_cast <ParameterList*>(pv [0])->appendParameters ( p2 );
				}
				else {
					string link_search_type = pv [0]->getStringValue ( "link_search_type" );
					bool crosslinking = link_search_type != "No Link" && link_search_type != "";
					if ( crosslinking ) {
						ParameterList p2 ( getDefaultResultsPath ( "searchCompare", "default_xl.xml" ), false, false, false, false );
						const_cast <ParameterList*>(pv [0])->appendParameters ( p2 );
					}
					else {
						ParameterList p2 ( getDefaultResultsPath ( "searchCompare" ), false, false, false, false );
						const_cast <ParameterList*>(pv [0])->appendParameters ( p2 );
					}
				}
				if ( paramList.getBoolValue ( "calibrate" ) )
					f = new SearchCompareFormCalibrationForm ( pv );
				else
					f = new SearchCompareForm ( pv );
			}
			else if ( form == "write_cal" ) {
				init_html ( os, "Create Calibrated Project" );
				string s = ProjectFile::updateProjectFile ( &paramList );
				std::cout << s << "<br />" << std::endl;
			}
			else if ( form == "results_management" ) {
				string action = paramList.getStringValue ( "action" );
				if ( action == "Delete" ) {
					deleteData ( paramList );
					f = new DeleteSelForm ( &paramList, false );
				}
				else if ( action == "Import" ) {
					init_html ( os, "Imported Projects" );
					try {
						ImportProject iProject ( paramList );
						ErrorHandler::genError ()->message ( "Project initialized.\n" );
						iProject.write ( os );
					}
					catch ( runtime_error e ) {
						ErrorHandler::genError ()->error ( e );
					}
					printHTMLFooter ( os, "2012" );
				}
				else if ( action == "Compress" || action == "Uncompress" ) {
					StringVector pr = paramList.getStringVectorValue ( "project_name" );
					if ( !pr.empty () ) {
						if ( action == "Compress" ) {
							init_html ( os, "Compress Projects" );
							compressProjects ( paramList );
						}
						else {
							init_html ( os, "Uncompress Projects" );
							uncompressProjects ( paramList );
						}
						ParameterList pL ( "msform", false, false, false, false, false );
						pL.appendParameters ( paramList );
						pL.addName ( "select_project", "1" );
						pL.removeName ( "project_name" );
						refreshJavascript ( os, 500, pL.getURL (), true );
						printHTMLFooter ( os, "2012" );
					}
					else {
						f = new DeleteSelForm ( &paramList, false );
					}
				}
				else if ( action == "Check" ) {
					StringVector pr = paramList.getStringVectorValue ( "project_name" );
					if ( !pr.empty () ) {
						init_html ( os, "Check Projects" );
						checkProjects ( paramList );
						printHTMLFooter ( os, "2012" );
					}
					else {
						f = new DeleteSelForm ( &paramList, false );
					}
				}
				else {										// Export
					StringVector pr = paramList.getStringVectorValue ( "project_name" );
					if ( !pr.empty () ) {
						init_html ( os, "Exported Projects" );
						exportData ( paramList );
						printHTMLFooter ( os, "2012" );
					}
					else {
						f = new DeleteSelForm ( &paramList, false );
					}
				}
			}
			else if ( form == "batchtagweb" ) {
				bool setCookie = false;
				if ( checkUser ( paramList, setCookie ) ) {
					pv.push_back ( getDefaultParameters ( "batchtag" ) );
					const_cast <ParameterList*>(pv [0])->addName ( "user", paramList.getStringValue ( "user" ) );
					f = new BatchTagWebForm ( pv, setCookie );
				}
				else
					f = new LoginForm ( &paramList );
			}
#endif
			else if ( form == "msseq" ) {
				pv.push_back ( getDefaultParameters ( "msseq" ) );
				f = new MSSeqForm ( pv );
			}
#ifdef MYSQL_DATABASE
			else if ( form == "add_user" ) {
				if ( paramList.getBoolValue ( "submit" ) ) {
					addUserToDatabase ( os, paramList );
				}
				else
					f = new AddUserForm ();
			}
#endif
			else if ( form == "mstagstandard" ) {
				pv.push_back ( getDefaultParameters ( "mstag" ) );
				f = new MSTagStandardForm ( pv );
			}
			else if ( form == "msfitstandard" ) {
				pv.push_back ( getDefaultParameters ( "msfit" ) );
				f = new MSFitStandardForm ( pv );
			}
			else if ( form == "msfitupload" ) {
				pv.push_back ( getDefaultParameters ( "msfitupload" ) );
				f = new MSFitUploadForm ( pv );
			}
			else if ( form == "dbstat" ) {
				pv.push_back ( getDefaultParameters ( "dbstat" ) );
				f = new DBStatForm ( pv );
			}
			else if ( form == "mspattern" ) {
				pv.push_back ( getDefaultParameters ( "mspattern" ) );
				f = new MSPatternForm ( pv );
			}
			else if ( form == "mshomology" ) {
				pv.push_back ( getDefaultParameters ( "mshomology" ) );
				f = new MSHomologyForm ( pv );
			}
			else if ( form == "msisotope" ) {
				pv.push_back ( getDefaultParameters ( "msisotope" ) );
				f = new MSIsotopeForm ( pv );
			}
			else if ( form == "msbridgestandard" ) {
				pv.push_back ( getDefaultParameters ( "msbridge" ) );
				f = new MSBridgeStandardForm ( pv );
			}
			else if ( form == "msdigest" ) {
				pv.push_back ( getDefaultParameters ( "msdigest" ) );
				f = new MSDigestForm ( pv );
			}
			else if ( form == "msviewer" ) {
				pv.push_back ( getDefaultParameters ( "msviewer" ) );
				f = new MSViewerForm ( pv );
			}
			else if ( form == "msfilter" ) {
				pv.push_back ( getDefaultParameters ( "msfilter" ) );
				f = new MSFilterForm ( pv );
			}
			else if ( form == "msproduct" ) {
				pv.push_back ( getDefaultParameters ( "msproduct" ) );
				f = new MSProductPasteForm ( pv );
			}
			//else if ( form == "msproductupload" ) {
			//	f = new MSProductUploadForm ( pv );
			//}
			else if ( form == "msnonspecific" ) {
				pv.push_back ( getDefaultParameters ( "msnonspecific" ) );
				f = new MSNonSpecificForm ( pv );
			}
			else if ( form == "mscomp" ) {
				pv.push_back ( getDefaultParameters ( "mscomp" ) );
				f = new MSCompForm ( pv );
			}
			else if ( form == "faindex" ) {
				f = new FAIndexForm ( pv );
			}
			else if ( form == "makeproject" ) {
				if ( paramList.getBoolValue ( "write_file" ) ) {
					makeProject ( os, paramList );
				}
#ifdef MYSQL_DATABASE
				else
					f = new MakeProjectForm ( &paramList );
#endif
			}
#ifdef MYSQL_DATABASE
			else if ( form == "delete_search_compare_cookie" ) {
				init_html ( os, "Delete Search Compare Cookie" );
				string p = ParameterList::getCookieValue ( "search_compare_params" );
				if ( !p.empty () ) deleteCookie ( os, "search_compare_params" );
				refreshJavascript ( os, 0, "../mshome.htm", true );
				return 0;
			}
			else if ( form == "delete_batchtag_cookie" ) {
				init_html ( os, "Delete Batch-Tag Cookie" );
				string p = ParameterList::getCookieValue ( "batchtag_params" );
				if ( !p.empty () ) deleteCookie ( os, "batchtag_params" );
				refreshJavascript ( os, 0, "../mshome.htm", true );
				return 0;
			}
			else if ( form == "logout" ) {
				init_html ( os, "Logout" );
				string key = ParameterList::getCookieValue ( "key" );
				if ( !key.empty () ) {
					deleteCookie ( os, "key" );
					MySQLPPSDDBase::instance ().deleteSession ( key );
				}
				refreshJavascript ( os, 0, "../mshome.htm", true );
				return 0;
			}
			else if ( form == "importproject" ) {
				importProject ( os, paramList );
			}
#endif
		}
		if ( f ) f->printHTML ( os );
	}
	catch ( runtime_error e ) {
		paramList.writeLogError ( e.what () );
	}
#ifdef MYSQL_DATABASE
	if ( !mySQLLoginFault ) {
		MySQLPPSDDBase::instance ( false, true );
	}
#endif
	return 0;
}
namespace {
#ifdef MYSQL_DATABASE
void importProject ( ostream& os, const ParameterList& paramList )
{
	init_html_premature_stop ( "Import Project", true );
	string projectName = paramList.getStringValue ( "project_name" );
	string user = paramList.getStringValue ( "user" );
	string uploadFpath = paramList.getStringValue ( "upload_project_filepath" );	// This is the full path name of the uploaded file as named by PP
	string uploadFname = paramList.getStringValue ( "upload_project_filename" );	// This is the original file name.
	if ( !isFileType ( uploadFname, XML ) ) {		// File type not xml
		genUnlink ( uploadFpath );
		ErrorHandler::genError ()->error ( "Not a project file.\n" );
	}
	if ( uploadFname.rfind ( ".cal." ) != string::npos ) {
		genUnlink ( uploadFpath );
		ErrorHandler::genError ()->error ( "This is a calibrated project file. Please select the uncalibrated project\n" );
	}
	if ( user.empty () ) {
		genUnlink ( uploadFpath );
		ErrorHandler::genError ()->error ( "User not specified.\n" );
	}
	UserInfo* userInfo = MySQLPPSDDBase::instance ().getUserInfo ( user );
	if ( !userInfo ) {													// Get the user information
		genUnlink ( uploadFpath );
		ErrorHandler::genError ()->error ( "Unknown user.\n" );
	}
	string userID = userInfo->getUserID ();
	if ( projectName.empty () ) {
		projectName = uploadFname.substr ( 0, uploadFname.length () - 4 );
	}
	if ( MySQLPPSDDBase::instance ().checkProject ( userID, projectName ) ) {
		genUnlink ( uploadFpath );
		ErrorHandler::genError ()->error ( "Project already exists.\n" );
	}
	StringVector sv;
	char* info = getFileInfo ( uploadFpath, '\n', 1, false );
	StringVector files = XMLParser::getStringVectorValue ( info, "centroid" );
	GenNameValueStream nvs ( MsparamsDir::instance ().getParamPath ( "inst_dir.txt" ) );
	string instrument;
	string oldInst;
	for ( StringVectorSizeType i = 0 ; i < files.size () ; i++ ) {
		if ( files [i][0] == '$' ) {
			sv.push_back ( files [i].substr ( 1 ) );
			int end = sv.back ().find ( "/" );
			string instDir = sv.back ().substr ( 0, end );
			if ( nvs.getValue ( instDir, instrument ) == false ) {
				genUnlink ( uploadFpath );
				ErrorHandler::genError ()->error ( "Instrument directory " + instDir + " not in inst_dir.txt parameter file.\n" );
			}
			if ( i != 0 ) {
				if ( instrument != oldInst ) {
					genUnlink ( uploadFpath );
					ErrorHandler::genError ()->error ( "Files from different instrument types can't be used in the same project.\n" );
				}
			}
			else oldInst = instrument;
		}
	}
	delete [] info;
	if ( files.empty () || ( sv.size () != files.size () ) ) {
		genUnlink ( uploadFpath );
		ErrorHandler::genError ()->error ( "Project file cannot be imported.\n" );
	}
	PPProject pf ( userInfo->getMaxMSMSSpectra (), userInfo->getAllowRaw (), projectName, sv );
	if ( pf.initialised () ) {
		UserRepository ur ( "batchtag", userInfo->getDirectoryName () );
		pf.createProjectFile ( ur.getFullProjectPath () + SLASH + projectName + ".xml" );
		MySQLPPSDDBase::instance ().submitProject ( userInfo->getUserID (), projectName, projectName + ".xml", ur.getProjectPath (), instrument );
	}
	genUnlink ( uploadFpath );
	printAbortFunctions ( os );
	if ( !pf.initialised () ) {
		ErrorHandler::genError ()->error ( "One of the selected files has an invalid format.\n" );
	}
	ParameterList pList ( "msform", false, false, false, false, false );
	pList.addName ( "form", "batchtag" );
	pList.addName ( "load_from", "Parameter File" );
	pList.addName ( "params_file", "default" );
	pList.addName ( "project_name", projectName );
	pList.addName ( "instrument_name", instrument );
	pList.addName ( "user", user );
	refreshJavascript ( os, 1000, pList.getURL (), true );

	os << "<hr />" << endl;
	printHTMLFooter ( os, "2008" );
}
void addUserToDatabase ( ostream& os, const ParameterList& paramList )
{
	init_html ( os, "Add User" );
	string user = paramList.getStringValue ( "user" );
	string email = paramList.getStringValue ( "email" );
	string firstName = paramList.getStringValue ( "first_name" );
	string lastName = paramList.getStringValue ( "last_name" );
	string password = paramList.getStringValue ( "password" );
	string retypePassword = paramList.getStringValue ( "retype_password" );
	os << "<br />" << endl;
	string directoryName = genToLower ( genRandomString ( 10 ) );	// Selects a unique directory name
	genSleep ( 1500 );								// Sleep a bit so a different random string produced
	if ( ParameterList::getMethod () != "post" ) {
		ErrorHandler::genError ()->error ( "User not added.\n" );
	}
	if ( !isValidUser ( user, 3, 10 ) ) {
		ErrorHandler::genError ()->error ( "Invalid user.\nUser names must be between 3 and 10 characters in length.\nThe first character must be a lower case letter.\nSubsequent characters may be lower case letters or numbers.\n" );
	}
	if ( !isValidEmailAddress ( email.c_str () ) ) {
		ErrorHandler::genError ()->error ( "Invalid email address.\n" );
	}
	if ( password != retypePassword ) {
		ErrorHandler::genError ()->error ( "The passwords don't match.\n" );
	}
	unsigned int minPasswordLength = InfoParams::instance ().getUIntValue ( "min_password_length" );
	if ( password.length () < minPasswordLength ) {
		ErrorHandler::genError ()->error ( "The password is too short. The minimum length is " + gen_itoa ( minPasswordLength ) + ".\n" );
	}
	if ( firstName.length () != 0 && !isValidName ( firstName ) ) {
		ErrorHandler::genError ()->error ( "Invalid first name.\n" );
	}
	if ( lastName.length () != 0 && !isValidName ( lastName ) ) {
		ErrorHandler::genError ()->error ( "Invalid last name.\n" );
	}
	os << "<p>" << endl;
	if ( MySQLPPSDDBase::instance ().submitUser ( user, password, email, directoryName, firstName, lastName ) ) {
		os << "User successfully added." << endl;
	}
	else {
		os << "User not added. Try using a different user name." << endl;
	}
	os << "</p>" << endl;
	printHTMLFooter ( os, "2007" );
}
string getUserFromCookieKey ()
{
	string key = ParameterList::getCookieValue ( "key" );
	string user;
	if ( !key.empty () ) {
		user = MySQLPPSDDBase::instance ().getUserName ( key );
	}
	return user;
}
bool checkGuest ( const string& form, const UserInfo* userInfo )
{
	if ( userInfo->getIsGuest () ) {
		if ( form == "batchtag" )		return false;	// guest user can't do Batch-Tag searches
		if ( form == "batchtagweb" )	return false;
		if ( form == "write_cal" )		return false;
		if ( form == "results_management" )	return false;
		if ( form == "makeproject" )		return false;
		if ( form == "importproject" )		return false;
	}
	return true;
}
bool checkUser ( ParameterList& paramList, bool& setCookie )
{
	string user = getUserFromCookieKey ();
	string form = paramList.getStringValue ( "form" );
	if ( user.empty () ) {	// Covers the case of user being set by login form
		user = paramList.getStringValue ( "user" );
		string password = paramList.getStringValue ( "password" );
		string form = paramList.getStringValue ( "form" );
		UserInfo* userInfo = MySQLPPSDDBase::instance ().getUserInfo ( user );
		if ( !userInfo ) return false;
		if ( userInfo->getPassword () != password ) return false;
		if ( user == "root" ) {
			if ( form == "batchtag" ) return false;			// root user can't do Batch-Tag searches
			if ( form == "batchtagweb" ) return false;
		}
		if ( !checkGuest ( form, userInfo ) ) return false;
		string newPassword = paramList.getStringValue ( "new_password" );
		string retypeNewPassword = paramList.getStringValue ( "retype_new_password" );
		if ( !newPassword.empty () ) {
			if ( newPassword == retypeNewPassword ) {
				if ( userInfo->getIsGuest () )
					ErrorHandler::genError ()->error ( "Guest user can't set password.\n" );
				else {
					unsigned int minPasswordLength = InfoParams::instance ().getUIntValue ( "min_password_length" );
					if ( newPassword.length () < minPasswordLength ) {
						ErrorHandler::genError ()->error ( "The new password is too short. The minimum length is " + gen_itoa ( minPasswordLength ) + ".\n" );
					}
					else {
						MySQLPPSDDBase::instance ().setPassword ( user, newPassword );
					}
				}
			}
		}
		setCookie = true;
	}
	else {
		paramList.addName ( "user", user );
		UserInfo* userInfo = MySQLPPSDDBase::instance ().getUserInfo ( user );
		if ( !checkGuest ( form, userInfo ) ) return false;
	}
	return true;
}
void deleteData ( const ParameterList& paramList )
{
	string user = paramList.getStringValue ( "user" );
	StringVector results = paramList.getStringVectorValue ( "results_file" );
	for ( StringVectorSizeType i = 0 ; i < results.size () ; i++ ) {
		StringVector sv = genGetSubstrings ( results [i], '/' );
		if ( sv.size () == 3 )
			deleteResults ( sv [0], sv [1], sv [2] );
		else
			deleteResults ( user, sv [0], sv [1] );
	}
	StringVector projects = paramList.getStringVectorValue ( "project_name" );
	for ( StringVectorSizeType j = 0 ; j < projects.size () ; j++ ) {
		StringVector sv = genGetSubstrings ( projects [j], '/' );
		bool flag;
		if ( sv.size () == 2 )
			flag = deleteProject ( sv [0], sv [1] );
		else
			flag = deleteProject ( user, sv [0] );
		if ( !flag ) {
			ErrorHandler::genError ()->error ( "The project " + sv.back () + " has submitted or running jobs so cannot be deleted.\n" );
		}
	}
}
#endif
void makeProject ( ostream& os, const ParameterList& paramList )
{
	string projectName = paramList.getStringValue ( "project_name" );
	string user = paramList.getStringValue ( "user" );
	if ( user.empty () ) {
		string peakListFilepath = paramList.getStringValue ( "peak_list_filepath" );
		PPProject* pf = 0;
		if ( !peakListFilepath.empty () ) {
			pf = new PPProject ( projectName, peakListFilepath );
		}
		else {

			StringVector dataFileList = paramList.getStringVectorValue ( "data_file_list" );
			if ( !dataFileList.empty () ) {
				pf = new PPProject ( projectName, dataFileList );
			}
		}
		if ( pf ) {
			if ( pf->initialised () ) {
				string projectFilepath = paramList.getStringValue ( "project_filepath" );
				pf->createProjectFile ( projectFilepath + SLASH + projectName + ".xml", true );
			}
			delete pf;
		}
	}
#ifdef MYSQL_DATABASE
	else {
		init_html_premature_stop ( "Make Project", true );
		StringVector sv = paramList.getStringVectorValue ( "data_file_list" );
		os << "<br />" << endl;
		if ( sv.empty () ) {
			ErrorHandler::genError ()->error ( "No centroid files selected.\n" );
		}
		if ( projectName.empty () ) {
			ErrorHandler::genError ()->error ( "Project name not specified.\n" );
		}
		if ( user.empty () ) {
			ErrorHandler::genError ()->error ( "User not specified.\n" );
		}
		UserInfo* userInfo = MySQLPPSDDBase::instance ().getUserInfo ( user );
		if ( !userInfo ) {													// Get the user information
			ErrorHandler::genError ()->error ( "Unknown user.\n" );
		}
		string userID = userInfo->getUserID ();
		if ( MySQLPPSDDBase::instance ().checkProject ( userID, projectName ) ) {
			ErrorHandler::genError ()->error ( "Project already exists.\n" );
		}
		string instrument;
		MapStringToStringVector reposParams;
		try {
			instrument = RepositoryInfo::exists ()	? RepositoryInfo::instance ().getInstrumentType ( sv )
													: InstrumentDir::instance ().getInstrumentType ( sv );
			if ( RepositoryInfo::exists () ) {
				reposParams = RepositoryInfo::instance ().getInstrumentParams ( sv );
			}
		}
		catch ( runtime_error e ) {
			ErrorHandler::genError ()->error ( e );
		}
		PPProject pf ( userInfo->getMaxMSMSSpectra (), userInfo->getAllowRaw (), projectName, sv );
		if ( pf.initialised () ) {
			UserRepository ur ( "batchtag", userInfo->getDirectoryName () );
			pf.createProjectFile ( ur.getFullProjectPath () + SLASH + projectName + ".xml" );
			MySQLPPSDDBase::instance ().submitProject ( userInfo->getUserID (), projectName, projectName + ".xml", ur.getProjectPath (), instrument );
		}
		printAbortFunctions ( os );
		if ( !pf.initialised () ) {
			ErrorHandler::genError ()->error ( "One of the selected files has an invalid format.\n" );
		}

		ParameterList pList ( "msform", false, false, false, false, false );
		pList.addName ( "form", "batchtag" );
		pList.addName ( "load_from", "Parameter File" );
		pList.addName ( "params_file", "default" );
		pList.addName ( "project_name", projectName );
		pList.addName ( "instrument_name", instrument );
		pList.addName ( "user", user );
		for ( MapStringToStringVectorConstIterator i = reposParams.begin () ; i != reposParams.end () ; i++ ) {
			for ( StringVectorSizeType j = 0 ; j < (*i).second.size () ; j++ ) {
				pList.addName ( (*i).first, (*i).second [j] );
			}
		}
		refreshJavascript ( os, 2000, pList.getURL (), true );

		os << "<hr />" << endl;
		printHTMLFooter ( os, "2005" );
	}
#endif
}
}
