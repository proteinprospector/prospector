/******************************************************************************
*                                                                             *
*  Program    : faindex                                                       *
*                                                                             *
*  Filename   : faindex_main.cpp                                              *
*                                                                             *
*  Created    : May 11th 1995                                                 *
*                                                                             *
*  Purpose    : A tool for managing FASTA formatted sequence databases and    *
*               indexing them for use with programs in the ProteinProspector  *
*               package.                                                      *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1995-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifndef VIS_C
#include <stdexcept>
#endif
#include <lu_html.h>
#include <lu_html_form.h>
#include <lu_acc_link.h>
#include <lu_faindex.h>
#include <lu_version.h>
#include <lu_faind_par.h>
#include <lu_dig_par.h>
#include <lu_hit.h>
#include <lu_fasta.h>
#include <lu_getfil.h>
#include <lu_param_list.h>
#include <lu_table.h>
using std::runtime_error;
using std::string;
using std::ostream;
using std::cout;
using std::endl;

static void faindexWebBrowser ( int argc, char** argv );
static void faindexCommandLine ( const string& database );
static void displayDatabase ( ostream& os, const FAIndexDatabaseSummaryParameters& params, const string& database );

int main ( int argc, char *argv [] )
{
	if ( argc != 2 ) faindexWebBrowser ( argc, argv );
	else faindexCommandLine ( argv [1] );
	return 0;
}
static void faindexWebBrowser ( int argc, char** argv )
{
	initialiseProspector ();
	ParameterList paramList ( argc, argv );
	try {
		if ( paramList.empty () ) {
			ErrorHandler::genError ()->error ( "No parameters passed to FA-Index.\n" );
		}
		//string ver = paramList.getStringValue ( "version" );
		//if ( !ver.empty () && ver != Version::instance ().getVersion () ) {
		//	ErrorHandler::genError ()->error ( "Form version mismatch. Try reloading the form.\n" );
		//}
		ProgramLink::setParams ( &paramList );

		ostream& os = cout;
		init_html ( os, "FA-Index" );

		if ( !paramList.getBoolValue ( "create_database_indicies" ) ) {
			FAIndexDatabaseSummaryParameters params ( &paramList );
			displayDatabase ( os, params, params.getDatabase () );
		}
		else {
			string createdDatabaseName;
			if ( paramList.getBoolValue ( "create_sub_database" ) ) {
				FAIndexSubsetDatabaseParameters params ( &paramList );
				string database = params.getDatabase ();
				FastaServer fs ( database );
				string subDatabaseID = "." + params.getSubDatabaseID ();
				createdDatabaseName = database + subDatabaseID;
				create_new_database_from_indicies_list ( &fs, database, subDatabaseID, params.getIndicies ( &fs ) );
			}
			else if ( paramList.getBoolValue ( "create_user_database" ) ) {
				FAIndexCreateOrAppendParameters params ( &paramList );
				createdDatabaseName = params.getDatabase ();
				add_single_database_entry ( createdDatabaseName, params.getUserProteinSequence (), params.getNameField (), params.getSpecies (), params.getAccessionNum () );
			}
			else {
				FAIndexNormalParameters params ( &paramList );
				if ( params.getDNAToProtein () )
					createdDatabaseName = dna_database_to_protein_database ( params.getDatabase (), params.getDeleteDNADatabase () );
				else if ( params.getRandomDatabase () )
					createdDatabaseName = convertDatabaseToRandomDatabase ( params.getDatabase (), false, params.getConcatDatabase () );
				else if ( params.getReverseDatabase () )
					createdDatabaseName = convertDatabaseToRandomDatabase ( params.getDatabase (), true, params.getConcatDatabase () );
				else
					createdDatabaseName = params.getDatabase ();
			}
			create_database_files ( createdDatabaseName, InfoParams::instance ().getBoolValue ( "faindex_parallel", false ) );
			os << "<br />" << endl;
			os << "<b>Indexed Database:</b> " << createdDatabaseName << "<br />" << endl;
			os << "<br />" << endl;
			StringVector databaseList = SeqdbDir::instance ().getDatabaseList ( true );
			os << "<b>Updated List of Available Databases</b><p />";
			for ( StringVectorSizeType i = 0 ; i < databaseList.size () ; i++ ) {
				os << databaseList [i] << "<br />" << endl;
			}
			os << "<p />" << endl;
			os << "<b>You must Reload Form to see the new database choices.</b><br />" << endl;
			os << "<b>If the reload button doesn't work, place the cursor in the URL location box of the browser and press return.</b><br />" << endl;
		}
		printProgramInformationHTML ( os, "FA-Index" );
	}
	catch ( runtime_error e ) {
		paramList.writeLogError ( e.what () );
	}
}
static void faindexCommandLine ( const string& database )
{
	create_database_files ( database, InfoParams::instance ().getBoolValue ( "faindex_parallel", false ) );
	cout << endl;
	cout << "Indexed New Database: " << database << endl << endl;
	printProgramInformation ( "FA-Index" );
}
static void displayDatabase ( ostream& os, const FAIndexDatabaseSummaryParameters& params, const string& database )
{
	int dnaReadingFrame = params.getDNAReadingFrame ();
	int startIndex = params.getStartIndexNumber ();
	int endIndex = params.getEndIndexNumber ();
	bool allIndicies = params.getAllIndicies ();
	bool showProteinSequences = !params.getHideProteinSequence ();
	FastaServer fs ( database );
	ProteinHit::addFS ( &fs, 0 );
	int numEntries = fs.getNumEntries ();

	if ( allIndicies ) {
		startIndex = 1;
		endIndex = numEntries;
	}
	else {
		startIndex = genMin ( startIndex, numEntries );
		endIndex = genMin ( endIndex, numEntries );
		endIndex = genMax ( endIndex, 1 );
	}

	startJavascript ( os );
	MSDigestLink digestLink ( "faindex" );
	digestLink.printHTML ( os );
	AccessionNumberLinkInfo anli;
	anli.printHTML ( os, database );
	endJavascript ( os );

	os << "\n\n<div class=\"header1\">Database Contents of " << database << "</div>";
	if ( fs.getDNADatabase () ) os << "\n\n<b>DNA Reading Frame " << dnaReadingFrame << "</b><p />";
	os << "\n<table>\n";

	for ( int i = startIndex ; i <= endIndex ; i++ ) {

		ProteinHit hit ( &fs, i, dnaReadingFrame, 0 );
		if ( i == startIndex ) {
			tableRowStart ( os );
				hit.printHTMLHeader3 ( os );
			tableRowEnd ( os );
		}
		hit.printHTMLMultiHit ( os );

		char* protein = fs.get_fasta_protein ( i, dnaReadingFrame );
		if ( showProteinSequences ) {
			int len = strlen ( protein );

			tableRowStart ( os );
				tableDataStart ( os, "", "left", false, 4 );
					os << "<tt>" << endl;
					for ( int j = 0 ; j < len ; j += 80 ) {
						for ( int k = 0 ; k < 80 ; k++ ) {
							if ( j + k < len ) os << protein [j + k];
						}
						os << "<br />" << endl;
					}
					os << "</tt>" << endl;
				tableDataEnd ( os );
			tableRowEnd ( os );
		}
	}
	os << "</table>" << endl;
	os << "<p />" << endl << endl;

	printHTMLFORMStart ( os, "post", "faindex" );

	printHTMLFORMSubmit ( os, "Next 20 Entries" );
	os << "<br />" << endl;

	printHTMLFORMHidden ( os, "database", database );
	printHTMLFORMHidden ( os, "start_index_number", endIndex + 1 );
	printHTMLFORMHidden ( os, "end_index_number", endIndex + 20 );
	printHTMLFORMHidden ( os, "all_indicies", 0 );
	printHTMLFORMHidden ( os, "dna_reading_frame", dnaReadingFrame );
	printHTMLFORMHidden ( os, "hide_protein_sequence", ( showProteinSequences ? 0 : 1 ) );

	os << "</form>" << endl << endl;
}
