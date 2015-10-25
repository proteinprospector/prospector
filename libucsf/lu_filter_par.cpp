/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_filter_par.cpp                                             *
*                                                                             *
*  Created    : September 3rd 2012                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2012-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
/*
Parameters ( <> required; [] optional; S == string ):
  [--oglcnac]     Use --loss=204.08 --tolerance=0.03 --tolerance-type='Da'
  [--hexnac]      Use --loss=203.08 --tolerance=0.6 --tolerance-type='Da'
  [--sialicacid]  Use --loss=291.09 --tolerance=0.6 --tolerance-type='Da'
  [--loss=F]      Neutral-loss mass to look for in the peaklist.
  [--tolerance=F] Mass tolerance to use when matching peaks.
  [--tolerance-type=S]
                  Type of mass tolerance to use, can be 'Da', 'mmu', '%',
                  or 'ppm'. Defaults to 'Da'.
  [--threshold=F] Threshold filter for each MSMS peaklist when performing
                  peak matches.
  [--threshold-type=S]
                  The type of thresholding to perform. Possible values are
                  'peaks', 'basepeak' (as percent of basepeak), 'tic' (as
                  percent of tic). Defaults to 'peaks'.
  [--print]       Output the matched peaklists.
  [--print-filter]
                  Output a list of the MS where the program found a neutral
                  loss. Output a tab-delimited list of 'file', 'pepmass', 
                  'charge', and 'rt'.
  [--print-etd]   Output a list of the corresponding ETD MSMS for the neutral
                  loss peaklists found in the CID file. The ETD file must 
                  have the same name as the CID file with 'etd' in place of
                  the 'cid' in the filename. e.g. file_cid.mgf -> file_etd.mgf
  [--output=S]    Output data to a named file instead of STDOUT.
  [--log=S]       Output log messages to this file instead of STDERR.
  [--verbose]     Increase verbosity.
  [--version]     Print out version number, then quit.
  [--help]        Print out this message, then quit.
*/
#ifndef VIS_C
#include <stdexcept>
#endif
#include <lg_string.h>
#include <lgen_error.h>
#include <lgen_file.h>
#include <lgen_uncompress.h>
#include <lu_file_type.h>
#include <lu_filter_par.h>
#include <lu_getfil.h>
#include <lu_inst.h>
#include <lu_param_list.h>
#include <lu_pk_filter.h>
#include <lu_pq_vector.h>

using std::string;
using std::runtime_error;

MSFilterParameters::MSFilterParameters ( const ParameterList* params ) :
	MSProgramParameters	( params ),
	aaInitInfo		( params ),
	lowMPlusH		( params->getDoubleValue ( "low_m_plus_h", 600.0 ) ),
	highMPlusH		( params->getDoubleValue ( "high_m_plus_h", 4000.0 ) ),
	fullMPlusHRange ( params->getBoolValue ( "full_m_plus_h_range" ) ),
	chargeFilter	( params->getStringVectorValue ( "charge_filter" ) ),
	allCharges		( params->getBoolValue ( "all_charges" ) ),
	lossFormula		( params->getStringValue ( "loss_composition" ) ),
	fragmentMZ		( params->getDoubleValue ( "fragment_mz" ) ),
	minMatches		( params->getIntValue ( "min_frag_matches" ) ),
	parentMassTolerance	( "msms_parent_mass", params ),
	productMassTolerance( "fragment_masses", params ),
	systematicError	( params->getDoubleValue ( "msms_parent_mass_systematic_error" ) ),
	instrumentName	( params->getStringValue ( "instrument_name" ) ),
	keepOrRemove	( initKeepOrRemove ( params->getStringValue ( "keep_or_remove" ) ) )
{
	initialise_amino_acid_weights ( MapStringConstModPtr (), ElementalFormulaVector (), true );
	initPeakList ( params );
	massInfo = new MassInfo ( params );
	msmsPeakFilterOptions = new MSMSPeakFilterOptions ( params );
	initialiseInstrumentName ( instrumentName );
	try {
		const char* value;
		if ( params->getValue ( "fragment_mzs", value ) )	getPostQueryVector ( value, fragmentMZs, '\n' );
	}
	catch ( runtime_error ) {
		ErrorHandler::genError ()->error ( "The fragment m/z field is incorrectly formatted. There should be a one m/z value per line.\n" );
	}
}
MSFilterParameters::~MSFilterParameters ()
{
	delete msmsPeakFilterOptions;
	delete massInfo;
}
string MSFilterParameters::initKeepOrRemove ( const string& s )
{
	if		( s == "Keep Spectra Matching Criteria" )	return "keep";
	else if ( s == "Remove Spectra Matching Criteria" )	return "remove";
	else												return "both";
}
void MSFilterParameters::initPeakList ( const ParameterList* params )
{
	string uploadPeakListFname = params->getStringValue ( "upload_temp_peak_list_filename" );
	string uploadPeakListFpath = params->getStringValue ( "upload_temp_peak_list_filepath" );
	archiveName = genShortFilenameFromPath ( uploadPeakListFname );
	if ( uploadPeakListFname.empty () ) {
		return;
	}
	uploadPeakListFname = genFilenameFromPath ( uploadPeakListFname );	// IE gives the full path whereas Mozilla give the filename (what we want)
	string uploadName = genPreprocessFile ( uploadPeakListFpath );
	if ( uploadName == uploadPeakListFpath || isCompressedUpload ( uploadName ) ) {		// If the file hasn't been preprocessed or it has just been uncompressed it must be a single file
		string shortName = genShortFilenameFromPath ( uploadName );
		string newDir = genDirectoryFromPath ( uploadName ) + SLASH + shortName;
		if ( newDir == uploadName ) {
			newDir += "_1";
		}
		bool ret = genCreateDirectory ( newDir );
		string newUploadName;
		if ( isCompressedUpload ( uploadName ) )
			newUploadName = newDir + SLASH + genShortFilenameFromPath ( uploadPeakListFname );
		else
			newUploadName = newDir + SLASH + uploadPeakListFname;
		genRename ( uploadName, newUploadName );
		uploadName = newDir;
	}
	PPTempFile pptf ( "", "" );					// Move the directory to the temporary directory
	string tempFileFullPathInput = pptf.getFullPathDir () + SLASH + genToLower ( genRandomString ( 10 ) );
	genCreateNewDirectory ( tempFileFullPathInput );
	inPeakListFpath = tempFileFullPathInput + SLASH + genShortFilenameFromPath ( uploadPeakListFname );
	genRename ( uploadName, inPeakListFpath );

	string randomDir = genToLower ( genRandomString ( 10 ) );
	string projName = genShortFilenameFromPath ( uploadPeakListFname );
	string tempFileFullPathOutput = pptf.getFullPathDir () + SLASH + randomDir;
	genCreateNewDirectory ( tempFileFullPathOutput );
	tempFileFullPathOutput += SLASH + projName;
	if ( keepOrRemove == "both" ) {
		genCreateNewDirectory ( tempFileFullPathOutput + "-matching" );
		genCreateNewDirectory ( tempFileFullPathOutput + "-non-matching" );
	}
	else {
		genCreateNewDirectory ( tempFileFullPathOutput );
	}
	outPeakListFpath = tempFileFullPathOutput;

	outPeakListURL = genDirectoryFromPath ( pptf.getURL () ) + "/" + randomDir + "/" + projName;
}
