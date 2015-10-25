/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_file_type.cpp                                              *
*                                                                             *
*  Created    : April 26th 2007                                               *
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
#include <lg_io.h>
#include <lg_string.h>
#include <lgen_file.h>
#include <lu_file_type.h>
#include <lu_mgf.h>
using std::ios;
using std::string;
using std::istringstream;
using std::getline;
using namespace FileTypes;

namespace {
unsigned char microsoftCompoundDocumentFormatPattern [] = { 0320, 0317, 021, 0340, 0241, 0261, 032, 0341 };
unsigned char pkZIPFile [] = { 'P', 'K', 003, 004 };

bool compareMagicNumber ( const string& fullPath, const unsigned char* pattern, int num )	// Used if some bytes have zero
{
	GenIFStream ist ( fullPath, ios::binary );
	unsigned char* fpattern = new unsigned char [num];
	ist.read ( (char*)fpattern, num );
	if ( ist.gcount () == num ) {
		return !memcmp ( (void*) pattern, (void*) fpattern, num );
	}
	return false;
}
bool compareMicrosoftCompoundDocumentFormatPattern ( const string& fullPath )
{
	return compareMagicNumber ( fullPath, microsoftCompoundDocumentFormatPattern, 8 ); 
}
bool comparePKZIPPattern ( const string& fullPath )
{
	return compareMagicNumber ( fullPath, pkZIPFile, 4 ); 
}
bool isPrecursorMPlusHChargeLine ( const string& line )
{
	string value1;
	istringstream istLine ( line );
	if ( istLine >> value1 ) {
		if ( !genStringIsPositiveNonZeroFloat ( value1 ) ) return false;
		double mPlusH = atof ( value1.c_str () );
		if ( mPlusH < 50.0 || mPlusH > 1000000.0 ) return false;
		string value2;
		if ( istLine >> value2 ) {
			if ( !genStringIsInteger ( value2 ) ) return false;
			int charge = atoi ( value2.c_str () );
			if ( charge < -100 || charge > 100 ) return false;
			string temp;
			if ( istLine >> temp ) return false;	// Too many fields
			else return true;
		}
		else return false;
	}
	else return false;
}
bool isPrecursorMOverZIntensityChargeLine ( const string& line )
{
	string value1;
	istringstream istLine ( line );
	if ( istLine >> value1 ) {
		if ( !genStringIsPositiveNonZeroFloat ( value1 ) ) return false;
		double mOverZ = atof ( value1.c_str () );
		if ( mOverZ < 50.0 || mOverZ > 10000.0 ) return false;
		string value2;
		if ( istLine >> value2 ) {
			if ( !genStringIsFloat ( value2 ) ) return false;
			string value3;
			if ( istLine >> value3 ) {		// Sometimes charges are written as floating point numbers
				int charge;
				if ( genStringIsInteger ( value3 ) )	charge = atoi ( value3.c_str () );
				else if ( genStringIsFloat ( value3 ) )	charge = static_cast <int> ( atof ( value3.c_str () ) );
				else return false;
				if ( charge < -100 || charge > 100 ) return false;
				string temp;
				if ( istLine >> temp ) return false;	// Too many fields
				else return true;
			}
			else return false;
		}
		else return false;
	}
	else return false;
}
bool isFragmentMOverZIntensityLine ( const string& line )
{
	string value1;
	istringstream istLine ( line );
	if ( istLine >> value1 ) {
		if ( !genStringIsPositiveNonZeroFloat ( value1 ) ) return false;
		double mOverZ = atof ( value1.c_str () );
		if ( mOverZ > 1000000.0 ) return false;
		string value2;
		if ( istLine >> value2 ) {
			if ( !genStringIsFloat ( value2 ) ) return false;
			string temp;
			if ( istLine >> temp ) return false;	// Too many fields
			else return true;
		}
		else return false;
	}
	else return false;
}
}
bool isMGFSpottingPlateFile ( const string& fullPath )
{
	MGFInfo::instance ().reset ();
	GenIFStream ist ( fullPath );
	string line;
	int i = 0;
	bool spottingPlate = false;
	while ( getline ( ist, line ) ) {
		if ( MGFInfo::instance ().isSpottingPlate ( line, spottingPlate ) ) {
			break;
		}
		if ( i++ == 100 ) break;
	}
	return spottingPlate;
}
bool isMGFFile ( const string& fullPath )
{
	GenIFStream ist ( fullPath );
	string line;
	int i = 0;
	while ( getline ( ist, line ) ) {
		if ( isPrefix ( line, "BEGIN IONS" ) ) {				// isPrefix deals with the case of a DOS file on a UNIX system
			return true;
		}
		if ( i++ == 100 ) break;
	}
	return false;
}
bool isMGFFile2 ( const string& fullPath )	// This used for files with an mgf suffix. It allows files with no scans (just comments).
{
	GenIFStream ist ( fullPath );
	string line;
	int i = 0;
	while ( getline ( ist, line ) ) {
		if ( line.length () == 0 || line [0] == '#' ) continue;
		if ( isPrefix ( line, "BEGIN IONS" ) ) {				// isPrefix deals with the case of a DOS file on a UNIX system
			return true;
		}
		if ( i++ == 100 ) break;
	}
	if ( i == 0 ) return true;		// Empty file
	return false;
}
bool mgfFileHasZeros ( const string& fullPath )
{
	GenIFStream ist ( fullPath );
	string line;
	int i = 0;
	while ( getline ( ist, line ) ) {
		if ( line.length () == 0 || line [0] == '#' ) continue;
		if ( isdigit ( line [0] ) ) {							// This line has data
			char* ptr = const_cast <char*> (line.c_str ());
			double mass = atof ( ptr );
			while ( *ptr && isspace ( *ptr ) ) ptr++;	// Skip leading white space
			while ( *ptr && !isspace ( *ptr ) ) ptr++;	// Skip mass
			if ( atof ( ptr ) == 0.0 ) {
				return true;
			}
		}
		if ( i++ == 100 ) break;						// No zeros found after 100 lines. Assume there are none.
	}
	return false;
}
bool isMS2File ( const string& fullPath )
{
	GenIFStream ist ( fullPath );
	string line;
	int i = 0;
	while ( getline ( ist, line ) ) {
		if ( isPrefix ( line, "S\t" ) ) {				// isPrefix deals with the case of a DOS file on a UNIX system
			return true;
		}
		if ( i++ == 100 ) break;
	}
	return false;
}
bool isMS2File2 ( const string& fullPath )	// This used for files with an ms2 suffix. It allows files with no scans (just comments).
{
	GenIFStream ist ( fullPath );
	string line;
	int i = 0;
	while ( getline ( ist, line ) ) {
		if ( line.length () == 0 || line [0] == '#' ) continue;
		if ( isPrefix ( line, "S\t" ) ) {				// isPrefix deals with the case of a DOS file on a UNIX system
			return true;
		}
		if ( i++ == 100 ) break;
	}
	if ( i == 0 ) return true;		// Empty file
	return false;
}
bool isAPLFile ( const string& fullPath )
{
	GenIFStream ist ( fullPath );
	string line;
	int i = 0;
	while ( getline ( ist, line ) ) {
		if ( isPrefix ( line, "peaklist start" ) ) {			// isPrefix deals with the case of a DOS file on a UNIX system
			return true;
		}
		if ( i++ == 100 ) break;
	}
	return false;
}
bool isAPLFile2 ( const string& fullPath )	// This used for files with an apl suffix. It allows files with no scans (just comments).
{
	GenIFStream ist ( fullPath );
	string line;
	int i = 0;
	while ( getline ( ist, line ) ) {
		if ( line.length () == 0 || line [0] == '#' ) continue;
		if ( isPrefix ( line, "peaklist start" ) ) {				// isPrefix deals with the case of a DOS file on a UNIX system
			return true;
		}
		if ( i++ == 100 ) break;
	}
	if ( i == 0 ) return true;		// Empty file
	return false;
}
bool isDTAFile ( const string& fullPath )
{
/*
File format
precursor_m_plus_h charge
m_over_z_1 intensity_1
m_over_z_2 intensity_2
etc
separated by at least 1 blank line
may not be any fragment ions
Eg:
1131.05 2
242.9 4
273.299 4
305.999 8
313.001 4
316.8 6
*/
	GenIFStream ist ( fullPath );
	string line;
	int i = 0;
	while ( getline ( ist, line ) ) {			// Look for precursor line
		if ( !genEmptyString ( line ) ) {		// This should be a spectrum
			if ( isPrecursorMPlusHChargeLine ( line ) ) {	// First line OK
				if ( getline ( ist, line ) ) {
					if ( isFragmentMOverZIntensityLine ( line ) ) return true;	// 2nd line is a fragment line
					else {
						if ( genEmptyString ( line ) ) {		// Only other valid thing here is a blank line
							while ( getline ( ist, line ) ) {	// Look for next non-blank line
								if ( !genEmptyString ( line ) ) {
									if ( isPrecursorMPlusHChargeLine ( line ) ) return true; // First line of next spectrum
									else return false;
								}
							}
							return true;		// File contains single precursor line followed by blank lines
						}
						else return false;
					}
				}
				else return true;				// File just contains a single precursor line
			}
			else return false;
		}
	}
	return false;	// Only blank line
}
bool isPKLFile ( const string& fullPath )
{
/*
File format
precursor_m_over_z intensity charge
m_over_z_1 intensity_1
m_over_z_2 intensity_2
etc
separated by at least 1 blank line
may not be any fragment ions
Eg:
446.1455 869.1702 1
90.9821 5.2086
149.0615 7.0952
151.0510 3.2608
158.9667 9.1224
172.9683 3.2086
*/
	GenIFStream ist ( fullPath );
	string line;
	int i = 0;
	while ( getline ( ist, line ) ) {			// Look for precursor line
		if ( !genEmptyString ( line ) ) {		// This should be a spectrum
			if ( isPrecursorMOverZIntensityChargeLine ( line ) ) {	// First line OK
				if ( getline ( ist, line ) ) {
					if ( isFragmentMOverZIntensityLine ( line ) ) return true;	// 2nd line is a fragment line
					else {
						if ( genEmptyString ( line ) ) {		// Only other valid thing here is a blank line
							while ( getline ( ist, line ) ) {	// Look for next non-blank line
								if ( !genEmptyString ( line ) ) {
									if ( isPrecursorMOverZIntensityChargeLine ( line ) ) return true; // First line of next spectrum
									else return false;
								}
							}
							return true;		// File contains single precursor line followed by blank lines
						}
						else return false;
					}
				}
				else return true;				// File just contains a single precursor line
			}
			else return false;
		}
	}
	return false;	// Only blank line
}
bool isPPSingleFile ( const string& fullPath )
{
/*
>M1$Spot:002$Run:1$JobRunItem:178368$
>M2$Spot:097$Run:2$JobRunItem:178595$
*/
	const char* contains [] = { "$Spot:", "$Run:", "$JobRunItem:", 0 };
	GenIFStream ist ( fullPath );
	string line;
	while ( getline ( ist, line ) ) {
		string start = ">M";
		if ( line.compare ( 0, start.length (), start ) ) return false;
		string::size_type idx = 0;
		for ( int i = 0 ; contains [i] = 0 ; i++ ) {
			idx = line.find ( contains [i], idx );
			if ( idx == string::npos ) return false;
			idx++;
		}
		if ( line.compare ( line.length () - 1, 1, "$" ) ) return false;
		return true;
	}
	return false;
}
bool isXMLFile ( const string& fullPath, const string& type )
{
	string checkStr = "<" + type;
	GenIFStream ist ( fullPath );
	string line;
	if ( getline ( ist, line ) ) {
		if ( line.find ( "<?xml" ) != string::npos ) {	// first line states this is an xml file
			while ( getline ( ist, line ) ) {
				if ( line.find ( "<!" ) == string::npos && line.find ( "<?" ) == string::npos ) {
					if ( line.find ( type ) != string::npos ) return true;
					else return false;
				}
			}
		}
	}
	return false;
}
bool isMZDataFile ( const string& fullPath )
{
	return isXMLFile ( fullPath, MZDATA );
}
bool isMZMLFile ( const string& fullPath )
{
	return isXMLFile ( fullPath, "indexedmzML" ) || isXMLFile ( fullPath, "mzML" );
}
bool isMZXMLFile ( const string& fullPath )
{
	return isXMLFile ( fullPath, MZXML );
}
bool isMascotSearchResultsFile ( const string& fullPath )
{
	return isXMLFile ( fullPath, "mascot_search_results" );
}
bool isPepXMLFile ( const string& fullPath )
{
	return isXMLFile ( fullPath, "msms_pipeline_analysis" );
}
bool isPrideXMLFile ( const string& fullPath )
{
	return isXMLFile ( fullPath, "ExperimentCollection" );
}

bool isWiffFile ( const string& fullPath )
{
	return compareMicrosoftCompoundDocumentFormatPattern ( fullPath );
}
bool isFinniganRawFile ( const string& fullPath )
{
	unsigned char pattern [] = { 1, 0241, 'F', 0, 'i', 0, 'n', 0 };
	return compareMagicNumber ( fullPath, pattern, 8 );
}
bool isBMPFile ( const string& fullPath )
{
	unsigned char pattern [] = { 'B', 'M' };
	return compareMagicNumber ( fullPath, pattern, 2 );
}
bool isRTFFile ( const string& fullPath )
{
	unsigned char pattern [] = { '{', '\\', 'r', 't', 'f' };
	return compareMagicNumber ( fullPath, pattern, 5 );
}
bool isExcelFile ( const string& fullPath )
{
	return compareMicrosoftCompoundDocumentFormatPattern ( fullPath );
}
bool isExcelXFile ( const string& fullPath )
{
	return comparePKZIPPattern ( fullPath );
}
bool isWordFile ( const string& fullPath )
{
	return compareMicrosoftCompoundDocumentFormatPattern ( fullPath );
}
bool isWordXFile ( const string& fullPath )
{
	return comparePKZIPPattern ( fullPath );
}
bool isPDFFile ( const string& fullPath )
{
	unsigned char pattern [] = { '%', 'P', 'D', 'F' };
	return compareMagicNumber ( fullPath, pattern, 4 );
}
bool isSQLiteFile ( const string& fullPath )
{
	unsigned char pattern [] = { 'S','Q','L','i','t','e',' ','f','o','r','m','a','t',' ','3' };
	return compareMagicNumber ( fullPath, pattern, 4 );
}
bool isMSPFile ( const string& fullPath )
{
	return isFileType ( fullPath, MSP );
}
bool isSPTXTFile ( const string& fullPath )
{
	return isFileType ( fullPath, SPTXT );
}
bool isBSCFile ( const string& fullPath )
{
	GenIFStream ist ( fullPath );
	string line;
	int i = 0;
	while ( getline ( ist, line ) ) {
		if ( line.length () == 0 || line [0] == '#' ) continue;
		if ( isPrefix ( line, "Profile Spectrum" ) ) {	// isPrefix deals with the case of a DOS file on a UNIX system
			return true;
		}
		if ( i++ == 10 ) break;
	}
	if ( i == 0 ) return true;		// Empty file
	return false;
}
bool isCompressedUpload ( const string& f )
{
	return isFileType ( f, GZ ) || isFileType ( f, BZ2 ) || isFileType ( f, Z ) || isFileType ( f, CMN );
}
