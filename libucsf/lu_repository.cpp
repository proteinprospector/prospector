/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_repository.cpp                                             *
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
#ifdef VIS_C
#pragma warning( disable : 4503 )	// Disable warnings about decorated name length being too long
#endif
#ifndef VIS_C
#include <stdexcept>
#endif
#include <lg_stdlib.h>
#include <lg_string.h>
#include <lg_time.h>
#include <lgen_error.h>
#include <lgen_file.h>
#include <lu_getfil.h>
#include <lu_repository.h>
#include <lu_spec_id.h>
#include <lu_count_scan.h>
#include <lu_xml.h>
#include <lu_xml_data.h>
#include <lu_mass_conv.h>
#include <lu_file_type.h>
#include <lu_param_list.h>
#include <lu_repos_info.h>
using std::string;
using std::endl;
using std::getline;
using std::istream;
using std::istringstream;
using std::find_if;
using std::cout;
using std::ostream;
using std::runtime_error;
using namespace FileTypes;

namespace {
string getFractionName ( const string& fName, const string& type )
{
	return fName.substr ( 0, fName.length () - ( type.length () + 1 ));
}
}

Repository::Repository ( const string& baseDir, const string& searchName ) :
	baseDir ( baseDir ),
	searchName ( searchName )
{
}
string Repository::getProjectFileDataPath ( bool upload ) const
{
	string projType = upload ? "#" : "$";
	return projType + genTranslateSlash ( dataPath );
}

string UserRepository::uploadRepository = InfoParams::instance ().getUserRepository ();

UserRepository::UserRepository ( const string& searchName, const string& userDir ) :
	Repository ( uploadRepository, searchName )
{
	string monthString = genCurrentYearAndIntMonth ( '_' );
	string homePath;
	homePath += userDir [0];
	homePath += SLASH;
	homePath += userDir [1];
	homePath += SLASH;
	homePath += userDir;
	homePath += SLASH;
	homePath += searchName;
	dataPath	= homePath + SLASH + "data" + SLASH + monthString;
	projectPath	= homePath + SLASH + "project" + SLASH + monthString;
	resultsPath	= homePath + SLASH + "results" + SLASH + monthString;
	genCreateDirectoryPath ( getFullDataPath () );
	genCreateDirectoryPath ( getFullProjectPath () );
	genCreateDirectoryPath ( getFullResultsPath () );
}
class FractionFile {
	string centroidName;
	string fractionName;
	string rawName;
protected:
	int numMSSpectra;
	int numMSMSSpectra;
public:
	FractionFile ( const std::string& centroidName, const std::string& fractionName, const std::string& rawName = "" ) :
		centroidName ( centroidName ),
		fractionName ( fractionName ),
		rawName ( rawName ),
		numMSSpectra ( 0 ) {}
	string getCentroidName () const { return centroidName; }
	string getFractionName () const { return fractionName; }
	string getRawName () const { return rawName; }
	int getNumMSSpectra () const { return numMSSpectra; }
	int getNumMSMSSpectra () const { return numMSMSSpectra; }
	void printNumSpectraHTML ( ostream& os ) const;
};
void FractionFile::printNumSpectraHTML ( ostream& os ) const
{
	if ( numMSSpectra ) {
		if ( numMSSpectra == 1 ) os << fractionName << " 1 MS spectrum<br />" << endl;
		else os << fractionName << " " << numMSSpectra << " MS spectra<br />" << endl;
	}
	if ( numMSMSSpectra ) {
		if ( numMSMSSpectra == 1 ) os << fractionName << " 1 MSMS spectrum<br />" << endl;
		else os << fractionName << " " << numMSMSSpectra << " MSMS spectra<br />" << endl;
	}
}
class MGFFractionFile : public FractionFile {
public:
	MGFFractionFile ( const string& fName, const string& fractionName, const string& centroidName, const string& rawName = "" ) :
		FractionFile ( centroidName, fractionName, rawName )
	{
		MascotCountScans mcs ( fName );
		numMSMSSpectra = mcs.getCount ();
	}
};

class MS2FractionFile : public FractionFile {
public:
	MS2FractionFile ( const string& fName, const string& fractionName, const string& centroidName, const string& rawName = "" ) :
		FractionFile ( centroidName, fractionName, rawName )
	{
		MS2CountScans ms2cs ( fName );
		numMSMSSpectra = ms2cs.getCount ();
	}
};

class APLFractionFile : public FractionFile {
public:
	APLFractionFile ( const string& fName, const string& fractionName, const string& centroidName, const string& rawName = "" ) :
		FractionFile ( centroidName, fractionName, rawName )
	{
		APLCountScans aplcs ( fName );
		numMSMSSpectra = aplcs.getCount ();
	}
};

class PPSFractionFile : public FractionFile {
public:
	PPSFractionFile ( const string& fName, const string& fractionName, const string& centroidName, const string& rawName = "" ) :
		FractionFile ( centroidName, fractionName, rawName )
	{
		PPCountMSScans ppms ( fName );
		numMSSpectra = ppms.getCount ();
		PPCountMSMSScans ppmsms ( fName );
		numMSMSSpectra = ppmsms.getCount ();
	}
};

class XMLFractionFile : public FractionFile {
public:
	XMLFractionFile ( const string& fName, const string& fractionName, const string& centroidName, const string& rawName = "" ) :
		FractionFile ( centroidName, fractionName, rawName )
	{
		XMLCountScans* xcs = 0;
		if ( isMZXMLFile ( fName ) )
			xcs = new XMLCountScans ( new PPExpatMZXMLCountScans, fName );
		else if ( isMZMLFile ( fName ) )
			xcs = new XMLCountScans  ( new PPExpatMZMLCountScans, fName );
		else if ( isMZDataFile ( fName ) )
			xcs = new XMLCountScans  ( new PPExpatMZDataCountScans, fName );
		else if ( isMascotSearchResultsFile ( fName ) )
			xcs = new XMLCountScans  ( new PPExpatMascotSearchResultsCountScans, fName );
		if ( xcs ) {
			numMSMSSpectra = xcs->getCount ();
			delete xcs;
		}
	}
};
PPProject::PPProject ( const string& name, const StringVector& fileList ) :
	name ( name ),															// Command line constructor
	init ( false ),
	totalMSSpectra ( 0 ),
	totalMSMSSpectra ( 0 ),
	maxMSMSSpectra ( 0 ),
	rawDataAllowed ( false ),
	duplicateScans ( false ),
	spottingPlate ( false ),
	deleteFlag ( false ),
	upload ( false )
{
	nFiles = 0;
	for ( StringVectorSizeType i = 0 ; i < fileList.size () ; i++ ) {
		string fullPath = fileList [i];
		string f = genFilenameFromPath ( fullPath );
		if ( !parseCentroidFile ( fullPath, f ) ) return;
		nFiles++;
	}
	if ( !parseCentroidFiles ( fileList ) ) return;
	if ( !createFractionFiles ( "" ) ) return;
	init = true;
}
PPProject::PPProject ( const string& name, const string& dir ) :
	name ( name ),															// Command line constructor
	init ( false ),
	totalMSSpectra ( 0 ),
	totalMSMSSpectra ( 0 ),
	maxMSMSSpectra ( 0 ),
	rawDataAllowed ( false ),
	duplicateScans ( false ),
	spottingPlate ( false ),
	deleteFlag ( false ),
	upload ( false )
{
	FileList fileList ( dir, "", "", false );
	nFiles = 0;
	for ( int i = 0 ; i < fileList.size () ; i++ ) {
		string f = fileList [i];
		if ( f [0] == '.' ) continue;					// Ignore dot files
		if ( f.find ( "__MACOSX" ) != std::string::npos )  continue;
		if ( !parseFile ( dir + SLASH + f, f ) ) return;
		nFiles++;
	}
	if ( !dtas.empty () ) {		// bunch dtas by fraction
		if ( !parseDTAs ( dir ) ) {
			if ( deleteFlag ) return;	// Delete flag has been set - known error
			if ( !parseDTAs2 ( dir ) ) return;	// Try parsing the dtas as single files
		}
	}
	if ( !pkls.empty () ) {		// convert pkls to mgf
		if ( !parsePKLs ( dir ) ) {
			if ( deleteFlag ) return;	// Delete flag has been set - known error
			if ( !parsePKLs2 ( dir ) ) return;	// Try parsing the pkls as single files
		}
	}
	if ( !apls.empty () ) {		// convert apls to mgf
		if ( !parseAPLs ( dir ) ) return;
	}
	if ( !xmls.empty () ) {	// convert xmls to mgf
		if ( !parseXMLs ( dir, xmls ) ) return;
	}
	if ( !wiffs.empty () ) {
		if ( !parseWiffFiles ( dir ) ) return;
	}
	else if ( !finnRaws.empty () ) {
		if ( !parseThermoFiles () ) return;
	}
	else if ( !t2ds.empty () ) {
		if ( !parseT2DFiles () ) return;
	}
	else {
		if ( !parseCentroidFiles () ) return;
	}
	if ( !createFractionFiles ( dir ) ) return;
	init = true;
}
PPProject::PPProject ( unsigned int maxMSMSSpectra, bool rawDataAllowed, const string& name, const StringVector& files ) :
	name ( name ),															// Constructor for repository
	init ( false ),
	totalMSSpectra ( 0 ),
	totalMSMSSpectra ( 0 ),
	maxMSMSSpectra ( maxMSMSSpectra ),
	rawDataAllowed ( rawDataAllowed ),
	duplicateScans ( false ),
	spottingPlate ( false ),
	deleteFlag ( false ),
	upload ( false )
{
	string centBaseDir = InfoParams::instance ().getCentroidDir ();
	nFiles = 0;
	for ( StringVectorSizeType i = 0 ; i < files.size () ; i++ ) {
		string f = genFilenameFromPath ( files [i] );
		if ( f [0] == '.' ) continue;					// Ignore dot files
		string fullPath = centBaseDir + SLASH + files [i];
		if ( !parseCentroidFile ( fullPath, files [i] ) ) return;
		nFiles++;
	}
	if ( !parseCentroidFiles () ) return;
	if ( !createFractionFiles () ) return;
	init = true;
}
PPProject::PPProject ( unsigned int maxMSMSSpectra, bool rawDataAllowed, const string& name, const string& uploadName, const string& searchKey, const Repository* reposit, const IntVector& chargeRange ) :
	name ( name ),
	init ( false ),
	totalMSSpectra ( 0 ),
	totalMSMSSpectra ( 0 ),
	maxMSMSSpectra ( maxMSMSSpectra ),
	rawDataAllowed ( rawDataAllowed ),
	duplicateScans ( InfoParams::instance ().getBoolValue ( "duplicate_scans" ) ),
	chargeRange ( chargeRange ),
	spottingPlate ( false ),
	deleteFlag ( false ),
	upload ( true )
{
	FileList fileList ( uploadName, "", "", false );
	nFiles = 0;
	for ( int i = 0 ; i < fileList.size () ; i++ ) {
		string f = fileList [i];
		if ( f [0] == '.' ) continue;					// Ignore dot files
		if ( f.find ( "__MACOSX" ) != std::string::npos )  continue;
		if ( !parseFile ( uploadName + SLASH + f, f ) ) return;
		nFiles++;
	}
	if ( !dtas.empty () ) {		// bunch dtas by fraction
		if ( !parseDTAs ( uploadName ) ) {
			if ( deleteFlag ) return;	// Delete flag has been set - known error
			if ( !parseDTAs2 ( uploadName ) ) return;	// Try parsing the dtas as single files
		}
	}
	if ( !pkls.empty () ) {		// convert pkls to mgf
		if ( !parsePKLs ( uploadName ) ) {
			if ( deleteFlag ) return;	// Delete flag has been set - known error
			if ( !parsePKLs2 ( uploadName ) ) return;	// Try parsing the pkls as single files
		}
	}
	if ( !apls.empty () ) {		// convert apls to mgf
		if ( !parseAPLs ( uploadName ) ) return;
	}
	if ( !xmls.empty () ) {	// convert xmls to mgf
		if ( !parseXMLs ( uploadName, xmls ) ) return;
	}
	if ( !wiffs.empty () ) {
		if ( !parseWiffFiles ( uploadName ) ) return;
	}
	else if ( !finnRaws.empty () ) {
		if ( !parseThermoFiles () ) return;
	}
	else if ( !t2ds.empty () ) {
		if ( !parseT2DFiles () ) return;
	}
	else {
		if ( !parseCentroidFiles () ) return;
	}
	if ( !createFractionFiles ( uploadName, searchKey, reposit ) ) return;
	init = true;
}
PPProject::PPProject ( const string& uploadName ) :			// MS-Viewer
	name ( "" ),
	init ( false ),
	totalMSSpectra ( 0 ),
	totalMSMSSpectra ( 0 ),
	maxMSMSSpectra ( 0 ),
	rawDataAllowed ( false ),
	duplicateScans ( false ),
	spottingPlate ( false ),
	deleteFlag ( false ),
	upload ( true )
{
	FileList fileList ( uploadName, "", "", false );
	nFiles = 0;
	for ( int i = 0 ; i < fileList.size () ; i++ ) {
		string f = fileList [i];
		if ( f [0] == '.' ) continue;					// Ignore dot files
		if ( f.find ( "__MACOSX" ) != std::string::npos )  continue;
		if ( !parseFile ( uploadName + SLASH + f, f ) ) return;
		nFiles++;
	}
	if ( !dtas.empty () ) {		// bunch dtas by fraction
		if ( !parseDTAs ( uploadName ) ) {
			if ( deleteFlag ) return;	// Delete flag has been set - known error
			if ( !parseDTAs2 ( uploadName ) ) return;	// Try parsing the dtas as single files
		}
	}
	if ( !pkls.empty () ) {		// convert pkls to mgf
		if ( !parsePKLs ( uploadName ) ) {
			if ( deleteFlag ) return;	// Delete flag has been set - known error
			if ( !parsePKLs2 ( uploadName ) ) return;	// Try parsing the pkls as single files
		}
	}
	if ( !apls.empty () ) {		// convert apls to mgf
		if ( !parseAPLs ( uploadName ) ) return;
	}
	if ( !xmls.empty () ) {	// convert xmls to mgf
		if ( !parseXMLs ( uploadName, xmls ) ) return;
	}
	if ( !parseCentroidFiles () ) return;
	init = true;
}
void PPProject::createProjectFile ( const string& filename, bool fullPath )
{
	string projType;
	if ( !fullPath ) {
		projType = upload ? "#" : "$";
	}
	GenOFStream os ( filename );
	printXMLHeader ( os );
	printXMLVersion ( os );
	os << "<project>" << endl;
		ParameterList::printXML ( os, "project_name", name );
		if ( spottingPlate ) ParameterList::printXML ( os, "spotting_plate", true );
		for ( VectorFractionFilePtrSizeType i = 0 ; i < fractionFiles.size () ; i++ ) {
			os << "<file>" << endl;
				ParameterList::printXML ( os, "centroid", projType + genTranslateSlash ( fractionFiles [i]->getCentroidName () ) );
				if ( fractionFiles [i]->getNumMSSpectra () != 0 ) ParameterList::printXML ( os, "num_ms_spectra", fractionFiles [i]->getNumMSSpectra () );
				ParameterList::printXML ( os, "num_msms_spectra", fractionFiles [i]->getNumMSMSSpectra () );
				ParameterList::printXML ( os, "centroid_name", fractionFiles [i]->getFractionName () );
				if ( fractionFiles [i]->getRawName () != "" ) {
					ParameterList::printXML ( os, "raw", projType + genTranslateSlash ( fractionFiles [i]->getRawName () ) );
				}
			os << "</file>" << endl;
		}
	os << "</project>" << endl;
}
bool PPProject::parseFile ( const string& fullPath, const string& f )
{
	if ( genIsDirectory ( fullPath ) ) {			// Sub directories only allowed if they contain t2d files
		if ( !rawDataAllowed ) {
			deleteFlag = true;
			errMessage = "RAW data uploads are not currently enabled.";
			return false;
		}
		FileList cal ( fullPath, "", "cal", false );
		FileList t2d ( fullPath, "", "t2d", false );
		FileList all ( fullPath, "", "", false );
		int nall = all.size () - 2; // Remove . and ..
		int nt2d = t2d.size ();
		int ncal = cal.size ();
		if ( nall && ( nt2d == ncal ) && ( nt2d + ncal == nall ) ) {
			t2ds.push_back ( f );
#ifndef VIS_C
			for ( int i = 0 ; i < ncal ; i++ ) {
				convertToUnixText ( fullPath + SLASH + cal [i] );
			}
#endif
		}
		else {
			if ( nall == 0 ) {
				errMessage = "The uploaded archive contains an empty subdirectory.";
			}
			else if ( nt2d == 0 ) {
				errMessage = "Subdirectories are only allowed in uploaded archives of ABI TOFTOF data.";
			}
			else if ( nt2d != ncal ) {
				errMessage = "There is an unmatched number of t2d and cal files in one of the subdirectories.";
			}
			else if ( nt2d + ncal != nall ) {
				errMessage = "No other files should be uploaded along with the t2d and cal files.";
			}
			deleteFlag = true;
			return false;
		}
	}
	else if ( isFileType ( f, WIFF_SCAN ) ) {
		if ( !rawDataAllowed ) {
			deleteFlag = true;
			errMessage = "RAW data uploads are not currently enabled.";
			return false;
		}
		wiffScans.push_back ( f );
	}
	else if ( isFileType ( f, WIFF ) ) {
		if ( !rawDataAllowed ) {
			deleteFlag = true;
			errMessage = "RAW data uploads are not currently enabled.";
			return false;
		}
		if ( isWiffFile ( fullPath ) ) wiffs.push_back ( f );
		else {
			errMessage = "One or more files with a wiff extension is not a wiff file.";
			deleteFlag = true;
			return false;
		}
	}
	else if ( isFileType ( f, RAW ) ) {
		if ( !rawDataAllowed ) {
			deleteFlag = true;
			errMessage = "RAW data uploads are not currently enabled.";
			return false;
		}
		if ( isFinniganRawFile ( fullPath ) ) finnRaws.push_back ( f );
		else {
			errMessage = "One or more files with a RAW extension is not a Thermo RAW file.";
			deleteFlag = true;
			return false;
		}
	}
	else if ( isFileType ( f, BMP ) && isBMPFile ( fullPath ) ) {
		deleteFlag = true;
		errMessage = "BMP files cannot be processed.";
		return false;
	}
	else if ( isFileType ( f, RTF ) && isRTFFile ( fullPath ) ) {
		deleteFlag = true;
		errMessage = "RTF files cannot be processed.";
		return false;
	}
	else if ( isFileType ( f, PDF ) && isPDFFile ( fullPath ) ) {
		deleteFlag = true;
		errMessage = "PDF files cannot be processed.";
		return false;
	}
	else if ( (isFileType ( f, XLS ) && isExcelFile ( fullPath )) || (isFileType ( f, XLSX ) && isExcelXFile ( fullPath )) ) {
		deleteFlag = true;
		errMessage = "EXCEL files cannot be processed.";
		return false;
	}
	else if ( (isFileType ( f, DOC ) && isWordFile ( fullPath )) || (isFileType ( f, DOCX ) && isWordXFile ( fullPath )) ) {
		deleteFlag = true;
		errMessage = "Word files cannot be processed.";
		return false;
	}
	else if ( isFileType ( f, HTM ) || isFileType ( f, HTML ) ) {
		deleteFlag = true;
		errMessage = "HTML files cannot be processed.";
		return false;
	}
	else if ( isFileType ( f, BSC ) && isBSCFile ( fullPath ) ) {
		deleteFlag = true;
		errMessage = "BSC files cannot be processed as the precursor m/z is not specified.";
		return false;
	}
	else
		return parseCentroidFile ( fullPath, f );
	return true;
}
bool PPProject::parseCentroidFile ( const string& fullPath, const string& f )
{
	bool mgfFlag = isFileType ( f, MGF );
	bool dtaFlag = isFileType ( f, DTA );
	bool pklFlag = isFileType ( f, PKL );
	bool ms2Flag = isFileType ( f, MS2 );
	bool aplFlag = isFileType ( f, APL );
#ifndef VIS_C
	// All currently supported centroid files are text files. On a UNIX platform we need to check
	// for dos or mac text files and do the necessary conversion.
	if ( upload ) {
		convertToUnixTextWithMessage ( fullPath );
	}
#endif
	if ( mgfFlag ) {
		if ( isMGFFile2 ( fullPath ) ) {
			if ( mgfFileHasZeros ( fullPath ) ) {
				removeMGFZeroIntensities ( fullPath );
			}
			mgfs.push_back ( f );
		}
		else {
			deleteFlag = true;
			errMessage = "One or more files with an mgf extension is not an mgf file.";
			return false;
		}
	}
	else if ( dtaFlag ) {
		static bool dtaChecked = false;
		if ( dtaChecked || isDTAFile ( fullPath ) ) {
			dtas.push_back ( f );
			dtaChecked = true;
		}
		else {
			deleteFlag = true;
			errMessage = "One or more files with a dta extension is not a dta file.";
			return false;
		}
	}
	else if ( ms2Flag ) {
		if ( isMS2File2 ( fullPath ) ) ms2s.push_back ( f );
		else {
			deleteFlag = true;
			errMessage = "One or more files with an ms2 extension is not an ms2 file.";
			return false;
		}
	}
	else if ( aplFlag ) {
		if ( isAPLFile2 ( fullPath ) ) apls.push_back ( f );
		else {
			deleteFlag = true;
			errMessage = "One or more files with an apl extension is not an apl file.";
			return false;
		}
	}
	else if ( pklFlag ) {
		if ( isPKLFile ( fullPath ) ) pkls.push_back ( f );
		else {
			deleteFlag = true;
			errMessage = "One or more files with a pkl extension is not a pkl file.";
			return false;
		}
	}
	else if ( isFileType ( f, MZDATA ) || isFileType ( f, MZML ) || isFileType ( f, MZXML ) ) {
		xmls.push_back ( f );
	}
	else if ( isFileType ( f, XML ) ) {
		if ( isMZDataFile ( fullPath ) || isMZMLFile ( fullPath ) || isMZXMLFile ( fullPath ) || isMascotSearchResultsFile ( fullPath ) ) {
			xmls.push_back ( f );
		}
		else {
			deleteFlag = true;
			errMessage = "One or more files with an xml extension does not have a recognized format.";
			return false;
		}
	}
	else {
		if ( isTextFile ( fullPath, 256 ) ) {		// Check first 256 bytes to see if it is a text file
			if ( isPPSingleFile ( fullPath ) ) ppsfs.push_back ( f );
			else if ( isMGFFile ( fullPath ) ) mgfs.push_back ( f );
			else if ( isPKLFile ( fullPath ) ) pkls.push_back ( f );
			else if ( isDTAFile ( fullPath ) ) dtas.push_back ( f );
			else if ( isAPLFile ( fullPath ) ) apls.push_back ( f );
			else {
				deleteFlag = true;
				errMessage = "The file " + f + " does not have a recognized format.";
				return false;
			}
		}
		else {
			deleteFlag = true;
			errMessage = "The file " + f + " does not have a recognized format.";
			return false;
		}
	}
	return true;
}
bool PPProject::getDTAPKLFractionMap ( MapStringToStringVector& fNameMap, const StringVector& fNames )
{
	for ( StringVectorSizeType i = 0 ; i < fNames.size () ; i++ ) {
		string d = fNames [i];
		string::size_type idx = d.length ();
		for ( int j = 0 ; j < 4 ; j++ ) {
			idx = d.rfind ( '.', idx-1 );
			if ( idx == string::npos ) {
				return false;
			}
		}
		string fraction = d.substr ( 0, idx );
		fNameMap [fraction].push_back ( d );
	}
	return true;
}
bool PPProject::parseDTAs ( const string& uploadName )
{
	MapStringToStringVector fNameMap;
	if ( !getDTAPKLFractionMap ( fNameMap, dtas ) ) return false;

	for ( MapStringToStringVectorConstIterator j = fNameMap.begin () ; j != fNameMap.end () ; j++ ) {
		StringVector files = (*j).second;
		string mgfFile = (*j).first + "." + MGF;
		ErrorHandler::genError ()->message ( "Creating file " + mgfFile + ".\n" );
		GenOFStream osf ( uploadName + SLASH + mgfFile );
		for ( StringVectorSizeType k = 0 ; k < files.size () ; k++ ) {
			string f = uploadName + SLASH + files [k];
			GenIFStream ifs ( f );
			string line;
			while ( getline ( ifs, line ) ) {
				if ( !genEmptyString ( line ) ) {		// Skip any more blank lines
					istringstream ist ( line );
					double mhPlus;
					int charge;
					if ( ist >> mhPlus ) {
						if ( ( ist >> charge ) == 0 ) charge = 1;
						double mOZ = mPlusHToMOverZ ( mhPlus, charge, true );
						processDTAPKLHeader ( osf, files [k], mOZ, charge, 0.0 );
						break;
					}
					else {
						deleteFlag = true;
						errMessage = "The file " + f + " does not have a recognized format.";
						return false;
					}
				}
			}
			while ( getline ( ifs, line ) ) {
				osf << line << endl;
			}
			osf << "END IONS" << endl;
			if ( k < files.size () - 1 ) osf << endl;
			ifs.close ();
			genUnlink ( f );
		}
		mgfs.push_back ( mgfFile );
	}
	nFiles += fNameMap.size () - dtas.size ();
	return true;
}
bool PPProject::parseDTAs2 ( const string& uploadName )
{
	for ( StringVectorSizeType i = 0 ; i < dtas.size () ; i++ ) {
		string f = dtas [i];
		string mgfFile = f.substr ( 0, f.length () - 3 ) + "mgf";
		GenOFStream osf ( uploadName + SLASH + mgfFile );
		GenIFStream ifs ( uploadName + SLASH + f );
		int scan = 1;
		for ( ; ; ) {		// Loop through the spectra
			string line;
			bool flag = false;
			while ( getline ( ifs, line ) ) {
				if ( !genEmptyString ( line ) ) {		// Skip any more blank lines
					istringstream ist ( line );
					double mhPlus;
					int charge;
					if ( ist >> mhPlus ) {
						if ( ( ist >> charge ) == 0 ) charge = 1;
						double mOZ = mPlusHToMOverZ ( mhPlus, charge, true );
						processDTAPKLHeader ( osf, "Scan " + gen_itoa ( scan ), mOZ, charge, 0.0 );
						flag = true;
						scan++;
						break;
					}
					else {
						deleteFlag = true;
						errMessage = "The file " + f + " is not a valid dta file.";
						return false;
					}
				}
			}
			if ( processDTAPKLFragmentIons ( osf, ifs, flag ) ) break;
		}
		mgfs.push_back ( mgfFile );
		ifs.close ();
		genUnlink ( uploadName + SLASH + f );
	}
	dtas.resize ( 0 );
	return true;
}
bool PPProject::parsePKLs ( const string& uploadName )
{
	MapStringToStringVector fNameMap;
	if ( !getDTAPKLFractionMap ( fNameMap, pkls ) ) return false;

	for ( MapStringToStringVectorConstIterator j = fNameMap.begin () ; j != fNameMap.end () ; j++ ) {
		StringVector files = (*j).second;
		string mgfFile = (*j).first + "." + MGF;
		ErrorHandler::genError ()->message ( "Creating file " + mgfFile + ".\n" );
		GenOFStream osf ( uploadName + SLASH + mgfFile );
		for ( StringVectorSizeType k = 0 ; k < files.size () ; k++ ) {
			string f = uploadName + SLASH + files [k];
			GenIFStream ifs ( f );
			string line;
			while ( getline ( ifs, line ) ) {
				if ( !genEmptyString ( line ) ) {		// Skip any more blank lines
					istringstream ist ( line );
					double mOZ;
					double intensity;
					int charge;
					if ( ist >> mOZ ) {
						if ( ist >> intensity ) {
							if ( ( ist >> charge ) == 0 ) charge = 1;
							processDTAPKLHeader ( osf, files [k], mOZ, charge, intensity );
							break;
						}
					}
				}
			}
			while ( getline ( ifs, line ) ) {
				osf << line << endl;
			}
			osf << "END IONS" << endl;
			if ( k < files.size () - 1 ) osf << endl;
			ifs.close ();
			genUnlink ( f );
		}
		mgfs.push_back ( mgfFile );
	}
	nFiles += fNameMap.size () - pkls.size ();
	return true;
}
bool PPProject::parsePKLs2 ( const string& uploadName )
{
	for ( StringVectorSizeType i = 0 ; i < pkls.size () ; i++ ) {
		string f = pkls [i];
		string mgfFile = f.substr ( 0, f.length () - 3 ) + "mgf";
		GenOFStream osf ( uploadName + SLASH + mgfFile );
		GenIFStream ifs ( uploadName + SLASH + f );
		int scan = 1;
		for ( ; ; ) {		// Loop through the spectra
			string line;
			bool flag = false;
			while ( getline ( ifs, line ) ) {
				if ( !genEmptyString ( line ) ) {		// Skip any more blank lines
					istringstream ist ( line );
					double mOZ;
					double intensity;
					int charge;
					if ( ist >> mOZ ) {
						if ( ist >> intensity ) {
							if ( ( ist >> charge ) == 0 ) charge = 1;
							processDTAPKLHeader ( osf, "Scan " + gen_itoa ( scan ), mOZ, charge, intensity );
							flag = true;
							scan++;
							break;
						}
					}
				}
			}
			if ( processDTAPKLFragmentIons ( osf, ifs, flag ) ) break;
		}
		mgfs.push_back ( mgfFile );
		ifs.close ();
		genUnlink ( uploadName + SLASH + f );
	}
	pkls.resize ( 0 );
	return true;
}
void PPProject::processDTAPKLHeader ( ostream& osf, const string& titleInfo, double mOZ, int charge, double intensity )
{
	osf << "BEGIN IONS" << endl;
	osf << "TITLE=" << titleInfo << endl;
	osf << "PEPMASS=";
	genPrint ( osf, mOZ, 4 );
	if ( intensity ) {
		osf << " ";
		genPrintSigFig ( osf, intensity, 4 );
	}
	osf << endl;
	if ( charge ) osf << "CHARGE=" << charge << "+" << endl;
}
bool PPProject::processDTAPKLFragmentIons ( ostream& osf, istream& ifs, bool flag )
{
	int c;
	string line;
	while ( ( c = ifs.peek () ) != EOF ) {	// Keep going until you come to a line with only blanks or the end of file.
		if ( !getline ( ifs, line ) ) break;
		if ( genEmptyString ( line ) ) break;
		if ( isdigit ( line [0] ) ) osf << line << endl;
	}
	if ( flag ) osf << "END IONS" << endl;
	if ( c == EOF ) return true;
	return false;
}
bool PPProject::parseAPLs ( const string& uploadName )
{
	for ( StringVectorSizeType i = 0 ; i < apls.size () ; i++ ) {
		string f = apls [i];
		ErrorHandler::genError ()->message ( "Processing input file " + f + ".\n" );
		parseAPLFile ( uploadName, f );
	}
	FileList fl ( uploadName, "", "mgf", false );
	mgfs = fl.getNameList ();
	apls.resize ( 0 );
	nFiles = mgfs.size ();
	return true;
}
void PPProject::parseAPLFile ( const string& uploadName, const string& f )
{
	GenIFStream ifs ( uploadName + SLASH + f );
	for ( ; ; ) {
		string line;
		double mz;
		int z;
		string title;
		string file;
		parseAPLHeaderLines ( ifs, mz, z, title, file );
		string curPath = uploadName + SLASH + file + ".mgf";
		//if ( !genFileExists ( curPath ) ) mgfs.push_back ( file + ".mgf" );
		GenOFStream ost ( curPath, std::ios_base::out | std::ios_base::app );
		ost << "BEGIN IONS" << endl;
		ost << "TITLE=" << title << endl;
		ost << "PEPMASS=";
		genPrint ( ost, mz, 4 );
		ost << endl;
		ost << "CHARGE=" << z << "+" << endl;

		int c;
		int numP = 0;
		while ( ( c = ifs.peek () ) != EOF ) {
			if ( c == 'p' ) numP++;
			if ( numP == 2 ) break;
			getline ( ifs, line );
			if ( line.length () != 0 ) {
				if ( c != 'p' ) {
					if ( line [line.length ()-1] == '\r' ) line = line.substr ( 0, line.length ()-1 );
					if ( line.length () != 0 ) ost << line << endl;
				}
			}
		}
		ost << "END IONS" << endl;
		if ( numP == 1 ) break;
	}
	ifs.close ();
	genUnlink ( uploadName + SLASH + f );
}
void PPProject::parseAPLHeaderLines ( istream& istr, double& mz, int& z, string& title, string& file )
{
	int c;
	string line;
	int numS = 0;
	while ( ( c = istr.peek () ) != EOF ) {
		if ( isdigit ( c ) ) {
			break;			// data section
		}
		if ( numS == 1 && c == 'p' ) {  //		2nd line, no data
			break;			// data section
		}
		getline ( istr, line );
		if ( line.length () != 0 ) {
			c = line [0];
			if ( c == 'm' ) {					// eg mz=1065.41113255839
				mz = atof ( line.substr ( 3 ).c_str () );
			}
			else if ( c == 'c' ) {				// charge=4
				z = atoi ( line.substr ( 7 ).c_str () );
			}
			else if ( c == 'h' ) {
				title = line.substr ( 7 );
				if ( title [title.length ()-1] == '\r' ) title = title.substr ( 0, title.length ()-1 );
				int start = line.find ( "RawFile: " ) + 9;
				int end = line.find ( "Index: " ) - 1;
				file = line.substr ( start, end - start );
			}
		}
	}
}
bool PPProject::parseXMLs ( const string& uploadName, StringVector& xmlFiles )
{
	for ( StringVectorSizeType i = 0 ; i < xmlFiles.size () ; i++ ) {
		string f = xmlFiles [i];
		string fullPath = uploadName + SLASH + f;
		MSMSDataPointVector msmsDataPointList;
		int fraction = -1;
		SpecID specID;
		specID.setFraction ( -1 );
		PPExpat* ppe = getPPExpatPtr ( fullPath, msmsDataPointList, fraction, specID );
		if ( ppe == 0 ) {
			errMessage = "Format problem with the file " + f;
			deleteFlag = true;
			return false;
		}
		try {
			ppe->parseXMLFromFile ( fullPath );
		}
		catch ( runtime_error e ) {
			deleteFlag = true;
			errMessage = string ( e.what () ) + "\nFilename: " + f + "\n";
			return false;
		}
		string mgfFile = getXMLFractionName ( f ) + ".mgf";
		GenOFStream osf ( uploadName + SLASH + mgfFile );
		for ( MSMSDataPointVectorSizeType j = 0 ; j < msmsDataPointList.size () ; j++ ) {
			if ( duplicateScans ) {
				MSMSDataPoint& msms = msmsDataPointList [j];
				int z = msms.getPrecursorCharge ();
				if ( z == -9999 ) {
					for ( IntVectorSizeType k = 0 ; k < chargeRange.size () ; k++ ) {
						writeMGFScan ( osf, msmsDataPointList [j], chargeRange [k] );
					}
				}
				else
					writeMGFScan ( osf, msmsDataPointList [j], z );
			}
			else
				writeMGFScan ( osf, msmsDataPointList [j], 0 );
		}
		mgfs.push_back ( mgfFile );
		delete ppe;
		genUnlink ( fullPath );
	}
	xmlFiles.resize ( 0 );
	return true;
}
string PPProject::getXMLFractionName ( const string& f )
{
	int len = 0;
	if ( isFileType ( f, XML ) ) {
		len = 4;
		if ( isFileType ( f, "mzxml.xml" ) )		len = 10;
		else if ( isFileType ( f, "mzml.xml" ) )	len = 9;
		else if ( isFileType ( f, "mzdata.xml" ) )	len = 11;
	}
	else if ( isFileType ( f, MZDATA ) )	len = 7;
	else if ( isFileType ( f, MZXML ) )		len = 6;
	else if ( isFileType ( f, MZML ) )		len = 5;
	return f.substr ( 0, f.length () - len );
}
PPExpat* PPProject::getPPExpatPtr ( const string& fullPath, MSMSDataPointVector& msmsDataPointList, int fraction, const SpecID& specID ) const
{
	PPExpat* ppe = 0;
	if ( isMZXMLFile ( fullPath ) )
		ppe = new PPExpatMZXMLData ( msmsDataPointList, fraction, SpectrumRange (), specID, -9999, IntVector () );
	else if ( isMZMLFile ( fullPath ) )
		ppe = new PPExpatMZMLData ( msmsDataPointList, fraction, SpectrumRange (), specID, -9999, IntVector () );
	else if ( isMZDataFile ( fullPath ) )
		ppe = new PPExpatMZDataData ( msmsDataPointList, fraction, SpectrumRange (), specID, -9999, IntVector () );
	else if ( isMascotSearchResultsFile ( fullPath ) )
		ppe = new PPExpatMascotSearchResults ( msmsDataPointList, fraction, SpectrumRange (), specID, -9999, IntVector () );
	return ppe;
}
void PPProject::writeMGFScan ( ostream& os, MSMSDataPoint& msms, int charge )
{
	os << "BEGIN IONS" << endl;

	os << "TITLE=Scan ";
	os << msms.getMSMSInfo ();
	os << " ";
	string spot = msms.getSpot ();
	if ( !spot.empty () ) {
		os << "(rt=";
		os << spot;
		os << ")";
		os << " ";
	}
	os << "[Prospector Created]";
	os << endl;

	os << "PEPMASS=";
	genPrint ( os, msms.getPrecursorMZ (), 4 );
	os << " ";
	genPrintSigFig ( os, msms.getPrecursorIntensity (), 3 );
	os << endl;
	if ( charge == 0 ) {
		int z = msms.getPrecursorCharge ();
		if ( z != -9999 ) os << "CHARGE=" << z << "+" << endl;
	}
	else {
		os << "CHARGE=" << charge << "+" << endl;
	}
	DataFilePeakVector& pks = msms.getDataPeaks ();
	for ( DataFilePeakVectorSizeType k = 0 ; k < pks.size () ; k++ ) {
		genPrint ( os, pks [k].getMOverZ (), 4 );
		os << " ";
		genPrintSigFig ( os, pks [k].getIntensity (), 3 );
		if ( pks [k].getCharge () != 1 )	os << " " << pks [k].getCharge ();
		os << endl;
	}
	os << "END IONS" << endl;
}
bool PPProject::parseWiffFiles ( const string& uploadName )
{
	int nWiffs = wiffs.size ();
	int nWiffScans = wiffScans.size ();
	if ( mgfs.empty () ) {	// Create mgf files
		if ( nWiffs != nFiles && nWiffs != nWiffScans ) {
			deleteFlag = true;
			errMessage = "All files must be wiff files or wiff files and wiff scan files.";
			return false;	// All files must be wiff files or wiff files and wiff scan files
		}
#ifdef VIS_C 
		for ( int i = 0 ; i < nWiffs ; i++ ) {
			string w = wiffs [i];
			string fractionName = getFractionName ( w, WIFF );
			string c = fractionName + "." + MGF;
			string command ( getSystemCall ( "wiffToCentroid.exe" ) );
			command += " ";
			command += "\"" + uploadName + SLASH + w +"\"";
			command += " ";
			command += "\"" + uploadName + SLASH + c +"\"";
			ErrorHandler::genError ()->message ( "Centroiding " + w + ".\n" );
			genSystem ( command, "", true );
			if ( !genFileExists ( uploadName + SLASH + c ) ) {
				deleteFlag = true;
				errMessage = "Centroid file " + c + " not created.";
				return false;
			}
			mgfs.push_back ( c );
		}
		ErrorHandler::genError ()->message ( "Centroiding complete.\n" );
		nFiles += nWiffs;
#else
		deleteFlag = true;
		errMessage = "You must upload peak list files along with wiff files.";
		return false;
#endif
	}
	if ( nWiffs != mgfs.size () ) {
		deleteFlag = true;
		errMessage = "There have to be equal numbers of wiff and peak list files.";
		return false; // Unbalanced number of wiff files/mgf files
	}
	else {
		if ( nWiffScans ) {
			if ( nWiffs != nFiles / 3 ) {
				deleteFlag = true;
				errMessage = "There have to be equal numbers of wiff, wiff.scan and peak list files.";
				return false;	// Some left over files
			}
		}
		else {
			if ( nWiffs != nFiles / 2 ) {
				deleteFlag = true;
				errMessage = "There have to be equal numbers of wiff and peak list files.";
				return false;	// Some left over files
			}
		}
		if ( nWiffs == 1 ) {						// In this case the centroid file name can be different to the raw file name
			rawFiles.push_back ( wiffs [0] );
			fractionNames.push_back ( getFractionName ( wiffs [0], WIFF ) );
			centroidFiles.push_back ( mgfs [0] );
		}
		else {
			for ( int i = 0 ; i < nWiffs ; i++ ) {
				rawFiles.push_back ( wiffs [i] );
				fractionNames.push_back ( getFractionName ( wiffs [i], WIFF ) );
				string centroidName = fractionNames.back () + "." + MGF;
				StringVectorIterator iter = find_if ( mgfs.begin (), mgfs.end (), CheckStrcasecmpEquality ( centroidName ) );
				if ( iter == mgfs.end () ) {
					deleteFlag = true;
					errMessage = "mgf file for " + wiffs [i] + " not found.";
					return false;
				}
				centroidFiles.push_back ( mgfs [iter-mgfs.begin()] );
			}
		}
	}
	return true;
}
bool PPProject::parseThermoFiles ()
{
	StringVector c = mgfs;		// The centroid files can either be mgf or mzxml
	string type = MGF;
	if ( c.empty () ) {
		c = xmls;
		type = MZXML;
	}
	int nFinnRaws = finnRaws.size ();
	if ( c.empty () ) {	// Create mgf files
		if ( nFinnRaws != nFiles ) return false;	// All files must be raw files to continue
		deleteFlag = true;
		errMessage = "You must upload peak list files along with RAW files.";
		return false;								// Not implemented yet
	}
	if ( nFinnRaws != c.size () ) {
		errMessage = "There have to be equal numbers of RAW and peak list files.";
		return false;
	}
	else {
		if ( nFinnRaws != nFiles / 2 ) {
			errMessage = "No other files should be uploaded along with the RAW and peak list files.";
			return false;// Some left over files
		}
		if ( nFinnRaws == 1 ) {						// In this case the centroid file name can be different to the raw file name
			rawFiles.push_back ( finnRaws [0] );
			fractionNames.push_back ( getFractionName ( finnRaws [0], RAW ) );
			centroidFiles.push_back ( c [0] );
		}
		else {
			for ( int i = 0 ; i < nFinnRaws ; i++ ) {
				rawFiles.push_back ( finnRaws [i] );
				fractionNames.push_back ( getFractionName ( finnRaws [i], RAW ) );
				string centroidName = fractionNames.back () + "." + type;
				StringVectorIterator iter = find_if ( c.begin (), c.end (), CheckStrcasecmpEquality ( centroidName ) );
				if ( iter == c.end () ) return false;
				centroidFiles.push_back ( c [iter-c.begin()] );
			}
		}
	}
	return true;
}
bool PPProject::parseT2DFiles ()
{
	int nt2ds = t2ds.size ();
	if ( nt2ds != ppsfs.size () ) return false;
	else {
		if ( nt2ds != nFiles / 2 ) return false;// Some left over files
		if ( nt2ds == 1 ) {						// In this case the centroid file name can be different to the raw file name
			rawFiles.push_back ( t2ds [0] );
			fractionNames.push_back ( t2ds [0] );
			centroidFiles.push_back ( ppsfs [0] );
		}
		else {
			for ( int i = 0 ; i < nt2ds ; i++ ) {
				rawFiles.push_back ( t2ds [i] );
				fractionNames.push_back ( t2ds [i] );
				string centroidName = fractionNames.back () + "." + TXT;
				StringVectorIterator iter = find_if ( ppsfs.begin (), ppsfs.end (), CheckStrcasecmpEquality ( centroidName ) );
				if ( iter == ppsfs.end () ) return false;
				centroidFiles.push_back ( ppsfs [iter-ppsfs.begin()] );
			}
		}
	}
	return true;
}
bool PPProject::parseCentroidFiles ()
{
	if ( !mgfs.empty () ) {					// Just mgf files
		if ( !parseCentroidFiles ( mgfs ) ) return false;
	}
	else if ( !xmls.empty () ) {			// Just xmls files
		if ( !parseCentroidFiles ( xmls ) ) return false;
	}
	else if ( !ppsfs.empty () ) {			// Just ppsf files
		if ( !parseCentroidFiles ( ppsfs ) ) return false;
	}
	else if ( !ms2s.empty () ) {			// Just ms2 files
		if ( !parseCentroidFiles ( ms2s ) ) return false;
	}
	else if ( !apls.empty () ) {			// Just apl files
		if ( !parseCentroidFiles ( apls ) ) return false;
	}
	else return false;
	return true;
}
bool PPProject::parseCentroidFiles ( const StringVector& files )
{
	int num = files.size ();
	if ( num != nFiles ) return false;
	for ( int i = 0 ; i < num ; i++ ) {
		fractionNames.push_back ( genShortFilenameFromPath ( files [i] ) );
		centroidFiles.push_back ( files [i] );
	}
	return true;
}
bool PPProject::createFractionFiles ( const string& uploadName, const string& searchKey, const Repository* reposit ) // upload
{
	for ( StringVectorSizeType i = 0 ; i < fractionNames.size () ; i++ ) {
		string fName = uploadName + SLASH + centroidFiles [i];
		string basePath = reposit->getDataPath () + SLASH + searchKey + SLASH;
		string centroidPath = basePath + centroidFiles [i];
		string rawPath;
		if ( !rawFiles.empty () ) {
			rawPath = basePath + rawFiles [i];
		}
		string frac = fractionNames [i];
		if ( !mgfs.empty () ) {
			if ( i == 0 ) spottingPlate = isMGFSpottingPlateFile ( fName );
			addFractionFile ( new MGFFractionFile ( fName, frac, centroidPath, rawPath ) );
		}
		if ( !ms2s.empty () )	addFractionFile ( new MS2FractionFile ( fName, frac, centroidPath, rawPath ) );
		if ( !apls.empty () )	addFractionFile ( new APLFractionFile ( fName, frac, centroidPath, rawPath ) );
		if ( !xmls.empty () )	addFractionFile ( new XMLFractionFile ( fName, frac, centroidPath, rawPath ) );
		if ( !ppsfs.empty () ) {
			if ( i == 0 ) spottingPlate = true;
			addFractionFile ( new PPSFractionFile ( fName, frac, centroidPath, rawPath ) );
		}
		totalMSSpectra += fractionFiles.back ()->getNumMSSpectra ();
		totalMSMSSpectra += fractionFiles.back ()->getNumMSMSSpectra ();
		if ( maxMSMSSpectra && ( totalMSMSSpectra > maxMSMSSpectra ) ) {
			deleteFlag = true;
			errMessage = "The maximum number of MSMS spectra have been exceeded.";
			return false;
		}
		fractionFiles.back ()->printNumSpectraHTML ( cout );
	}
	printTotalSpectraHTML ( cout );
	return true;
}
bool PPProject::createFractionFiles () // repository
{
	string centBaseDir = InfoParams::instance ().getCentroidDir ();
	string rawBaseDir = InfoParams::instance ().getRawDir ();
	for ( StringVectorSizeType i = 0 ; i < fractionNames.size () ; i++ ) {
		string cent = centroidFiles [i];
		string centroidFullPath = centBaseDir + SLASH + cent;
		string frac = fractionNames [i];
		if ( RepositoryInfo::exists () ) {
			try {
				string suffix = RepositoryInfo::instance ().getCentroidSuffix ( cent );
				if ( !suffix.empty () ) {
					frac = frac.substr ( 0, frac.length () - suffix.length () );
				}
			}
			catch ( runtime_error e ) {
				ErrorHandler::genError ()->error ( e );
			}
		}
		string rPath = centroidFiles [i];
		if ( RepositoryInfo::exists () ) rPath = RepositoryInfo::instance ().getAdjustedRawPath ( rPath );
		string dir = genDirectoryFromPath ( rawBaseDir + SLASH + rPath );
		FileList fileList ( dir, frac, "", false );
		string rawPath;
		for ( int j = 0 ; j < fileList.size () ; j++ ) {
			string f = fileList [j];
			string fFrac;
			if ( genIsDirectory ( dir + SLASH + f ) ) fFrac = f;
			else if ( isFileType ( f, WIFF ) ) fFrac = getFractionName ( f, WIFF );
			else if ( isFileType ( f, RAW ) ) fFrac = getFractionName ( f, RAW );
			if ( frac == fFrac ) {
				rawPath = genDirectoryFromPath ( rPath ) + "/" + f;
				break;
			}
		}
		if ( RepositoryInfo::exists () && rawPath.empty () ) {
			string fileType = RepositoryInfo::instance ().getRawType ( centroidFiles [i] );
			if ( !fileType.empty () ) {
				rawPath = genDirectoryFromPath ( rPath ) + "/" + frac + "." + fileType;
			}
		}
		if ( !mgfs.empty () ) {
			if ( i == 0 ) spottingPlate = isMGFSpottingPlateFile ( centroidFullPath );
			addFractionFile ( new MGFFractionFile ( centroidFullPath, frac, cent, rawPath ) );
		}
		if ( !ms2s.empty () )	addFractionFile ( new MS2FractionFile ( centroidFullPath, frac, cent, rawPath ) );
		if ( !apls.empty () )	addFractionFile ( new APLFractionFile ( centroidFullPath, frac, cent, rawPath ) );
		if ( !xmls.empty () )	addFractionFile ( new XMLFractionFile ( centroidFullPath, frac, cent, rawPath ) );
		if ( !ppsfs.empty () ) {
			if ( i == 0 ) spottingPlate = true;
			addFractionFile ( new PPSFractionFile ( centroidFullPath, frac, cent, rawPath ) );
		}
		totalMSSpectra += fractionFiles.back ()->getNumMSSpectra ();
		totalMSMSSpectra += fractionFiles.back ()->getNumMSMSSpectra ();
		if ( maxMSMSSpectra && ( totalMSMSSpectra > maxMSMSSpectra ) ) return false;
		fractionFiles.back ()->printNumSpectraHTML ( cout );
	}
	printTotalSpectraHTML ( cout );
	return true;
}
bool PPProject::createFractionFiles ( const string& dir ) // command line
{
	for ( StringVectorSizeType i = 0 ; i < fractionNames.size () ; i++ ) {
		string fName;
		if ( !dir.empty () ) fName += dir + SLASH;
		fName += centroidFiles [i];
		string frac = fractionNames [i];
		if ( !mgfs.empty () ) {
			if ( i == 0 ) spottingPlate = isMGFSpottingPlateFile ( fName );
			addFractionFile ( new MGFFractionFile ( fName, frac, fName, "" ) );
		}
		if ( !ms2s.empty () )	addFractionFile ( new MS2FractionFile ( fName, frac, fName, "" ) );
		if ( !apls.empty () )	addFractionFile ( new APLFractionFile ( fName, frac, fName, "" ) );
		if ( !xmls.empty () )	addFractionFile ( new XMLFractionFile ( fName, frac, fName, "" ) );
		if ( !ppsfs.empty () ) {
			if ( i == 0 ) spottingPlate = true;
			addFractionFile ( new PPSFractionFile ( fName, frac, fName, "" ) );
		}
		totalMSSpectra += fractionFiles.back ()->getNumMSSpectra ();
		totalMSMSSpectra += fractionFiles.back ()->getNumMSMSSpectra ();
		fractionFiles.back ()->printNumSpectraHTML ( cout );
	}
	printTotalSpectraHTML ( cout );
	return true;
}
void PPProject::printTotalSpectraHTML ( ostream& os ) const
{
	os << "<br />" << endl;
	if ( totalMSSpectra != 0 ) {
		if ( totalMSSpectra == 1 ) os << "Total " << "1 MS spectrum<br />" << endl;
		else os << "Total " << totalMSSpectra << " MS spectra<br />" << endl;
	}
	if ( totalMSMSSpectra != 0 ) {
		if ( totalMSMSSpectra == 1 ) os << "Total " << "1 MSMS spectrum<br />" << endl;
		else os << "Total " << totalMSMSSpectra << " MSMS spectra<br />" << endl;
	}
}
void PPProject::removeMGFZeroIntensities ( const string& name )
{
	ErrorHandler::genError ()->message ( "Removing zero intensities from mgf file " + name + ".\n" );
	string name2 = name + ".new";
	removeMGFZeroIntensities ( name, name2 );
	genUnlink ( name );
	genRename ( name2, name );
}
void PPProject::removeMGFZeroIntensities ( const string& name1, const string& name2 )
{
	GenIFStream ist ( name1 );
	GenOFStream osf ( name2 );
	string line;
	int i = 0;
	while ( getline ( ist, line ) ) {
		if ( line.length () != 0 && isdigit ( line [0] ) ) {							// This line has data
			char* ptr = const_cast <char*> (line.c_str ());
			double mass = atof ( ptr );
			while ( *ptr && isspace ( *ptr ) ) ptr++;	// Skip leading white space
			while ( *ptr && !isspace ( *ptr ) ) ptr++;	// Skip mass
			if ( atof ( ptr ) == 0.0 ) {
				continue;								// Don't output this line
			}
		}
		osf << line << endl;
	}
}
