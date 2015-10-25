/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_file_type.h                                                *
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

#ifndef __lu_file_type_h
#define __lu_file_type_h

#include <string>

namespace FileTypes {
	const std::string WIFF = "wiff";
	const std::string WIFF_SCAN = "wiff.scan";
	const std::string MGF = "mgf";
	const std::string MS2 = "ms2";
	const std::string APL = "apl";
	const std::string DTA = "dta";
	const std::string PKL = "pkl";
	const std::string BSC = "bsc";
	const std::string RAW = "raw";
	const std::string MZXML = "mzXML";
	const std::string MZDATA = "mzData";
	const std::string MZML = "mzML";
	const std::string XML = "xml";
	const std::string TXT = "txt";
	const std::string T2D = "t2d";
	const std::string BMP = "bmp";
	const std::string RTF = "rtf";
	const std::string PDF = "pdf";
	const std::string XLS = "xls";
	const std::string XLSX = "xlsx";
	const std::string DOC = "doc";
	const std::string DOCX = "docx";
	const std::string GZ = "gz";
	const std::string Z = "z";
	const std::string BZ2 = "bz2";
	const std::string Z7 = "7z";
	const std::string CMN = "cmn";
	const std::string HTM = "htm";
	const std::string HTML = "html";
	const std::string MSP = "msp";
	const std::string SPTXT = "sptxt";
	const std::string MSF = "msf";
	const std::string PEPXML = "pep.xml";
}

bool isMGFSpottingPlateFile ( const std::string& fullPath );
bool isMGFFile ( const std::string& fullPath );
bool isMS2File ( const std::string& fullPath );
bool isAPLFile ( const std::string& fullPath );
bool isMGFFile2 ( const std::string& fullPath );
bool mgfFileHasZeros ( const std::string& fullPath );
bool isMS2File2 ( const std::string& fullPath );
bool isAPLFile2 ( const std::string& fullPath );
bool isDTAFile ( const std::string& fullPath );
bool isPKLFile ( const std::string& fullPath );
bool isPPSingleFile ( const std::string& fullPath );
bool isXMLFile ( const std::string& fullPath, const std::string& type );
bool isMZDataFile ( const std::string& fullPath );
bool isMZMLFile ( const std::string& fullPath );
bool isMZXMLFile ( const std::string& fullPath );
bool isMascotSearchResultsFile ( const std::string& fullPath );
bool isPepXMLFile ( const std::string& fullPath );
bool isPrideXMLFile ( const std::string& fullPath );
bool isWiffFile ( const std::string& fullPath );
bool isFinniganRawFile ( const std::string& fullPath );
bool isBMPFile ( const std::string& fullPath );
bool isRTFFile ( const std::string& fullPath );
bool isExcelFile ( const std::string& fullPath );
bool isExcelXFile ( const std::string& fullPath );
bool isWordFile ( const std::string& fullPath );
bool isWordXFile ( const std::string& fullPath );
bool isSQLiteFile ( const std::string& fullPath );
bool isMSPFile ( const std::string& fullPath );
bool isSPTXTFile ( const std::string& fullPath );
bool isPDFFile ( const std::string& fullPath );
bool isBSCFile ( const std::string& fullPath );
bool isCompressedUpload ( const std::string& fullPath );

#endif /* ! __lu_file_type_h */
