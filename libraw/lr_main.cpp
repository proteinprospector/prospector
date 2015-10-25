/******************************************************************************
*                                                                             *
*  Program    : libraw                                                        *
*                                                                             *
*  Filename   : lr_main.cpp                                                   *
*                                                                             *
*  Created    : April 23rd 2007                                               *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2007-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifndef VIS_C
#include <stdexcept>
#endif
#include <lg_string.h>
#include <lgen_math.h>
#include <lgen_file.h>
#include <lgen_error.h>
#include <lu_df_info.h>
#include <lu_proj_file.h>
#include <lu_param_list.h>
#include <lu_spec_id.h>
#include <lr_main.h>
#include <lu_t2d.h>
#include <lu_file_type.h>
#ifdef XCALIBUR
#include <lx_raw.h>
#endif
#ifdef ANALYST
#include <la_wiff.h>
#endif
using std::string;
using std::vector;
using std::ostringstream;
using std::runtime_error;
using namespace FileTypes;

namespace {
StringVector getTimeStringsFromDouble ( const DoubleVector& t )
{
	StringVector sv;
	for ( DoubleVectorSizeType i = 0 ; i < t.size () ; i++ ) {
		sv.push_back ( gen_ftoa ( t [i], "%.3f" ) );
	}
	return sv;
}
StringVector getTimeStringsInMinutesFromDouble ( const DoubleVector& t )
{
	StringVector sv;
	for ( DoubleVectorSizeType i = 0 ; i < t.size () ; i++ ) {
		sv.push_back ( gen_ftoa ( t [i] / 60.0, "%.3f" ) );
	}
	return sv;
}
PairStringVectorString getStringTimesInMinutes ( const DoubleVector& t, const string& spot )
{
	StringVector sv = getTimeStringsInMinutesFromDouble ( t );
	int chosenIndex = genNearestIndex ( t, atof ( spot.c_str () ) * 60.0 );
	if ( sv [chosenIndex] != spot ) {
		sv.insert ( sv.begin () + chosenIndex, spot );
	}
	PairStringVectorString psvs;
	psvs.first = sv;
	psvs.second = sv [chosenIndex];
	return psvs;
}
}

class RawInstance {
public:
	virtual ~RawInstance ();
	virtual void getXYData ( vector <XYData>& vXYData, const SpecID& specID, double mOverZ, double startMass = 0.0, double endMass = 0.0 ) = 0;
	virtual void getXYData ( vector <XYData>& vXYData, PairStringVectorString& sTimes, bool msFullScan, const SpecID& specID, double mOverZ, double startMass = 0.0, double endMass = 0.0 ) = 0;
	virtual void getQuantitationXYData ( vector <XYData>& vXYData, PairStringVectorString& sTimes, bool ITRAQ, const SpecID& specID, double mOverZ, const string& version, double startMass = 0.0, double endMass = 0.0 ) = 0;
	virtual string getRawType () const = 0;
};
RawInstance::~RawInstance () {}

RawFile::~RawFile ()
{
	delete rawInstance;
}
void RawFile::getXYData ( vector <XYData>& vXYData, const SpecID& specID, double mOverZ, double startMass, double endMass )
{
	rawInstance->getXYData ( vXYData, specID, mOverZ, startMass, endMass );
}
void RawFile::getXYData ( vector <XYData>& vXYData, PairStringVectorString& sTimes, bool msFullScan, const SpecID& specID, double mOverZ, double startMass, double endMass )
{
	rawInstance->getXYData ( vXYData, sTimes, msFullScan, specID, mOverZ, startMass, endMass );
}
void RawFile::getQuantitationXYData ( vector <XYData>& vXYData, PairStringVectorString& sTimes, bool ITRAQ, const SpecID& specID, double mOverZ, const string& version, double startMass, double endMass )
{
	rawInstance->getQuantitationXYData ( vXYData, sTimes, ITRAQ, specID, mOverZ, version, startMass, endMass );
}
string RawFile::getRawType () const
{
	return rawInstance->getRawType ();
}

class CentroidInstance : public RawInstance {
	MSMSQuanDataSetInfo dsi;
public:
	CentroidInstance ( const string& fName );
	void getXYData ( vector <XYData>& vXYData, const SpecID& specID, double mOverZ, double startMass = 0.0, double endMass = 0.0 );
	void getXYData ( vector <XYData>& vXYData, PairStringVectorString& sTimes, bool msFullScan, const SpecID& specID, double mOverZ, double startMass = 0.0, double endMass = 0.0 );
	void getQuantitationXYData ( vector <XYData>& vXYData, PairStringVectorString& sTimes, bool ITRAQ, const SpecID& specID, double mOverZ, const string& version, double startMass = 0.0, double endMass = 0.0 );
	string getRawType () const;
};
CentroidInstance::CentroidInstance ( const string& fName ) :
	dsi ( fName )
{
}
void CentroidInstance::getXYData ( vector <XYData>& vXYData, const SpecID& specID, double mOverZ, double startMass, double endMass )
{
	// stub - not implemented yet
}
void CentroidInstance::getXYData ( vector <XYData>& vXYData, PairStringVectorString& sTimes, bool msFullScan, const SpecID& specID, double mOverZ, double startMass, double endMass )
{
	// stub - not implemented yet
}
void CentroidInstance::getQuantitationXYData ( vector <XYData>& vXYData, PairStringVectorString& sTimes, bool ITRAQ, const SpecID& specID, double mOverZ, const string& version, double startMass, double endMass )
{
	if ( ITRAQ ) {		// Only relevant if iTRAQ
		dsi.getMSMSQuanData ( vXYData [0], specID, startMass, endMass, version );
	}
}
string CentroidInstance::getRawType () const
{
	return "";		// return empty string - not a raw type
}

class T2DInstance : public RawInstance {
	string fName;
	static string getTOFTOFRawFilename ( const string& path, const SpecID& specID );
public:
	T2DInstance ( const string& fName ) :
		fName ( fName ) {}
	void getXYData ( vector <XYData>& vXYData, const SpecID& specID, double mOverZ, double startMass = 0.0, double endMass = 0.0 );
	void getXYData ( vector <XYData>& vXYData, PairStringVectorString& sTimes, bool msFullScan, const SpecID& specID, double mOverZ, double startMass = 0.0, double endMass = 0.0 );
	void getQuantitationXYData ( vector <XYData>& vXYData, PairStringVectorString& sTimes, bool ITRAQ, const SpecID& specID, double mOverZ, const string& version, double startMass = 0.0, double endMass = 0.0 );
	string getRawType () const;
};

string T2DInstance::getTOFTOFRawFilename ( const string& path, const SpecID& specID )
{
	int jobRunItemID = atoi ( specID.getMSMSInfo ().c_str () );
	int msmsRun = specID.getRun ();
	if ( jobRunItemID == 0 ) {		// ms
		ostringstream os;
		os << specID.getID ();
		os << "_";
		os << "MS";
		os << "_";
		string prefix = os.str ();
		FileList fList ( path, prefix, ".t2d", true );
		StringVector sv = fList.getNameList ();
		IntVector iv;
		for ( StringVectorSizeType i = 0 ; i < sv.size () ; i++ ) {
			iv.push_back ( atoi ( sv [i].substr ( prefix.length () ).c_str () ) );
		}
		int idx = getT2DMSRunNumberIndex ( msmsRun, iv );
		if ( idx == -1 ) {
			ErrorHandler::genError ()->error ( "No matching MS file.\n" );
			return "";
		}
		else {
			return path + SLASH + sv [idx];
		}
	}
	else {
		ostringstream os;
		os << path;
		os << SLASH;
		os << specID.getID ();
		os << "_";
		os << "MSMS";
		os << "_";
		os << jobRunItemID;
		os << "_";
		os << msmsRun;
		return os.str ();
	}
}
void T2DInstance::getXYData ( vector <XYData>& vXYData, const SpecID& specID, double mOverZ, double startMass, double endMass )
{
	T2DFile t2d ( getTOFTOFRawFilename ( fName, specID ), startMass, endMass );
	vXYData [0] = t2d.getData ();
}
void T2DInstance::getXYData ( vector <XYData>& vXYData, PairStringVectorString& sTimes, bool msFullScan, const SpecID& specID, double mOverZ, double startMass, double endMass )
{
	SpecID s = specID;
	if ( msFullScan ) {
		s.setMSMSInfo ( "0" );
	}
	T2DFile t2d ( getTOFTOFRawFilename ( fName, s ), startMass, endMass );
	vXYData [0] = t2d.getData ();
}
void T2DInstance::getQuantitationXYData ( vector <XYData>& vXYData, PairStringVectorString& sTimes, bool ITRAQ, const SpecID& specID, double mOverZ, const string& version, double startMass, double endMass )
{
	getXYData ( vXYData, sTimes, !ITRAQ, specID, mOverZ, startMass, endMass );
}
string T2DInstance::getRawType () const
{
	return T2D;
}
#ifdef ANALYST
class WiffInstance : public RawInstance {
	WiffFile wf;
	DoubleVector times;
public:
	WiffInstance ( const string& fName, double timeWindowStart, double timeWindowEnd ) :
		wf ( fName, timeWindowStart, timeWindowEnd ),
		times ( wf.getTimes () ) {}
	void getXYData ( vector <XYData>& vXYData, const SpecID& specID, double mOverZ, double startMass = 0.0, double endMass = 0.0 );
	void getXYData ( vector <XYData>& vXYData, PairStringVectorString& sTimes, bool msFullScan, const SpecID& specID, double mOverZ, double startMass = 0.0, double endMass = 0.0 );
	void getQuantitationXYData ( vector <XYData>& vXYData, PairStringVectorString& sTimes, bool ITRAQ, const SpecID& specID, double mOverZ, const string& version, double startMass = 0.0, double endMass = 0.0 );
	string getRawType () const;
};
void WiffInstance::getXYData ( vector <XYData>& vXYData, const SpecID& specID, double mOverZ, double startMass, double endMass )
{
	wf.getDataFromSpecList ( vXYData [0], specID.getMSMSInfo (), startMass, endMass );
}
void WiffInstance::getXYData ( vector <XYData>& vXYData, PairStringVectorString& sTimes, bool msFullScan, const SpecID& specID, double mOverZ, double startMass, double endMass )
{
	sTimes = getStringTimesInMinutes ( times, specID.getID () );
	if ( msFullScan )
		wf.getDataFromRT ( vXYData [0], specID.getID (), startMass, endMass );
	else
		wf.getDataFromSpecList ( vXYData [0], specID.getMSMSInfo (), startMass, endMass );
}
void WiffInstance::getQuantitationXYData ( vector <XYData>& vXYData, PairStringVectorString& sTimes, bool ITRAQ, const SpecID& specID, double mOverZ, const string& version, double startMass, double endMass )
{
	getXYData ( vXYData, sTimes, !ITRAQ, specID, mOverZ, startMass, endMass );
}
string WiffInstance::getRawType () const
{
	return WIFF;
}
#endif
#ifdef XCALIBUR
class ThermoRawInstance : public RawInstance {
	ThermoRawFile ixr;
	StringVector stringTimes;
	int startScan;
	int endScan;
public:
	ThermoRawInstance ( const string& fName, double timeWindowStart, double timeWindowEnd ) :
		ixr ( fName, timeWindowStart, timeWindowEnd ),
		stringTimes ( getTimeStringsFromDouble ( ixr.getTimes () ) ),
		startScan ( ixr.getStartScan () ),
		endScan ( ixr.getEndScan () )
	{
		if ( startScan == -1 ) {
			ostringstream err;
			err << "The start scan cannot be read from the file " << genFilenameFromPath ( fName ) << ".";
			throw runtime_error ( err.str () );
		}
	}
	void getXYData ( vector <XYData>& vXYData, const SpecID& specID, double mOverZ, double startMass, double endMass );
	void getXYData ( vector <XYData>& vXYData, PairStringVectorString& sTimes, bool msFullScan, const SpecID& specID, double mOverZ, double startMass, double endMass );
	void getQuantitationXYData ( vector <XYData>& vXYData, PairStringVectorString& sTimes, bool ITRAQ, const SpecID& specID, double mOverZ, const string& version, double startMass = 0.0, double endMass = 0.0 );
	string getRawType () const;
};
void ThermoRawInstance::getXYData ( vector <XYData>& vXYData, const SpecID& specID, double mOverZ, double startMass, double endMass )
{
	ixr.getXYData ( vXYData, specID.getMSMSInfo (), mOverZ, false );
}
void ThermoRawInstance::getXYData ( vector <XYData>& vXYData, PairStringVectorString& sTimes, bool msFullScan, const SpecID& specID, double mOverZ, double startMass, double endMass )
{
	sTimes.first = stringTimes;
	string msmsInfo = specID.getMSMSInfo ();
	if ( msmsInfo.empty () ) {
		ostringstream err;
		err << "Scan information not present in the centroid file.";
		throw runtime_error ( err.str () );
	}
	int index = atoi ( msmsInfo.c_str() );
	if ( index > stringTimes.size () ) {
		ostringstream err;
		err << "Scan number " << index << " from centroid file is greater than the number of scans in the run " << genFilenameFromPath ( ixr.getFilename () ) << ".";
		throw runtime_error ( err.str () );
	}
	sTimes.second = stringTimes [index-1];
	ixr.getXYData ( vXYData, msmsInfo, mOverZ, msFullScan );
}
void ThermoRawInstance::getQuantitationXYData ( vector <XYData>& vXYData, PairStringVectorString& sTimes, bool ITRAQ, const SpecID& specID, double mOverZ, const string& version, double startMass, double endMass )
{
	sTimes.first = stringTimes;
	string msmsInfo = specID.getMSMSInfo ();
	if ( msmsInfo.empty () ) {
		ostringstream err;
		err << "Scan information not present in the centroid file.";
		throw runtime_error ( err.str () );
	}
	int index = atoi ( msmsInfo.c_str() );
	if ( index > stringTimes.size () ) {
		ostringstream err;
		err << "Scan number " << index << " from centroid file is greater than the number of scans in the run " << genFilenameFromPath ( ixr.getFilename () ) << ".";
		throw runtime_error ( err.str () );
	}
	sTimes.second = stringTimes [index-1];
	if ( ITRAQ )
		ixr.getXYData ( vXYData, msmsInfo, mOverZ, false, false, startMass, endMass );
	else
		ixr.getXYData ( vXYData, msmsInfo, mOverZ, true, true, startMass, endMass );
}
string ThermoRawInstance::getRawType () const
{
	return RAW;
}
#endif
RawFile::RawFile ( const ParameterList* params, int fraction, double rtIntervalStart, double rtIntervalEnd )
{
	string filePath = getRawDataFilename ( params, fraction );
	if ( filePath.empty () ) {					// The raw file isn't specified in the project file.
		filePath = getCentroidDataFilename ( params, fraction );
		rawInstance = new CentroidInstance ( filePath );
		return;
	}
	if ( !genFileExists ( filePath ) ) {		// The raw file is specified in the project file but missing from disk.
		throw runtime_error ( "Raw data file " + filePath + " not present.\n" );
	}
	if ( genIsDirectory ( filePath ) )				// TOF-TOF
		rawInstance = new T2DInstance ( filePath );
#ifdef ANALYST
	else if ( isFileType ( filePath, WIFF ) )		// WIFF
		rawInstance = new WiffInstance ( filePath, rtIntervalStart, rtIntervalEnd );
#endif
#ifdef XCALIBUR
	else if ( isFileType ( filePath, RAW ) )		// RAW
		rawInstance = new ThermoRawInstance ( filePath, rtIntervalStart, rtIntervalEnd );
#endif
	else											// The raw file is specified in the project file but can't be read.
		throw runtime_error ( "Raw data display is not currently possible on the chosen instrument type." );
}
